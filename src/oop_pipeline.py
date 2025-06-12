# Object-oriented repeat phasing pipeline
import gzip
from dataclasses import dataclass
from typing import Iterator, Tuple, List
from time import perf_counter

import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq


# ----------------------- Encoding helpers -----------------------

# ASCII base to index (A,C,G,T,N -> 0..4)
BASE2IDX = np.full(256, 4, dtype=np.int8)
BASE2IDX[ord('A')] = 0
BASE2IDX[ord('C')] = 1
BASE2IDX[ord('G')] = 2
BASE2IDX[ord('T')] = 3
BASE2IDX[ord('N')] = 4

IDX2BASE = np.array(list("ACGT"), dtype="<U1")


def encode_seq(seq: str) -> np.ndarray:
    ascii_vals = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return BASE2IDX[ascii_vals]


def encode_qual(qual: str) -> np.ndarray:
    return (np.frombuffer(qual.encode("ascii"), dtype=np.uint8) - 33).astype(int)


# ----------------------- FastqStream -----------------------------

class FastqStream:
    """Stream paired FASTQ records."""

    def __init__(self, r1_path: str, r2_path: str):
        self.r1_path = r1_path
        self.r2_path = r2_path

    def _open(self, path: str):
        return gzip.open(path, "rt") if path.endswith('.gz') else open(path, 'r')

    def __iter__(self) -> Iterator[Tuple[SeqIO.SeqRecord, SeqIO.SeqRecord]]:
        with self._open(self.r1_path) as f1, self._open(self.r2_path) as f2:
            it1 = SeqIO.parse(f1, "fastq")
            it2 = SeqIO.parse(f2, "fastq")
            for r1, r2 in zip(it1, it2):
                yield r1, r2

    def sample(self, n: int) -> List[Tuple[SeqIO.SeqRecord, SeqIO.SeqRecord]]:
        pairs = []
        with self._open(self.r1_path) as f1, self._open(self.r2_path) as f2:
            it1 = SeqIO.parse(f1, "fastq")
            it2 = SeqIO.parse(f2, "fastq")
            for idx, pair in enumerate(zip(it1, it2)):
                if idx >= n:
                    break
                pairs.append(pair)
        return pairs


# ----------------------- ConsensusMatrix ------------------------

class ConsensusMatrix:
    """Maintain a 3-D base/quality matrix for consensus building."""

    def __init__(self, d: int):
        self.d = d
        self.mat = np.zeros((4, d, 1), dtype=int)

    def _ensure_cycles(self, cycle: int):
        if cycle >= self.mat.shape[2]:
            pad = cycle + 1 - self.mat.shape[2]
            self.mat = np.pad(self.mat, ((0, 0), (0, 0), (0, pad)))

    def update(self, seq: str, qual: str, shift: int = 0):
        idx = encode_seq(seq)
        q = encode_qual(qual)
        positions = shift + np.arange(len(idx), dtype=np.int64)
        cols = positions % self.d
        cycles = positions // self.d
        max_cyc = int(cycles.max())
        self._ensure_cycles(max_cyc)
        for b in range(4):
            mask = idx == b
            if not mask.any():
                continue
            c = cycles[mask]
            col = cols[mask]
            w = q[mask]
            for cyc, col_id, weight in zip(c, col, w):
                self.mat[b, col_id, cyc] += int(weight)

    def to_consensus(self) -> Tuple[str, List[int], np.ndarray, int]:
        sum4 = self.mat.sum(axis=2)
        best_idx = np.argmax(sum4, axis=0)
        total = np.sum(sum4, axis=0)
        best = sum4[best_idx, np.arange(self.d)]
        cons_q = best - (total - best)
        cons_seq = "".join(IDX2BASE[best_idx].tolist())
        return cons_seq, cons_q.tolist(), self.mat, self.mat.shape[2]


# ----------------------- RepeatDetector -------------------------

class NoRepeats(Exception):
    pass


class RepeatDetector:
    def __init__(self, k: int = 20, max_errors: int = 2):
        self.k = k
        self.max_errors = max_errors

    def _compare(self, a: str, b: str, errors: int) -> bool:
        err = 0
        for x, y in zip(a, b):
            if x != y:
                err += 1
                if err > errors:
                    return False
        return True

    def find_repeat_distance(self, seq: str) -> int:
        k = self.k
        max_e = self.max_errors
        anchor = seq[0:k]
        max_range = len(seq) - k
        for i in range(max_range):
            cand = seq[k + i : 2 * k + i]
            if self._compare(anchor, cand, max_e):
                return k + i
        raise NoRepeats


# ----------------------- OrientationDecider ---------------------

class OrientationDecider:
    """Decide global orientation of R2 reads."""

    def __init__(self, detector: RepeatDetector):
        self.detector = detector

    def _score_orientation(self, cons_idx: np.ndarray, seq: str, qual: str, d: int):
        idx = encode_seq(seq)
        phred = encode_qual(qual)
        positions = np.arange(len(idx), dtype=np.int64)
        best_score = -1
        best_phi = 0
        for phi in range(d):
            cols = (phi + positions) % d
            mask = idx == cons_idx[cols]
            score = int(phred[mask].sum()) if mask.any() else 0
            if score > best_score:
                best_score, best_phi = score, phi
        return best_score, best_phi

    def decide_orientation(self, stream: FastqStream, sample_size: int = 1000) -> str:
        pairs = stream.sample(sample_size)
        forward = 0
        rc = 0
        for r1, r2 in pairs:
            seq1 = str(r1.seq)
            try:
                d = self.detector.find_repeat_distance(seq1)
            except NoRepeats:
                continue
            qual1 = "".join(chr(q + 33) for q in r1.letter_annotations["phred_quality"])
            cm = ConsensusMatrix(d)
            cm.update(seq1, qual1)
            cons_seq, cons_q, mat, _ = cm.to_consensus()
            cons_idx = encode_seq(cons_seq)

            seq2 = str(r2.seq)
            qual2 = "".join(chr(q + 33) for q in r2.letter_annotations["phred_quality"])
            s_f, _ = self._score_orientation(cons_idx, seq2, qual2, d)
            rc_seq = str(Seq(seq2).reverse_complement())
            rc_qual = qual2[::-1]
            s_rc, _ = self._score_orientation(cons_idx, rc_seq, rc_qual, d)
            if s_rc > s_f:
                rc += 1
            else:
                forward += 1
        return "RC" if rc > forward else "forward"


# ----------------------- PhaseAligner ---------------------------

class PhaseAligner:
    def best_shift(self, cons_idx: np.ndarray, read_idx: np.ndarray, read_q: np.ndarray, d: int) -> int:
        positions = np.arange(len(read_idx), dtype=np.int64)
        best_score = -1
        best_phi = 0
        for phi in range(d):
            cols = (phi + positions) % d
            mask = read_idx == cons_idx[cols]
            score = int(read_q[mask].sum()) if mask.any() else 0
            if score > best_score:
                best_score, best_phi = score, phi
        return best_phi

    def merge(self, matrix: ConsensusMatrix, seq: str, qual: str, shift: int):
        matrix.update(seq, qual, shift=shift)


# ----------------------- RepeatPhasingPipeline ------------------

@dataclass
class PhasingResult:
    read_id: str
    phase_shift: int
    consensus_len: int
    consensus: str
    qualities: List[int]
    matrix: np.ndarray
    elapsed: float


class RepeatPhasingPipeline:
    def __init__(self, r1_path: str, r2_path: str, sample_size: int = 1000, k: int = 40, max_errors: int = 2):
        self.stream = FastqStream(r1_path, r2_path)
        self.detector = RepeatDetector(k=k, max_errors=max_errors)
        self.decider = OrientationDecider(self.detector)
        self.aligner = PhaseAligner()
        self.sample_size = sample_size

    def run(self) -> Iterator[PhasingResult]:
        orientation = self.decider.decide_orientation(self.stream, self.sample_size)
        for r1, r2 in self.stream:
            start = perf_counter()
            seq1 = str(r1.seq)
            try:
                d = self.detector.find_repeat_distance(seq1)
            except NoRepeats:
                continue
            qual1 = "".join(chr(q + 33) for q in r1.letter_annotations["phred_quality"])
            cm = ConsensusMatrix(d)
            cm.update(seq1, qual1)
            cons_seq, _, _, _ = cm.to_consensus()
            cons_idx = encode_seq(cons_seq)

            seq2 = str(r2.seq)
            qual2 = "".join(chr(q + 33) for q in r2.letter_annotations["phred_quality"])
            if orientation == "RC":
                seq2 = str(Seq(seq2).reverse_complement())
                qual2 = qual2[::-1]
            read_idx = encode_seq(seq2)
            read_q = encode_qual(qual2)
            phi = self.aligner.best_shift(cons_idx, read_idx, read_q, d)
            self.aligner.merge(cm, seq2, qual2, phi)
            final_seq, final_q, _, _ = cm.to_consensus()
            elapsed = perf_counter() - start
            yield PhasingResult(
                read_id=r1.id,
                phase_shift=phi,
                consensus_len=d,
                consensus=final_seq,
                qualities=final_q,
                matrix=cm.mat.copy(),
                elapsed=elapsed,
            )

