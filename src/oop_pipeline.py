# Object-oriented repeat phasing pipeline
import gzip
from dataclasses import dataclass
from typing import Iterator, Tuple, List
from time import perf_counter
import numpy as np
from numba import njit, uint8, uint64
from Bio import SeqIO
from Bio.Seq import Seq


@njit
def myers_dist(text: uint8[:], PM: uint64[:], pattern_length: int, max_errors: int) -> int:
    """Return index of first match within `max_errors` edits using Myers's algorithm."""
    m = pattern_length
    Pv = uint64(-1)
    Mv = uint64(0)
    score = m
    for i in range(len(text)):
        b = text[i]
        Eq = PM[b] if b < 4 else uint64(0)
        Xv = Eq | Mv
        Xh = (((Eq & Pv) + Pv) ^ Pv) | Eq
        Ph = Mv | ~(Xh | Pv)
        Mh = Pv & Xh
        if (Ph >> (m - 1)) & 1:
            score += 1
        elif (Mh >> (m - 1)) & 1:
            score -= 1
        if score <= max_errors:
            return i
        Ph <<= 1
        Mh <<= 1
        Pv = Mh | ~(Xv | Ph)
        Mv = Ph & Xv
    return -1


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


# Cache 64-bit pattern masks so each distinct pattern is processed only once
PATTERN_CACHE = {}

def build_pattern_mask(pattern: bytes) -> np.ndarray:
    cached = PATTERN_CACHE.get(pattern)
    if cached is not None:
        return cached
    pm = [0, 0, 0, 0]
    for j, b in enumerate(pattern):
        idx = BASE2IDX[b]
        if idx < 4:
            pm[idx] |= 1 << j
    arr = np.array(pm, dtype=np.uint64)
    PATTERN_CACHE[pattern] = arr
    return arr


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
    """Quality-weighted consensus matrix with per-position coverage."""

    def __init__(self, d: int):
        self.d = d
        # Aggregate quality for each base (A,C,G,T) at each position 0..d-1
        self.quality_matrix = np.zeros((4, d), dtype=int)
        # Number of bases contributing to each position
        self.coverage_vector = np.zeros(d, dtype=int)

    def update(self, seq: str, qual: str, shift: int = 0):
        idx = encode_seq(seq)
        q = encode_qual(qual)
        positions = shift + np.arange(len(idx), dtype=np.int64)
        cols = positions % self.d

        valid = idx < 4
        if not np.any(valid):
            return

        np.add.at(self.quality_matrix, (idx[valid], cols[valid]), q[valid])
        np.add.at(self.coverage_vector, cols[valid], 1)

    def to_consensus(self) -> Tuple[str, List[int], np.ndarray, np.ndarray]:
        sum4 = self.quality_matrix
        best_idx = np.argmax(sum4, axis=0)
        total = np.sum(sum4, axis=0)
        best = sum4[best_idx, np.arange(self.d)]
        cons_q = best - (total - best)
        cons_seq = "".join(IDX2BASE[best_idx].tolist())
        return (
            cons_seq,
            cons_q.tolist(),
            self.quality_matrix,
            self.coverage_vector,
        )


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
        pattern_bytes = seq[:k].encode("ascii")
        text_idx = encode_seq(seq[k:]).astype(np.uint8)

        PM = build_pattern_mask(pattern_bytes)

        idx = myers_dist(text_idx, PM, k, max_e)
        if idx >= 0:
            return k + idx
        raise NoRepeats


# ----------------------- OrientationDecider ---------------------

class OrientationDecider:
    """Decide global orientation of R2 reads."""

    def __init__(self, detector: RepeatDetector):
        self.detector = detector

    def _score_orientation(self, cons_idx: np.ndarray, seq: str, qual: str, d: int):
        idx = encode_seq(seq)
        phred = encode_qual(qual)
        pos = np.arange(len(idx))
        cols = (pos[np.newaxis, :] + np.arange(d)[:, np.newaxis]) % d
        matches = (idx[np.newaxis, :] == cons_idx[cols])
        scores = (phred[np.newaxis, :] * matches).sum(axis=1)
        phi = int(np.argmax(scores))
        return int(scores[phi]), phi

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
    def best_shift(self, cons_idx: np.ndarray, read_idx: np.ndarray, read_q: np.ndarray, d: int) -> Tuple[int, int]:
        pos = np.arange(len(read_idx))
        cols = (pos[np.newaxis, :] + np.arange(d)[:, np.newaxis]) % d
        matches = (read_idx[np.newaxis, :] == cons_idx[cols])
        scores = (read_q[np.newaxis, :] * matches).sum(axis=1)
        phi = int(np.argmax(scores))
        return phi, int(scores[phi])

    def merge(self, matrix: ConsensusMatrix, seq: str, qual: str, shift: int):
        matrix.update(seq, qual, shift=shift)


# ----------------------- RepeatPhasingPipeline ------------------

@dataclass
class PhasingResult:
    read_id: str
    phase_shift: int
    consensus_len: int
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
            cons_idx = np.argmax(cm.quality_matrix, axis=0)

            seq2 = str(r2.seq)
            qual2 = "".join(chr(q + 33) for q in r2.letter_annotations["phred_quality"])
            if orientation == "RC":
                seq2 = str(Seq(seq2).reverse_complement())
                qual2 = qual2[::-1]
            read_idx = encode_seq(seq2)
            read_q = encode_qual(qual2)
            phi, _ = self.aligner.best_shift(cons_idx, read_idx, read_q, d)
            self.aligner.merge(cm, seq2, qual2, phi)
            elapsed = perf_counter() - start
            result = PhasingResult(
                read_id=r1.id,
                phase_shift=phi,
                consensus_len=d,
                elapsed=elapsed,
            )
            del cm
            yield result

