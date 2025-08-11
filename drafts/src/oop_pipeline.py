"""Object-oriented repeat phasing pipeline.

The classes defined here break the repeat phasing logic into small reusable
components.  They are used by ``main.py`` to process paired-end FASTQ files and
produce phased consensus reads.
"""

import gzip
from dataclasses import dataclass
from typing import Iterator, Tuple, List
from time import perf_counter
from collections import defaultdict
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from numba import njit
import os

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
    """Convert an A/C/G/T string to numeric indices used by the matrix."""
    ascii_vals = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
    return BASE2IDX[ascii_vals]
    
def encode_qual(qual: str) -> np.ndarray:
    """Turn ASCII-encoded Phred qualities into integer scores."""
    return (np.frombuffer(qual.encode("ascii"), dtype=np.uint8) - 33).astype(int) # Return a 4-element mask array for Myers's algorithm, caching results.

# ----------------------- FastqStream -----------------------------

class FastqStream:
    """Stream paired FASTQ records."""

    def __init__(self, r1_path: str, r2_path: str):
        """Store file paths for later streaming."""
        self.r1_path = r1_path
        self.r2_path = r2_path

    def _open(self, path: str):
        """Open plain or gzipped FASTQ files transparently."""
        return gzip.open(path, "rt") if path.endswith('.gz') else open(path, 'r')

    def __iter__(self) -> Iterator[Tuple[SeqIO.SeqRecord, SeqIO.SeqRecord]]:
        """Yield pairs of R1/R2 records from the FASTQ streams."""
        with self._open(self.r1_path) as f1, self._open(self.r2_path) as f2:
            it1 = SeqIO.parse(f1, "fastq")
            it2 = SeqIO.parse(f2, "fastq")
            for r1, r2 in zip(it1, it2):
                yield r1, r2

    def sample(self, n: int) -> List[Tuple[SeqIO.SeqRecord, SeqIO.SeqRecord]]:
        """Return the first ``n`` pairs without consuming the iterator."""
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
        """Initialize empty matrices for a repeat unit of length ``d``."""
        self.d = d
        # Aggregate quality for each base (A,C,G,T) at each position 0..d-1
        self.quality_matrix = np.zeros((4, d), dtype=int)
        # Number of bases contributing to each position
        self.coverage_vector = np.zeros(d, dtype=int)

    def update(self, seq: str, qual: str, shift: int = 0):
        """Add a read to the matrix, optionally shifted by ``shift`` bases."""
        idx = encode_seq(seq)
        q = encode_qual(qual)
        self.positions = shift + np.arange(len(idx), dtype=np.int64)
        cols = self.positions % self.d

        valid = idx < 4
        if not np.any(valid):
            return

        np.add.at(self.quality_matrix, (idx[valid], cols[valid]), q[valid])
        np.add.at(self.coverage_vector, cols[valid], 1)
        

    def to_consensus(self) -> Tuple[str, List[int], np.ndarray, np.ndarray]:
        """Return consensus string, qualities, and underlying matrices."""
        sum4 = self.quality_matrix

        best_idx = np.argmax(sum4, axis=0)
        total = np.sum(sum4, axis=0)
        best = sum4[best_idx, np.arange(self.d)]
        cons_q = best - (total - best)
        cons_bases = IDX2BASE[best_idx]
        cons_seq = "".join(cons_bases.tolist())
        return cons_seq, cons_q.tolist(), self.quality_matrix, self.coverage_vector

    def consensus_length(self) -> int:
        """Return the fixed consensus length ``d``."""
        return self.d


# ----------------------- RepeatDetector -------------------------

class NoRepeats(Exception):
    pass


class RepeatDetector:
    """Find the repeat unit length using a k-mer search."""

    def __init__(self, k: int = 20, max_errors: int = 2):
        """Configure search window size and allowed mismatches."""
        self.k = k
        self.max_errors = max_errors

    def _compare(self, a: str, b: str, max_errors: int) -> bool:
        """Return ``True`` if ``a`` and ``b`` differ by at most ``errors`` bases."""
        err = 0
        for x, y in zip(a, b):
            if x != y:
                err += 1
                if err > max_errors:
                    return None
        return err

    def find_repeat_distance(self, seq: str) -> int:
        """Return the repeat unit length using a rolling k-mer search."""
        k = self.k
        anchor = seq[:k]
        n = len(seq)
        max_start = n - k

        for start in range(k, max_start + 1):
            window = seq[start:start+k]
            mis = self._compare(anchor, window, self.max_errors)
            if mis is not None:
                # we found a match within budget
                return start

        raise NoRepeats #(f"No repeat of the first {k}-mer within {self.max_errors} errors")


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

    def decide_orientation(self, stream: FastqStream, sample_size: int = 10000) -> str:
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
        print(f"RC Count: {rc}\nForward Count: {forward}")
        return "RC" if rc > forward else "forward"


# ----------------------- PhaseAligner ---------------------------
@njit
def lazy_best_shift(cons_idx, read_idx, d, max_errors, window: int = 25):
    search_limit = len(read_idx) - window + 1
    for phi in range(search_limit):  # read shift (phase)
        mismatches = 0
        for i in range(window):
            r = read_idx[phi + i]
            if r > 3:
                continue
            c = cons_idx[i % d]
            if r != c:
                mismatches += 1
                if mismatches > max_errors:
                    break
        if mismatches <= max_errors:
            return phi
    return -1

class PhaseAligner:
    # def best_shift(self, cons_idx: np.ndarray, read_idx: np.ndarray, read_q: np.ndarray, d: int) -> Tuple[int, int]:
    #     pos = np.arange(len(read_idx))
    #     cols = (pos[np.newaxis, :] + np.arange(d)[:, np.newaxis]) % d
    #     matches = (read_idx[np.newaxis, :] == cons_idx[cols])
    #     scores = (read_q[np.newaxis, :] * matches).sum(axis=1)
    #     phi = int(np.argmax(scores))
    #     return phi, int(scores[phi])

    def __init__(self, max_errors, k):
        self.fallback_count = 0
        self.max_errors = max_errors
        self.k = k

    def best_shift(self, cons_idx, read_idx, read_q, d):
        # 1) Try the fast, thresholded scan first
        window = self.k
        max_errors = self.max_errors
        phi = lazy_best_shift(cons_idx, read_idx, d, max_errors, window)
        if phi >= 0:
            return phi, None      # only phi, no score

        # 2) Fallback: do a window-restricted search 
        #    (mismatch-minimizing, quality tie-breaker)
        self.fallback_count += 1
        best_phi   = None
        best_mis   = window+1
        best_score = -1

        max_phi = len(read_idx) - window + 1
        for p in range(max_phi):
            mis, score = 0, 0
            for i in range(window):
                r = read_idx[p+i]
                if r > 3: 
                    continue
                if r != cons_idx[i]:
                    mis += 1
                    if mis > max_errors:
                        break
                else:
                    score += read_q[p+i]
            else:
                # we never broke → mis ≤ max_errors
                if mis < best_mis or (mis == best_mis and score > best_score):
                    best_mis, best_score, best_phi = mis, score, p

        # If best_phi stayed None, you could decide to pick some default (e.g. 0),
        # but whatever you pick, **only** return that phi:
        if best_phi is not None:
            while best_phi >= d:
                best_phi -= d
        return best_phi if best_phi is not None else 0, None


    # def best_shift(self, cons_idx: np.ndarray, read_idx: np.ndarray, read_q: np.ndarray, d: int, max_errors: int = 2, window: int = 40) -> Tuple[int, int | None]:
    #     # Attempt fast mismatch-limited scan
    #     max_errors = 2
    #     phi = lazy_best_shift(cons_idx, read_idx, d, max_errors, window)
    #     if phi >= 0:
    #         return phi, None

    #     # Fallback: quality-weighted scoring
    #     self.fallback_count += 1
    #     pos = np.arange(len(read_idx))
    #     cols = (pos[np.newaxis, :] + np.arange(d)[:, np.newaxis]) % d
    #     matches = (read_idx[np.newaxis, :] == cons_idx[cols])
    #     scores = (read_q[np.newaxis, :] * matches).sum(axis=1)
    #     phi = int(np.argmax(scores))
    #     return phi, int(scores[phi])

    def merge(self, matrix: ConsensusMatrix, seq: str, qual: str, shift: int):
        matrix.update(seq, qual, shift=shift)

# ----------------------- ResultWriter Class ------------------

class ResultWriter:
    """Helper for writing output FASTQ and metadata files."""
    def __init__(self, fastq_path: str, meta_path: str):
        # open FASTQ output (gzipped) and metadata output
        meta_dir = os.path.dirname(meta_path)
        if meta_dir:
            os.makedirs(meta_dir, exist_ok=True)

        fastq_dir = os.path.dirname(fastq_path)
        if fastq_dir:
            os.makedirs(fastq_dir, exist_ok=True)

        self.fastq_fh = gzip.open(fastq_path, "wt")
        self.meta_fh  = gzip.open(meta_path, "wt")

    def write(self, read_id: str, seq: str, quals: str, length: int, phase: int):
        # write a FASTQ record
        print(f"@{read_id}", file=self.fastq_fh)
        print(seq,       file=self.fastq_fh)
        print("+",       file=self.fastq_fh)
        print(quals,     file=self.fastq_fh)
        # write metadata: read_id, consensus length, phase shift
        self.meta_fh.write(f"{read_id}\t{length}\t{phase}\n")

    def close(self):
        self.fastq_fh.close()
        self.meta_fh.close()

# ----------------------- Helper Functions -----------------------

def rotate_back(seq, q, phi):
    """Rotate ``seq`` and ``q`` backwards by ``phi`` bases."""
    return seq[-phi:] + seq[:-phi], q[-phi:] + q[:-phi]

# ----------------------- RepeatPhasingPipeline ------------------

@dataclass
class PhasingResult:
    """Container for results yielded by :class:`RepeatPhasingPipeline`."""

    read_id: str
    phase_shift: int
    consensus_len: int
    elapsed: float


class RepeatPhasingPipeline:
    """Final Outputs: Read ID, Length of Consensus, Phase Shift, Double Consensus (with Quality Scores)"""
    def __init__(self, r1_path: str, r2_path: str, sample_size: int = 1000, k: int = 25, max_errors: int = 2, max_reads: int | None = None, fastq_out: str = "output.fastq.gz", meta_out: str  = "metadata.txt.gz"):
        self.stream = FastqStream(r1_path, r2_path)
        self.detector = RepeatDetector(k=k, max_errors=max_errors)
        self.decider = OrientationDecider(self.detector)
        self.aligner = PhaseAligner(max_errors=max_errors, k=k)
        self.sample_size = sample_size
        self.max_reads = max_reads
        self.total_reads = 0
        self.no_repeat_count = 0
        self.writer = ResultWriter(fastq_out, meta_out)

    def run(self) -> Iterator[PhasingResult]:
        timings = defaultdict(float)
        orientation = self.decider.decide_orientation(self.stream, self.sample_size)
        processed = 0
        for r1, r2 in self.stream:
            self.total_reads += 1
            if self.max_reads is not None and processed >= self.max_reads:
                break

            pair_start = perf_counter()

            t = perf_counter()
            seq1 = str(r1.seq)
            try:
                d = self.detector.find_repeat_distance(seq1)
            except NoRepeats:
                self.no_repeat_count += 1
                continue
            timings["repeat detection"] += perf_counter() - t

            t = perf_counter()
            qual1 = "".join(chr(q + 33) for q in r1.letter_annotations["phred_quality"])
            cm = ConsensusMatrix(d)
            cm.update(seq1, qual1)
            timings["consensus build"] += perf_counter() - t

            t = perf_counter()
            cons_seq, _, _, _ = cm.to_consensus()
            cons_idx = encode_seq(cons_seq)
            timings["consensus extract"] += perf_counter() - t

            t = perf_counter()
            seq2 = str(r2.seq)
            qual2 = "".join(chr(q + 33) for q in r2.letter_annotations["phred_quality"])
            if orientation == "RC":
                seq2 = str(Seq(seq2).reverse_complement())
                qual2 = qual2[::-1]
            read_idx = encode_seq(seq2)
            read_q = encode_qual(qual2)
            timings["R2 prep"] += perf_counter() - t

            t = perf_counter()
            phi, _ = self.aligner.best_shift(cons_idx, read_idx, read_q, d)
            timings["phase alignment"] += perf_counter() - t

            # t = perf_counter()
            # self.aligner.merge(cm, seq2, qual2, phi)
            # timings["merge"] += perf_counter() - t

            t = perf_counter()
            rot_seq2 = seq2[phi:] + seq2[:phi]
            rot_qual2 = qual2[phi:] + qual2[:phi]
            cm.update(rot_seq2, rot_qual2, shift=0)
            timings["merge"] += perf_counter() - t

            # derive final consensus length (single)
            single_len = cm.consensus_length()

            # extract final consensus string and quality list
            single_seq, single_q, _, _ = cm.to_consensus()
            # single_seq, single_q = rotate_back(single_seq, single_q, phi)
            # build ASCII-encoded quality string for single consensus, clamp values to [0,42]
            qstr_single = "".join(
                chr(max(0, min(int(q), 42)) + 33)
                for q in single_q
            )

            # create double consensus (sequence and qualities concatenated twice)
            double_seq = single_seq + single_seq
            qstr_double = qstr_single + qstr_single

            # write double-consensus FASTQ and metadata (metadata uses single consensus length)
            self.writer.write(r1.id, double_seq, qstr_double, single_len, phi)

            elapsed = perf_counter() - pair_start
            result = PhasingResult(
                read_id=r1.id,
                phase_shift=phi,
                consensus_len=single_len,
                elapsed=elapsed,
            )
            del cm
            processed += 1
            yield result

        total = sum(timings.values())
        if total:
            print("Timing summary:")
            for name, secs in sorted(timings.items(), key=lambda x: x[1], reverse=True):
                pct = secs / total * 100
                print(f"  {name:16s}{secs:.3f}s ({pct:.1f}%)")
        
        # close output files
        self.writer.close()

