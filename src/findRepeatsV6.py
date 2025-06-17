from Bio import SeqIO
from Bio.Seq import Seq
import gzip
import sys
import numpy as np
from tqdm import tqdm
from time import perf_counter
from datetime import datetime

class NoRepeats(Exception):
    pass

class RepeatFinder:
    def __init__(self):
        # Build once: ASCII→{0,1,2,3,4}
        self.base2idx = np.full(256, 4, dtype=np.int8)
        self.base2idx[ord('A')] = 0
        self.base2idx[ord('C')] = 1
        self.base2idx[ord('G')] = 2
        self.base2idx[ord('T')] = 3
        self.base2idx[ord('N')] = 4

        # Reverse: {0,1,2,3}→"ACGT"
        self.idx2base = np.array(list("ACGT"), dtype="<U1")

    def compare_string(self, str1, str2, errors=None):
        # Edge case: check length difference
        delta_len = len(str1) - len(str2)
        if errors is not None and abs(delta_len) > errors:
            return False

        # Default error margin is 1% of string length
        if errors is None:
            errors = len(str1) * 0.01

        error_count = 0
        for a, b in zip(str1, str2):
            if a != b:
                error_count += 1
                if error_count > errors:
                    return False
        return True

    def find_repeat_dist(self, sequence: str, k=20, max_errors=2):
        seq_len = len(sequence)
        anchor = sequence[0:k]
        max_range = seq_len - k

        for i in range(max_range):
            candidate = sequence[k + i : 2*k + i]
            if self.compare_string(anchor, candidate, max_errors):
                return k + i

        raise NoRepeats

    def generate_stack(self, sequence: str, delta_d: int) -> list:
        segments = []
        for i in range(0, len(sequence), delta_d):
            segments.append(sequence[i : i + delta_d])
        return segments

    def phreds_to_ascii_string(self, phreds):
        shift = 33
        def _single_phred_to_ascii(q):
            val = int(q)
            val = max(0, min(val, 42))
            return chr(val + shift)

        return "".join(_single_phred_to_ascii(q) for q in phreds)

    def consensus_from_long_read_numpy(self, long_seq: str, long_qual: str, d: int):
        M = len(long_seq)
        if len(long_qual) != M:
            raise ValueError("Sequence and quality must be same length")

        # 1) ASCII→0..4
        ascii_vals = np.frombuffer(long_seq.encode("ascii"), dtype=np.uint8)
        seq_idx = self.base2idx[ascii_vals]  # shape (M,)

        # 2) ASCII→Phred (0..)
        phreds = (np.frombuffer(long_qual.encode("ascii"), dtype=np.uint8) - 33).astype(np.int16)

        # 3) column indices
        cols = np.arange(M, dtype=np.int64) % d

        # 4) build 4×d
        sum4 = np.zeros((4, d), dtype=np.int32)
        for b in range(4):
            mask = (seq_idx == b)
            if not mask.any():
                continue
            col_ids = cols[mask]
            wgt     = phreds[mask]
            sum4[b, :] = np.bincount(col_ids, weights=wgt, minlength=d)

        # 5) best base per column
        best_idx = np.argmax(sum4, axis=0)         # shape (d,)
        total    = np.sum(sum4, axis=0)            # shape (d,)
        best     = sum4[best_idx, np.arange(d)]    # shape (d,)
        cons_q   = best - (total - best)           # shape (d,)

        # 6) turn best_idx → letters
        cons_bases = self.idx2base[best_idx]       # e.g. array(['A','C','G',...], dtype='<U1')
        cons_seq = "".join(cons_bases.tolist())

        # 7) Return (string, list_of_phreds, index_array)
        return cons_seq, cons_q.tolist(), best_idx

    def build_sum4_matrix(self, seq: str, phred_ints: list, d: int) -> np.ndarray:
        """
        Build the same 4×d matrix that consensus_from_long_read_numpy constructed,
        but just return the raw sum4 (not yet turned into a consensus).
        """
        M = len(seq)
        ascii_vals = np.frombuffer(seq.encode("ascii"), dtype=np.uint8)
        idx_arr    = self.base2idx[ascii_vals]          # shape (M,)

        phreds = np.array(phred_ints, dtype=int)        # shape (M,)
        cols   = np.arange(M, dtype=np.int64) % d

        sum4 = np.zeros((4, d), dtype=int)
        for b in range(4):
            mask = (idx_arr == b)
            if not mask.any():
                continue
            col_ids = cols[mask]
            wgt     = phreds[mask]
            sum4[b, :] = np.bincount(col_ids, weights=wgt, minlength=d)
        return sum4

    def find_best_shift(self, cons_idx: np.ndarray, read_idx: np.ndarray,
                        read_phred: np.ndarray, d: int) -> (int, int):
        """
        For a given read (read_idx, read_phred), and a consensus index array (length d),
        slide read over all φ ∈ [0..d−1] to find:
          - best_score = sum of Phreds where read matches consensus
          - best_phi   = argmax φ
        """
        M = len(read_idx)
        best_score = -1
        best_phi = 0
        positions = np.arange(M, dtype=np.int64)

        for φ in range(d):
            cols = (φ + positions) % d
            mask_match = (read_idx == cons_idx[cols])
            if not mask_match.any():
                score = 0
            else:
                score = int(np.sum(read_phred[mask_match]))
            if score > best_score:
                best_score, best_phi = score, φ
        return best_score, best_phi

    def score_orientation(self, cons_idx: np.ndarray, r2_seq: str, r2_qual: str, d: int):
        """
        Return (score_forward, phi_forward, score_rc, phi_rc) for this single R2 read.
        """
        # (1) Forward
        ascii_f   = np.frombuffer(r2_seq.encode("ascii"), dtype=np.uint8)
        r2_f_idx  = self.base2idx[ascii_f]
        r2_f_phred = (np.frombuffer(r2_qual.encode("ascii"), dtype=np.uint8) - 33).astype(int)
        score_f, phi_f = self.find_best_shift(cons_idx, r2_f_idx, r2_f_phred, d)

        # (2) Reverse complement
        rc_seq   = str(Seq(r2_seq).reverse_complement())
        rc_qual  = r2_qual[::-1]
        ascii_rc = np.frombuffer(rc_seq.encode("ascii"), dtype=np.uint8)
        r2_rc_idx   = self.base2idx[ascii_rc]
        r2_rc_phred = (np.frombuffer(rc_qual.encode("ascii"), dtype=np.uint8) - 33).astype(int)
        score_rc, phi_rc = self.find_best_shift(cons_idx, r2_rc_idx, r2_rc_phred, d)

        return score_f, phi_f, score_rc, phi_rc

    def find_R2_orientation(self, cons_idx: np.ndarray, r2_seq: str, r2_qual: str, d: int) -> str:
        """
        Return just "forward" or "RC" for this single R2 read, by comparing scores.
        """
        score_f, phi_f, score_rc, phi_rc = \
            self.score_orientation(cons_idx, r2_seq, r2_qual, d)

        if score_rc > score_f:
            return "RC"
        else:
            return "forward"

    def decide_global_orientation(self, r1_path: str, r2_path: str,
                                  sample_size: int = 10000,
                                  k: int = 40, max_errors: int = 2) -> str:
        """
        Look at the first `sample_size` read‐pairs, build R1 consensus for each,
        then see whether R2 aligns better “forward” or “RC.”  Whichever majority wins
        is returned as the global orientation.
        """
        def open_fastq(path):
            return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

        # File handles for reading fastq files
        fh1 = open_fastq(r1_path) 
        fh2 = open_fastq(r2_path)

        # Iterators for reading records in the FASTQ files based on the above file handles
        it1 = SeqIO.parse(fh1, "fastq") 
        it2 = SeqIO.parse(fh2, "fastq")

        count_forward = 0
        count_rc = 0
        n = 0

        for rec1, rec2 in zip(it1, it2):
            if n >= sample_size:
                break
            n += 1

            
            r1_seq = str(rec1.seq) # Maybe changing operations to directly use Seq types instead of strings would speed things up?
            try:
                d = self.find_repeat_dist(r1_seq, k=k, max_errors=max_errors)
            except NoRepeats:
                continue # Shouldn't we handle the lack of repeats? This may cause downstream errors

            # Build R1 consensus (3 outputs now)
            r1_qual_ascii = "".join(chr(q + 33) for q in rec1.letter_annotations["phred_quality"])
            cons_seq, cons_qual_list, cons_idx = self.consensus_from_long_read_numpy(r1_seq, r1_qual_ascii, d) # Why isnt the cons_seq or 

            # Score R2 (just need raw R2 in forward + RC)
            r2_seq_raw  = str(rec2.seq)
            r2_qual_raw = "".join(chr(q + 33) for q in rec2.letter_annotations["phred_quality"])
            orientation = self.find_R2_orientation(cons_idx, r2_seq_raw, r2_qual_raw, d)
            if orientation == "RC":
                count_rc += 1
            else:
                count_forward += 1

        fh1.close()
        fh2.close()

        return "forward" if count_forward >= count_rc else "RC"

    def process_all_pairs(self, r1_path: str, r2_path: str, global_orientation: str,
                          k: int = 40, max_errors: int = 2, report_every=100000):
        """
        Now that we have a single fixed `global_orientation`, process every read‐pair:
          1) compute R1’s repeat‐unit consensus (string, qual_list, idx_array),
          2) apply the fixed `global_orientation` to R2,
          3) do one find_best_shift(...) to get φ,
          4) merge R2 into R1’s 4×d via build_sum4_matrix + R2 binning,
          5) recompute and print final consensus.
        """
        def open_fastq(path):
            return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "r")

        fh1 = open_fastq(r1_path)
        fh2 = open_fastq(r2_path)

        with gzip.open(r1_path, "rt") as f1, gzip.open(r2_path, "rt") as f2:
            it = zip(SeqIO.parse(f1, "fastq"), SeqIO.parse(f2, "fastq"))
            for idx, (rec1, rec2) in enumerate(tqdm(it,
                                                # total=total,
                                                desc="Processing read-pairs",
                                                unit="pairs"), start=1):
                r1_seq = str(rec1.seq)
                try:
                    d = self.find_repeat_dist(r1_seq, k=k, max_errors=max_errors)
                except NoRepeats:
                    continue

                # (a) Build R1 consensus
                r1_qual_ascii = "".join(chr(q + 33) for q in rec1.letter_annotations["phred_quality"])
                cons_seq, cons_qual_list, cons_idx = \
                    self.consensus_from_long_read_numpy(r1_seq, r1_qual_ascii, d)

                # (b) Apply fixed orientation to R2
                r2_seq_raw = str(rec2.seq)
                r2_qual_raw = "".join(chr(q + 33) for q in rec2.letter_annotations["phred_quality"])
                if global_orientation == "RC":
                    r2_seq_used = str(Seq(r2_seq_raw).reverse_complement())
                    r2_qual_used = r2_qual_raw[::-1]
                else:
                    r2_seq_used = r2_seq_raw
                    r2_qual_used = r2_qual_raw

                # (c) Encode r2 → idx, phred arrays
                ascii_u     = np.frombuffer(r2_seq_used.encode("ascii"), dtype=np.uint8)
                r2_u_idx    = self.base2idx[ascii_u]
                r2_u_phred  = (np.frombuffer(r2_qual_used.encode("ascii"), dtype=np.uint8) - 33).astype(int)

                # (d) Find best φ
                score, phi = self.find_best_shift(cons_idx, r2_u_idx, r2_u_phred, d)

                # (e) Build R1’s 4×d
                sum4_R1 = self.build_sum4_matrix(r1_seq, rec1.letter_annotations["phred_quality"], d)

                # (f) Build R2’s 4×d at that φ
                M2 = len(r2_u_idx)
                cols2 = (phi + np.arange(M2, dtype=np.int64)) % d
                sum4_R2 = np.zeros((4, d), dtype=int)
                for b in range(4):
                    mask = (r2_u_idx == b)
                    if mask.any():
                        sum4_R2[b] = np.bincount(cols2[mask], weights=r2_u_phred[mask], minlength=d)

                # (g) Merge and recompute unit consensus
                sum4_total = sum4_R1 + sum4_R2
                best_idx   = np.argmax(sum4_total, axis=0)
                total      = np.sum(sum4_total, axis=0)
                best       = sum4_total[best_idx, np.arange(d)]
                cons_q     = best - (total - best)
                cons_b     = self.idx2base[best_idx]
                final_cons_seq  = "".join(cons_b.tolist())
                final_cons_qual = cons_q.tolist()

                # (h) Example output—I just print ID, orientation, φ, and final consensus
                # print(f"{rec1.id}\t{global_orientation}\tphi={phi}\t{final_cons_seq}")

        fh1.close()
        fh2.close()


if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python findRepeatsV5.py <R1.fastq(.gz)> <R2.fastq(.gz)>")
        sys.exit(1)

    r1_path = sys.argv[1]
    r2_path = sys.argv[2]
    rf = RepeatFinder()

    # 1) Sample the first 10 000 (or fewer) read‐pairs to decide "forward" vs "RC"
    global_orient = rf.decide_global_orientation(r1_path, r2_path, sample_size=10000)
    print(f">> Decided global orientation: {global_orient!r}")

    # 2) Now process ALL pairs with that orientation
    rf.process_all_pairs(r1_path, r2_path, global_orient)
