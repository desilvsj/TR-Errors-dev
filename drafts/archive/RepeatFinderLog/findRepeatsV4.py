"""
Assumptions:
- No indels
- RCA-like input
- Valid Bases "AGCTN"
"""

"""Virtual Environment at ~/venvs/tr-errors/bin/activate"""

import subprocess
from Bio import SeqIO
from collections import defaultdict    
from collections import Counter
import gzip
import math
import random
from tqdm import tqdm
import numpy as np

class NoRepeats(Exception):
    pass

class RepeatFinder():
    def compare_string(self, str1, str2, errors=None):

        #Edge case: check length of strings
        delta_len = len(str1) - len(str2)
        if abs(delta_len) > errors:
            return False
        
        # Default error margin is 1% of string
        if errors is None:
            errors = len(str1)*0.01

        # Calibrate error count
        error_count = 0

        for a, b in zip(str1, str2): # Generates tuples one at a time
            if a != b:
                error_count += 1
                if error_count > errors:
                    return False
        return True
    
    def compare_string_v2(self, str1, str2, errors=None):

        #Edge case: check length of strings
        delta_len = len(str1) - len(str2)
        if abs(delta_len) > errors:
            return False
        
        # Default error margin is 1% of string
        if errors is None:
            errors = len(str1)*0.01

        # Calibrate error count
        error_count = 0
        max_len = len(str1) + delta_len

        for i in range(max_len):
            if str1[i] != str2[i]:
                error_count += 1
                if error_count > errors:
                    return False
        return True

    def hamming_distance(self, anchor, candidate):
        return sum(c1 != c2 for c1, c2 in zip(anchor, candidate))
    
    def compare_string_hd(self, str1, str2, errors=None):

        # Default error margin is 1% of string
        if errors is None:
            errors = len(str1)*0.01

        error_count = self.hamming_distance(str1, str2)

        if error_count > errors:
            return False
        else:
            return True
        
    def find_repeat_dist(self, sequence: str, k=20, max_errors=2):
        
        """Algorithm starts with the anchor k-mer of seq[0:k] and candidate k-mer of seq[k:2k]. Intial delta_d = k and sequence length of l
        If comparison passes, returns a distance of delta_d (equivalent to k here)
        If comparison fails, repeats comparison with the same anchor k-mer but a candidate k-met of seq[k+i:2k+i] with d +=i; where i=1
        If this fails, repeat for all values of i until k+i = l"""

        seq_len = len(sequence)
        delta_d = k
        anchor = sequence[0:k]
        max_range = seq_len-k

        for i in range(max_range):
            candidate = sequence[k+i:2*k+i]
            if self.compare_string(anchor, candidate, max_errors):
                return delta_d + i
        
        raise NoRepeats
    
    def generate_stack(self, sequence: str, delta_d: int) -> list:
        """
        Splits the input sequence into substrings of length delta_d.
        Includes the final segment even if it's shorter than delta_d.
        """
        segments = []
        for i in range(0, len(sequence), delta_d):
            segments.append(sequence[i:i + delta_d])
        return segments
    
    def find_fastq_repeats(self, filepath, k=40, max_errors=2):

        records = list(SeqIO.parse(filepath, "fastq"))

        for record in tqdm(records, desc="Processing records"):
            print(f"Record ID: {record.id}")
            try:
                repeat_distance = self.find_repeat_dist(str(record.seq), k, max_errors)
                print(f"Detected repeat distance: {repeat_distance}")
                stack = self.generate_stack(str(record.seq), repeat_distance)
                print("Stack:")
                for line in stack:
                    print(line)
            except NoRepeats:
                print("No Repeats Found")
            print()
        
    def find_gzip_repeats(self, filepath, k=40, max_errors=2):

        nbDone = 0

        with gzip.open(filepath, "rt") as handle:
            repeat_count = 0
            no_repeat_count = 0
            for record in tqdm(SeqIO.parse(handle, "fastq"), desc="Processing records"):
                nbDone += 1
                if nbDone % 100000 == 0:
                    print(f"Number Done: {nbDone}")
                try:
                    repeat_distance = self.find_repeat_dist(str(record.seq), k, max_errors)
                    stack = self.generate_stack(str(record.seq), repeat_distance)
                    repeat_count += 1
                except NoRepeats:
                    no_repeat_count += 1
        
        print(f"Repeat Count: {repeat_count}\nNo Repeat Count: {no_repeat_count}")

    def gen_consensus_seq(self, stack):
        """This function should take in the stack produced using the above function and then generate a consensus sequence.
        This should be built as follows:
        - for i in range(len(stack[0])), generate a cumulative quality index score for each potential base there.
        - if one base has a higher cumulative quality score (in most cases), add that to the consensus sequence.
        - Calculate the quality score for this base in the consensus sequence as the difference in scores for the chosen base and the sum of those not chosen
        - if the scores equal each other, do a coin flip (expression) and choose one base"""
        pass

    def phreds_to_ascii_string(self, phreds):
        """
        Convert a sequence of numeric Phred‐quality scores into a single ASCII‐encoded string,
        capping each score at 42 (so any q > 42 becomes 42, which maps to 'K').

        Parameters
        ----------
        phreds : iterable of int
            A sequence (e.g. list or tuple) of nonnegative integer Phred scores.

        Returns
        -------
        str
            A string where each character is chr(min(q,42) + 33).

        ref: https://www.drive5.com/usearch/manual/quality_score.html
        """
        shift = 33 

        # Helper for a single Phred score:
        def _single_phred_to_ascii(q):
            capped = q if q <= 42 else 42
            return chr(capped + shift)

        # Build the output by mapping each integer in `phreds` through the helper:
        return "".join(_single_phred_to_ascii(q) for q in phreds)


    def consensus_from_long_read(self, long_seq: str, long_qual: str, d: int):
        """
        Build a consensus of length d from a single long RCA read (long_seq,long_qual),
        by treating every d-mer as one “row” without explicitly cutting out N rows.
        Returns: (consensus_str, consensus_qualities_list_of_length_d)
        """
        M = len(long_seq)
        if len(long_qual) != M:
            raise ValueError("Sequence and quality must be same length")

        # Prepare a list of d dictionaries, each mapping base→cumulative Phred.
        # e.g. col_scores[i] might end up as {'A':142, 'G':98, 'T':12}
        col_scores = [defaultdict(int) for _ in range(d)]

        # One pass through the entire read:
        for k in range(M):
            col = k % d
            b   = long_seq[k]
            q   = ord(long_qual[k]) - 33
            col_scores[col][b] += q

        # Now build consensus and "consensus_quals" (best_score - sum(other_scores)) for each column i
        consensus_bases = []
        consensus_quals = []

        for i in range(d):
            score_dict = col_scores[i]  # e.g. {'A':142, 'G':98, 'C':30}
            if not score_dict:
                # no bases ever landed in this column (possible if M < d). 
                # You can decide to pick 'N' or skip entirely. Here we pick 'N' with quality 0:
                # consensus_bases.append('N')
                # consensus_quals.append(0)
                continue

            # Find the maximum cumulative score (and detect ties):
            max_score = max(score_dict.values())
            candidates = [base for base, s in score_dict.items() if s == max_score]

            # If tie, pick one at random:
            chosen = random.choice(candidates) if len(candidates) > 1 else candidates[0]
            consensus_bases.append(chosen)

            # Sum of all other bases’ scores:
            other_sum = sum(score_dict.values()) - score_dict[chosen]
            consensus_quals.append(max_score - other_sum)

        return ''.join(consensus_bases), consensus_quals
    
    def consensus_from_long_read_fast(self, long_seq: str, long_qual: str, d: int):
        M = len(long_seq)
        if len(long_qual) != M:
            raise ValueError("Sequence and quality must be same length")

        # 1) Map bases → indices 0..4
        base2idx = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4}

        # 2) Pre‐allocate an (d × 5) integer table of zeros
        scores = [[0]*5 for _ in range(d)]
        #    scores[i][b_idx] will hold the sum of Phred's for base b_idx at column i.

        # 3) Single pass over M bases
        for k in range(M):
            col = k % d
            b   = long_seq[k]
            q   = ord(long_qual[k]) - 33        # Phred score
            idx = base2idx.get(b, 4)            # treat any unexpected char as 'N'→idx=4
            scores[col][idx] += q

        # 4) Build consensus string + per‐column “quality”
        cons_bases = []
        cons_quals  = []
        for i in range(d):
            col_scores = scores[i]             # a list [sumA, sumC, sumG, sumT, sumN]
            # We only want to pick among A/C/G/T for the “true” consensus; ignore row‐4 ('N') unless all others are zero.
            # Find the max among indices 0..3:
            max_score = -1
            best_idx  = 0
            for b_idx in range(4):
                if col_scores[b_idx] > max_score:
                    max_score = col_scores[b_idx]
                    best_idx  = b_idx

            # Check for ties among A/C/G/T
            ties = [b_idx for b_idx in range(4) if col_scores[b_idx] == max_score]
            if len(ties) > 1:
                chosen_idx = random.choice(ties)
                chosen_score = col_scores[chosen_idx]
            else:
                chosen_idx   = best_idx
                chosen_score = max_score

            # If all four are zero (max_score == 0) and col_scores[4] > 0 (i.e. only 'N' or padding),
            # you can choose to emit 'N' or force a tie‐break among A/C/G/T = 0. Here we’ll emit 'N' if
            # A/C/G/T are all zero but ‘N’ has a positive sum:
            if max_score == 0 and col_scores[4] > 0:
                # cons_bases.append('N')
                # cons_quals.append(0)
                continue

            # Otherwise map chosen_idx → actual base letter:
            idx2base = ['A','C','G','T','N']
            cons_bases.append(idx2base[chosen_idx])

            # “Consensus quality” = chosen_score − sum(others among A/C/G/T)
            other_sum = sum(col_scores[b] for b in range(4) if b != chosen_idx)
            cons_quals.append(chosen_score - other_sum)

        return "".join(cons_bases), cons_quals


    def consensus_from_long_read_numpy(self, long_seq: str, long_qual: str, d: int):
        """
        One‐pass NumPy version (C-accelerated) that never builds Python dicts or per‐column dict lookups.

        Returns (consensus_str, consensus_qualities_list_of_length_d).
        """

        M = len(long_seq)
        if len(long_qual) != M:
            raise ValueError("Sequence and quality must be same length")

        # 1) Build a 256→base_index map in a NumPy array so we can vectorize ASCII→{0..4} in one shot.
        #    We will use:
        #       0 = 'A'
        #       1 = 'C'
        #       2 = 'G'
        #       3 = 'T'
        #       4 = 'N'  (or anything else)
        base2idx = np.full(256, 4, dtype=np.int8)
        base2idx[ord('A')] = 0
        base2idx[ord('C')] = 1
        base2idx[ord('G')] = 2
        base2idx[ord('T')] = 3
        base2idx[ord('N')] = 4

        # 2) Convert long_seq → array of uint8 (ASCII codes), then map → 0..4:
        ascii_vals = np.frombuffer(long_seq.encode("ascii"), dtype=np.uint8)  # shape (M,)
        seq_idx = base2idx[ascii_vals]                                    # shape (M,), values ∈ {0..4}

        # 3) Convert long_qual → array of Phred scores (0..~40+):
        phreds = np.frombuffer(long_qual.encode("ascii"), dtype=np.uint8) - 33
        phreds = phreds.astype(np.int16)                                  # shape (M,)

        # 4) Build column indices: each position k maps to column (k % d)
        cols = np.arange(M, dtype=np.int64) % d                                # shape (M,)

        # 5) For each base b in {0,1,2,3}, do a weighted bincount on cols:
        sum4 = np.zeros((4, d), dtype=np.int32)  # row b, column i = sum of phreds where seq_idx==b and (k % d)==i

        for b in range(4):
            mask = (seq_idx == b)         # boolean mask of length M
            if not mask.any():
                continue
            # Weighted sum of phreds[mask] grouped by column indices cols[mask].
            # This returns a length‐d array of sums.
            col_ids = cols[mask]          # which column each b‐base falls into
            wgt     = phreds[mask]        # its Phred score
            sum4[b, :] = np.bincount(col_ids, weights=wgt, minlength=d)

        # 6) Now cum4.shape = (4, d).  cum4[b, i] = total Phred of base b at column i.
        #    Find the best base (0..3) per column:
        best_idx = np.argmax(sum4, axis=0)                   # shape (d,), values ∈ {0,1,2,3}
        total    = np.sum(sum4, axis=0)                      # shape (d,)
        best     = sum4[best_idx, np.arange(d)]              # shape (d,)
        cons_q   = best - (total - best)                     # shape (d,)

        # 7) Turn best_idx → actual letter
        idx2base = np.array(list("ACGT"), dtype="<U1")       # array(['A','C','G','T'])
        cons_bases = idx2base[best_idx]                       # shape (d,), dtype='<U1'
        cons_seq = "".join(cons_bases.tolist())

        # 8) In case you need a Python list of ints for qualities:
        return cons_seq, cons_q.tolist()
            


if __name__ == "__main__":

    """
    Refer codes:
    
    File Type      : Code
    String (1 Seq) : 0
    String (6 Seq) : 1
    Short FASTQ    : 2
    Long FASTQ     : 3
    GZIP           : 4

    """

    CODE = 0
    rf = RepeatFinder()

    if CODE == 0:
        sequence = "GGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGCCAAAGACCCGTGGTGTTGCCATGAATCCAGTTGATCACACTCACGGTGGTGGTAACCATCAACATATTGTTAAGGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGC"
        quality = "FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFF,FF,:,F,,F:F,:F:FF,F,FFFFFFFFFF,,FF:F,,FFFF,:FFFFF,F,,F,FFFFF:FFFFFF,,FFFF::FFFFF:::F:FF:,FFFF:FFFFF:F,F:,FFF,,FF,FFF,:F:F::FFFFFFF,F"
        try:
            repeat_distance = rf.find_repeat_dist(sequence, k=20, max_errors=2)
            print(f"Detected repeat distance: {repeat_distance}")
            stack = rf.generate_stack(sequence, repeat_distance)
            print("Stack:")
            for line in stack:
                print(line)
            print(f"{quality[:repeat_distance]}\n{quality[repeat_distance:]}")
        except NoRepeats:
            print("No Repeats Found")
        
        print("Consensus")
        consensus_seq, consensus_qual = rf.consensus_from_long_read_numpy(sequence, quality, repeat_distance)
        print(consensus_seq)
        print(rf.phreds_to_ascii_string(consensus_qual))
    elif CODE == 1:
        sequences = ["GGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGCCAAAGACCCGTGGTGTTGCCATGAATCCAGTTGATCACACTCACGGTGGTGGTAACCATCAACATATTGTTAAGGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGC",
                 "CTTACAGACCAAGTTGGAGTAACGACGACCAACACCCTTAATAGTGGTCAAAGCGTAAAGGACAATTCTTTCCAATTCTTCTTGGGTCAATTCACCAGCTCTCTTGTGCAAATCAACATCAGCCTTCTTACAGACCAAGTTGGAGTAACGACGACCAACACCCTAAATAGTGGTCAAAGCGTAAAGGACAATTCTTTCGAATTCTTCTTGGGTCAATTCACCAGCTCTCTTGTGCAAATCAACATCAGCCT",
                 "TAGCACCTAACAAATGAACGGTTTCTGGAACTTGAACGTCTAAGCCAAATTCATTATTATCGTGAGAACGTTCAAGTCCATCTTCTTCTTCATTTCTGTCTTCGTCACTGTAGTAGCACCTAACAAATGAACGGTTTCTGGAACTTGAACGTCTAAGCCAAATTCATTATTATCGTGAGAACGTTCAAGTCCATCTTCTTCTTCATTTCTGTCTTCGTCACTGTAGTCGCACCTAACAAATGAACGGTTTC",
                 "GGCATAGTTTATGGTTAAGACTACGACGGTATCTGATCATCTTCGATCCCCTAACTTTCGTTCTTCGGAATCATAGATGGGTGTCGTGCCGAGTGGGTCATTAAAAAAACACCACCCGATCCCTAGTCGGCATAGTTTATGGTTAAGACTACGACGGTATCTGATCATCT",
                 "AGATAAGTGACCAAGAGACTGAGAAAAGGATTATACACCGTTCACCTTCAATCGATGCTTCTCAACCGGAAAATAATTCATCCCTAAGACCGTTTATATTAGGATTGTCAAGACACTCCGGTATTAATGAAGCGGAGCTGGAATTCATTTTCCACGTTCTAGCATTCAAGGTCCCATTCGGGGCTGATCCGGGTTGAAGACATTGTCAGGGGGGGAGTTTGGCTATTAACTTCTCATCTAGCACACAAAGA",
                 "CAAAGCAAACTAACTCCCCCGTCCTAGTACCATTCAACAACAAAGCAAACTAACTCCCCCGTCCTAGTACCATTCAACAACAAAGCAAACTAACTCCCCCGTCCTAGTACCATTCAACAACAAAGCAAACGAACTCCCC"]
        for sequence in sequences:
            try:
                repeat_distance = rf.find_repeat_dist(sequence, k=20, max_errors=2)
                print(f"Detected repeat distance: {repeat_distance}")
                stack = rf.generate_stack(sequence, repeat_distance)
                print("Stack:")
                for line in stack:
                    print(line)
            except NoRepeats:
                print("No Repeats Found")
            print()
    elif CODE == 2:
        rf.find_fastq_repeats("src/example_r1.fastq", k=40, max_errors=2)
    elif CODE == 3:
        pass
    elif CODE == 4:
        rf.find_gzip_repeats("trimmedReads/T1_S30_L002_R1_001.fastq.gz", k=40, max_errors=2)

    
    

    


        
        