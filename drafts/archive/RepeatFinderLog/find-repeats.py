from collections import defaultdict
from Bio import SeqIO
from tqdm import tqdm

class RepeatFinder():
    def __init__(self):
        pass
    
    def hamming_distance(self, s1, s2):
        return sum(c1 != c2 for c1, c2 in zip(s1, s2))

    def extend_repeat(self, seq, i1, i2, max_len=500, tolerance=0.1):
        max_possible = min(len(seq) - i1, len(seq) - i2, max_len)
        mismatches = 0
        match_len = 0

        for offset in range(max_possible):
            if seq[i1 + offset] != seq[i2 + offset]:
                mismatches += 1
            match_len += 1
            if mismatches / match_len > tolerance:
                break

        if match_len >= 20 and mismatches / match_len <= tolerance:
            return {
                "start1": i1,
                "start2": i2,
                "length": match_len,
                "mismatches": mismatches,
                "repeat1": seq[i1:i1+match_len],
                "repeat2": seq[i2:i2+match_len]
            }

        return None

    def find_first_long_repeat(self, seq, k=20, k_mismatch=2, extend_tol=0.1, stride=None):
        seen = defaultdict(list)
        stride = stride or k
        seq = str(seq)

        for i in range(0, len(seq) - k + 1, stride):
            kmer = seq[i:i+k]
            bucket = kmer[:5]
            for j, prev_kmer in seen[bucket]:
                if self.hamming_distance(kmer, prev_kmer) <= k_mismatch:
                    result = self.extend_repeat(seq, j, i, max_len=300, tolerance=extend_tol)
                    if result:
                        return result
            seen[bucket].append((i, kmer))
        return None

if __name__ == "__main__":
    fastq_path = "src/example_r1.fastq"
    records = list(SeqIO.parse(fastq_path, "fastq"))
    repeat_finder = RepeatFinder()

    for idx, record in enumerate(records):
        print(f"\n🔬 Record {idx + 1}: {record.id}")
        result = repeat_finder.find_first_long_repeat(record.seq, k=20, k_mismatch=2, extend_tol=0.08, stride=1)
        if result:
            print(f"✅ Repeat from {result['start1']} to {result['start2']} (length = {result['length']}, mismatches = {result['mismatches']})")
            print(result['repeat1'])
        else:
            print("❌ No repeat found")