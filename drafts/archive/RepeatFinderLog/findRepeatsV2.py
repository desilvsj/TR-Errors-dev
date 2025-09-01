"""
Assumptions:
- No indels
- RCA-like input
"""
class RepeatFinder():
    def recursive_repeat_finder(self, seq, i=0, k=10, z=None, max_depth=10, error_margin=2):
        seq = str(seq)
        if z is None:
            z = k

        if z < 1 or max_depth == 0 or i + k >= len(seq):
            return None

        anchor = seq[i:i+k]
        j = i + z
        while j + k <= len(seq):
            candidate = seq[j:j+k]
            # print(f"Anchor: {anchor}\nCandidate: {candidate}")
            mismatches = sum(c1 != c2 for c1, c2 in zip(anchor, candidate)) # Stop checking as soon as it reaches the threshold value 
            if mismatches <= error_margin:
                d = j - i
                sub_result = self.recursive_repeat_finder(seq[i:], i=0, k=k, z=z // 2, max_depth=max_depth - 1, error_margin=error_margin)
                return sub_result if sub_result is not None else d
            j += z

        return None

    def find_repeat_distance(self, seq, k=20, error_margin=2, max_depth=10):
        seq = str(seq)
        # count = 0
        for i in range(len(seq) - k):
            # count +=1
            result = self.recursive_repeat_finder(seq, i=i, k=k, z=k, max_depth=max_depth, error_margin=error_margin)
            if result is not None:
                return result
        # print(f"Count: {count}")
        return None



if __name__ == "__main__":
    str1 = "CAAAGCAAACTAACTCCCCCGTCCTAGTACCATTCAACAA"
    str2 = "CAAAGCAAACGAACTCCCCCGTCCTATTACCATTCAACAA"

    sequence = "GGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGCCAAAGACCCGTGGTGTTGCCATGAATCCAGTTGATCACACTCACGGTGGTGGTAACCATCAACATATTGTTAAGGCTTCTACCATCTCTAGAGGTGCTGTTTCTCAAGTACAGATTGAAGAGAAACTCTTGGC"

    rf = RepeatFinder()
    repeat_distance = rf.find_repeat_distance(sequence, k=20, error_margin=2)

    if repeat_distance:
        print(f"🧬 Repeat unit distance: {repeat_distance}")
    else:
        print("❌ No repeat found.")
