import subprocess
from Bio import SeqIO
from collections import defaultdict

def get_k_tuple_map(sequence, k):
    """
    Given a sequence and a value k, returns a dictionary mapping each unique k-tuple
    (substring of length k) to a list of indices where it occurs in the sequence.
    """
    k_tuple_map = {}

    for i in range(len(sequence) - k + 1):
        k_tuple = sequence[i:i + k]

        if k_tuple not in k_tuple_map:
            k_tuple_map[k_tuple] = []

        k_tuple_map[k_tuple].append(i)

    return k_tuple_map

def count_delta_occurrences(k_tuple_map):
    delta_counts = defaultdict(lambda: [0, None])  # value = [count, first_index]

    for positions in k_tuple_map.values():
        if len(positions) < 2:
            continue
        for i in range(len(positions) - 1):
            delta = positions[i+1] - positions[i]
            if delta_counts[delta][1] is None:
                delta_counts[delta][1] = positions[i]  # first index seen
            delta_counts[delta][0] += 1

    return dict(delta_counts)

# This algorithm searches for the first repeat, records the difference in index and uses this to generate the stack of sequences

# This function takes in the sequence and the min/max of k for the k-mer range. Returns the 
def get_repeat_delta(sequence, k):
    seen = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i + k]
        if kmer in seen:
            return i - seen[kmer]  # Return distance between first and second occurrence
        seen[kmer] = i

    return None  # No repeat found

# def 






# Apply to your file
fastq_path = "src/example_r1.fastq"

records = list(SeqIO.parse(fastq_path, "fastq"))
for record in records:
    print(f"Record ID: {record.id}")
    k_tuple_map = get_k_tuple_map(record.seq, 10)
    print(k_tuple_map)
    D_d = count_delta_occurrences(k_tuple_map)

    if D_d:
        for delta, (count, start_index) in D_d.items():
            print(f"Delta: {delta}, Count: {count}, Start index: {start_index}")
            seq = str(record.seq)  # Ensure it's a string
            i = start_index
            while i < len(seq):
                print(seq[i:i + delta])
                i += delta
    print("\n")

