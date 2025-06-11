import subprocess
from Bio import SeqIO
from collections import defaultdict    
from collections import Counter
import gzip
from tqdm import tqdm

"""Virtual Environment at ~/venvs/tr-errors/bin/activate"""

"""This algorithm searches for the first repeat, records the difference in index and uses this to generate the stack of sequences"""

class NoRepeats(Exception):
    pass

# This function takes in the sequence and the min/max of k for the k-mer range. Returns the 
def get_repeat_delta(sequence, k=20):
    seen = {}
    # count = 0
    for i in range(len(sequence) - k + 1):
        # count += 1
        kmer = sequence[i:i + k]
        if kmer in seen:
            # print("get_repear_delta count: ", count)
            return (i - seen[kmer], seen[kmer])  # Return distance between first and second occurrence
        seen[kmer] = i

    # print("get_repear_delta count: ", count)
    return None  # No repeat found



# This function should take in a sequence, a delta, and an index. It should output the sequence starting and <index> continue for <delta> bases, followed by a newline character and repeat until it hits the character before index <index> (i.e., looping around)
def generate_stack(sequence, delta_tuple):
    if delta_tuple:
        sequence_list = list()
        delta = delta_tuple[0]
        index = delta_tuple[1]
        n = len(sequence)
        pos = index
        end_index = (index - 1) % n

        while True:
            chunk = ''
            for _ in range(delta):
                chunk += sequence[pos]
                pos = (pos + 1) % n
                if pos == end_index:
                    break
            sequence_list.append(chunk)
            if pos == end_index:
                break
        return sequence_list
    else:
        raise NoRepeats

def generate_consensus(sequences, delta):
    """
    Generate a consensus sequence of length `delta` and per-base confidence scores.

    Args:
        sequences (list of str): List of repeated sequences.
        delta (int): Desired consensus sequence length.

    Returns:
        tuple: (consensus_sequence: str, confidence_scores: list of float)
    """
    # Truncate all sequences to length delta
    # Pad shorter sequences with 'N' (unknown base)
    truncated = [seq[:delta].ljust(delta, 'N') for seq in sequences]

    if not truncated:
        return "", []

    consensus = []
    confidences = []

    # Iterate over each base position
    for i in range(delta):
        column = [seq[i] for seq in truncated]
        base_counts = Counter(column)

        if base_counts:
            most_common_base, count = base_counts.most_common(1)[0]

            # Exclude 'N' bases from confidence calculation
            total_valid = sum(v for b, v in base_counts.items() if b != 'N')
            if total_valid > 0:
                confidence = count / total_valid if most_common_base != 'N' else 0.0
            else:
                most_common_base = 'N'
                confidence = 0.0
        else:
            most_common_base = 'N'
            confidence = 0.0

        consensus.append(most_common_base)
        confidences.append(confidence)

    return ''.join(consensus), confidences

def confidence_to_symbol(confidence):
    if confidence == 1.0:
        return '#'
    elif confidence >= 0.9:
        return '@'
    elif confidence >= 0.75:
        return '+'
    elif confidence >= 0.5:
        return '='
    elif confidence >= 0.25:
        return '-'
    elif confidence > 0.0:
        return '.'
    else:
        return ' '


# Apply to your file
if __name__ == "__main__":
    GZIP = False

    if GZIP:
        fastq_path = "reads/T1_S30_L002_R1_001.fastq.gz"
        with gzip.open(fastq_path, "rt") as handle:
            output_path = "results/consensus_output.txt"
            with open(output_path, "w") as output_file:
                for record in tqdm(SeqIO.parse(handle, "fastq"), desc="Processing records", dynamic_ncols=True, ascii=True, mininterval=0.5):
                    # print(f"Processing: {record.id}")  # stays in terminal
                    try:
                        delta_tuple = get_repeat_delta(record.seq, 20)
                        sequence_list = generate_stack(record.seq, delta_tuple)
                        consensus, confidences = generate_consensus(sequence_list, delta_tuple[0])
                        
                        # Write results to output file
                        output_file.write(f"Record ID: {record.id}\n")
                        output_file.write("Consensus & Confidence:\n")
                        output_file.write(f"{consensus}\n")
                        output_file.write(f"{''.join(confidence_to_symbol(c) for c in confidences)}\n\n")

                    except NoRepeats:
                        output_file.write(f"Record ID: {record.id}\nNo Repeats found\n\n")

    else:
        fastq_path = "src/example_r1.fastq"
        records = list(SeqIO.parse(fastq_path, "fastq"))

    for record in tqdm(records, desc="Processing records"):
        print(f"Record ID: {record.id}")
        try:
            delta_tuple = get_repeat_delta(record.seq, 20)
            sequence_list = generate_stack(record.seq, delta_tuple)
            consensus, confidences = generate_consensus(sequence_list, delta_tuple[0])
            # print("    Consensus    :", consensus)  
            # print("Confidence Visual:", ''.join(confidence_to_symbol(c) for c in confidences))
            # print("Consensus & Confidence:")
            # print(f"{consensus}\n{''.join(confidence_to_symbol(c) for c in confidences)}")
            for sequence in sequence_list:
                print(sequence)
            print()
        except NoRepeats:
            print("No Repeats found\n")

