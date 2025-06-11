import subprocess

"""
Assumes the prep used Illumina TruSeq (Paired-end, Read 2) with adapted sequence = AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT
"""

def trim_reads_cutadapt(input_fastq, output_fastq, adapter_sequence, quality_threshold=20):
    command = [
        "cutadapt",
        "-a", adapter_sequence,
        "-q", str(quality_threshold),
        "-o", output_fastq,
        input_fastq
    ]
    subprocess.run(command, check=True)

# Example usage
input_file = "reads/T1_S30_L002_R1_001.fastq.gz"
output_file = "trimmed_reads.fastq.gz"
adapter = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT"
trim_reads_cutadapt(input_file, output_file, adapter)
