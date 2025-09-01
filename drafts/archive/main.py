import subprocess
import os

# Paths to your input files
read1 = "trimmedReads/T1_S30_L002_R1_001.fastq.gz"
read2 = "trimmedReads/T1_S30_L002_R2_001.fastq.gz"

# Output files
merged = "merged/merged.fastq.gz"
unmerged_R1 = "merged/unmerged_R1.fastq.gz"
unmerged_R2 = "merged/unmerged_R2.fastq.gz"

# Path to bbmerge.sh (use full path if not in PATH)
bbmerge_cmd = [
    "bbmap/bbmerge.sh",
    f"in1={read1}",
    f"in2={read2}",
    f"out={merged}",
    f"outu1={unmerged_R1}",
    f"outu2={unmerged_R2}",
    "strict=t"  # Optional: stricter matching, you can remove or tune it
]

# Run BBMerge
# subprocess.run(bbmerge_cmd, check=True)

# Define input/output paths
# input_output_map = {
#     "merged/merged.fastq.gz": "fasta/merged.fasta",
#     "merged/unmerged_R1.fastq.gz": "fasta/unmerged_R1.fasta",
#     "merged/unmerged_R2.fastq.gz": "fasta/unmerged_R2.fasta"
# }

# # Ensure output folder exists
# os.makedirs("fasta", exist_ok=True)

# # Run seqtk for each file
# for fastq, fasta in input_output_map.items():
#     print(f"Converting {fastq} → {fasta}")
#     subprocess.run(["seqtk", "seq", "-A", fastq], stdout=open(fasta, "w"), check=True)

# # Optionally concatenate
# with open("fasta/all_reads.fasta", "w") as outfile:
#     for fasta_file in input_output_map.values():
#         with open(fasta_file) as infile:
#             outfile.write(infile.read())

# Path to TRF binary
trf_path = "./trf"  # or full path like "/home/username/bin/trf"

# Input and output
fasta_path = "fasta/all_reads.fasta"
output_prefix = "fasta/all_reads.fasta"  # TRF adds .2.7.7.80.10.50.500.dat to this

# TRF parameters
params = ["2", "7", "7", "80", "10", "50", "500", "-d", "-h"]

# Build command
trf_cmd = [trf_path, fasta_path] + params

# Run TRF
print("Running Tandem Repeats Finder...")
subprocess.run(trf_cmd, check=True)
print("TRF run complete.")

# Output will be: fasta/all_reads.fasta.2.7.7.80.10.50.500.dat and others

