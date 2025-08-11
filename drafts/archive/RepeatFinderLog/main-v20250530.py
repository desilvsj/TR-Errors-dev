from findRepeatsV2 import RepeatFinder
from Bio import SeqIO

repeat_finder = RepeatFinder()
record = next(SeqIO.parse("src/example_r1.fastq", "fastq"))
result = repeat_finder.find_repeat_unit(record.seq, k_min=20, k_max=150, error_margin=3)

if result:
    k, i = result
    print(f"🧬 Repeat of length {k} starts at index {i}")
else:
    print("❌ No repeat found")
