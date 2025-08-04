import argparse
import pysam
# from Bio import SeqIO, pairwise2
# from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.Seq import Seq
from collections import defaultdict
import edlib
import re

def fuzzy_aligner():
    # Example DNA sequences
    seq1 = "ATGAATCCTCAAGTCAGTAACATCATCATCATGTTGGTCATGATGCAACTCTCCCGTCGCATTGACATGGAGGACCCAACCATCATCATGTACATTAGAATTTTATACTGTTCTTCCATCGGTATCTCTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACTACCGTCAGAGACTACGATTTGAAGGAAATAGACAGTGCTATCAAGTCTATCTACACTGGTATGGCTATGATGGGTTTCATGCATTTGTACTTGAAATACACCAACCCATTGTTCATGCAATCCATTTCTCCAGTGAAAAGCGCTTTGGAACACAACGAAGTGAAAATTCACCTCTTCGGTAAGCCTGCAACCGGCGATTTGAAGAGACCATTCAAGGCTCCATCTTTGTTTGGTGGTATGGGTCAAACTGGTCCAAAGACCGACAAGAAATCTATCGAAGAAGCTGAAAGAGCCGGTAACGCTGGTGTTAAGGCTGAATGA"
    seq2 = "CTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACTACCGTCAGAGACTACGATTTGAAGGAAATAGACAGTGATATCAAGTCTATCTACACTGGTATGGCTATGATGGGTTTCATAATGCAGCTCCGGTATCT"

    # Run local alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"

    # Edit scoring
    aligner.match_score = 2       # score for matching bases
    aligner.mismatch_score = -1   # penalty for mismatches
    aligner.open_gap_score = -2   # penalty for opening a gap
    aligner.extend_gap_score = -0.5  # penalty for extending a gap

    # Perform alignment
    alignments = aligner.align(seq1, seq2)

    # Print the best alignment
    for i in alignments:
        print(f"Score: {i.score}")
        print(f"Starting mRNA index: {i.aligned[0][0][0]}")
        print(f"Starting double consensus index: {i.aligned[1][0][0]}")
        print(i)

def main():
    fuzzy_aligner()

if __name__ == "__main__":
    main()