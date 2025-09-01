import argparse
import pysam
# from Bio import SeqIO, pairwise2
# from Bio.pairwise2 import format_alignment
from Bio import Align
from Bio.Seq import Seq
from collections import defaultdict
import edlib
import re
import parasail

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
        print(int(i.aligned[1][-1][1]) - int(i.aligned[1][0][0]))
        print(i)

def test():
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"
    aligner.match_score = 2
    aligner.mismatch_score = -10
    aligner.open_gap_score = -100
    aligner.extend_gap_score = -20

    ref_seq = "ATGAATCCTCAAGTCAGTAACATCATCATCATGTTGGTCATGATGCAACTCTCCCGTCGCATTGACATGGAGGACCCAACCATCATCATGTACATTAGAATTTTATACTGTTCTTCCATCGGTATCTCTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACTACCGTCAGAGACTACGATTTGAAGGAAATAGACAGTGCTATCAAGTCTATCTACACTGGTATGGCTATGATGGGTTTCATGCATTTGTACTTGAAATACACCAACCCATTGTTCATGCAATCCATTTCTCCAGTGAAAAGCGCTTTGGAACACAACGAAGTGAAAATTCACCTCTTCGGTAAGCCTGCAACCGGCGATTTGAAGAGACCATTCAAGGCTCCATCTTTGTTTGGTGGTATGGGTCAAACTGGTCCAAAGACCGACAAGAAATCTATCGAAGAAGCTGAAAGAGCCGGTAACGCTGGTGTTAAGGCTGAATGA"
    double_consensus = "CTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACCAACGTCAGGCACTACGAATTAAAGGAAATAGACACCCGTAAGAAATCTCTCTACAATGGAAAGACAATGAGCGGTCTCATAATGCAGCACCGGTATCGCTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACCAACGTCAGGCACTACGAATTAAAGGAAATAGACACCCGTAAGAAATCTCTCTACAATGGAAAGACAATGAGCGGTCTCATAATGCAGCACCGGTATCG"

    alignments = aligner.align(ref_seq, double_consensus)
    consensus_length = int(len(double_consensus)/2)
    best = alignments[0]
    print(best)
    print(best.aligned)
    ref_start = best.aligned[0][0][0]
    query_start = int(best.aligned[1][0][0])
    query_end = int(best.aligned[1][-1][1])
    matches = query_end - query_start

    print(f"phase: {2}\nref_start: {ref_start}\nis_reverse: {False}\nphi: {query_start}\nconsensus_length: {consensus_length}\nsingle_consensus: {double_consensus[query_start: query_start + consensus_length]}\nnum_matches: {matches}")

def test_parasail():
    ref = "ATGAATCCTCAAGTCAGTAACATCATCATCATGTTGGTCATGATGCAACTCTCCCGTCGCATTGACATGGAGGACCCAACCATCATCATGTACATTAGAATTTTATACTGTTCTTCCATCGGTATCTCTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACTACCGTCAGAGACTACGATTTGAAGGAAATAGACAGTGCTATCAAGTCTATCTACACTGGTATGGCTATGATGGGTTTCATGCATTTGTACTTGAAATACACCAACCCATTGTTCATGCAATCCATTTCTCCAGTGAAAAGCGCTTTGGAACACAACGAAGTGAAAATTCACCTCTTCGGTAAGCCTGCAACCGGCGATTTGAAGAGACCATTCAAGGCTCCATCTTTGTTTGGTGGTATGGGTCAAACTGGTCCAAAGACCGACAAGAAATCTATCGAAGAAGCTGAAAGAGCCGGTAACGCTGGTGTTAAGGCTGAATGA"
    dc = "CTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACCAACGTCAGGCACTACGAATTAAAGGAAATAGACACCCGTAAGAAATCTCTCTACAATGGAAAGACAATGAGCGGTCTCATAATGCAGCACCGGTATCGCTTGGATCATCTACCAAATGGCCAGAAAGAGAATTGTTGCTAAAAACGACATGACTACCATGAAGTACGTCGAACCTGGTAATGCTATGTCCGGCGAAGGTGAGAAGCTGCAAGTTACCAACGTCAGGCACTACGAATTAAAGGAAATAGACACCCGTAAGAAATCTCTCTACAATGGAAAGACAATGAGCGGTCTCATAATGCAGCACCGGTATCG"
    
    matrix = parasail.matrix_create("ACGT", 10, -20)

    # Perform local (Smith–Waterman) with striped vectors, 16-bit
    result = parasail.sw_trace_striped_16(dc, ref, 200, 10, matrix)

    # Extract the aligned strings
    aligned_ref   = result.traceback.ref
    aligned_query = result.traceback.query

    # Build the visual match line
    match_line = ''.join(
        '|' if r == q else ' '
        for r, q in zip(aligned_ref, aligned_query)
    )

    # Pull out start/end positions
    start_ref   = result.cigar.beg_ref    # CIGAR start on reference :contentReference[oaicite:2]{index=2}
    end_ref     = result.end_ref          # CIGAR end on reference :contentReference[oaicite:3]{index=3}
    start_query = result.cigar.beg_query  # CIGAR start on query :contentReference[oaicite:4]{index=4}
    end_query   = result.end_query        # CIGAR end on query :contentReference[oaicite:5]{index=5}

    # Print summary
    print(f"Alignment score: {result.score}")
    print(f"Reference aligned from index {start_ref} to {end_ref}")
    print(f"Query     aligned from index {start_query} to {end_query}\n")

    # Print the alignment
    print("Alignment:")
    print(aligned_ref)
    print(match_line)
    print(aligned_query)


def main():
    test_parasail()

if __name__ == "__main__":
    main()



"""I need you to refactor this code. Can you extract the workings of phase 1? I want to create a master function (like a refiner pipeline?) that does the following:

- takes in the bam file from kallisto and the fasta fill with all the reference transcripts (allReads.fa). 
- If the reads from allReads.fa are read into memory, it may be more efficient to sort (using a helper function) the bam file by mRNA so that once the specific mRNA is loaded into memory from allReads.fa all the bam reads can be processed until that mRNA is no longer needed, and can move on to the next set of reads that use the next mRNA (let me know if this doesnt make sense either logically, semantically, or efficiently)
- pass the relevant information (i think this would be reference sequence (from the allReads.fa), double consensus (bam), potential starting index (bam), length of consensus (bam). to the phase 1 function
- if it passes the phase 1 function skip this step. If it fails, pass the information [reference sequence (from the allReads.fa), double consensus (bam), length of consensus (bam)] to the phase 2 function.
- return the consensus, reference starting index, single consensus (which needs to be cropped using a helper function - see below), double consensus starting index, consensus length, and any other information that may be relevant

sorting function:
- takes in the bam file. sorts it according to the mRNA, grouping those with the same mRNA together.
- it may make more sense to do this before entering this pipeline, if so suggest it. Then we can do this before the refiner pipeline and produce a new sorted bam file to begin the pipeline

phase 1 function
- Extract this logic from the existing refiner_v3.py

phase 2 function
- Keep a place holder for now but define the function domain etc.

consensus crop function
- this should crop the double consensus based on the proposed alignment. Lets say the reference and double consensus aligned at indexes 35 and 3 respectively - and the consensus length was 70 - then the double consensus should be cropped from 3 to 73 (or so)
- if the logic exists in the refiner_v3 code go ahead and try to implement it. If not keep a placeholder

You can develop this in a collaborative doc/canvas"""