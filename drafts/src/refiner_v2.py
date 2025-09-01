#!/usr/bin/env python3

import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq


def find_exact_consensus_match(ref_seq: str, double_consensus: str, kallisto_index: int):
    """
    Search for any phase-shifted single-consensus substring within the doubled-consensus
    by trying all rotations (phi = 0..d). Returns (new_start, is_reverse, phi) if found.
    """
    d = len(double_consensus) // 2
    # Try every rotation of the single-consensus within the doubled sequence
    for phi in range(d + 1):
        candidate = double_consensus[phi : phi + d]
        rc_candidate = str(Seq(candidate).reverse_complement())
        # search forward
        pos = ref_seq.find(candidate, kallisto_index)
        if pos != -1:
            return pos, False, phi
        # search reverse
        pos_rc = ref_seq.find(rc_candidate, kallisto_index)
        if pos_rc != -1:
            return pos_rc, True, phi
    return None, None, None


def phase1_align_reads(bam_path, fasta_path, output_path="phase1.bam", max_reads=None):
    """
    Phase-1: slide single-consensus into its exact match, update flags/tags, and write out.
    Unmatched/unmapped reads are flagged PH=2 and passed through.
    """
    # Load reference sequences
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    outfile = pysam.AlignmentFile(output_path, "wb", template=bamfile)

    total = phase2 = 0
    for read in bamfile.fetch(until_eof=True):
        if max_reads and total >= max_reads:
            break
        total += 1

        # skip unmapped or missing-reference reads
        if read.is_unmapped or read.reference_name not in ref_dict:
            read.set_tag("PH", 2)
            read.is_unmapped = True
            outfile.write(read)
            continue

        ref_seq = str(ref_dict[read.reference_name].seq)
        kidx = read.reference_start
        dc = read.query_sequence  # doubled-consensus string

        match_pos, is_rev, phi = find_exact_consensus_match(ref_seq, dc, kidx)
        if match_pos is not None:
            # update reference start
            read.reference_start = match_pos
            # extract the rotated consensus
            d = len(dc) // 2
            seq = dc[phi : phi + d]
            if is_rev:
                seq = str(Seq(seq).reverse_complement())

            # overwrite sequence & qualities
            read.query_sequence = seq
            # assign dummy quality if needed
            try:
                original_quals = read.query_qualities
                read.query_qualities = original_quals[phi : phi + d]
            except Exception:
                read.query_qualities = [40] * len(seq)

            # set CIGAR to dM
            read.cigartuples = [(0, len(seq))]

            # update flags & tags
            read.is_unmapped = False
            read.is_reverse = is_rev
            read.set_tag("PH", 1)
            read.set_tag("BP", phi)
            read.set_tag("ST", "-" if is_rev else "A")
            outfile.write(read)
        else:
            # phase2 fallback
            phase2 += 1
            read.is_unmapped = True
            read.set_tag("PH", 2)
            outfile.write(read)

    bamfile.close()
    outfile.close()
    print(f"Processed {total} reads, {phase2} need phase2")


def main():
    parser = argparse.ArgumentParser(description="Phase 1 alignment for TR-Errors Pipeline")
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("fasta", help="Reference FASTA file")
    parser.add_argument("-o", "--output", default="phase1.bam", help="Output BAM path")
    parser.add_argument("-n", "--max_reads", type=int, help="Maximum reads to process")
    args = parser.parse_args()
    phase1_align_reads(
        bam_path=args.bam,
        fasta_path=args.fasta,
        output_path=args.output,
        max_reads=args.max_reads,
    )

if __name__ == "__main__":
    main()
