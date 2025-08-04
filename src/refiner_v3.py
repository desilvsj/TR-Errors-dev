#!/usr/bin/env python3

import argparse
import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
import edlib
import re


def parse_edlib_cigar(cigar: str):
    """
    Convert an Edlib CIGAR string (e.g., '10M1I5M2D') into pysam cigar tuples.
    0=M, 1=I, 2=D. If cigar is None, return empty list.
    """
    if not cigar:
        return []
    ops = {'M': 0, 'I': 1, 'D': 2}
    return [(ops[op], int(length)) for length, op in re.findall(r'(\d+)([MID])', cigar)]


def find_exact_consensus_match(ref_seq: str, double_consensus: str, kallisto_index: int):
    """
    Phase 1: exact-match search anywhere downstream of kallisto_index.
    Returns (new_start, is_reverse, phi) if a perfect d-long match is found.
    """
    d = len(double_consensus) // 2
    for phi in range(d + 1):
        candidate = double_consensus[phi:phi + d]
        rc_candidate = str(Seq(candidate).reverse_complement())
        pos = ref_seq.find(candidate, kallisto_index)
        if pos != -1:
            return pos, False, phi
        pos_rc = ref_seq.find(rc_candidate, kallisto_index)
        if pos_rc != -1:
            return pos_rc, True, phi
    return None, None, None


def phase1_align_reads(bam_path, fasta_path, output_path="phase1.bam", max_reads=None):
    """
    Phase 1: align reads by exact consensus match and tag failures for Phase 2.
    Returns (total_reads, phase2_count).
    """
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    outfile = pysam.AlignmentFile(output_path, "wb", template=bamfile)

    total = 0
    phase2_count = 0
    for read in bamfile.fetch(until_eof=True):
        if max_reads and total >= max_reads:
            break
        total += 1

        if read.is_unmapped or read.reference_name not in ref_dict:
            read.set_tag("PH", 2)
            read.is_unmapped = True
            outfile.write(read)
            phase2_count += 1
            continue

        ref_seq = str(ref_dict[read.reference_name].seq)
        kidx = read.reference_start
        dc = read.query_sequence

        match_pos, is_rev, phi = find_exact_consensus_match(ref_seq, dc, kidx)
        if match_pos is not None:
            read.reference_start = match_pos
            d = len(dc) // 2
            seq = dc[phi:phi + d]
            if is_rev:
                seq = str(Seq(seq).reverse_complement())
            read.query_sequence = seq
            try:
                quals = read.query_qualities
                read.query_qualities = quals[phi:phi + d]
            except Exception:
                read.query_qualities = [40] * len(seq)
            read.cigartuples = [(0, len(seq))]
            read.is_unmapped = False
            read.is_reverse = is_rev
            read.set_tag("PH", 1)
            read.set_tag("BP", phi)
            read.set_tag("ST", "-" if is_rev else "A")
            outfile.write(read)
        else:
            phase2_count += 1
            read.is_unmapped = True
            read.set_tag("PH", 2)
            outfile.write(read)

    bamfile.close()
    outfile.close()
    return total, phase2_count


def seed_index(window: str, k_len: int) -> dict:
    """
    Build a k-mer to positions index for the given string window.
    """
    idx = defaultdict(list)
    for i in range(len(window) - k_len + 1):
        idx[window[i:i + k_len]].append(i)
    return idx


def find_fuzzy_consensus_match(ref_seq: str,
                               double_consensus: str,
                               kallisto_index: int,
                               k_len: int = 13,
                               max_edits: int = 2,
                               padding: int = 100):
    """
    Phase 2: fuzzy-match search within a padded window around kallisto_index.
    Uses k-mer seeding to limit full edit-distance alignments.
    Returns (new_start, is_reverse, phi, edit_distance) or (None, None, None, None).
    """
    d = len(double_consensus) // 2
    start = max(0, kallisto_index - padding)
    end = min(len(ref_seq), kallisto_index + padding + d)
    window = ref_seq[start:end]

    idx = seed_index(window, k_len)
    best = (None, None, None, None)
    best_dist = max_edits + 1

    for phi in range(d + 1):
        candidate = double_consensus[phi:phi + d]
        seed = candidate[d//2 : d//2 + k_len]
        for pos in idx.get(seed, []):
            res = edlib.align(candidate, window, mode="distance", k=max_edits)
            dist = res.get("editDistance")
            if dist is not None and dist < best_dist:
                best = (start + pos, False, phi, dist)
                best_dist = dist
                if best_dist == 0:
                    break
        if best_dist == 0:
            break

    if best[0] is None:
        rc_double = str(Seq(double_consensus).reverse_complement())
        for phi in range(d + 1):
            candidate = rc_double[phi:phi + d]
            seed = candidate[d//2 : d//2 + k_len]
            for pos in idx.get(seed, []):
                res = edlib.align(candidate, window, mode="distance", k=max_edits)
                dist = res.get("editDistance")
                if dist is not None and dist < best_dist:
                    best = (start + pos, True, phi, dist)
                    best_dist = dist
                    if best_dist == 0:
                        break
            if best_dist == 0:
                break

    return best


def phase2_rescue_reads(bam_path, fasta_path, output_path="phase2.bam",
                        k_len=13, max_edits=2, padding=100):
    """
    Phase 2: rescue reads tagged PH=2 by Phase 1 using fuzzy matching.
    Returns (rescued_count, final_fail_count).
    """
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    outfile = pysam.AlignmentFile(output_path, "wb", template=bamfile)

    rescued_count = 0
    final_fail_count = 0
    for read in bamfile.fetch(until_eof=True):
        if read.get_tag("PH") != 2:
            outfile.write(read)
            continue

        if read.reference_name not in ref_dict:
            read.is_unmapped = True
            read.set_tag("PH", 4)
            outfile.write(read)
            final_fail_count += 1
            continue

        ref_seq = str(ref_dict[read.reference_name].seq)
        kidx = read.reference_start
        dc = read.query_sequence

        new_start, is_rev, phi, dist = find_fuzzy_consensus_match(
            ref_seq, dc, kidx, k_len, max_edits, padding
        )
        if new_start is not None:
            read.reference_start = new_start
            d = len(dc) // 2
            seq = dc[phi:phi + d]
            if is_rev:
                seq = str(Seq(seq).reverse_complement())
            read.query_sequence = seq
            try:
                quals = read.query_qualities
                read.query_qualities = quals[phi:phi + d]
            except Exception:
                read.query_qualities = [40] * len(seq)
            window_seq = ref_seq[new_start:new_start + d + padding]
            path = edlib.align(seq, window_seq, mode="path", k=dist).get("cigar")
            cigartuples = parse_edlib_cigar(path)
            if cigartuples:
                read.cigartuples = cigartuples
            else:
                # fallback to simple match
                read.cigartuples = [(0, len(seq))]
            read.is_reverse = is_rev
            read.set_tag("PH", 3)
            read.set_tag("BP", phi)
            read.set_tag("NH", dist)
            outfile.write(read)
            rescued_count += 1
        else:
            read.is_unmapped = True
            read.set_tag("PH", 4)
            outfile.write(read)
            final_fail_count += 1

    bamfile.close()
    outfile.close()
    return rescued_count, final_fail_count


def main():
    parser = argparse.ArgumentParser(description="TR-Errors Refinement Pipeline")
    parser.add_argument("bam", help="Input BAM file from kallisto")
    parser.add_argument("fasta", help="Reference FASTA file of transcripts")
    parser.add_argument("-n", "--max_reads", type=int, help="Max reads for Phase 1")
    parser.add_argument("-k", "--k_len", type=int, default=13, help="k-mer seed length for Phase 2")
    parser.add_argument("-e", "--max_edits", type=int, default=2, help="Max edit distance for Phase 2")
    parser.add_argument("-p", "--padding", type=int, default=100, help="Window padding size for Phase 2")
    parser.add_argument("-o", "--output_prefix", default="refined", help="Output prefix for BAM files")
    args = parser.parse_args()

    phase1_out = f"{args.output_prefix}_phase1.bam"
    phase2_out = f"{args.output_prefix}_phase2.bam"

    total, to_phase2 = phase1_align_reads(args.bam, args.fasta, phase1_out, args.max_reads)
    rescued, failed = phase2_rescue_reads(phase1_out, args.fasta, phase2_out,
                                         args.k_len, args.max_edits, args.padding)

    print(f"Phase 1: {total} reads processed; {to_phase2} sent to Phase 2")
    print(f"Phase 2: {rescued} rescued; {failed} failed both phases")
    print(f"Final output BAM: {phase2_out}")

if __name__ == "__main__":
    main()
