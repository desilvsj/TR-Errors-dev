import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from typing import Tuple, Optional
import argparse
from Bio import Align


def parse_edlib_cigar(cigar: str):
    import re
    if not cigar:
        return []
    ops = {'M': 0, 'I': 1, 'D': 2}
    return [(ops[op], int(length)) for length, op in re.findall(r'(\d+)([MID])', cigar)]


def sort_bam_by_reference_name(bam_path: str, sorted_bam_path: str):
    pysam.sort("-o", sorted_bam_path, "-n", bam_path)
    return sorted_bam_path


def crop_double_consensus(dc: str, phi: int, d: int, reverse: bool) -> str:
    seq = dc[phi:phi + d]
    return str(Seq(seq).reverse_complement()) if reverse else seq


def find_exact_consensus_match(ref_seq: str, double_consensus: str, kallisto_index: int) -> Tuple[Optional[int], Optional[bool], Optional[int]]:
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


def run_phase1(ref_seq: str, double_consensus: str, kallisto_index: int):
    match_pos, is_rev, phi = find_exact_consensus_match(ref_seq, double_consensus, kallisto_index)
    if match_pos is not None:
        d = len(double_consensus) // 2
        cropped = crop_double_consensus(double_consensus, phi, d, is_rev)
        return {
            "phase": 1,
            "ref_start": match_pos,
            "is_reverse": is_rev,
            "phi": phi,
            "consensus_length": d,
            "single_consensus": cropped
        }
    return None

def run_phase2(ref_seq: str, double_consensus: str, consensus_length: int):
    # Run local alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = "local"

    # Edit scoring
    aligner.match_score = 2       # score for matching bases
    aligner.mismatch_score = -1   # penalty for mismatches
    aligner.open_gap_score = -2   # penalty for opening a gap
    aligner.extend_gap_score = -0.5  # penalty for extending a gap

    # Perform alignment
    alignments = aligner.align(ref_seq, double_consensus)
    ref_start = alignments[0].aligned[0][0][0]

    return {
        "phase": 2,
        "ref_start": ref_start,
        "is_reverse": False,
        "phi": int(alignments[0].aligned[1][0][0]),
        "consensus_length": consensus_length,
        "single_consensus": double_consensus[alignments[0].aligned[1][0][0]: alignments[0].aligned[1][0][0] + consensus_length]
    }


def run_phase2_placeholder(ref_seq: str, double_consensus: str, consensus_length: int):
    # Placeholder implementation
    return {
        "phase": 2,
        "ref_start": 150,
        "is_reverse": False,
        "phi": 3,
        "consensus_length": consensus_length,
        "single_consensus": double_consensus[3:3 + consensus_length]
    }


def refiner_pipeline(bam_path: str, fasta_path: str, output_bam_path: str, max_reads: Optional[int] = None):
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    outfile = pysam.AlignmentFile(output_bam_path, "wb", template=bamfile)

    results = []
    count = 0
    for read in bamfile.fetch(until_eof=True):
        if max_reads is not None and count >= max_reads:
            break

        if read.is_unmapped or read.reference_name not in ref_dict:
            continue

        ref_seq = str(ref_dict[read.reference_name].seq)
        dc = read.query_sequence
        kidx = read.reference_start

        phase1_result = run_phase1(ref_seq, dc, kidx)
        if phase1_result:
            d = phase1_result["consensus_length"]
            phi = phase1_result["phi"]
            is_rev = phase1_result["is_reverse"]
            new_seq = phase1_result["single_consensus"]

            read.reference_start = phase1_result["ref_start"]
            read.query_sequence = new_seq
            try:
                quals = read.query_qualities
                read.query_qualities = quals[phi:phi + d]
            except Exception:
                read.query_qualities = [40] * len(new_seq)
            read.cigartuples = [(0, len(new_seq))]
            read.is_unmapped = False
            read.is_reverse = is_rev

            # BAM tags:
            # PH: Phase assignment (1 = Phase 1 success, 2 = failed Phase 1)
            # BP: Phi shift (start position offset used within the double consensus)
            # ST: Strand tag ('A' = forward, '-' = reverse)
            # CL: Consensus length (length of aligned single consensus)
            read.set_tag("PH", 1)
            read.set_tag("BP", phi)
            read.set_tag("ST", "-" if is_rev else "A")
            read.set_tag("CL", d)
            results.append(phase1_result)
            outfile.write(read)

        else:
            d = len(dc) // 2
            phase2_result = run_phase2(ref_seq, dc, d)

            if phase2_result["ref_start"] is not None:
                d = phase2_result["consensus_length"]
                phi = phase2_result["phi"]
                is_rev = phase2_result["is_reverse"]
                new_seq = phase2_result["single_consensus"]

                read.reference_start = phase2_result["ref_start"]
                read.query_sequence = new_seq
                try:
                    quals = read.query_qualities
                    read.query_qualities = quals[phi:phi + d]
                except Exception:
                    read.query_qualities = [40] * len(new_seq)
                read.cigartuples = [(0, len(new_seq))]
                read.is_unmapped = False
                read.is_reverse = is_rev
                read.set_tag("PH", 3)
                read.set_tag("BP", phi)
                read.set_tag("ST", "-" if is_rev else "A")
                read.set_tag("CL", d)
            else:
                read.is_unmapped = True
                read.set_tag("PH", 4)  # Phase 2 failed

            results.append(phase2_result)
            outfile.write(read)

        count += 1

    bamfile.close()
    outfile.close()
    return results


def main():
    parser = argparse.ArgumentParser(description="TR-Errors Refinement Pipeline")
    parser.add_argument("bam", help="Input BAM file from kallisto")
    parser.add_argument("fasta", help="Reference FASTA file of transcripts")
    parser.add_argument("-o", "--output_prefix", default="refined", help="Output prefix for logging/results")
    parser.add_argument("-n", "--max_reads", type=int, default=None, help="Maximum number of reads to process")
    args = parser.parse_args()

    output_bam_path = f"{args.output_prefix}.bam"
    results = refiner_pipeline(args.bam, args.fasta, output_bam_path, max_reads=args.max_reads)

    # Simple summary print
    phase1_count = sum(1 for r in results if r["phase"] == 1)
    phase2_count = sum(1 for r in results if r["phase"] == 2)

    print(f"Pipeline complete: {len(results)} reads processed.")
    print(f" - Phase 1 successful: {phase1_count}")
    print(f" - Phase 2 (fallback): {phase2_count}")




if __name__ == "__main__":
    main()