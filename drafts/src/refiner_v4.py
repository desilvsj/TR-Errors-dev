import pysam
from Bio import SeqIO
from Bio.Seq import Seq
from collections import defaultdict
from typing import Tuple, Optional
import argparse
from Bio import Align
import time
import re
"""
Parasail Library
Daily, Jeff. (2016). Parasail: SIMD C library for global, semi-global, and local pairwise sequence alignments. BMC Bioinformatics, 17(1), 1-11. doi:10.1186/s12859-016-0930-z
http://dx.doi.org/10.1186/s12859-016-0930-z
"""
import parasail

class Read:
    def __init__(self):
        self.phase = 0
        self.ref_start = None
        self.phi = None
        self.consensus_length = None
        self.consensus_type = "double" #either single or double
        self.match_count = None
        self.sequence = None
        self.match_len = None


def parse_edlib_cigar(cigar: str):
    if not cigar:
        return []
    ops = {'M': 0, 'I': 1, 'D': 2}
    return [(ops[op], int(length)) for length, op in re.findall(r'(\d+)([MID])', cigar)]

def sort_bam_by_reference_name(bam_path: str, sorted_bam_path: str):
    pysam.sort("-o", sorted_bam_path, "-n", bam_path)
    return sorted_bam_path

def crop_double_consensus(dc: str, phi: int, d: int) -> str:
    seq = dc[phi:phi + d]
    return seq

def find_exact_consensus_match(ref_seq: str, double_consensus: str, kallisto_index: int) -> Tuple[Optional[int], Optional[bool], Optional[int]]:
    d = len(double_consensus) // 2
    for phi in range(d + 1):
        candidate = double_consensus[phi:phi + d]
        pos = ref_seq.find(candidate, kallisto_index)
        if pos != -1:
            return pos, phi
    return None, None

def run_phase1(ref_seq: str, double_consensus: str, kallisto_index: int):
    match_pos, phi = find_exact_consensus_match(ref_seq, double_consensus, kallisto_index)
    if match_pos is not None:
        d = len(double_consensus) // 2
        cropped = crop_double_consensus(double_consensus, phi, d)
        return {
            "phase": 1,
            "ref_start": match_pos,
            "phi": phi,
            "consensus_length": d,
            "single_consensus": cropped,
            "num_matches": d
        }
    return None


def run_phase2_parasail(ref_seq: str, double_consensus: str, consensus_length: int):
    matrix = parasail.matrix_create("ACGT", 10, -20)

    # Perform local (Smith–Waterman) with striped vectors, 16-bit
    result = parasail.sw_trace_striped_16(double_consensus, ref_seq, 200, 10, matrix)

    # Extract the aligned strings
    # aligned_ref   = result.traceback.ref
    # aligned_query = result.traceback.query

    start_ref   = result.cigar.beg_ref    # CIGAR start on reference :contentReference[oaicite:2]{index=2}
    end_ref     = result.end_ref          # CIGAR end on reference :contentReference[oaicite:3]{index=3}
    start_query = result.cigar.beg_query  # CIGAR start on query :contentReference[oaicite:4]{index=4}
    end_query   = result.end_query        # CIGAR end on query :contentReference[oaicite:5]{index=5}

    matches = end_query - start_query

        # --- NEW: compute matches/mismatches/indels from traceback ---
    aln_ref  = result.traceback.ref
    aln_qry  = result.traceback.query
    mt = mm = ins = dels = 0
    for r, q in zip(aln_ref, aln_qry):
        if r != '-' and q != '-':
            if r == q: mt += 1
            else:      mm += 1
        elif r == '-' and q != '-':
            ins += 1
        elif r != '-' and q == '-':
            dels += 1
    nm = mm + ins + dels  # SAM-standard edit distance

    return {
        "phase": 2,
        "ref_start": start_ref,
        "is_reverse": False,
        "phi": start_query,
        "consensus_length": consensus_length,
        "single_consensus": double_consensus[start_query:end_query],
        "num_matches": matches,
        "num_mismatches": nm
    }

def refiner_pipeline(bam_path: str, fasta_path: str, output_bam_path: str, max_reads: Optional[int] = None):
    t_load_start = time.time()
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    t_load_end = time.time()

    bamfile = pysam.AlignmentFile(bam_path, "rb")
    outfile = pysam.AlignmentFile(output_bam_path, "wb", template=bamfile)

    results = []
    processed_count = 0
    chimera_count = 0
    discard_count = 0
    total_touched = 0
    reversed_count = 0

    t0 = time.time()
    t_phase1 = 0
    t_phase2 = 0

    for read in bamfile.fetch(until_eof=True):
        total_touched += 1
        # Filter to primary alignments only
        if read.is_secondary or read.is_supplementary:
            discard_count += 1
            continue

        if max_reads is not None and processed_count >= max_reads:
            break
        if read.is_unmapped or read.reference_name not in ref_dict:
            continue

        ref_seq = str(ref_dict[read.reference_name].seq)
        dc = read.query_sequence
        kidx = read.reference_start
        reversal_stat = False

        if read.is_reverse:
            # ref_seq = Seq(ref_seq).reverse_complement()
            reversal_stat = True
            reversed_count+=1


        t1_start = time.time()
        phase1_result = run_phase1(ref_seq, dc, kidx)
        t_phase1 += time.time() - t1_start

        if phase1_result:
            d = phase1_result["consensus_length"]
            phi = phase1_result["phi"]
            new_seq = phase1_result["single_consensus"]
            is_rev = reversal_stat
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

            read.set_tag("PH", 1)
            read.set_tag("BP", phi)
            # read.set_tag("ST", "-" if is_rev else "A")
            read.set_tag("CL", d)
            read.set_tag("RC", reversal_stat)
            results.append(phase1_result)
            outfile.write(read)

        else:
            d = len(dc) // 2
            t2_start = time.time()
            phase2_result = run_phase2_parasail(ref_seq, dc, d)
            t_phase2 += time.time() - t2_start

            if phase2_result["ref_start"] is not None:
                d = phase2_result["consensus_length"]
                phi = phase2_result["phi"]
                is_rev = phase2_result["is_reverse"]
                new_seq = phase2_result["single_consensus"]
                matches = phase2_result["num_matches"]
                nm = phase2_result["num_mismatches"]

                read.reference_start = phase2_result["ref_start"]
                read.query_sequence = new_seq

                chimeric_ratio = (d - matches) / d
                if chimeric_ratio > 0.05:
                    chimera = 1
                    chimera_count += 1
                else:
                    chimera = 0
                try:
                    quals = read.query_qualities
                    read.query_qualities = quals[phi:phi + d]
                except Exception:
                    read.query_qualities = [40] * len(new_seq)
                read.cigartuples = [(0, len(new_seq))]
                read.is_unmapped = False
                read.is_reverse = is_rev
                read.set_tag("PH", 2)
                read.set_tag("BP", phi)
                # read.set_tag("ST", "-" if is_rev else "A")
                read.set_tag("CL", d)
                read.set_tag("MT", matches)
                read.set_tag("NM", nm)
                read.set_tag("CH", chimera)
                read.set_tag("RC", reversal_stat)

            else:
                read.is_unmapped = True
                read.set_tag("PH", 3)

            results.append(phase2_result)
            outfile.write(read)

        processed_count += 1

    bamfile.close()
    outfile.close()

    total_time = time.time() - t0
    rps = total_touched / total_time if total_time > 0 else 0

    print(f"Time to load fasta: {t_load_end - t_load_start:.2f}s")
    print(f"Chimeras detected: {chimera_count}")
    print(f"Reversed Reads: {reversed_count}")
    print(f"Reads discarded (non-primary): {discard_count}")
    print(f"Total reads touched: {total_touched}")
    print(f"Total time: {total_time:.2f}s")
    print(f" - Phase 1 total time: {t_phase1:.2f}s")
    print(f" - Phase 2 total time: {t_phase2:.2f}s")
    print(f" - Reads per second (all reads): {rps:.2f}")

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

    phase1_count = sum(1 for r in results if r and r.get("phase") == 1)
    phase2_count = sum(1 for r in results if r and r.get("phase") == 2)

    print(f"Pipeline complete: {len(results)} reads processed.")
    print(f" - Phase 1 successful: {phase1_count}")
    print(f" - Phase 2 (fallback): {phase2_count}")

if __name__ == "__main__":
    main()
