import pysam
from Bio import SeqIO

def count_mismatches_limited(a, b, max_errors):
    mismatches = 0
    for x, y in zip(a, b):
        if x != y:
            mismatches += 1
            if mismatches > max_errors:
                return mismatches
    return mismatches

def find_breakpoint(ref_seq, double_consensus, k=20, max_errors=5):
    for i in range(len(double_consensus) - k + 1):
        window = double_consensus[i:i + k]
        target = ref_seq[i:i + k]
        if len(target) < k:
            break
        errors = count_mismatches_limited(window, target, max_errors)
        if errors <= max_errors:
            return i
    return None

def lazy_align_reads(bam_path, fasta_path, output_path="refined.bam", k=20, max_errors=5, max_reads=None):
    ref_dict = SeqIO.to_dict(SeqIO.parse(fasta_path, "fasta"))
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    outfile = pysam.AlignmentFile(output_path, "wb", template=bamfile)

    count = 0
    for read in bamfile.fetch(until_eof=True):
        if max_reads is not None and count >= max_reads:
            break
        count += 1

        ref_name = read.reference_name
        pos = read.reference_start
        read_id = read.query_name

        double_consensus = str(read.query_sequence)
        consensus = double_consensus[:len(double_consensus) // 2]

        if ref_name not in ref_dict:
            print(f"[WARN] {read_id} – reference {ref_name} not found in FASTA.")
            read.is_unmapped = True
            outfile.write(read)
            continue

        ref_seq = str(ref_dict[ref_name].seq)
        ref_slice = ref_seq[pos:pos + len(double_consensus)*2]

        if len(ref_slice) < k:
            print(f"[WARN] {read_id} – reference slice too short.")
            read.is_unmapped = True
            outfile.write(read)
            continue

        offset = find_breakpoint(ref_slice, double_consensus, k, max_errors)

        if offset is not None:
            refined_pos = pos + offset
            refined_seq = double_consensus[offset:offset + len(consensus)//2]
            read.reference_start = refined_pos
            read.query_sequence = refined_seq
            read.query_qualities = pysam.qualitystring_to_array("K" * len(refined_seq))
            read.cigar = [(0, len(refined_seq))]
            read.set_tag("NH", 0)
            read.set_tag("BP", offset)
            read.set_tag("ST", "-", value_type='A')
            read.set_tag("MS", ref_name)
            outfile.write(read)
            print(f"[MATCH] {read_id}: Refined POS = {refined_pos} (original POS: {pos})")
        else:
            read.is_unmapped = True
            outfile.write(read)
            print(f"[FLAG] {read_id}: No confident match found – written as unmapped")

    bamfile.close()
    outfile.close()

if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Lazy breakpoint finder from BAM using consensus matching.")
    parser.add_argument("bam", help="Path to input BAM file")
    parser.add_argument("fasta", help="Path to reference FASTA file")
    parser.add_argument("-o", "--output", default="refined.bam", help="Path to output BAM file (default: refined.bam)")
    parser.add_argument("-k", type=int, default=5, help="K-mer size (default: 5)")
    parser.add_argument("-e", "--max-errors", type=int, default=1, help="Max allowed errors in k-mer match (default: 1)")
    parser.add_argument("-n", "--max-reads", type=int, default=None, help="Maximum number of reads to process (default: all)")

    args = parser.parse_args()

    lazy_align_reads(
        bam_path=args.bam,
        fasta_path=args.fasta,
        output_path=args.output,
        k=args.k,
        max_errors=args.max_errors,
        max_reads=args.max_reads
    )
