#!/usr/bin/env python3
"""
Align R2 to R1 reads using phase-shift data from an output file.

For each record in the phase data file (ID, PHI, LENGTH),
read the next FASTQ record from R1 and R2 in lock-step.
Print the full R1 sequence and then the RC of R2 shifted by PHI bases.

Usage:
    python test.py \
        --output output.txt \
        --r1 file_R1.fastq.gz \
        --r2 file_R2.fastq.gz \
        --out alignment.txt \
        [--max-reads N] [--no-rc]
"""
import argparse
import gzip
from Bio.Seq import Seq


def read_fastq_record(handle):
    """
    Read one FASTQ entry: header and sequence. Skip + and quality lines.
    Returns (header, seq) or None if EOF.
    """
    header = handle.readline().rstrip()
    if not header:
        return None
    seq = handle.readline().rstrip()
    handle.readline()  # skip '+' line
    handle.readline()  # skip quality line
    return header, seq


def main():
    """Command line entry point for aligning R2 to R1 using phase data."""
    parser = argparse.ArgumentParser(description="Align R2 to R1 using phase-shift data.")
    parser.add_argument('--output', '-o', default='output.txt',
                        help='Phase data file (ID, PHI, LENGTH)')
    parser.add_argument('--r1', default="reads/trimmedReads/T1_S30_L002_R1_001.fastq.gz",
                        help='Path to R1 FASTQ.gz')
    parser.add_argument('--r2', default="reads/trimmedReads/T1_S30_L002_R2_001.fastq.gz",
                        help='Path to R2 FASTQ.gz')
    parser.add_argument('--out', '-O', default='alignment.txt',
                        help='Output file for aligned reads')
    parser.add_argument('--max-reads', type=int, default=50,
                        help='Maximum number of records to process')
    parser.add_argument('--no-rc', dest='rc', action='store_false',
                        help='Disable reverse-complement of R2 (default: RC on)')
    args = parser.parse_args()

    with open(args.output) as phase_fh, \
         gzip.open(args.r1, 'rt') as r1_fh, \
         gzip.open(args.r2, 'rt') as r2_fh, \
         open(args.out, 'w') as out_fh:

        count = 0
        for line in phase_fh:
            if count >= args.max_reads:
                break

            parts = line.strip().split()
            if len(parts) < 2:
                continue
            read_id = parts[0]
            phi = int(parts[1])

            # skip FASTQ until ID matches
            rec1 = rec2 = None
            while True:
                rec1 = read_fastq_record(r1_fh)
                rec2 = read_fastq_record(r2_fh)
                if rec1 is None or rec2 is None:
                    # reached end of FASTQ without match
                    return
                header1, seq1 = rec1
                actual_id = header1.lstrip('@').split()[0]
                if actual_id == read_id:
                    break
                # else: skip and continue looping

            # now rec1/rec2 correspond to read_id
            _, seq2 = rec2

            # reverse-complement R2 (unless disabled)
            r2_seq = Seq(seq2)
            if args.rc:
                r2_seq = r2_seq.reverse_complement()

            # shift by phi bases
            shifted_r2 = str(r2_seq)[phi:]

            # write output
            out_fh.write(f"@{read_id}\n")
            out_fh.write(f"R1: {seq1}\n")
            out_fh.write(f"R2: {shifted_r2}\n\n")

            count += 1


if __name__ == '__main__':  # noqa
    main()
