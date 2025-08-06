import pysam
from collections import defaultdict
import argparse
import csv

def tally_chimeras_by_rname(bam_path: str, output_csv_path: str):
    bamfile = pysam.AlignmentFile(bam_path, "rb")
    chimera_tally = defaultdict(int)

    for read in bamfile.fetch(until_eof=True):
        if read.has_tag("CH") and read.get_tag("CH") == 1:
            rname = read.reference_name
            chimera_tally[rname] += 1

    bamfile.close()

    # Write results to CSV
    with open(output_csv_path, mode='w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["RNAME", "Chimeric_Read_Count"])
        for rname, count in chimera_tally.items():
            writer.writerow([rname, count])

    print(f"Chimera tally written to: {output_csv_path}")

def main():
    parser = argparse.ArgumentParser(description="Count chimeric reads by RNAME and write to CSV")
    parser.add_argument("bam", help="Input BAM file")
    parser.add_argument("output_csv", help="Output CSV file path")
    args = parser.parse_args()

    tally_chimeras_by_rname(args.bam, args.output_csv)

if __name__ == "__main__":
    main()
