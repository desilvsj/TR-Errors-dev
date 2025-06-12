import sys
import argparse
from time import perf_counter
from oop_pipeline import RepeatPhasingPipeline


def phreds_to_ascii(phreds):
    shift = 33
    return ''.join(chr(min(q, 42) + shift) for q in phreds)


def main():
    parser = argparse.ArgumentParser(description="Run OOP repeat phasing pipeline")
    parser.add_argument("r1", help="R1 FASTQ (.gz accepted)")
    parser.add_argument("r2", help="R2 FASTQ (.gz accepted)")
    parser.add_argument("-o", "--output", help="Optional output file")
    parser.add_argument("-s", "--sample-size", type=int, default=1000,
                        help="Pairs to sample for orientation")
    args = parser.parse_args()

    pipeline = RepeatPhasingPipeline(args.r1, args.r2, sample_size=args.sample_size)
    out_fh = open(args.output, "w") if args.output else sys.stdout

    t0 = perf_counter()
    count = 0
    for result in pipeline.run():
        count += 1
        qual_str = phreds_to_ascii(result.qualities)
        line = f"{result.phase_shift}\t{result.cycles}\t{result.consensus}\t{qual_str}"
        print(line, file=out_fh)

    if args.output:
        out_fh.close()

    duration = perf_counter() - t0
    rate = count / duration if duration else 0
    print(f"Processed {count} results in {duration:.2f}s ({rate:.1f} pairs/s)")


if __name__ == "__main__":
    main()

