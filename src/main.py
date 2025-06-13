import sys
import argparse
from time import perf_counter
from oop_pipeline import RepeatPhasingPipeline
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(description="Run OOP repeat phasing pipeline")
    parser.add_argument("r1", help="R1 FASTQ (.gz accepted)")
    parser.add_argument("r2", help="R2 FASTQ (.gz accepted)")
    parser.add_argument("-o", "--output", help="Optional output file")
    parser.add_argument("-s", "--sample-size", type=int, default=1000,
                        help="Pairs to sample for orientation")
    parser.add_argument("--quiet", action="store_true",
                        help="Disable per-read output")
    parser.add_argument("--progress", action="store_true",
                        help="Print live progress (pairs/s)")
    args = parser.parse_args()

    pipeline = RepeatPhasingPipeline(args.r1, args.r2, sample_size=args.sample_size)
    out_fh = None
    if not args.quiet:
        out_fh = open(args.output, "w") if args.output else sys.stdout

    t0 = perf_counter()
    count = 0
    total_pair_time = 0.0
    progress = tqdm(disable=not args.progress, unit="pairs", leave=False)

    for result in pipeline.run():
        count += 1
        total_pair_time += result.elapsed
        if out_fh:
            line = (
                f"{result.read_id}\t{result.phase_shift}\t{result.consensus_len}\t"
                f"{result.elapsed:.6f}"
            )
            print(line, file=out_fh)

        progress.update(1)
    progress.close()

    if out_fh and out_fh is not sys.stdout:
        out_fh.close()

    duration = perf_counter() - t0
    rate = count / duration if duration else 0
    avg_pair = total_pair_time / count if count else 0.0
    print(
        f"Processed {count} results in {duration:.2f}s ({rate:.1f} pairs/s,"
        f" avg {avg_pair:.4f}s per pair)"
    )


if __name__ == "__main__":
    main()

