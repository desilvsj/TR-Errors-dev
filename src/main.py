import sys
import argparse
from time import perf_counter
from oop_pipeline import RepeatPhasingPipeline
import numpy as np


def phreds_to_ascii(phreds):
    shift = 33
    chars = []
    for q in phreds:
        # q can be numpy type, ensure python int
        val = int(q)
        # cap within 0..42 before converting to ASCII
        val = max(0, min(val, 42))
        chars.append(chr(val + shift))
    return ''.join(chars)


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
    last_report = t0
    count = 0
    total_pair_time = 0.0
    for result in pipeline.run():
        count += 1
        total_pair_time += result.elapsed
        if out_fh:
            mat_str = np.array2string(result.matrix, separator=",", max_line_width=1000000)
            qual_str = phreds_to_ascii(result.qualities)
            line = (
                f"{result.read_id}\t{result.phase_shift}\t{result.consensus_len}\t"
                f"{result.elapsed:.6f}\t{result.consensus}\t{qual_str}\t{mat_str}"
            )
            print(line, file=out_fh)

        now = perf_counter()
        if args.progress and (now - last_report) >= 1.0:
            rate = count / (now - t0)
            print(f"Processed {count} pairs \u2014 {rate:.1f} pairs/s", file=sys.stderr)
            last_report = now

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

