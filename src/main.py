import sys
import argparse
from time import perf_counter
from oop_pipeline import RepeatPhasingPipeline
from tqdm import tqdm


def main():
    parser = argparse.ArgumentParser(description="Run OOP repeat phasing pipeline")
    parser.add_argument("r1", help="R1 FASTQ (.gz accepted)")
    parser.add_argument("r2", help="R2 FASTQ (.gz accepted)")
    parser.add_argument(
        "--fastq-out",
        default="output.fastq.gz",
        help="Path to output gzipped FASTQ file containing the double-consensus sequences (default: output.fastq.gz)",
    )
    parser.add_argument(
        "--meta-out",
        default="metadata.txt.gz",
        help="Path to output gzipped metadata file with read ID, consensus length, and phase shift (default: metadata.txt.gz)",
    )
    parser.add_argument("-o", "--output", help="Optional output file")
    parser.add_argument("-n", "--max-reads", type=int, default=None,
                        help="Maximum number of read pairs to process")
    parser.add_argument("-s", "--sample-size", type=int, default=10000,
                        help="Pairs to sample for orientation")
    parser.add_argument("--quiet", action="store_true",
                        help="Disable per-read output")
    parser.add_argument("--progress", action="store_true",
                        help="Print live progress (pairs/s)")
    parser.add_argument(
        "--max-errors", type=int, default=2,
        help="Maximum allowed errors in repeat detection"
    )
    parser.add_argument(
        "--k", type=int, default=25,
        help="k-mer length for repeat detection and alignment"
    )
    args = parser.parse_args()

    pipeline = RepeatPhasingPipeline(
        args.r1,
        args.r2,
        sample_size=args.sample_size,
        k=args.k,
        max_errors=args.max_errors,
        max_reads=args.max_reads,
        fastq_out=args.fastq_out,
        meta_out=args.meta_out,
    )

    # out_fh = None
    # if not args.quiet:
    #     out_fh = open(args.output, "w") if args.output else sys.stdout

    t0 = perf_counter()
    count = 0
    total_pair_time = 0.0
    progress = tqdm(total=args.max_reads, disable=not args.progress, unit="pairs", leave=False)

    # for result in pipeline.run():
    #     count += 1
    #     total_pair_time += result.elapsed
    #     if out_fh:
    #         line = (
    #             f"{result.read_id}\t{result.phase_shift}\t{result.consensus_len}\t"
    #             f"{result.elapsed:.6f}"
    #         )
    #         print(line, file=out_fh)

    #     progress.update(1)
    # progress.close()

    # if out_fh and out_fh is not sys.stdout:
    #     out_fh.close()

    for result in pipeline.run():
        count += 1
        total_pair_time += result.elapsed
        progress.update(1)

    duration = perf_counter() - t0
    rate = count / duration if duration else 0
    avg_pair = total_pair_time / count if count else 0.0
    print(
        f"Processed {count} results in {duration:.2f}s ({rate:.1f} pairs/s,"
        f" avg {avg_pair:.4f}s per pair)"
        f"\nPercentage of Repeat Reads: {(pipeline.total_reads-pipeline.no_repeat_count)*100/pipeline.total_reads}"
        f"\nNumber of reads that used Aligner Backup: {pipeline.aligner.fallback_count}"
    )

if __name__ == "__main__":
    main()
