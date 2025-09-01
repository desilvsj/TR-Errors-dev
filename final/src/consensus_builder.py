"""
Command Line Interface (CLI) for the Repeat Phasing Pipeline.

Purpose:
    Provides a command-line interface to run the repeat phasing pipeline defined
    in `consensus_pipeline.py`. Handles argument parsing, pipeline execution, 
    and reporting of performance statistics.

Usage:
    python consensus_builder.py R1.fastq(.gz) R2.fastq(.gz) [options]

    Options allow controlling input/output paths, sample size for orientation
    detection, max reads processed, k-mer length, and max allowed errors.
"""

import sys
import argparse
from time import perf_counter
from oop_pipeline import RepeatPhasingPipeline
from tqdm import tqdm


def main():
    """
    Parse arguments and execute the repeat phasing pipeline.
    - Purpose: Defines CLI interface for running the phasing pipeline.
    - Defaults: Output paths pre-set to outputs/.
    - Progress Flag: Allows live performance monitoring.
    """
    parser = argparse.ArgumentParser(description="Run OOP repeat phasing pipeline")
    parser.add_argument("r1", help="Path to R1 FASTQ (.gz supported)")
    parser.add_argument("r2", help="Path to R2 FASTQ (.gz supported)")
    parser.add_argument("--fastq-out", default="outputs/output.fastq.gz",
        help="Gzipped FASTQ file for double-consensus sequences")
    parser.add_argument("--meta-out", default="outputs/metadata.txt.gz",
        help="Gzipped metadata file: read_id, consensus length, phase shift")
    parser.add_argument("-n", "--max-reads", type=int, help="Limit number of pairs")
    parser.add_argument("-s", "--sample-size", type=int, default=10000,
        help="Read pairs sampled for deciding R2 orientation")
    parser.add_argument("--quiet", action="store_true",
        help="Disable per-read output printing")
    parser.add_argument("--progress", action="store_true",
        help="Show live progress bar with pairs/sec")
    parser.add_argument("--max-errors", type=int, default=2,
        help="Max allowed mismatches in repeat detection")
    parser.add_argument("--k", type=int, default=25,
    help="k-mer length for repeat detection/alignment")
    args = parser.parse_args()

    # ----------------------- Pipeline Initialization ---------------------------
    """
    - Purpose: Creates pipeline instance with parameters from CLI.
    - Link: This calls the same processing engine annotated in consensus_pipeline.py.
    """
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

    # -------------------------- Excecution & Timing  -----------------------------
    """
    - Purpose: Iterates over results from pipeline, updating count and timing stats.
    - result.elapsed: Time to process one read pair.
    - tqdm: Displays progress if enabled.
    """
    t0 = perf_counter()
    count = 0
    total_pair_time = 0.0
    progress = tqdm(total=args.max_reads, disable=not args.progress, unit="pairs", leave=False)

    for result in pipeline.run():
        count += 1
        total_pair_time += result.elapsed
        progress.update(1)


    # -------------------------- Performance Summary -----------------------------
    """
    Outputs
        1. Total processed pairs and runtime.
        2. Processing rate (pairs/sec).
        3. Average time per pair.
        4. Percentage of reads containing repeats.
        5. Count of reads using PhaseAligner fallback method.
    """
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
