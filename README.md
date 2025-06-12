# findRepeatsV5.py

## Overview

`findRepeatsV5.py` implements a high-throughput pipeline for detecting and phasing tandem repeats in paired-end FastQ data (including `.gz`). It builds a per-cycle consensus from the long RCA reads (R1), determines global orientation of the partner reads (R2), and merges both strands into a final high-quality consensus.

An object oriented refactor of this logic is available in `src/oop_pipeline.py`.
This module exposes small classes such as `FastqStream`, `ConsensusMatrix` and
`RepeatPhasingPipeline` for easier extension and testing.

## Features

* **Automatic orientation**: samples the first N read-pairs to decide if R2 needs reverse-complementing
* **Vectorized consensus**: uses NumPy to build and merge 4×d quality-weighted matrices for speed
* **Flexible repeat detection**: sliding k-mer search (`find_repeat_dist`) with bounded mismatches
* **Streaming processing**: handles millions of reads via Biopython’s `SeqIO.parse`
* **Progress & timing hooks**: optional `--progress` shows a tqdm progress bar with
  current rate; use `--quiet` to suppress per-read output

## Requirements

* Python 3.7+
* [NumPy](https://numpy.org/)
* [Biopython](https://biopython.org/)
* [tqdm](https://github.com/tqdm/tqdm) *(optional, for live progress bar)*
* Standard library: `gzip`, `time`, `datetime`, `sys`

## Installation

```bash
pip install numpy biopython tqdm
```

## Usage

```bash
python findRepeatsV5.py <reads_R1.fastq(.gz)> <reads_R2.fastq(.gz)>
```

* **R1/R2 paths**: can be plain `.fastq` or gzipped (`.gz`)
* Output is printed to `stdout` (ID, repeat distance `d`, phase shift φ, consensus sequence, quality string)

### Running the OOP pipeline

```
python -m src.main <R1.fastq(.gz)> <R2.fastq(.gz)> [-o results.txt] [--sample-size N] [--quiet] [--progress]
```

The script iterates over all read pairs using `RepeatPhasingPipeline` and writes one line per processed pair unless
`--quiet` is given. Each line contains the read ID, the inferred phase shift, the consensus length, the processing time,
the consensus sequence and its quality string. When `--progress` is enabled a `tqdm` progress bar shows the current
processing rate. If `-o/--output` is omitted, results are printed to the console. When finished a summary reports total
runtime and the processing rate in pairs per second.

## Configuration

* **Sample size** (orientation): change `sample_size` in `decide_global_orientation()`
* **Report interval**: pass `report_every=<int>` to `process_all_pairs()` (e.g. every 100 000 pairs)
* **TQDM**: wrap the read-pair iterator in `tqdm(...)` for ETA, rate, elapsed

## High-Level Flow

1. **Initialization** (`__init__`)

   * Build `base2idx` (ASCII→{0\:A,1\:C,2\:G,3\:T,4\:N}) and `idx2base` vectors
2. **Orientation Sampling** (`decide_global_orientation`)

   * Parse first N pairs
   * For each:

     * Detect repeat distance `d` via `find_repeat_dist(seq, k_min, k_max)`
     * Build R1 consensus (`consensus_from_long_read_numpy`)
     * Score R2 in forward vs. RC (`score_orientation`→`find_best_shift`)
   * Majority vote → global orientation
3. **Full Processing** (`process_all_pairs`)

   * Iterate all pairs (with optional `tqdm`)
   * Per pair:

     1. `d = find_repeat_dist(r1.seq, k_min, k_max)`
     2. `cons_idx, cons_qual = consensus_from_long_read_numpy(r1.seq, r1.qual, d)`
     3. Orient R2 according to global setting
     4. φ = `find_best_shift(cons_idx, r2_idx, d)`
     5. Merge sum4 matrices & derive final consensus + quality
     6. Print or collect results
4. **Helpers & Logging**

   * `phreds_to_ascii_string(qualities)`: caps at 42 → ASCII (shift 33)
   * Timing: wrap calls in `time.perf_counter()` or inject into loop with `datetime.now()` logs

## Core Methods

### `decide_global_orientation(r1_path, r2_path, sample_size=10000)`

* Samples first `sample_size` read-pairs
* Uses `find_repeat_dist` + consensus + orientation scoring
* Returns `"forward"` or `"RC"`

### `process_all_pairs(r1_path, r2_path, orient, report_every=100000)`

* Main loop over all read-pairs
* Accepts `orient` from orientation step
* Optional `report_every` to print `"[YYYY-MM-DD HH:MM:SS] processed N pairs — elapsed XXs"`

### `consensus_from_long_read_numpy(long_seq: str, long_qual: str, d: int)`

* Vectorizes bases → indices, qualities → weights
* Builds 4×d matrix of summed weights per base per cycle
* `best_idx = np.argmax(sum4, axis=0)` → numeric consensus
* Maps `best_idx` → bases (string) and quality scores

### `find_repeat_dist(sequence: str, k_min: int, k_max: int=0)`

* Slides k-mers to detect repeat unit length `d`
* Allows bounded mismatches via `compare_string(seq1, seq2, max_diff)`

### `find_best_shift(cons_idx: ndarray, read_idx: ndarray, d: int)`

* Scores all cyclic shifts of the read against the consensus index
* Returns phase shift φ maximizing quality-weighted matches

## Timing & Progress

* **Entire run**: wrap `process_all_pairs()` in:

  ```python
  t0 = perf_counter()
  rf.process_all_pairs(...)
  print(f"Total runtime: {perf_counter() - t0:.1f}s")
  ```
* **With `tqdm`**: replace loop with:

  ```python
  for idx, (r1, r2) in enumerate(tqdm(it, desc="Processing", unit="pairs"), start=1):
      ...
      if idx % report_every == 0:
          tqdm.write(...)
  ```

## Customization & Extension

* **Adjust k-mer bounds**: modify `k_min`, `k_max` in `find_repeat_dist`
* **Change mismatch tolerance** in `compare_string`
* **Enable logging**: swap `print` → Python `logging` for levels/handlers
* **Output formats**: adapt print statements to CSV, JSON, or write to file

## Troubleshooting

* **Slow startup**: sample orientation on smaller `sample_size`
* **Memory spikes**: ensure NumPy arrays are reused or process in chunks
* **Gzip hangs**: verify file integrity; consider using buffered reads

---

**Author**: Janesh De Silva
**License**: MIT
**Contact**: [janeshdsilva@gmail.com](mailto:janeshdsilva@gmail.com)
