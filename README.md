# TR‑Errors Pipeline

Refine Rolling Circle Amplification (RCA) read alignments by cropping a **double** consensus into a **single** consensus and placing it on the reference transcript. The pipeline runs a fast exact search (Phase‑1) and falls back to a SIMD‑accelerated Smith–Waterman (Phase‑2 via Parasail). Output is a BAM with lightweight, interpretable tags.

---

## Table of contents

* [Features](#features)
* [Install](#install)
* [Quick start](#quick-start)
* [Inputs](#inputs)
* [How it works](#how-it-works)
* [CLI](#cli)
* [Output & tags](#output--tags)
* [Design choices](#design-choices)
* [Performance](#performance)
* [Repo structure](#repo-structure)
* [Testing](#testing)
* [Troubleshooting](#troubleshooting)
* [Citations](#citations)
* [License](#license)

---

## Features

* **Two‑phase placement**

  * **Phase‑1**: exact match near the upstream aligner’s index (cheap substring search).
  * **Phase‑2**: local alignment with Parasail’s striped SW (robust fallback).
* **Reference‑oriented consensus**: output consensus is reported 5′→3′ in the same orientation as the reference transcript.
* **Minimal yet rich BAM tags**: `PH`, `BP`, `CL`, `MT`, `CH`, `RC` (and optional `NM`, `MM`). See `/docs/TAGS.md`.
* **Simple CIGAR** by default (all `M`) to keep IO fast; detailed alignment strength captured in tags.
* **Run metrics** printed: total time, per‑phase time, reads/second (all reads touched), counts of chimeras/discards.

---

## Install

Requirements: Python ≥ 3.9. Tested on Linux/macOS; Windows via WSL is recommended.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# or
pip install pysam biopython parasail numpy tqdm
```

> Parasail: `pip install parasail` uses prebuilt wheels on many platforms; otherwise build from source.

---

## Quick start

```bash
# Refine an upstream BAM against a transcript FASTA
python refiner_v4.py input.bam transcripts.fasta -o refined -n 200000

# Inspect a few records
samtools view refined.bam | head -n 3
```

---

## Inputs

* **BAM**: upstream alignments (e.g., kallisto output). Secondary/supplementary reads are skipped.
* **FASTA**: transcriptome/reference; sequence IDs must match `read.reference_name` values from the BAM.

Assumptions:

* Each read’s query sequence in the input BAM contains a **double consensus** (single consensus repeated twice).
* The upstream aligner provides a starting index that’s near the true placement (used to anchor Phase‑1 search).

---

## How it works

1. **Load references** (`SeqIO.to_dict`).
2. **Iterate reads** (primary alignments only). Compute counters and timing.
3. **Phase‑1 (exact)**: for `phi ∈ [0..d]`, try `ref.find(dc[phi:phi+d], kallisto_index)`. On success, crop to single consensus, slice qualities, and write BAM with `PH=1`.
4. **Phase‑2 (SW fallback)**: run `parasail.sw_trace_striped_16(dc, ref)`. Use `beg_ref` and `beg_query:end_query` to crop single consensus. Optionally compute matches/mismatches from traceback and write `PH=2`.
5. **Failure**: mark read unmapped (`PH=4`).

See `/docs/ARCHITECTURE.md` for data flow and invariants.

---

## CLI

```bash
python refiner_v4.py <input.bam> <reference.fasta> [options]

Options:
  -o, --output_prefix  Output prefix (default: refined) → writes <prefix>.bam
  -n, --max_reads      Limit number of reads processed (for testing)
```

**Example**

```bash
python refiner_v4.py kallisto.bam transcripts.fasta -o yeast_refined -n 100000
```

---

## Output & tags

The pipeline writes at most one refined alignment per primary input read.

* **Placed via Phase‑1**: `PH=1`, `BP` set, `CL=d`, `RC` from original FLAG.
* **Placed via Phase‑2**: `PH=2`, `BP` from `beg_query`; `MT` provides alignment strength (see below). Optional: `NM` (edit distance), `MM` (mismatches only).
* **Unplaced**: `PH=4`, record marked unmapped.

**Default CIGAR policy**: we emit `M`‑only with length `len(single_consensus)`. Indels are not represented in CIGAR; use `MT`/`NM` for accuracy.

**Key tags** (see `/docs/TAGS.md` for the full legend):

* `PH` — phase (1 exact, 2 SW, 4 unmapped)
* `BP` — phi offset used to crop single consensus from the double
* `CL` — intended single‑consensus length (`d`)
* `MT` — matches for the placed segment (query span by default; can be true matches if using traceback)
* `CH` — chimera heuristic: `(CL − MT)/CL > 0.05` → `1` else `0`
* `RC` — original orientation from input FLAG (1 reverse, 0 forward)
* `NM` (optional) — SAM edit distance (`mismatches + insertions + deletions`)
* `MM` (optional) — mismatches only (from traceback)

---

## Design choices

* **Consensus orientation policy**: consensus is always written forward (reference orientation). Original orientation is preserved in `RC`.
* **CIGAR simplification**: IO and compatibility first; tags carry finer alignment details.
* **Chimera flag**: cheap but tunable threshold; adjust per dataset in code if needed.

Rationales and trade‑offs are captured in `/docs/DECISIONS.md`.

---

## Performance

* Phase‑1 is O(d) on the double consensus and typically dominates throughput when successful.
* Phase‑2 uses Parasail’s striped vectors (SIMD). Consider reducing Phase‑2 frequency by improving Phase‑1 anchoring.
* The CLI prints: load time, per‑phase totals, total runtime, reads/second (counts **all reads touched**, including discarded), and chimera/reversal counts.

Tips:

* Ensure Parasail installs with vectorization (SSE/AVX) on your machine.
* Use `-n` during development to iterate quickly.

---

## Repo structure

```
refiner_v4.py         # Main pipeline (Phase‑1/Phase‑2, tag writing, timing)
oop_pipeline.py       # OOP utilities / prior components
/docs/                # Human docs
  TAGS.md            # BAM tag legend
  ARCHITECTURE.md    # Data flow, invariants, conventions
  USAGE.md           # CLI and examples
  DECISIONS.md       # ADRs (design choices & trade‑offs)
  LOGGING.md         # Logging conventions
  CONTRIBUTING.md    # How to set up, style, and PR
  STYLE.md           # Code/docstyle conventions
  CHANGELOG.md       # Human‑readable history
  TESTING.md         # Test strategy and fixtures (optional)
```

---

## Testing

We recommend **pytest** with a tiny fixture FASTA/BAM pair.

```bash
pytest -q
```

See `/docs/TESTING.md` for unit and golden end‑to‑end guidance.

---

## Troubleshooting

* **`ModuleNotFoundError: parasail`** — install with `pip install parasail`; if wheels aren’t available, build from source.
* **No placements / many PH=4** — check that `reference_name` values in BAM match FASTA IDs.
* **Windows shell issues** — prefer WSL or Git Bash; for grepping in `.gz` files on Windows PowerShell use: `gzip -dc file.fastq.gz | findstr "pattern"`.
* **Unexpected CIGAR/CL** — by design CIGAR is `M`‑only; rely on `MT`/`NM` for alignment strength.

---

## Citations

* Parasail: Jeff Daily (2016). *Parasail: SIMD C library for global, semi‑global, and local pairwise sequence alignments.* BMC Bioinformatics 17(1):81. doi:10.1186/s12859‑016‑0930‑z

---

## License

MIT (or your preferred license). Add a `LICENSE` file at repo root.
