# TR‑Errors Pipeline

Refine Rolling Circle Amplification (RCA) read alignments in **two stages**. Stage 1 builds a **double‑consensus** FASTQ from paired R1/R2 reads (OOP pipeline). Stage 2 consumes an external aligner’s BAM (e.g., **Kalipso/Kallisto**) and produces a refined BAM where a **single consensus** is placed on the reference using a fast exact pass (Phase‑1) with a SIMD‑accelerated Smith–Waterman fallback (Phase‑2 via Parasail).

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

* **Two‑stage pipeline**

  * **Stage 1 — OOP phasing (**`** / **`**)**: consumes R1/R2 FASTQ (gz OK), detects repeats, and writes a gzipped FASTQ of **double‑consensus** sequences plus a gzipped **metadata** file.
  * **Stage 2 — Refinement (**\`\`**)**: consumes the external aligner’s BAM and places a **single consensus** onto the reference with Phase‑1 exact and Phase‑2 Parasail fallback.
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

### Stage 1 — Build double‑consensus FASTQ (OOP pipeline)

```bash
# R1/R2 can be .fastq or .fastq.gz
python main.py R1.fastq.gz R2.fastq.gz \
  --fastq-out outputs/output.fastq.gz \
  --meta-out  outputs/metadata.txt.gz \
  -n 200000 --progress
```

This writes:

* `outputs/output.fastq.gz` — FASTQ of **double‑consensus** sequences
* `outputs/metadata.txt.gz` — per‑read metadata (read id, consensus length `d`, phase shift `phi`)

### Stage 2 — Align externally, then refine

1. **Align** the Stage‑1 FASTQ with your tool of choice (e.g., **Kalipso/Kallisto** — not bundled here). You are responsible for installing/running it.
2. **Refine** with this repo using the aligner’s BAM plus your transcript FASTA:

```bash
python refiner_v4.py aligner_output.bam transcripts.fasta -o refined -n 200000
samtools view refined.bam | head -n 3
```

---

## Inputs

### Stage 1 (OOP phasing)

* **R1 FASTQ** and **R2 FASTQ** (gzip accepted). These are the left/right RCA tails per read id.
* Outputs a gz FASTQ of **double consensus** and a gz metadata file.

### Stage 2 (Refinement)

* **BAM**: output from your external aligner (e.g., Kalipso/Kallisto) run **against the Stage‑1 FASTQ**. Secondary/supplementary reads are skipped.
* **FASTA**: transcriptome/reference; sequence IDs must match `reference_name` values in the BAM.

Assumptions:

* Each input read for Stage‑2 contains a **double consensus** in its query sequence.
* The aligner’s reported start index is near the true placement (used to anchor Phase‑1 search).

---

## How it works

### Stage 1 — OOP repeat phasing (see `main.py`, `oop_pipeline.py`)

1. Read paired FASTQs (R1/R2) and, per read id, detect periodic repeats.
2. Build a **double‑consensus** sequence (single consensus repeated twice) and compute metadata (`phi`, `d`).
3. Write `outputs/output.fastq.gz` (double consensus as FASTQ) and `outputs/metadata.txt.gz`.

### Stage 2 — Refinement (see `refiner_v4.py`)

1. Load the transcript FASTA and iterate primary alignments from the external aligner’s BAM.
2. **Phase‑1 (exact)**: for `phi ∈ [0..d]`, attempt `ref.find(dc[phi:phi+d], start_index)` near the aligner’s index. On success, crop single consensus, slice qualities, write BAM (`PH=1`).
3. **Phase‑2 (SW fallback)**: run `parasail.sw_trace_striped_16(dc, ref)`. Use `beg_ref` and `beg_query:end_query` to crop single consensus; optionally compute matches/mismatches (`MT`, `NM`, `MM`). Write BAM (`PH=2`).
4. On failure to place, mark unmapped (`PH=4`).

See `/docs/ARCHITECTURE.md` for data flow and invariants.

---

## CLI

### Stage 1 — `main.py`

```bash
python main.py <R1.fastq[.gz]> <R2.fastq[.gz]> [options]

Options:
  --fastq-out   Path to gz FASTQ with double‑consensus (default: outputs/output.fastq.gz)
  --meta-out    Path to gz metadata (default: outputs/metadata.txt.gz)
  -n, --max-reads    Limit pairs processed (testing)
  -s, --sample-size  Pairs to sample for orientation (default: 10000)
  --progress         Live progress (pairs/s)
  --quiet            Disable per‑read output
  --max-errors INT   Max allowed errors in repeat detection (default: 2)
  --k INT            k‑mer length for repeat detection/alignment (default: 25)
```

Refer to `main.py` for the full list and defaults.

### Stage 2 — `refiner_v4.py`

```bash
python refiner_v4.py <aligner_output.bam> <reference.fasta> [options]

Options:
  -o, --output_prefix  Output prefix (default: refined) → writes <prefix>.bam
  -n, --max_reads      Limit number of reads processed (testing)
```

---

## Output & tags

### Stage 1 (OOP phasing)

* `outputs/output.fastq.gz` — FASTQ with **double‑consensus** sequences.
* `outputs/metadata.txt.gz` — Gzipped TSV: `<read_id>	<phi>	<d>` (plus any additional fields you add).

### Stage 2 (Refinement)

The pipeline writes at most one refined alignment per primary input read.

* **Placed via Phase‑1**: `PH=1`, `BP` set, `CL=d`, `RC` from original FLAG.
* **Placed via Phase‑2**: `PH=2`, `BP` from `beg_query`; `MT` provides alignment strength (see `/docs/TAGS.md`). Optional: `NM` (edit distance), `MM` (mismatches only).
* **Unplaced**: `PH=4`, record marked unmapped.

**Default CIGAR policy**: we emit `M`‑only with length `len(single_consensus)`. Indels are not represented in CIGAR; use `MT`/`NM` for accuracy.

Key tags are summarized in `/docs/TAGS.md`.

---

## Design choices

* **Consensus orientation policy**: consensus is always written forward (reference orientation). Original orientation is preserved in `RC`.
* **CIGAR simplification**: IO and compatibility first; tags carry finer alignment details.
* **Chimera flag**: cheap but tunable threshold; adjust per dataset in code if needed.

Rationales and trade‑offs are captured in `/docs/DECISIONS.md`.

---

## Performance

* Stage 1 is linear in read length and fast; tune `k`, `sample_size`, and `max_errors` for your data.
* Phase‑1 in Stage 2 is O(d) on the double consensus and typically dominates throughput when successful.
* Phase‑2 uses Parasail’s striped vectors (SIMD). Consider reducing Phase‑2 frequency by improving Phase‑1 anchoring.
* The CLI prints: load time, per‑phase totals, total runtime, reads/second (counts **all reads touched**, including discarded), and chimera/reversal counts.

Tips:

* Ensure Parasail installs with vectorization (SSE/AVX) on your machine.
* Use `-n` during development to iterate quickly.

---

## Repo structure

```
main.py               # CLI entrypoint for Stage 1 (OOP phasing)
oop_pipeline.py       # OOP repeat phasing logic (consider rename to rca_phasing.py)
refiner_v4.py         # Stage 2 refinement (Phase‑1/Phase‑2, tag writing, timing)
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

* \`\` — install with `pip install parasail`; if wheels aren’t available, build from source.
* **No placements / many PH=4** — check that `reference_name` values in BAM match FASTA IDs.
* **Windows shell issues** — prefer WSL or Git Bash; for grepping in `.gz` files on Windows PowerShell use: `gzip -dc file.fastq.gz | findstr "pattern"`.
* **Unexpected CIGAR/CL** — by design CIGAR is `M`‑only; rely on `MT`/`NM` for alignment strength.

---

## Citations

* Parasail: Jeff Daily (2016). *Parasail: SIMD C library for global, semi‑global, and local pairwise sequence alignments.* BMC Bioinformatics 17(1):81. doi:10.1186/s12859‑016‑0930‑z

---

## License

MIT (or your preferred license). Add a `LICENSE` file at repo root.
