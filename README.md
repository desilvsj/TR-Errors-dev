# TR‑Errors-dev

## Overview
This project develops a modernized pipeline to detect transcription errors from Rolling Circle Amplification (RCA) sequencing data. Transcription errors occur far less frequently than sequencing errors, making them difficult to identify with standard methods. RCA generates tandem repeats of the same RNA molecule, allowing sequencing noise to be averaged out and the true consensus sequence to be recovered.

The pipeline takes paired-end RCA reads, builds a consensus verified across both R1 and R2, and produces a “double consensus” sequence for mapping. Using kallisto, these sequences are pseudoaligned to a reference, then refined to locate the exact consensus window and its true placement. The final goal is to return precise genomic or transcriptomic positions where transcription errors can be confidently identified. For further explanation, please refer [HOW-IT-WORKS.md](/docs/HOW-IT-WORKS.md).

This repository contains the development implementation of the Transcription (TR) Error Pipeline. Unlike the companion Nextflow repo [(TR-Errors-Pipeline)](https://github.com/desilvsj/TR-Errors-Pipeline), which provides a production-ready workflow, this repo hosts the raw Python scripts, prototypes, and experimental modules. This repository can be used to explore how the core algorithms work, run modules directly, and test out new ideas before integration into the formal pipeline.

---

Refine Rolling Circle Amplification (RCA) read alignments in **three stages**:

* **Stage 1 — Consensus Builder**: take paired R1/R2 FASTQ (\*.fastq or \*.fastq.gz), detect repeats, and emit a gzipped FASTQ of **double‑consensus** reads plus a gzipped **metadata** file.
* **Step 2 — Alignment (Kallisto)**: map the Step-1 FASTQ using Kallisto to produce a BAM of pseudoalignments.
* **Stage 3 — Alignment Refinement**: refine those placements by cropping to a **single consensus** via a fast exact search (Phase‑1) with a SIMD‑accelerated Smith–Waterman fallback (Phase‑2 via Parasail).

> The external alignment (Step 2) is **not** bundled here; you install/run it yourself. See `main.py` for Stage‑1 arguments.

---

## Features
* **Reference‑oriented output**: consensus is written 5′→3′ in the reference orientation; original input orientation is preserved in `RC`.
* **Informative BAM tags**: `PH`, `BP`, `CL`, `MT`, `CH`, `RC` (and optional `NM`, `MM`). See `/docs/HOW-IT-WORKS.md`.
* **Simple CIGAR** by default (all `M`); alignment strength captured in tags.
* **Run metrics**: total time, per‑phase time, reads/s (all reads touched), chimera/discard counts.

---

## How it works

Please see [HOW-IT-WORKS.md](/docs/HOW-IT-WORKS.md) for a complete breakdown of the process.

---

## Modules at a Glance

* **consensus\_pipeline.py** → Builds consensus from R1/R2, determines repeat unit, orientation, and produces a double consensus.
* **consensus\_builder.py** → Command-line interface for running the pipeline on FASTQ inputs.
* **refiner.py** → Refines double consensus alignments vs reference using exact-match and Smith–Waterman fallback.

Detailed descriptions of how each module works are in [`HOW-IT-WORKS.md`](/docs/HOW-IT-WORKS.md).

---

## Install

Requirements: Python ≥ 3.9. Linux/macOS recommended; Windows via WSL.

```bash
python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# or
pip install pysam biopython parasail numpy tqdm
```

> Parasail: `pip install parasail` uses prebuilt wheels where available; otherwise build from source.

---


## Quickstart

### Prerequisites

* Python 3.10+
* Conda or pip/uv for environment management
* Recommended packages listed in `requirements.txt`

### Setup

```bash
# Clone repo
git clone <repo-url>
cd <repo-name>

# Create environment (example with conda)
conda create -n rca-dev python=3.10
conda activate rca-dev
pip install -r requirements.txt
```

### Running the Code

* **Consensus Builder** (Step 1):

```bash
python consensus_builder.py R1.fastq.gz R2.fastq.gz \
    --fastq-out outputs/output.fastq.gz \
    --meta-out outputs/metadata.txt.gz \
    -n 1000
```

* **Refiner** (Step 3, after Kallisto):

```bash
python refiner.py pseudoalignments.bam reference.fa -o refined
```

Outputs include:

* Double consensus FASTQ (`output.fastq.gz`)
* Metadata with consensus length + phase shift (`metadata.txt.gz`)
* Refined BAM (`refined.bam`)

---

## CLI

### Step 1 — `consensus_builder.py`

```bash
python consensus_builder.py <R1.fastq[.gz]> <R2.fastq[.gz]> [options]

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

Refer to `consensus_builder.py` for the full list and defaults.

### Step 3 — `refiner.py`

```bash
python refiner.py <aligner_output.bam> <reference.fasta> [options]

Options:
  -o, --output_prefix  Output prefix (default: refined) → writes <prefix>.bam
  -n, --max_reads      Limit number of reads processed (testing)
```

---

## Inputs

### Step 1 — Consensus Builder

* **R1 FASTQ** and **R2 FASTQ** (gzip accepted). These are the left/right RCA tails per read id.
* Outputs a gz FASTQ of **double consensus** and a gz metadata file.

### Step 2 — Alignment (Kallisto)

* **Double consensus FASTQ**: output from Step‑1, aligned externally with **Kallisto** against your transcriptome/reference FASTA.
* Produces a BAM of pseudoalignments to be passed to Step‑3.

### Step 3 — Alignment Refinement

* **BAM**: output from Kallisto run **against the Step‑1 FASTQ**. Secondary/supplementary reads are skipped.
* **FASTA**: transcriptome/reference; sequence IDs must match `reference_name` values in the BAM.

Assumptions:

* Each Step‑3 input read contains a **double consensus** in its query sequence.
* The aligner’s reported start index is near the true placement (used to anchor Phase‑1 search).

---

## Outputs & Tags

### Step 1 — Consensus Builder

* `outputs/output.fastq.gz` — FASTQ with **double‑consensus** sequences.
* `outputs/metadata.txt.gz` — Gzipped TSV: `<read_id>\t<phi>\t<d>` (plus any additional fields you add).

### Step 3 — Alignment Refinement

The pipeline writes at most one refined alignment per primary input read.

* **Placed via Phase‑1**: `PH=1`, `BP` set, `CL=d`, `RC` from original FLAG.
* **Placed via Phase‑2**: `PH=2`, `BP` from `beg_query`; `MT` provides alignment strength. Optional: `NM` (edit distance), `MM` (mismatches only).
* **Unplaced**: `PH=4`, record marked unmapped.

**Default CIGAR policy**: we emit `M`‑only with length `len(single_consensus)`. Indels are not represented in CIGAR; use `MT`/`NM` for accuracy.

Key tags are summarized in `/docs/HOW-IT-WORKS`.

---

## Design choices & Notes

* **Consensus orientation policy**: consensus is always written forward (reference orientation). Original orientation is preserved in `RC`.
* **CIGAR simplification**: IO and compatibility first; tags carry finer alignment details.
* **Chimera flag**: cheap but tunable threshold; adjust per dataset in code if needed.
* Parameters like `k` (k-mer size), `max_errors`, and `sample_size` control sensitivity and performance.
* Assumes paired-end RCA reads; orientation determined once per sample.
* Refiner supports fallback alignment with Parasail when exact matches are not found.

---

## Performance

* Stage 1 is linear in read length and fast; tune `k`, `sample_size`, and `max_errors` for your data.
* Phase‑1 in Stage 3 is O(d) on the double consensus and typically dominates throughput when successful.
* Phase‑2 uses Parasail’s striped vectors (SIMD). Consider reducing Phase‑2 frequency by improving Phase‑1 anchoring.
* The CLI prints: load time, per‑phase totals, total runtime, reads/second (counts **all reads touched**, including discarded), and chimera/reversal counts.

Tips:

* Ensure Parasail installs with vectorization (SSE/AVX) on your machine.
* Use `-n` during development to iterate quickly.

---

## Troubleshooting

* **`ModuleNotFoundError: parasail`** — install with `pip install parasail`; if wheels aren’t available, build from source.
* **No placements / many PH=4** — ensure `reference_name` values in BAM match FASTA IDs.
* **Windows shell issues** — prefer WSL or Git Bash; on PowerShell: `gzip -dc file.fastq.gz | findstr "pattern"`.
* **Unexpected CIGAR/CL** — by design CIGAR is `M`‑only; rely on `MT`/`NM` for alignment strength.

---

## Citations

* Parasail: Jeff Daily (2016). *Parasail: SIMD C library for global, semi‑global, and local pairwise sequence alignments.* BMC Bioinformatics 17(1):81. doi:10.1186/s12859‑016‑0930‑z

---

## References

* Illumina Quality Scores: https://www.drive5.com/usearch/manual/quality_score.html
* Understanding Illumina Quality Scores: https://www.illumina.com/content/dam/illumina-marketing/documents/products/technotes/technote_understanding_quality_scores.pdf
* Kallisto: https://pachterlab.github.io/kallisto/about
