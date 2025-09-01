# HOW IT WORKS

This document explains how the development repository’s core modules work and why each step exists. It follows the pipeline flow from paired-end RCA reads → consensus (double consensus) → pseudoalignment → refinement.

---

## Purpose & Scope

**Goal:** Detect true **transcription errors** in RNA by leveraging **Rolling Circle Amplification (RCA)** to reduce sequencing noise. RCA produces tandem repeats of the same molecule within each read; by phasing and collapsing these repeats, we recover a high‑confidence **consensus** for each original transcript. We then anchor that consensus to a reference and refine its exact placement so the resulting coordinates can be used for downstream error analysis.

**Inputs (typical):**

* Paired-end FASTQ reads (R1/R2) generated from RCA libraries.
* A reference FASTA used for mapping (commonly transcriptome‑level; the workflow is agnostic as long as your mapper reference matches your goals).

**Outputs (current):**

* Double‑consensus FASTQ and a metadata file (per‑read consensus length `d` and phase/offset `phi`).
* Refined BAM where each record’s sequence is the **single** consensus window trimmed from the double consensus and placed on the reference (with informative tags).

**Planned:**

* Emission of **absolute references** (e.g., genomic/transcriptomic position, index, and strand orientation) consolidated for downstream transcription‑error calling.

---

## Step 1 — Consensus Builder (Double Consensus)

### 1. Repeat Length Detection (find `d`)

RCA reads contain tandem repeats of the same unit. We estimate the **repeat unit length** `d` by taking the first `k` bases from R1 and scanning forward for the next window within an **error budget** (`max_errors`) that matches this `k`‑mer. The distance from the anchor to that matching window is `d`.

* **Parameters:** `k` (default e.g. 20–25), `max_errors` (default 2).
* **Why this approach?** It is much faster than full alignment and robust to small sequencing errors within the anchor.
* **Failure mode:** If no repeat is detected, the read is skipped and counted as `NoRepeats`.

### 2. Quality‑Weighted Consensus Matrix

For a discovered repeat length `d`, we build a **4×d quality matrix** for bases **A,C,G,T** across positions `0 … d−1`.

* **Update rule:** For each base observed at position `i` of the repeat frame, add its **Phred** quality score to the matrix cell `[base, i]`. Maintain a coverage vector to count contributions per position.
* **Consensus call:** For each column `i`, choose the base with the largest accumulated quality. This yields the **single‑consensus** string and an accompanying list of per‑position quality scores (e.g., best minus non‑best totals).
* **Intuition:** Multiple repeat copies act like technical replicates. Aggregating **quality‑weighted** evidence per frame position suppresses random sequencing errors while preserving real discrepancies.

### 3. Orientation Decision (R2 forward vs reverse‑complement)

R2 may be the reverse complement of R1. We **sample** a subset of read pairs and score R2 in both orientations against the provisional consensus to decide a **global orientation** for the dataset (forward vs RC). This avoids per‑read orientation guessing and improves speed/consistency.

### 4. Phase Alignment (compute `phi`)

Even with correct orientation, R2 is offset relative to the consensus frame. We compute the **phase shift** `phi` by scanning for the best in‑frame alignment within a window of length `k` using a two‑tier strategy:

* **Fast scan:** A mismatch‑budgeted search (Numba‑accelerated) returns the first phase meeting the error threshold.
* **Fallback:** If the fast path fails, perform a mismatch‑minimizing search with a **quality‑aware** tie‑breaker, then normalize `phi` into `[0, d)`.

### 5. Merge & Build the Double Consensus

We rotate R2 by `phi` to align it to the same frame as R1 and update the quality matrix again. We then extract the final **single consensus** and **duplicate it back‑to‑back** to form the **double consensus (DC)**:

```
single_consensus  = S
double_consensus  = S + S
```

**Why double consensus?** Downstream pseudoalignment is local; duplicating the consensus ensures that an in‑frame window matching the true single consensus exists regardless of where the mapper anchors, simplifying later trimming.

### 6. Step‑1 Outputs

* **`output.fastq.gz`** — FASTQ of double‑consensus sequences with qualities (single consensus qualities duplicated).
* **`metadata.txt.gz`** — Tab‑separated: `read_id`, `consensus_length=d`, `phase_shift=phi`.

---

## Step 2 — Pseudoalignment (Kallisto)

We pseudoalign the **double consensus** against a chosen reference using **kallisto** to rapidly identify likely targets and coarse coordinates.

* **Reference choice:** Often a **transcriptome FASTA** (mapping within transcript space); genome‑level is possible but affects interpretation. Choose the reference that matches your downstream coordinate system.
* **Output:** A BAM (e.g., `pseudoalignments.bam`) that anchors each DC to a reference target with a **local** placement.

---

## Step 3 — Refiner (Trim to the Single Consensus & Re‑anchor)

Kallisto’s placement is local and the input is **double** consensus, so we must recover the **single** consensus window and its exact placement.

### Phase 1 — Exact‑Match Near the Anchor

Using kallisto’s reported reference start as a seed, scan the reference sequence for an **exact** `d`‑base match to **any in‑frame window** within the DC. If found, we:

* Crop the DC at the discovered `phi` to the single consensus.
* Update the BAM record to set `reference_start`, sequence, qualities (sliced), and CIGAR to `dM`.

### Phase 2 — Smith–Waterman (Fallback)

If no perfect window is found, perform a local alignment of DC vs reference using a SIMD‑accelerated Smith–Waterman (e.g., Parasail). From the highest‑scoring local segment we:

* Derive the consensus window (`start_query:end_query` on DC) and corresponding reference placement (`beg_ref:end_ref`).
* Compute basic edit metrics (matches, mismatches/indels) for QC.
* Rewrite the BAM with the trimmed **single consensus** and clean CIGAR.

### Output Tags & QC

Each refined read carries informative tags (illustrative):

* `PH`: phase used (1 = exact, 2 = SW fallback, 3 = unmapped)
* `BP`: phase/offset `phi` on DC (query start)
* `CL`: consensus length `d`
* `MT`: matches (Phase 2)
* `NM`: edit distance (Phase 2)
* `CH`: chimeric flag (e.g., >5% non‑matches in window)
* `RC`: whether the original input was marked reverse

**Refined output:** `refined.bam` containing single‑consensus sequences placed on the reference with normalized CIGARs and tags for downstream analysis.

---

## Parameters & Tuning (Highlights)

* `k` — k‑mer/window length for detection and phase alignment (speed/robustness trade‑off).
* `max_errors` — mismatch budget for detection/fast scan.
* `sample_size` — number of pairs used to decide global R2 orientation.
* `max_reads` — cap for test runs and benchmarking.

**Tips:**

* Short `d` (small repeat units) may require a smaller `k` to avoid over‑penalizing mismatches.
* If Phase‑2 fallback is frequent, consider adjusting `k`, `max_errors`, or input quality filters.

---

## Edge Cases & Failure Modes

* **NoRepeats:** Repeat unit not detected → skip read (tracked for run stats).
* **Orientation mis‑call:** Extremely imbalanced/low‑quality samples may confound the global decision; revisit `sample_size` or force orientation for debugging.
* **Low coverage per frame:** If too few repeat copies contribute, consensus quality drops; consider filtering short/low‑quality reads.
* **Chimeric signals:** Mixed templates within an RCA read inflate `CH`; inspect flagged reads downstream.

---

## Coordinate Systems & Future Work

* **Current:** Placement is on the reference used by the mapper (commonly transcriptome). The refiner normalizes to the **single consensus** with clean CIGARs and tags.
* **Planned:** Emit consolidated **absolute positions** (index and strand) consistent with your chosen reference. If genome‑level reporting is required after transcriptome mapping, add a transcript‑to‑genome lift‑over step using annotation.

---

## Glossary

* **`d` (repeat length):** Size of the repeated unit within an RCA read.
* **Consensus matrix:** 4×d accumulator of Phred‑weighted evidence for A/C/G/T at each position; the per‑column argmax yields the base call.
* **`phi` (phase shift):** Offset aligning R2 to the consensus frame.
* **Double consensus (DC):** Single consensus concatenated with itself (`S+S`) to facilitate robust local mapping and later trimming.
