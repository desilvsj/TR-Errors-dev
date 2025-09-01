Install

Python ≥ 3.9

Dependencies: pysam, biopython, parasail, numpy, tqdm (and numba if using other modules).

python -m venv .venv && source .venv/bin/activate
pip install -r requirements.txt
# or
pip install pysam biopython parasail numpy tqdm

Parasail: use pip install parasail for prebuilt wheels when available; otherwise build from source.

Run

python refiner_v4.py <input.bam> <reference.fasta> -o refined -n 200000

-o → output prefix; writes refined.bam.

-n → (optional) limit reads processed (for testing).

Output

Each primary input alignment yields 0 or 1 refined BAM records.

Placed via Phase‑1: PH=1, BP set, CL=d, RC reflects original orientation.

Placed via Phase‑2: PH=2, BP set from parasail query start. MT is matches (see below).

Unplaced: PH=4 and is_unmapped=True.

CIGAR policy (default): we emit M‑only with length = consensus length. Use MT/NM tags for detailed accuracy.

Interpreting tags quickly

Was this Phase‑1 or Phase‑2? → PH

Where did we crop the double consensus? → BP

How long is the single consensus? → CL

How many matches did Phase‑2 achieve? → MT (query span by default; see DECISIONS.md for options)

Is it chimeric by heuristic? → CH ((CL−MT)/CL > 0.05 → 1)

What was the original orientation? → RC

Example

samtools view refined.bam | head -n 1

Then cross‑reference the fields using /docs/TAGS.md.