FROM python:3.11-slim

RUN apt-get update && apt-get install -y --no-install-recommends \
      samtools kallisto pigz gcc g++ make curl ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir biopython numpy numba tqdm pysam edlib parasail

RUN useradd -m -u 1000 appuser
USER appuser
WORKDIR /work
