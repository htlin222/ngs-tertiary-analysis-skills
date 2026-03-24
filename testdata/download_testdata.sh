#!/usr/bin/env bash
# download_testdata.sh — Fetch test BAM files for pipeline validation
#
# Usage: bash testdata/download_testdata.sh [--small|--full]
#
# Options:
#   --small   Download small test BAM from hartwigmedical (default, ~50MB)
#   --full    Download SEQC2 TSO500 BAM from SRA (~5GB, requires sra-tools)

set -euo pipefail
cd "$(dirname "$0")"

MODE="${1:---small}"

case "$MODE" in
  --small)
    echo "==> Downloading small test BAM from hartwigmedical/testdata..."
    echo "    (100k random reads, suitable for smoke tests)"

    if ! command -v curl &>/dev/null; then
      echo "ERROR: curl is required" >&2
      exit 1
    fi

    # Clone minimal test data from Hartwig Medical Foundation
    if [ ! -d "hartwigmedical-testdata" ]; then
      git clone --depth 1 https://github.com/hartwigmedical/testdata.git hartwigmedical-testdata
    fi

    # Use the available BAM or create a synthetic one
    echo "==> Test data available in testdata/hartwigmedical-testdata/"
    echo "    Set BAM_PATH to a BAM file from this directory"
    ;;

  --full)
    echo "==> Downloading SEQC2 TSO500 BAM from SRA..."
    echo "    BioProject: PRJNA489865 (HCC1395 tumor cell line)"
    echo "    This will take a while (~5GB)"

    if ! command -v prefetch &>/dev/null; then
      echo "ERROR: sra-tools required. Install with: brew install sratoolkit" >&2
      exit 1
    fi

    if ! command -v samtools &>/dev/null; then
      echo "ERROR: samtools required. Install with: brew install samtools" >&2
      exit 1
    fi

    # SEQC2 TSO500 accession (HCC1395 tumor)
    ACCESSION="SRR8861483"

    echo "  Prefetching ${ACCESSION}..."
    prefetch "$ACCESSION" -O .

    echo "  Converting to BAM..."
    sam-dump "$ACCESSION" | samtools view -bS -o "seqc2_tso500_tumor.bam" -

    echo "  Indexing BAM..."
    samtools index "seqc2_tso500_tumor.bam"

    echo "==> Done! BAM at: testdata/seqc2_tso500_tumor.bam"
    echo "    Run pipeline with:"
    echo "    BAM_PATH=testdata/seqc2_tso500_tumor.bam SAMPLE_ID=SEQC2_HCC1395 Rscript -e 'targets::tar_make()'"
    ;;

  *)
    echo "Usage: $0 [--small|--full]"
    exit 1
    ;;
esac
