#!/usr/bin/env bash
# =============================================================================
# 01_qc.sh
# Run per-sample quality control with FastQC, then aggregate with MultiQC
#
# Usage:
#   bash 01_qc.sh
#
# Input:
#   ../data/raw/fastq/<GSE>/<SRR>/*.fastq.gz   — all FASTQ files (paired or single-end)
#
# Output:
#   ../results/qc/fastqc/<GSE>/       — per-sample FastQC HTML + zip reports
#   ../results/qc/multiqc/<GSE>/multiqc_report.html  — aggregated MultiQC report
#
# Notes:
#   - FastQC is run on ALL .fastq.gz files found recursively under ../fastq/
#   - MultiQC is run AFTER all FastQC jobs complete (not per-sample)
#   - Parallel execution: up to THREADS FastQC jobs at once
#   - Key QC thresholds to check in the MultiQC report before proceeding:
#       * Mapping rate (checked later in Salmon logs) >= 70%
#       * Per-base sequence quality: median Phred >= 28
#       * Adapter content: < 5% (if high, add Trimmomatic/fastp step)
#       * Duplicate rate: < 50%
#       * rRNA contamination: check overrepresented sequences
# =============================================================================

set -euo pipefail

FASTQ_DIR="../data/raw/fastq"
QC_DIR="../results/qc/fastqc"
MULTIQC_DIR="../results/qc/multiqc"
THREADS=8

mkdir -p "${QC_DIR}"

# Collect all FASTQ files recursively
FASTQ_FILES=( $(find "${FASTQ_DIR}" -name "*.fastq.gz" | sort) )

if [[ ${#FASTQ_FILES[@]} -eq 0 ]]; then
    echo "ERROR: No .fastq.gz files found under ${FASTQ_DIR}"
    exit 1
fi

echo "Found ${#FASTQ_FILES[@]} FASTQ files. Running FastQC..."

# Run FastQC in parallel (xargs -P for parallel jobs)
printf '%s\n' "${FASTQ_FILES[@]}" | \
    xargs -P "${THREADS}" -I{} \
    fastqc {} --outdir "${QC_DIR}" --threads 1 --quiet

echo "FastQC complete. Running MultiQC..."

multiqc "${QC_DIR}" \
    --outdir "${MULTIQC_DIR}" \
    --filename "multiqc_report" \
    --force

echo ""
echo "QC complete!"
echo "  FastQC reports : ${QC_DIR}/"
echo "  MultiQC report : ${MULTIQC_DIR}/multiqc_report.html"
echo ""
echo "Review the MultiQC report before proceeding to quantification."
echo "Flag any samples with:"
echo "  - Per-base quality median < 28"
echo "  - Adapter content > 5%"
echo "  - Duplicate rate > 50%"
echo "  - Overrepresented sequences suggesting rRNA contamination"
