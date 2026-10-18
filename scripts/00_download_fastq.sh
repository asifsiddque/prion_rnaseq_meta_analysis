#!/usr/bin/env bash
# =============================================================================
# 00_download_fastq.sh
# Download FASTQ files for all prion disease RNA-seq datasets
#
# Usage:
#   bash 00_download_fastq.sh
#
# Input:
#   ../config/sample_metadata.tsv — tab-separated file with columns:
#       SRR  GSE  LAYOUT  ORGANISM
#   (867 rows; LAYOUT is PAIRED or SINGLE per-sample)
#
# Output:
#   ../data/raw/fastq/<GSE>/<SRR>/<SRR>.fastq.gz        — single-end samples
#   ../data/raw/fastq/<GSE>/<SRR>/<SRR>_1.fastq.gz      — paired-end read 1
#   ../data/raw/fastq/<GSE>/<SRR>/<SRR>_2.fastq.gz      — paired-end read 2
#
# Notes:
#   - Uses --split-3 (NOT --split-files) so that:
#       SINGLE-end → <SRR>.fastq.gz          (no numeric suffix)
#       PAIRED-end → <SRR>_1.fastq.gz + <SRR>_2.fastq.gz
#     This matches the filename convention expected by 03_reads_quant.sh.
#     --split-files always appends _1 even for single-end reads, which would
#     cause 03_reads_quant.sh to fail to find <SRR>.fastq.gz.
#   - Layout is read per-SRR from sample_metadata.tsv (not assumed globally).
#     GSE144738 is mixed: 131 PAIRED + 6 SINGLE — handled correctly.
#   - Resume-safe: skips any SRR whose expected FASTQ(s) already exist.
#   - Requires: parallel-fastq-dump (sra-tools conda package)
#
# Dataset layout summary (from sample_metadata.tsv):
#   GSE149805 : SINGLE  (N=186)
#   GSE198063 : PAIRED  (N=159)
#   GSE144738 : MIXED   (PAIRED=131, SINGLE=6)   ← per-SRR layout used
#   GSE189527 : SINGLE  (N=98)
#   GSE201249 : PAIRED  (N=50)
#   GSE210128 : PAIRED  (N=173)
#   GSE214373 : PAIRED  (N=14)
#   GSE90977  : SINGLE  (N=16)
#   GSE124932 : PAIRED  (N=16)
#   GSE207251 : PAIRED  (N=6)
#   GSE202275 : PAIRED  (N=12)
# =============================================================================

set -euo pipefail

METADATA="../config/sample_metadata.tsv"
FASTQ_DIR="../data/raw/fastq"
THREADS=4

# --- Validate input ---
if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: sample_metadata.tsv not found at ${METADATA}"
    echo "Expected columns: SRR  GSE  LAYOUT  ORGANISM"
    exit 1
fi

echo "============================================"
echo "  00_download_fastq.sh — FASTQ Download"
echo "  Metadata : ${METADATA}"
echo "  Output   : ${FASTQ_DIR}/<GSE>/<SRR>/"
echo "  Tool     : parallel-fastq-dump --split-3"
echo "============================================"
echo ""

# Skip header line; process each SRR
# NOTE: process substitution (not pipe) keeps the loop in the current shell,
# so 'exit 1' on error correctly terminates the whole script.
while IFS=$'\t' read -r SRR GSE LAYOUT ORGANISM; do

    [[ -z "${SRR}" ]] && continue

    LAYOUT="${LAYOUT^^}"   # normalise to uppercase: PAIRED or SINGLE
    SAMPLE_DIR="${FASTQ_DIR}/${GSE}/${SRR}"

    # --- Resume check: skip if expected FASTQ(s) already present ---
    if [[ "${LAYOUT}" == "PAIRED" ]]; then
        if [[ -f "${SAMPLE_DIR}/${SRR}_1.fastq.gz" && \
              -f "${SAMPLE_DIR}/${SRR}_2.fastq.gz" ]]; then
            echo "[SKIP] ${SRR} (${GSE}, ${LAYOUT}) — FASTQs already present"
            continue
        fi
    elif [[ "${LAYOUT}" == "SINGLE" ]]; then
        if [[ -f "${SAMPLE_DIR}/${SRR}.fastq.gz" ]]; then
            echo "[SKIP] ${SRR} (${GSE}, ${LAYOUT}) — FASTQ already present"
            continue
        fi
    else
        echo "ERROR: Unknown LAYOUT '${LAYOUT}' for ${SRR}. Must be PAIRED or SINGLE."
        exit 1
    fi

    mkdir -p "${SAMPLE_DIR}"

    echo "[DL] ${SRR}  |  ${GSE}  |  ${LAYOUT}"

    # --- Download with --split-3 ---
    # --split-3 behaviour:
    #   PAIRED-end SRA: writes <SRR>_1.fastq.gz + <SRR>_2.fastq.gz
    #   SINGLE-end SRA: writes <SRR>.fastq.gz   (no numeric suffix)
    # This is the correct convention for 03_reads_quant.sh.
    parallel-fastq-dump \
        --sra-id  "${SRR}" \
        --threads "${THREADS}" \
        --outdir  "${SAMPLE_DIR}" \
        --split-3 \
        --gzip

    # --- Verify output ---
    if [[ "${LAYOUT}" == "PAIRED" ]]; then
        if [[ ! -f "${SAMPLE_DIR}/${SRR}_1.fastq.gz" || \
              ! -f "${SAMPLE_DIR}/${SRR}_2.fastq.gz" ]]; then
            echo "ERROR: Expected paired FASTQs not found after download for ${SRR}"
            echo "  Expected: ${SAMPLE_DIR}/${SRR}_1.fastq.gz"
            echo "  Expected: ${SAMPLE_DIR}/${SRR}_2.fastq.gz"
            exit 1
        fi
        echo "[OK]  ${SRR} → ${SAMPLE_DIR}/${SRR}_1.fastq.gz + _2.fastq.gz"
    else
        if [[ ! -f "${SAMPLE_DIR}/${SRR}.fastq.gz" ]]; then
            echo "ERROR: Expected single-end FASTQ not found after download for ${SRR}"
            echo "  Expected: ${SAMPLE_DIR}/${SRR}.fastq.gz"
            exit 1
        fi
        echo "[OK]  ${SRR} → ${SAMPLE_DIR}/${SRR}.fastq.gz"
    fi

done < <(tail -n +2 "${METADATA}")

echo ""
echo "============================================"
echo "  Download complete."
echo "  FASTQs in: ${FASTQ_DIR}/<GSE>/<SRR>/"
echo "  Next step: bash 01_qc.sh"
echo "============================================"
