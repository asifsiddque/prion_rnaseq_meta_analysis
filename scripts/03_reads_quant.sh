#!/usr/bin/env bash
# =============================================================================
# 03_reads_quant.sh
# Run Salmon quantification for all prion disease RNA-seq datasets
#
# Usage:
#   bash 03_reads_quant.sh
#
# Input:
#   ../data/raw/fastq/<GSE>/<SRR>/              — per-sample FASTQ directories from 00_download_fastq.sh
#   ../resources/genome_indices/mouse_salmon_index/ — built by 02_build_salmon_index.sh
#   ../resources/genome_indices/human_salmon_index/ — built by 02_build_salmon_index.sh
#   ../config/sample_metadata.tsv — tab-separated file mapping SRR -> GSE -> layout -> organism
#                                   (see format below)
#
# Output:
#   ../results/salmon_quant/<GSE>/<SRR>/quant.sf
#   ../results/salmon_quant/<GSE>/<SRR>/aux_info/
#   ../results/salmon_quant/<GSE>/<SRR>/logs/
#
# sample_metadata.tsv format (tab-separated, with header):
#   SRR         GSE         LAYOUT    ORGANISM
#   SRR19907647 GSE207251   PAIRED    mouse
#   SRR11684314 GSE149805   SINGLE    mouse
#   SRR...      GSE214373   PAIRED    human
#
# Notes:
#   - LAYOUT must be PAIRED or SINGLE (case-insensitive)
#   - ORGANISM must be mouse or human (case-insensitive)
#   - Mouse datasets -> mouse_salmon_index (10 of 11 studies)
#   - Human dataset  -> human_salmon_index (GSE214373 only)
#   - --gcBias is applied to all samples (paired and single-end)
#   - --seqBias is applied to PAIRED-END samples only; omitted for SINGLE-END
#     because Salmon's sequence bias model requires the fragment length
#     distribution, which is unavailable for single-end reads
#     (cohort spans HiSeq 2500/4000, NovaSeq 6000, NextSeq 500/550)
#   - --numBootstraps is NOT used (DESeq2 does not consume bootstrap replicates)
#   - --validateMappings enables selective alignment (more accurate than quasi-mapping)
#   - Skips samples already quantified (resume-safe)
#   - After quantification, check mapping rates in logs:
#       grep "Mapping rate" ../results/salmon_quant/*/*/logs/salmon_quant.log
#     Flag any sample with mapping rate < 70%
# =============================================================================

set -euo pipefail

# --- Locale fix (required for Salmon >= 1.10 in minimal environments) ---
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

# --- Paths ---
FASTQ_DIR="../data/raw/fastq"
MOUSE_INDEX="../resources/genome_indices/mouse_salmon_index"
HUMAN_INDEX="../resources/genome_indices/human_salmon_index"
OUTPUT_DIR="../results/salmon_quant"
METADATA="../config/sample_metadata.tsv"
THREADS=8

mkdir -p "${OUTPUT_DIR}"

# --- Validate inputs ---
if [[ ! -f "${METADATA}" ]]; then
    echo "ERROR: sample_metadata.tsv not found at ${METADATA}"
    echo "Create it with columns: SRR  GSE  LAYOUT  ORGANISM"
    exit 1
fi

if [[ ! -d "${MOUSE_INDEX}" ]]; then
    echo "ERROR: Mouse Salmon index not found at ${MOUSE_INDEX}"
    echo "Run 02_build_salmon_index.sh first."
    exit 1
fi

if [[ ! -d "${HUMAN_INDEX}" ]]; then
    echo "ERROR: Human Salmon index not found at ${HUMAN_INDEX}"
    echo "Run 02_build_salmon_index.sh first."
    exit 1
fi

# --- Process each sample ---
# Skip header line (tail -n +2)
# NOTE: process substitution (not pipe) keeps the loop in the current shell,
# so 'exit 1' on error correctly terminates the whole script.
while IFS=$'\t' read -r SRR GSE LAYOUT ORGANISM; do

    # Skip blank lines
    [[ -z "${SRR}" ]] && continue

    # Normalize case
    LAYOUT="${LAYOUT^^}"       # PAIRED or SINGLE
    ORGANISM="${ORGANISM,,}"   # mouse or human

    SAMPLE_FASTQ_DIR="${FASTQ_DIR}/${GSE}/${SRR}"
    SAMPLE_OUT_DIR="${OUTPUT_DIR}/${GSE}/${SRR}"

    # --- Resume check: skip if quant.sf already exists ---
    if [[ -f "${SAMPLE_OUT_DIR}/quant.sf" ]]; then
        echo "[SKIP] ${SRR} (${GSE}) — quant.sf already present"
        continue
    fi

    echo "============================================"
    echo "[QUANT] ${SRR}  |  ${GSE}  |  ${LAYOUT}  |  ${ORGANISM}"
    echo "============================================"

    # --- Select index based on organism ---
    if [[ "${ORGANISM}" == "human" ]]; then
        SALMON_INDEX="${HUMAN_INDEX}"
    elif [[ "${ORGANISM}" == "mouse" ]]; then
        SALMON_INDEX="${MOUSE_INDEX}"
    else
        echo "ERROR: Unknown organism '${ORGANISM}' for ${SRR}. Must be 'mouse' or 'human'."
        exit 1
    fi

    mkdir -p "${SAMPLE_OUT_DIR}"

    # --- Run Salmon based on library layout ---
    if [[ "${LAYOUT}" == "PAIRED" ]]; then
        R1="${SAMPLE_FASTQ_DIR}/${SRR}_1.fastq.gz"
        R2="${SAMPLE_FASTQ_DIR}/${SRR}_2.fastq.gz"

        if [[ ! -f "${R1}" || ! -f "${R2}" ]]; then
            echo "ERROR: FASTQ files not found for ${SRR}"
            echo "  Expected: ${R1}"
            echo "  Expected: ${R2}"
            exit 1
        fi

        salmon quant \
            --index "${SALMON_INDEX}" \
            --libType A \
            -1 "${R1}" \
            -2 "${R2}" \
            --threads "${THREADS}" \
            --validateMappings \
            --gcBias \
            --seqBias \
            --output "${SAMPLE_OUT_DIR}"

    elif [[ "${LAYOUT}" == "SINGLE" ]]; then
        R="${SAMPLE_FASTQ_DIR}/${SRR}.fastq.gz"

        if [[ ! -f "${R}" ]]; then
            echo "ERROR: FASTQ file not found for ${SRR}"
            echo "  Expected: ${R}"
            exit 1
        fi

        # NOTE: --seqBias is intentionally omitted for single-end reads.
        # Salmon's sequence bias correction relies on the fragment length
        # distribution, which is only available for paired-end data.
        # Applying --seqBias to single-end reads produces unreliable corrections
        # and would introduce inconsistent bias handling across the cohort.
        salmon quant \
            --index "${SALMON_INDEX}" \
            --libType A \
            -r "${R}" \
            --threads "${THREADS}" \
            --validateMappings \
            --gcBias \
            --output "${SAMPLE_OUT_DIR}"

    else
        echo "ERROR: Unknown layout '${LAYOUT}' for ${SRR}. Must be 'PAIRED' or 'SINGLE'."
        exit 1
    fi

    echo "[DONE] ${SRR} -> ${SAMPLE_OUT_DIR}/quant.sf"

done < <(tail -n +2 "${METADATA}")

# --- Summary: print mapping rates for all completed samples ---
echo ""
echo "============================================"
echo "Quantification complete. Mapping rate summary:"
echo "============================================"
find "${OUTPUT_DIR}" -name "salmon_quant.log" | sort | while read -r LOG; do
    RATE=$(grep "Mapping rate" "${LOG}" 2>/dev/null | tail -1 || echo "N/A")
    SAMPLE=$(echo "${LOG}" | awk -F'/' '{print $(NF-2)"/"$(NF-1)}')
    echo "  ${SAMPLE}: ${RATE}"
done

echo ""
echo "Flag any sample with mapping rate < 70% before proceeding to DESeq2."
