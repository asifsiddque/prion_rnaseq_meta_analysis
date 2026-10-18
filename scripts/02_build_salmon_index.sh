#!/usr/bin/env bash
# =============================================================================
# 02_build_salmon_index.sh
# Build decoy-aware Salmon indices for BOTH mouse (GRCm39) and human (GRCh38)
#
# Usage:
#   bash 02_build_salmon_index.sh
#
# Output:
#   ../resources/genome_indices/mouse_salmon_index/   — used for all mouse datasets (Strata A, B, D)
#   ../resources/genome_indices/human_salmon_index/   — used for human dataset only (Stratum C: GSE214373)
#
# Reference: Ensembl Release 111
#   Mouse : GRCm39 cDNA + ncRNA + genome decoy
#   Human : GRCh38 cDNA + ncRNA + genome decoy
#
# Notes:
#   - Both indices use decoy-aware selective alignment (gentrome approach)
#   - --keepDuplicates retains duplicate transcript sequences (important for
#     multi-copy gene families; prevents silent read loss)
#   - cDNA + ncRNA are combined so lncRNA/miRNA reads do not map spuriously
#     to protein-coding transcripts
#   - The locale fix (LC_ALL=en_US.UTF-8) is required for Salmon >= 1.10 in
#     minimal Docker/conda environments to prevent index write failure
#   - Requires ~20 GB disk space total (mouse ~6 GB, human ~7 GB, temp files)
#   - Requires ~16 GB RAM during index build
# =============================================================================

set -euo pipefail

# --- Locale fix (required for Salmon >= 1.10 in minimal environments) ---
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

THREADS=8
RESOURCES_DIR="../resources/genome_indices"
mkdir -p "${RESOURCES_DIR}"

# =============================================================================
# SECTION 1: MOUSE INDEX (GRCm39, Ensembl r111)
# Used for: GSE149805, GSE198063, GSE144738, GSE189527, GSE201249 (Stratum A)
#           GSE210128 (Stratum B)
#           GSE90977, GSE124932, GSE207251, GSE202275 (Stratum D)
# =============================================================================

echo "============================================"
echo "Building MOUSE index (GRCm39, Ensembl r111)"
echo "============================================"

MOUSE_CDNA_URL="https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/cdna/Mus_musculus.GRCm39.cdna.all.fa.gz"
MOUSE_NCRNA_URL="https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/ncrna/Mus_musculus.GRCm39.ncrna.fa.gz"
MOUSE_GENOME_URL="https://ftp.ensembl.org/pub/release-111/fasta/mus_musculus/dna/Mus_musculus.GRCm39.dna.primary_assembly.fa.gz"

MOUSE_CDNA="${RESOURCES_DIR}/mouse_cdna.fa.gz"
MOUSE_NCRNA="${RESOURCES_DIR}/mouse_ncrna.fa.gz"
MOUSE_GENOME="${RESOURCES_DIR}/mouse_genome.fa.gz"
MOUSE_GENTROME="${RESOURCES_DIR}/mouse_gentrome.fa.gz"
MOUSE_DECOYS="${RESOURCES_DIR}/mouse_decoys.txt"
MOUSE_INDEX="${RESOURCES_DIR}/mouse_salmon_index"

# Download (skip if already present)
[[ -f "${MOUSE_CDNA}" ]]   || wget -c "${MOUSE_CDNA_URL}"   -O "${MOUSE_CDNA}"
[[ -f "${MOUSE_NCRNA}" ]]  || wget -c "${MOUSE_NCRNA_URL}"  -O "${MOUSE_NCRNA}"
[[ -f "${MOUSE_GENOME}" ]] || wget -c "${MOUSE_GENOME_URL}" -O "${MOUSE_GENOME}"

# Build decoy list from genome chromosome names
echo "Building mouse decoy list..."
grep '^>' <(gunzip -c "${MOUSE_GENOME}") | cut -d ' ' -f 1 | sed 's/>//' > "${MOUSE_DECOYS}"

# Build gentrome: transcriptome (cDNA + ncRNA) concatenated with genome
echo "Building mouse gentrome..."
cat "${MOUSE_CDNA}" "${MOUSE_NCRNA}" "${MOUSE_GENOME}" > "${MOUSE_GENTROME}"

# Build Salmon index
echo "Running salmon index (mouse)..."
salmon index \
    --transcripts "${MOUSE_GENTROME}" \
    --decoys "${MOUSE_DECOYS}" \
    --index "${MOUSE_INDEX}" \
    --threads "${THREADS}" \
    --keepDuplicates

echo "Mouse index complete: ${MOUSE_INDEX}/"

# =============================================================================
# SECTION 2: HUMAN INDEX (GRCh38, Ensembl r111)
# Used for: GSE214373 ONLY (Stratum C — human sCJD bulk RNA-seq)
# =============================================================================

echo ""
echo "============================================"
echo "Building HUMAN index (GRCh38, Ensembl r111)"
echo "============================================"

HUMAN_CDNA_URL="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/cdna/Homo_sapiens.GRCh38.cdna.all.fa.gz"
HUMAN_NCRNA_URL="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/ncrna/Homo_sapiens.GRCh38.ncrna.fa.gz"
HUMAN_GENOME_URL="https://ftp.ensembl.org/pub/release-111/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz"

HUMAN_CDNA="${RESOURCES_DIR}/human_cdna.fa.gz"
HUMAN_NCRNA="${RESOURCES_DIR}/human_ncrna.fa.gz"
HUMAN_GENOME="${RESOURCES_DIR}/human_genome.fa.gz"
HUMAN_GENTROME="${RESOURCES_DIR}/human_gentrome.fa.gz"
HUMAN_DECOYS="${RESOURCES_DIR}/human_decoys.txt"
HUMAN_INDEX="${RESOURCES_DIR}/human_salmon_index"

# Download (skip if already present)
[[ -f "${HUMAN_CDNA}" ]]   || wget -c "${HUMAN_CDNA_URL}"   -O "${HUMAN_CDNA}"
[[ -f "${HUMAN_NCRNA}" ]]  || wget -c "${HUMAN_NCRNA_URL}"  -O "${HUMAN_NCRNA}"
[[ -f "${HUMAN_GENOME}" ]] || wget -c "${HUMAN_GENOME_URL}" -O "${HUMAN_GENOME}"

# Build decoy list from genome chromosome names
echo "Building human decoy list..."
grep '^>' <(gunzip -c "${HUMAN_GENOME}") | cut -d ' ' -f 1 | sed 's/>//' > "${HUMAN_DECOYS}"

# Build gentrome: transcriptome (cDNA + ncRNA) concatenated with genome
echo "Building human gentrome..."
cat "${HUMAN_CDNA}" "${HUMAN_NCRNA}" "${HUMAN_GENOME}" > "${HUMAN_GENTROME}"

# Build Salmon index
echo "Running salmon index (human)..."
salmon index \
    --transcripts "${HUMAN_GENTROME}" \
    --decoys "${HUMAN_DECOYS}" \
    --index "${HUMAN_INDEX}" \
    --threads "${THREADS}" \
    --keepDuplicates

echo "Human index complete: ${HUMAN_INDEX}/"

# =============================================================================
# Summary
# =============================================================================
echo ""
echo "============================================"
echo "Both indices built successfully."
echo "  Mouse index : ${MOUSE_INDEX}/"
echo "  Human index : ${HUMAN_INDEX}/"
echo ""
echo "Routing:"
echo "  Mouse index -> GSE149805, GSE198063, GSE144738, GSE189527, GSE201249"
echo "                 GSE210128, GSE90977, GSE124932, GSE207251, GSE202275"
echo "  Human index -> GSE214373 ONLY"
echo "============================================"
