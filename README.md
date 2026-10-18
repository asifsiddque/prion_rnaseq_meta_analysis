# Prion Disease Bulk RNA-Seq Meta-Analysis Pipeline

A five-stream computational pipeline for identifying and validating the
**Mouse Conserved Prion Signature (MCPS)** — a set of genes consistently
dysregulated across CNS prion disease models, validated in human sCJD, and
detectable in peripheral tissues before clinical onset.

---

## Table of Contents

1. [Overview](#overview)
2. [Study Design](#study-design)
3. [Directory Structure](#directory-structure)
4. [Software Requirements](#software-requirements)
5. [Run Order](#run-order)
6. [Script Reference](#script-reference)
   - [00 — Download FASTQs](#00--download-fastqs)
   - [01 — Quality Control](#01--quality-control)
   - [02 — Build Salmon Indices](#02--build-salmon-indices)
   - [03 — Salmon Quantification](#03--salmon-quantification)
   - [04 — Build Coldata](#04--build-coldata)
   - [05 — Stream 1: CNS Meta-Analysis](#05--stream-1-cns-meta-analysis)
   - [06 — Stream 2: Human Validation](#06--stream-2-human-validation)
   - [07 — Stream 3: Peripheral Analysis](#07--stream-3-peripheral-analysis)
   - [08 — Cross-Stream Integration](#08--cross-stream-integration)
   - [09 — Cell-Type Deconvolution](#09--cell-type-deconvolution)
7. [Input/Output Dependency Map](#inputoutput-dependency-map)
8. [Key Parameters](#key-parameters)
9. [Evidence Tier Definitions](#evidence-tier-definitions)
10. [Reproducibility Notes](#reproducibility-notes)
    - [R Package Reproducibility (renv)](#r-package-reproducibility-renv)
    - [Sample Exclusion Documentation](#sample-exclusion-documentation)
11. [Known Limitations](#known-limitations)

---

## Overview

The pipeline processes **11 GEO datasets** (867 samples) across four strata:

| Stratum | Role | Studies | N samples |
|---------|------|---------|-----------|
| A | CNS Discovery (meta-analysis) | GSE149805, GSE198063, GSE144738, GSE189527, GSE201249 | 630 → 516 analyzed (post-exclusion/QC) |
| B | Peripheral longitudinal | GSE210128 | 173 |
| C | Human sCJD validation | GSE214373 | 14 |
| D | CNS held-out replication | GSE90977, GSE124932, GSE207251, GSE202275 | 50 → 36 analyzed (post-exclusion) |

**Three analysis streams run in parallel after Stream 1:**

```
00 → 01 → 02 → 03 → 04 → 05 (Stream 1: CNS meta-analysis)
                              ├── 06 (Stream 2: Human sCJD validation)
                              ├── 07 (Stream 3: Peripheral longitudinal)
                              └── 08 (Cross-stream integration)
```

---

## Study Design

### Stream 1 — CNS Meta-Analysis (05)
- **5 Stratum A studies** processed independently with DESeq2 (study-specific designs)
- **Random-effects meta-analysis** (DerSimonian-Laird) across all 5 studies
- **MCPS defined** as: meta-padj ≤ 0.05, |meta-log2FC| ≥ 0.5, I² ≤ 75%
- **4 Stratum D studies** used for held-out replication (not included in meta-analysis)
- **Leave-one-out sensitivity** analysis confirms MCPS stability

### Stream 2 — Human sCJD Validation (06)
- DESeq2 on GSE214373 (N=14 frontal cortex; 5 sCJD, 9 control)
- Primary design: `~ age + sex + condition` (age confound documented; VIF = 2.06)
- Sensitivity analysis: `~ sex + condition` (quantifies age-driven false positives)
- GSEA tests whether MCPS is enriched in the human sCJD ranked gene list
- Ortholog mapping: human Ensembl → mouse MGI symbols (BioMart, Ensembl r111, cached)

### Stream 3 — Peripheral Longitudinal Analysis (07)
- **Arm A (LRT):** Interaction test `condition × time_point_wpi` per tissue
  - Terminal timepoint excluded from continuous LRT (exact wpi unknown in GEO)
  - Factor-based sensitivity LRT run in parallel
- **Arm B (Wald):** Terminal-timepoint DESeq2 per tissue (apeglm shrinkage)
- **Arm C (MCPS concordance):** Direction-of-effect test + binomial test per tissue
- **First-detection:** Per-timepoint Wald tests with dispersion priors from full model

### Stream 4 — Integration (08)
- Master table: one row per MCPS gene, all evidence columns joined
- Evidence tiers 1–5 (see [Evidence Tier Definitions](#evidence-tier-definitions))
- Cross-stream LFC correlations (7 Pearson r pairs)
- Replication cascade (6 stages)
- Pathway convergence analysis
- **Literature convergence table:** MCPS genes cross-referenced against published prion transcriptomics
- **Pathway-level human validation (GSE124571):** GSEA of 362 Stream 1 consensus pathways against human sCJD brain (N=10+10); 314/362 pathways significant (86.7%), 311/314 with NES > 0

### Stream 5 — Cell-Type Deconvolution (09)
- Marker-gene z-score deconvolution across all 5 Stratum A studies (N=171 samples for violin; N=135 for MCPS-microglia correlation, GSE144738 only)
- 8 cell types: Microglia, A1/A2/pan-Astrocytes, OPC, Neurons, Oligodendrocytes, Endothelial
- MCPS genes excluded from marker gene sets (prevents circular correlation)
- Microglia score: Δ = +0.842, p = 2.45×10⁻⁷ (Wilcoxon, prion vs. control)
- 5/6 MCPS genes correlated with microglia score (Spearman r = 0.47–0.64, all p < 0.001)
- Slc15a4/Tasl: r = 0.061, NS — confirms non-microglial mechanism

---

## Directory Structure

The project follows the scientific community standard layout used by
[Harvard HBC Training](https://hbctraining.github.io/Intro-to-rnaseq-fasrc-salmon-flipped/lessons/04b_data_organization.html),
[Cookiecutter Data Science](https://cookiecutter-data-science.drivendata.org/), and the
[Max Planck Bioinformatics template](https://github.com/maxplanck-ie/cookiecutter-bioinformatics-project):
`scripts/` and `config/` at the top level, raw data in `data/`, all pipeline outputs in `results/`,
publication figures in `reports/figures/`, reusable reference assets in `resources/`, and
manuscript writing work in `manuscript/`.

```
prion_rnaseq_meta_analysis/          ← project root (clone this repo)
│
├── README.md                         ← start here; run order + all instructions
├── environment.yml                   ← conda env (salmon-rnaseq); pin after first run
├── .gitignore                        ← excludes data/raw/, results/, resources/cache/
│
│  # ── SCRIPTS ─────────────────────────────────────────────────────────────
│
├── scripts/                          ← all executable scripts; numbered for run order
│   ├── 00_download_fastq.sh          ← download FASTQs from SRA (run on your system)
│   ├── 01_qc.sh                      ← FastQC + MultiQC (run on your system)
│   ├── 02_build_salmon_index.sh      ← build mouse + human Salmon indices
│   ├── 03_reads_quant.sh             ← Salmon quantification; FASTQs deleted after
│   ├── 04_build_coldata.R            ← fetch GEO metadata; write coldata CSVs
│   ├── 05_stream1_cns_meta_analysis.R← DESeq2 × 5 studies + meta-analysis (MCPS)
│   ├── 06_stream2_human_validation.R ← human sCJD validation (GSE214373)
│   ├── 07_stream3_peripheral_analysis.R← peripheral longitudinal (GSE210128)
│   ├── 08_integrate_streams.R        ← cross-stream integration + master table
│   ├── 09_cell_type_deconvolution.R   ← marker-gene deconvolution + MCPS-microglia correlation
│   └── notebooks/                    ← executed Colab notebooks (.ipynb) — the exact
│       ├── 04_build_coldata_v2.ipynb    analysis record behind the .R scripts
│       ├── 05_stream1_cns_meta_analysis_v2.ipynb   (kept for provenance; the .R
│       ├── 06_stream2_human_validation_v2.ipynb     scripts above are a 1:1
│       ├── 07_stream3_peripheral_analysis_v2.ipynb  conversion of these notebooks)
│       ├── 08_integrate_streams_v2.ipynb
│       └── 09_cell_type_deconvolution_v1.ipynb
│
│  # ── CONFIG ───────────────────────────────────────────────────────────────
│
├── config/                           ← curated metadata; read by scripts 00, 04–07
│   ├── sample_metadata.tsv           ← master SRR/GSE/LAYOUT/ORGANISM table (867 rows; input to scripts 00, 03, 04)
│   └── coldata/                      ← per-study coldata CSVs (written by script 04)
│       ├── GSE149805_coldata.csv     ← Stratum A: cell_type, time_point, condition
│       ├── GSE198063_coldata.csv     ← Stratum A: condition
│       ├── GSE144738_coldata.csv     ← Stratum A: time_point, condition
│       ├── GSE189527_coldata.csv     ← Stratum A: cell_type, time_point, condition
│       ├── GSE201249_coldata.csv     ← Stratum A: brain_region, cell_type, condition
│       ├── GSE210128_coldata.csv     ← Stratum B: tissue, time_point, cohort, condition
│       ├── GSE214373_coldata.csv     ← Stratum C: condition, sex, age
│       ├── GSE90977_coldata.csv      ← Stratum D: time_point, condition
│       ├── GSE124932_coldata.csv     ← Stratum D: condition (WT-only)
│       ├── GSE207251_coldata.csv     ← Stratum D: condition
│       └── GSE202275_coldata.csv     ← Stratum D: condition (vehicle-only)
│
│  # ── DATA ─────────────────────────────────────────────────────────────────
│  # FASTQs are TRANSIENT — downloaded by script 00, deleted by script 03
│  # after Salmon quantification. The folder exists only during processing.
│
├── data/
│   └── raw/
│       ├── fastq/                    ← transient; deleted after script 03 completes
│       │   ├── GSE149805/            ← one subfolder per GSE study
│       │   │   ├── SRR19907647_1.fastq.gz
│       │   │   ├── SRR19907647_2.fastq.gz
│       │   │   └── ...
│       │   ├── GSE198063/
│       │   ├── GSE144738/
│       │   ├── GSE189527/
│       │   ├── GSE201249/
│       │   ├── GSE210128/
│       │   ├── GSE214373/
│       │   ├── GSE90977/
│       │   ├── GSE124932/
│       │   ├── GSE207251/
│       │   └── GSE202275/
│       └── geo_soft/                 ← raw GEO SOFT files (cached by script 04)
│           ├── GSE149805_family.soft.gz
│           ├── GSE198063_family.soft.gz
│           └── .../ (11 studies total)
│
│  # ── RESOURCES ─────────────────────────────────────────────────────────────
│
├── resources/                        ← reusable reference assets; safe to delete + rebuild
│   ├── genome_indices/               ← Salmon indices (script 02); ~20 GB total
│   │   ├── mouse_salmon_index/       ← GRCm39, Ensembl r111 (all 10 mouse studies)
│   │   └── human_salmon_index/       ← GRCh38, Ensembl r111 (GSE214373 only)
│   ├── gene_sets/                    ← MSigDB gene sets (downloaded by script 05)
│   │   ├── msigdb_hallmark_mouse.rds
│   │   ├── msigdb_kegg_mouse.rds
│   │   └── msigdb_reactome_mouse.rds
│   └── cache/                        ← API caches; safe to delete (will re-download)
│       ├── tximeta/                  ← tximeta transcript-level reference metadata
│       │   └── BiocFileCache.sqlite
│       ├── biomart/                  ← BioMart ortholog + annotation caches
│       │   ├── GSE214373_orthologs_cache.rds
│       │   └── GSE210128_gene_annotation_cache.rds
│       └── geo/                      ← GEO metadata cache (script 04)
│
│  # ── RESULTS ───────────────────────────────────────────────────────────────
│  # All pipeline outputs. Sub-folders mirror the script that produced them.
│  # stream1 = script 05, stream2 = script 06, stream3 = script 07.
│
├── results/
│   ├── qc/                           ← script 01 outputs; organised by GSE
│   │   ├── fastqc/
│   │   │   ├── GSE149805/            ← FastQC HTML + zip for each sample in this study
│   │   │   ├── GSE198063/
│   │   │   ├── GSE144738/
│   │   │   ├── GSE189527/
│   │   │   ├── GSE201249/
│   │   │   ├── GSE210128/
│   │   │   ├── GSE214373/
│   │   │   ├── GSE90977/
│   │   │   ├── GSE124932/
│   │   │   ├── GSE207251/
│   │   │   └── GSE202275/
│   │   └── multiqc/
│   │       ├── GSE149805/            ← multiqc_report.html for this study
│   │       ├── GSE198063/
│   │       ├── GSE144738/
│   │       ├── GSE189527/
│   │       ├── GSE201249/
│   │       ├── GSE210128/
│   │       ├── GSE214373/
│   │       ├── GSE90977/
│   │       ├── GSE124932/
│   │       ├── GSE207251/
│   │       └── GSE202275/
│   ├── salmon_quant/                 ← script 03 outputs; GSE → SRR → quant.sf
│   │   ├── GSE149805/
│   │   │   ├── SRR19907647/
│   │   │   │   ├── quant.sf          ← main file read by tximeta
│   │   │   │   ├── aux_info/
│   │   │   │   └── logs/
│   │   │   │       └── salmon_quant.log  ← check mapping rate (flag < 70%)
│   │   │   └── .../
│   │   ├── GSE198063/
│   │   ├── GSE144738/
│   │   ├── GSE189527/
│   │   ├── GSE201249/
│   │   ├── GSE210128/
│   │   ├── GSE214373/
│   │   ├── GSE90977/
│   │   ├── GSE124932/
│   │   ├── GSE207251/
│   │   └── GSE202275/
│   │
│   │  # ── SCRIPT 05 OUTPUTS ─────────────────────────────────────────────
│   │
│   ├── deseq2/                       ← per-study DESeq2 results (script 05)
│   │   ├── GSE149805_deseq2_results.csv  ← all genes tested, with padj
│   │   ├── GSE149805_vst.rds             ← VST-normalised SummarizedExperiment
│   │   ├── GSE198063_deseq2_results.csv
│   │   ├── GSE198063_vst.rds
│   │   ├── GSE144738_deseq2_results.csv
│   │   ├── GSE144738_vst.rds
│   │   ├── GSE189527_deseq2_results.csv
│   │   ├── GSE189527_vst.rds
│   │   ├── GSE201249_deseq2_results.csv
│   │   └── GSE201249_vst.rds
│   ├── meta_analysis/                ← meta-analysis computation outputs (script 05)
│   │   ├── meta_analysis_full_results.csv  ← all genes: pooled LFC, SE, p, I²
│   │   ├── MCPS_genes.csv                  ← filtered MCPS (padj ≤ 0.05, |LFC| ≥ 0.5)
│   │   ├── MCPS_per_study_lfc_matrix.csv   ← per-study LFC for each MCPS gene
│   │   ├── stream1_summary.txt
│   │   └── session_info.txt
│   ├── stream1/                      ← canonical stream1 outputs read by scripts 06–08
│   │   ├── MCPS_genes.csv            ← *** KEY FILE — read by scripts 06, 07, 08 ***
│   │   ├── MCPS.csv                  ← alias with column names expected by script 06
│   │   ├── MCPS_per_study_lfc_matrix.csv
│   │   ├── meta_analysis_full_results.csv
│   │   ├── replication/              ← held-out Stratum D replication results
│   │   │   ├── GSE90977_replication.csv
│   │   │   ├── GSE124932_replication.csv
│   │   │   ├── GSE207251_replication.csv
│   │   │   ├── GSE202275_replication.csv
│   │   │   ├── replication_summary.csv
│   │   │   └── LOO_sensitivity_summary.csv
│   │   └── enrichment/               ← ORA + GSEA enrichment results
│   │       ├── GSE149805_ORA_results.csv
│   │       ├── GSE198063_ORA_results.csv
│   │       ├── GSE144738_ORA_results.csv
│   │       ├── GSE189527_ORA_results.csv
│   │       ├── GSE201249_ORA_results.csv
│   │       ├── consensus_pathways_3of5.csv
│   │       └── meta_GSEA_results.csv
│   │
│   │  # ── SCRIPT 06 OUTPUTS ─────────────────────────────────────────────
│   │
│   ├── stream2/                      ← human sCJD validation outputs (script 06)
│   │   ├── GSE214373_deseq2_primary.csv
│   │   ├── GSE214373_deseq2_sensitivity.csv
│   │   ├── GSE214373_sensitivity_comparison.csv
│   │   ├── GSE214373_orthologs.csv
│   │   ├── GSE214373_MCPS_concordance.csv
│   │   ├── GSE214373_GSEA_MCPS.csv   ← *** read by scripts 07 and 08 ***
│   │   ├── GSE214373_power_analysis.csv
│   │   ├── power_curves_multi_dispersion.csv
│   │   ├── min_detectable_lfc_by_dispersion.csv
│   │   └── stream2_summary.txt
│   │
│   │  # ── SCRIPT 07 OUTPUTS ─────────────────────────────────────────────
│   │
│   ├── stream3/                      ← peripheral longitudinal outputs (script 07)
│   │   ├── GSE210128_whole_blood_LRT_results.csv
│   │   ├── GSE210128_spleen_LRT_results.csv
│   │   ├── GSE210128_skeletal_muscle_LRT_results.csv
│   │   ├── GSE210128_whole_blood_terminal_results.csv
│   │   ├── GSE210128_spleen_terminal_results.csv
│   │   ├── GSE210128_skeletal_muscle_terminal_results.csv
│   │   ├── GSE210128_whole_blood_MCPS_concordance.csv
│   │   ├── GSE210128_spleen_MCPS_concordance.csv
│   │   ├── GSE210128_skeletal_muscle_MCPS_concordance.csv
│   │   ├── GSE210128_MCPS_concordance_summary.csv
│   │   ├── GSE210128_MCPS_first_detection.csv
│   │   ├── GSE210128_MCPS_first_detection_genes.csv
│   │   ├── GSE210128_per_timepoint_wald_all.csv
│   │   ├── GSE210128_cross_tissue_deg_overlap.csv
│   │   ├── GSE210128_triple_validated_candidates.csv  ← *** read by script 08 ***
│   │   ├── GSE210128_stream3_summary.txt
│   │   └── GSE210128_session_info.txt
│   │
│   │  # ── SCRIPT 08 OUTPUTS ─────────────────────────────────────────────
│   │
│   └── integration/                  ← cross-stream integration outputs (script 08)
│       ├── master_integration_table.csv  ← *** single source of truth for manuscript ***
│       ├── triple_validated_candidates.csv
│       ├── cross_stream_lfc_correlations.csv
│       ├── replication_cascade.csv
│       ├── pathway_convergence.csv
│       ├── integration_summary.txt
│       └── integration_session_info.txt
│
│  # ── REPORTS ───────────────────────────────────────────────────────────────
│  # Publication figures only — separated from CSV results so you can zip
│  # reports/figures/ and send directly to collaborators or journal.
│
├── reports/
│   └── figures/
│       ├── stream1/                  ← figures from script 05 (PDF + PNG pairs)
│       │   ├── Fig1A_meta_volcano.pdf / .png
│       │   ├── Fig1B_forest_plots_top20.pdf / .png
│       │   ├── Fig1C_I2_distribution.pdf / .png
│       │   ├── Fig2_per_study_DEG_counts.pdf / .png
│       │   ├── Fig3_MCPS_heatmap.pdf / .png
│       │   ├── Fig4_consensus_pathways.pdf / .png
│       │   ├── Fig5_GSEA_Hallmark.pdf / .png
│       │   ├── Fig6_replication_barplot.pdf / .png
│       │   └── Fig7_LOO_sensitivity.pdf / .png
│       ├── stream2/                  ← figures from script 06
│       │   ├── fig1a_pca.pdf / .png
│       │   ├── fig1b_sample_distance_heatmap.pdf / .png
│       │   ├── fig2_age_confound.pdf / .png
│       │   ├── fig3_sensitivity_overlap.pdf / .png
│       │   ├── fig4_volcano_primary.pdf / .png
│       │   ├── fig5_lfc_scatter.pdf / .png
│       │   ├── fig6_gsea_enrichment.pdf / .png
│       │   └── fig7_power_curve.pdf / .png
│       ├── stream3/                  ← figures from script 07
│       │   ├── fig1_sample_overview.pdf / .png
│       │   ├── fig2_qc_libsize_pca.pdf / .png
│       │   ├── fig3_*_trajectory_heatmap.pdf / .png  (3 tissues)
│       │   ├── fig4_mcps_trajectory_lines.pdf / .png
│       │   ├── fig5_terminal_volcanos.pdf / .png
│       │   ├── fig6_lfc_scatter_peripheral_vs_mcps.pdf / .png
│       │   ├── fig7_mcps_concordance_bar.pdf / .png
│       │   ├── fig8_mcps_first_detection.pdf / .png
│       │   ├── fig9_per_timepoint_deg_count.pdf / .png
│       │   ├── fig10_*_pca_longitudinal.pdf / .png  (3 tissues)
│       │   └── fig11_mcps_cross_tissue_heatmap.pdf / .png
│       ├── integration/              ← figures from script 08
│       │   ├── fig1_replication_cascade.pdf / .png
│       │   ├── fig2_lfc_scatter_matrix.pdf / .png
│       │   ├── fig3_triple_validated_dotplot.pdf / .png  ← absent (0 triple-validated candidates)
│       │   ├── fig4_detection_timeline.pdf / .png
│       │   ├── fig5_pathway_convergence.pdf / .png
│       │   └── fig6_human_pathway_validation_bubble.pdf / .png
│       └── deconvolution/            ← figures from script 09
│           ├── fig09a_celltype_violin.pdf / .png
│           └── fig09b_MCPS_microglia_corr.pdf / .png
│
│  # ── MANUSCRIPT ────────────────────────────────────────────────────────────
│  # All writing work lives here. Tracked in git (text files are small).
│  # Keep .docx / .tex drafts here; final submission goes to journal portal.
│
├── manuscript/
│   ├── literature_review/            ← reading notes, annotated paper summaries
│   │   ├── prion_transcriptomics_notes.md   ← key findings from reviewed papers
│   │   ├── methods_comparison.md            ← meta-analysis method comparisons
│   │   └── papers/                          ← PDFs of key references (optional)
│   ├── evidence_synthesis/           ← cross-study tables, evidence grading
│   │   ├── evidence_table.csv               ← gene × study evidence matrix
│   │   ├── pathway_evidence_summary.md      ← narrative pathway synthesis
│   │   └── MCPS_biological_annotation.md    ← biological interpretation of MCPS genes
│   └── draft/                        ← manuscript text files
│       ├── main_text.md              ← working draft (Markdown → Word/LaTeX)
│       ├── supplementary.md          ← supplementary methods + tables
│       └── cover_letter.md           ← journal cover letter
│
│  # ── LOGS ──────────────────────────────────────────────────────────────────
│
├── logs/                             ← shell script stdout/stderr (nohup output)
│   ├── 00_download.log
│   ├── 01_qc.log
│   ├── 02_index.log
│   └── 03_quant.log
│
│  # ── DOCS ──────────────────────────────────────────────────────────────────
│
└── docs/                             ← pipeline documentation
    ├── docs_sample_exclusion_rationale.md   ← per-study exclusion justification + Methods paragraph
    ├── pipeline_documentation.html          ← rendered end-to-end pipeline documentation
    ├── report_integration.html              ← cross-stream integration methods report
    └── report_prion_metaanalysis_study_design_v2.html  ← study design report
```

### What is excluded from git (`.gitignore`)

| Path | Size | Reason |
|------|------|--------|
| `data/raw/fastq/` | transient | FASTQs are deleted after Salmon quantification; re-downloadable from SRA via script 00 |
| `data/raw/geo_soft/` | — | Re-downloadable via script 04 |
| `results/salmon_quant/` | ~7 GB | Regenerable intermediate (script 03); raw data are public on SRA/GEO. Exceeds GitHub limits — see [Large data files](#large-data-files) |
| `results/qc/` | ~2 GB | Regenerable intermediate (script 01); per-sample FastQC HTML/zip. Exceeds practical repo size — see [Large data files](#large-data-files) |
| `resources/genome_indices/` | ~20 GB | Rebuilt by script 02 |
| `resources/cache/` | — | Transient; re-downloaded automatically |
| `manuscript/papers/` | — | PDFs are large; cite by DOI instead |

All *final* result tables (`results/deseq2/`, `results/stream1/`, `results/stream2/`,
`results/stream3/`, `results/integration/`, `results/meta_analysis/`,
`results/deconvolution/`) and all publication figures (`reports/figures/`) **are**
included in the repository, so every manuscript figure and table can be inspected
and regenerated without re-running the compute-heavy steps.

### Large data files

GitHub enforces a 100 MB per-file hard limit (and warns above 50 MB), and
recommends repositories stay well under ~1 GB. The two regenerable intermediates
(`results/salmon_quant/` ≈ 7 GB, `results/qc/` ≈ 2 GB) therefore cannot — and
need not — be hosted here:

- They are fully regenerated by scripts 00–03 from the public SRA/GEO accessions
  listed in `config/sample_metadata.tsv`.
- For long-term archival of these intermediates (and a citable DOI for the
  Data Availability statement), we recommend uploading them to
  [Zenodo](https://zenodo.org) (free, 50 GB per record) and linking the DOI here:
  `Zenodo DOI: [TO BE INSERTED AFTER ARCHIVAL]`.
- One result file, `results/stream3/GSE210128_per_timepoint_wald_all.csv`
  (52 MB), is stored **gzipped** (`.csv.gz`, 15 MB) to stay under GitHub's 50 MB
  warning threshold. Decompress with `gunzip` before use, or read directly in R
  with `readr::read_csv()` / `data.table::fread()` (both accept `.gz`).

---

## Software Requirements

### System tools (must be in PATH)
| Tool | Version tested | Purpose |
|------|---------------|---------|
| SRA Toolkit (`fasterq-dump`) | ≥ 3.0 | FASTQ download |
| pigz | any | Parallel gzip compression |
| FastQC | ≥ 0.12 | Per-sample QC |
| MultiQC | ≥ 1.14 | QC aggregation |
| Salmon | ≥ 1.10 | RNA-seq quantification |
| R | ≥ 4.3 | All downstream analysis |

### R packages (auto-installed by each script)
**Bioconductor:**
`tximeta`, `tximport`, `DESeq2`, `apeglm`, `BiocParallel`, `SummarizedExperiment`,
`clusterProfiler`, `fgsea`, `org.Mm.eg.db`, `AnnotationDbi`, `enrichplot`,
`msigdbr`, `biomaRt`, `ComplexHeatmap`, `circlize`

**CRAN:**
`metafor`, `ggplot2`, `ggrepel`, `dplyr`, `tidyr`, `readr`, `purrr`, `stringr`,
`scales`, `RColorBrewer`, `Matrix`, `tibble`, `patchwork`, `cowplot`

### Internet access required at runtime
- `04_build_coldata.R` — NCBI GEO FTP + Entrez API
- `06_stream2_human_validation.R` — Ensembl BioMart (first run only; cached thereafter)
- `05_stream1_cns_meta_analysis.R` — KEGG enrichment via `enrichKEGG`
- `02_build_salmon_index.sh` — Ensembl FTP (reference downloads)

---

## Script Formats (.R and .ipynb)

Each analysis step (04–09) ships in two equivalent forms:

- **`scripts/*.R`** — the primary, review-ready scripts. Run these.
- **`scripts/notebooks/*.ipynb`** — the executed Google Colab notebooks that
  produced the reported results, kept as the provenance record (cell outputs
  included). The .R scripts are a faithful 1:1 conversion of these notebooks;
  only Colab environment-setup cells (Drive mount, `rpy2` magic, package
  installers) were removed. Package dependencies are managed centrally via
  `environment.yml` + `renv.lock` instead.

One repair was necessary during conversion: the executed notebook for script 06
had been edited after its final run and was missing the guard line
`if ("gene_symbol" %in% colnames(res_primary_df)) {` opening the expanded
ortholog-recovery block (its `else` branch names the condition). The line was
restored in the .R script so it parses and runs; executed results are
unaffected. The original notebook is preserved unmodified in `scripts/notebooks/`.

---

## Run Order

Scripts must be run in the order shown. Scripts 06, 07 can run in parallel
after 05 completes. Script 08 requires both 06 and 07 to be complete.

```
Step 1:  bash scripts/00_download_fastq.sh
Step 2:  bash scripts/01_qc.sh
         → Review multiqc_report.html before proceeding
Step 3:  bash scripts/02_build_salmon_index.sh
Step 4:  bash scripts/03_reads_quant.sh
         → Check mapping rates: grep "Mapping rate" results/salmon_quant/*/*/logs/salmon_quant.log
         → Flag any sample with mapping rate < 70%
Step 5:  Rscript scripts/04_build_coldata.R
Step 6:  Rscript scripts/05_stream1_cns_meta_analysis.R
         → Produces results/stream1/MCPS_genes.csv (required by 06, 07, 08)
Step 7a: Rscript scripts/06_stream2_human_validation.R   # can run in parallel with 7b
Step 7b: Rscript scripts/07_stream3_peripheral_analysis.R
Step 8:  Rscript scripts/08_integrate_streams.R
Step 9:  Rscript scripts/09_cell_type_deconvolution.R
         → Requires VST .rds files from Step 6 (results/deseq2/*_vst.rds)
```

**Estimated runtimes** (8-core workstation, 32 GB RAM):

| Script | Estimated time |
|--------|---------------|
| 00 | 4–8 hours (network-bound; 867 SRR accessions) |
| 01 | 1–2 hours |
| 02 | 30–60 min |
| 03 | 3–6 hours |
| 04 | 15–20 min (NCBI API rate-limited) |
| 05 | 2–4 hours (5 DESeq2 runs + meta-analysis) |
| 06 | 20–40 min |
| 07 | 1–3 hours (3 tissues × LRT + Wald + per-timepoint) |
| 08 | 5–15 min |
| 09 | 10–20 min |

---

## Script Reference

### 00 — Download FASTQs
**File:** `00_download_fastq.sh`

Downloads all 867 SRR accessions from SRA using `fasterq-dump`.

**Input:**
- `config/sample_metadata.tsv` — one SRR accession per line

**Output:**
- `fastq/<SRR>/<SRR>_1.fastq.gz` (paired-end R1)
- `fastq/<SRR>/<SRR>_2.fastq.gz` (paired-end R2)
- `fastq/<SRR>/<SRR>.fastq.gz` (single-end)

**Notes:**
- Resume-safe: skips accessions already downloaded
- Temporary files written to `/tmp/fasterq_tmp`
- Uses `pigz` for parallel compression (falls back to `gzip`)

---

### 01 — Quality Control
**File:** `01_qc.sh`

Runs FastQC on all FASTQs, then aggregates with MultiQC.

**Input:** `fastq/<SRR>/*.fastq.gz`

**Output:**
- `results/qc/fastqc/` — per-sample FastQC HTML + zip
- `results/qc/multiqc_report.html`

**QC thresholds to check before proceeding:**
- Per-base quality: median Phred ≥ 28
- Adapter content: < 5%
- Duplicate rate: < 50%
- Overrepresented sequences: check for rRNA contamination

---

### 02 — Build Salmon Indices
**File:** `02_build_salmon_index.sh`

Builds decoy-aware Salmon indices for mouse (GRCm39) and human (GRCh38).

**Reference:** Ensembl Release 111

**Output:**
- `resources/genome_indices/mouse_salmon_index/` — used for all 10 mouse datasets
- `resources/genome_indices/human_salmon_index/` — used for GSE214373 only

**Notes:**
- Uses gentrome (cDNA + ncRNA + genome decoy) for selective alignment
- `--keepDuplicates` retains multi-copy gene family transcripts
- Requires ~20 GB disk, ~16 GB RAM during build

---

### 03 — Salmon Quantification
**File:** `03_reads_quant.sh`

Quantifies all samples with Salmon selective alignment.

**Input:**
- `fastq/<SRR>/` — per-sample FASTQ directories
- `resources/genome_indices/mouse_salmon_index/` or `resources/genome_indices/human_salmon_index/`
- `config/sample_metadata.tsv` — SRR / GSE / LAYOUT / ORGANISM table

**Output:** `results/salmon_quant/<GSE>/<SRR>/quant.sf`

**Notes:**
- `--gcBias` applied to all samples
- `--seqBias` applied to paired-end only (requires fragment length distribution)
- Resume-safe: skips already-quantified samples
- Check mapping rates after completion (flag < 70%)

---

### 04 — Build Coldata
**File:** `04_build_coldata.R`

Fetches GEO metadata for all 11 studies and writes per-study coldata CSVs.

**Input:**
- `config/sample_metadata.tsv` — SRR/GSE/LAYOUT/ORGANISM table (867 rows)
- Internet access (NCBI GEO FTP + Entrez API)

**Output:** `config/coldata/<GSE>_coldata.csv` for all 11 studies

**Coldata CSV format:**

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | character | SRR accession |
| `files` | character | Path to `quant.sf` |
| `condition` | factor | `"prion"` or `"control"` (reference = control) |
| `time_point` | factor | wpi (study-specific; where applicable) |
| `cell_type` | factor | Cell type (GSE149805, GSE189527, GSE201249) |
| `brain_region` | factor | Brain region (GSE201249) |
| `sex` | factor | Sex (GSE214373) |
| `tissue` | factor | Tissue (GSE210128) |
| `cohort` | factor | Main/validation cohort (GSE210128) |

**Study-specific decisions documented in script header:**
- GSE198063: Input fraction only (RiboTag IP excluded); FFI + CJD → `prion`
- GSE124932: WT arms only (KO arms excluded)
- GSE202275: Vehicle arms only (VU846-treated arms excluded)

---

### 05 — Stream 1: CNS Meta-Analysis
**File:** `05_stream1_cns_meta_analysis.R`

**Sections:**
| Section | Content |
|---------|---------|
| 0 | Package installation, global configuration |
| 1 | Directory paths, significance thresholds, study registries |
| 2 | Helper functions: `load_coldata()`, `import_with_tximeta()`, `filter_low_counts()`, `detect_pca_outliers()`, `run_deseq2()`, `annotate_mouse_genes()` |
| 3 | Per-study DESeq2 (5 Stratum A studies, sequential to manage RAM) |
| 4 | DerSimonian-Laird random-effects meta-analysis; MCPS definition |
| 5 | Stratum D held-out replication; leave-one-out sensitivity |
| 6 | Pathway enrichment: ORA (per study) + GSEA (meta-analysis ranking) |
| 7 | Publication-ready figures (forest plots, I² histogram, MCPS heatmap, etc.) |
| 8 | Output tables and session info |

**Key outputs:**
```
results/stream1/MCPS_genes.csv              ← MCPS gene list (required by 06, 07, 08)
results/meta_analysis/meta_analysis_full_results.csv
results/stream1/replication/<GSE>_replication.csv   ← one per Stratum D study
results/stream1/replication/LOO_sensitivity_summary.csv
results/stream1/enrichment/meta_GSEA_*.csv
```

**DESeq2 designs per study:**

| Study | Design | Notes |
|-------|--------|-------|
| GSE149805 | `~ time_point + condition` | Applied per cell_type stratum |
| GSE198063 | `~ condition` | Input fraction only |
| GSE144738 | `~ time_point + condition` | |
| GSE189527 | `~ time_point + condition` | |
| GSE201249 | `~ brain_region + cell_type + condition` | Rank-deficiency guard included |
| GSE90977 | `~ time_point + condition` | Stratum D |
| GSE124932 | `~ condition` | Stratum D, WT only |
| GSE207251 | `~ condition` | Stratum D |
| GSE202275 | `~ condition` | Stratum D, vehicle only |

**MCPS definition:**
```
meta-padj ≤ 0.05  AND  |meta-log2FC| ≥ 0.5  AND  I² ≤ 75%
(all thresholds inclusive; Benjamini-Hochberg FDR)
```

---

### 06 — Stream 2: Human Validation
**File:** `06_stream2_human_validation.R`

**Prerequisite:** `results/stream1/MCPS_genes.csv` from Script 05

**Sections:**
| Section | Content |
|---------|---------|
| 0 | Libraries, configuration, `theme_pub`, `save_fig()` |
| 1 | Load and validate coldata; age confound documentation |
| 2 | Load Salmon quantifications via tximeta |
| 3 | Sample-level QC: PCA + sample distance heatmap |
| 4 | Primary DESeq2: `~ age + sex + condition`; apeglm shrinkage |
| 5 | Sensitivity analysis: `~ sex + condition` (age omitted) |
| 6 | Human → Mouse ortholog mapping (BioMart, Ensembl r111, cached) |
| 7 | MCPS concordance: direction-of-effect test + binomial test |
| 8 | GSEA: MCPS as gene set on ranked human DEGs |
| 9 | Permutation-based power analysis (multi-dispersion) |
| 10 | Publication-ready figures (7 figures) |
| 11 | Output tables and summary report |

**Key outputs:**
```
results/stream2/GSE214373_deseq2_primary.csv
results/stream2/GSE214373_MCPS_concordance.csv
results/stream2/GSE214373_GSEA_MCPS.csv          ← read by 07 (Section 9.7) and 08
results/stream2/GSE214373_orthologs.csv
results/stream2/GSE214373_power_analysis.csv
results/stream2/stream2_summary.txt
```

**Age confound (documented):**
- sCJD mean age 71.6 ± 4.2 y vs control 83.3 ± 9.6 y (Welch t p = 0.024)
- Age-condition Pearson r = −0.60
- VIF(condition) = 2.06 in primary model — full rank, age included as covariate
- Sensitivity analysis quantifies age-driven false positives

**Ortholog mapping:**
- BioMart host: `jan2024.archive.ensembl.org` (Ensembl r111, matching Salmon index)
- Results cached at `resources/cache/GSE214373_orthologs_cache.rds`
- One-to-many mappings resolved by taking the ortholog with highest `pct_id_human2mouse`

---

### 07 — Stream 3: Peripheral Analysis
**File:** `07_stream3_peripheral_analysis.R`

**Prerequisite:** `results/stream1/MCPS_genes.csv` from Script 05

**Tissues:** `whole_blood`, `spleen`, `skeletal_muscle`

**Timepoints:** 4, 8, 12, 14, 16, 18, 20 wpi + terminal (encoded as 22)

**Sections:**
| Section | Content |
|---------|---------|
| 0 | Libraries, configuration, helper functions |
| 1 | Load and validate coldata; CONSORT sample flow |
| 2 | Load Salmon quantifications; summarizeToGene |
| 3 | Sample-level QC: library size, VST PCA per tissue |
| 4 | Load and validate MCPS gene list |
| 5 | Arm A: Longitudinal LRT per tissue (continuous + factor sensitivity) |
| 6 | Arm B: Terminal-timepoint Wald test per tissue (apeglm shrinkage) |
| 7 | Gene symbol annotation (biomaRt, cached) |
| 8 | Arm C: MCPS concordance at terminal timepoint |
| 9 | MCPS first-detection: per-timepoint Wald tests |
| 9.7 | Cross-stream: Stream 2 leading-edge cross-reference (triple-validated candidates) |
| 10 | Publication-ready figures (11 figures including Fig 11 cross-tissue heatmap) |
| 11 | Output tables, summary report, session info |

**Key outputs:**
```
results/stream3/GSE210128_<tissue>_LRT_results.csv          ← 3 files
results/stream3/GSE210128_<tissue>_terminal_results.csv     ← 3 files
results/stream3/GSE210128_<tissue>_MCPS_concordance.csv     ← 3 files
results/stream3/GSE210128_MCPS_concordance_summary.csv
results/stream3/GSE210128_MCPS_first_detection.csv
results/stream3/GSE210128_MCPS_first_detection_genes.csv
results/stream3/GSE210128_per_timepoint_wald_all.csv
results/stream3/GSE210128_triple_validated_candidates.csv   ← read by 08
results/stream3/GSE210128_stream3_summary.txt
```

**Statistical design notes:**
1. LRT tests ONLY the `condition × time_point_wpi` interaction (symmetric formulation)
2. Terminal timepoint excluded from continuous LRT; used only in Arm B
3. Size factors estimated from full tissue dataset, transferred to all subsets
4. Dispersion priors transferred from full-tissue dds to per-timepoint subsets (Arm C)
5. `cohort` included as blocking factor where both cohorts present (not in blood)
6. MCPS concordance denominator = all MCPS genes (conservative)

---

### 08 — Cross-Stream Integration
**File:** `08_integrate_streams.R`

**Prerequisites:** Outputs from Scripts 05, 06, and 07 (graceful fallback if any are missing)

**Sections:**
| Section | Content |
|---------|---------|
| 0 | Libraries, `BASE_DIR` resolution, `save_fig()`, `save_heatmap()`, `safe_read()` |
| 1 | Load 15 CSV files across all three streams; availability flags |
| 2 | Validate join keys; detect gene symbol case mismatches |
| 3 | Master integration table (one row per MCPS gene; evidence tiers 1–5) |
| 4 | Triple-validated candidates (Tier 1 + Tier 2) |
| 5 | Cross-stream LFC correlations (7 Pearson r pairs) |
| 6 | Replication cascade (6 stages) |
| 7 | Pathway convergence analysis |
| 8 | Publication-ready figures (5 figures) |
| 9 | Output tables (5 CSVs) |
| 10 | Plain-text summary report + session info |

**Key outputs:**
```
results/integration/master_integration_table.csv
results/integration/triple_validated_candidates.csv
results/integration/cross_stream_lfc_correlations.csv
results/integration/replication_cascade.csv
results/integration/pathway_convergence.csv
results/integration/fig1_replication_cascade.pdf/png
results/integration/fig2_lfc_scatter_matrix.pdf/png
results/integration/fig3_triple_validated_dotplot.pdf/png
results/integration/fig4_detection_timeline.pdf/png
results/integration/fig5_pathway_convergence.pdf/png
results/integration/integration_summary.txt
```

**Graceful fallback:** All 15 input files are loaded with `safe_read()`. If any file
is missing (e.g., Stream 2 not yet run), the corresponding columns are set to NA
and downstream analyses skip gracefully. The script always produces a master table
and summary report, even with partial inputs.

---

### 09 — Cell-Type Deconvolution
**File:** `09_cell_type_deconvolution_v1.ipynb`

**Prerequisite:** VST-normalised count matrices from Script 05 (Stratum A studies)

**Sections:**
| Section | Content |
|---------|---------|
| 0 | Libraries, configuration, `theme_pub`, `save_fig()` |
| 1 | Load VST matrices for all 5 Stratum A studies; merge into N=171 sample matrix |
| 2 | Define marker gene sets (8 cell types); exclude MCPS genes from all marker sets |
| 3 | Compute marker-gene z-scores per sample; Wilcoxon tests (prion vs. control) |
| 4 | MCPS-microglia correlation (Spearman r per MCPS gene vs. microglia score) |
| 5 | Publication-ready figures (violin plot + correlation bar chart) |
| 6 | Output tables and summary report |

**Key outputs:**
```
results/deconvolution/celltype_score_stats.csv
results/deconvolution/MCPS_microglia_correlation.csv
reports/figures/deconvolution/fig09a_celltype_violin.pdf/png
reports/figures/deconvolution/fig09b_MCPS_microglia_corr.pdf/png
```

**Key findings:**
- Microglia: Δ = +0.842, p = 2.45×10⁻⁷ (***); A1 Astrocytes: Δ = +0.787 (*); Neurons: Δ = −0.468 (***)
- 5/6 MCPS genes correlated with microglia score (Spearman r = 0.47–0.64, all p < 0.001)
- Slc15a4/Tasl: r = 0.061, NS — non-microglial mechanism confirmed

---

## Input/Output Dependency Map

```
config/sample_metadata.tsv
    └── 00_download_fastq.sh
            └── data/raw/fastq/<GSE>/<SRR>/
                    ├── 01_qc.sh → results/qc/<fastqc|multiqc>/<GSE>/
                    └── 03_reads_quant.sh ←── resources/genome_indices/mouse_salmon_index/  (from 02)
                                              resources/genome_indices/human_salmon_index/   (from 02)
                                              config/sample_metadata.tsv
                            └── results/salmon_quant/<GSE>/<SRR>/quant.sf
                                    └── 04_build_coldata.R ←── config/sample_metadata.tsv
                                                                NCBI GEO API
                                            └── config/coldata/<GSE>_coldata.csv
                                                    └── 05_stream1_cns_meta_analysis.R
                                                            └── results/stream1/MCPS_genes.csv
                                                                    ├── 06_stream2_human_validation.R
                                                                    │       └── results/stream2/GSE214373_GSEA_MCPS.csv
                                                                    │           results/stream2/GSE214373_MCPS_concordance.csv
                                                                    │           results/stream2/GSE214373_deseq2_primary.csv
                                                                    │
                                                                    ├── 07_stream3_peripheral_analysis.R
                                                                    │       └── results/stream3/GSE210128_*.csv
                                                                    │           results/stream3/GSE210128_triple_validated_candidates.csv
                                                                    │
                                                                    ├── 08_integrate_streams.R  ←── outputs from 05, 06, 07
                                                                    │       └── results/integration/
                                                                    └── 09_cell_type_deconvolution.R ←── results/deseq2/*_vst.rds (from 05)
                                                                            └── results/deconvolution/
```

---

## Key Parameters

All significance thresholds are defined as named constants at the top of each
script and applied consistently throughout. To change a threshold, edit the
constant in the relevant script — do not change individual filter calls.

| Parameter | Value | Scripts | Description |
|-----------|-------|---------|-------------|
| `PADJ_THRESH` | 0.05 | 05, 06, 07, 08 | BH-FDR threshold (inclusive) |
| `LFC_THRESH` | 0.5 | 05, 06, 07, 08 | \|log2FC\| threshold (inclusive) |
| `I2_THRESH` | 75 | 05 | I² exclusion threshold for MCPS (exclusive: I² > 75% excluded) |
| `MIN_STUDIES` | 3 | 05 | Minimum studies for meta-analysis inclusion |
| `GSEA_NPERM` | 100,000 | 06 | GSEA permutations (publication-quality; fgsea eps=0) |
| `POWER_ALPHA` | 0.05 | 06 | Type I error rate for power analysis |
| `POWER_TARGET` | 0.80 | 06 | Target power for power analysis |
| `TERMINAL_WPI` | 22L | 07 | Sentinel value encoding terminal timepoint |

---

## Evidence Tier Definitions

Applied in `08_integrate_streams.R` to classify MCPS genes by cross-stream support:

| Tier | Criteria | Interpretation | Genes (this run) |
|------|----------|---------------|-----------------|
| **Tier 1** | MCPS + Stream 2 GSEA leading edge + concordant peripheral detection | Strongest: CNS signal replicated in human sCJD AND peripheral tissue (same direction) | 0 |
| **Tier 2** | MCPS + Stream 2 GSEA leading edge + peripheral detection (direction may differ) | Strong: CNS + human, peripheral detected but direction uncertain | 0 |
| **Tier 3** | MCPS + Stream 2 GSEA leading edge only | Moderate: CNS + human; no peripheral signal | 0 |
| **Tier 4** | MCPS + peripheral detection only | Moderate: CNS + peripheral; not in human leading edge | 2 (Igsf6, Fcgr4) |
| **Tier 5** | MCPS only | Baseline: CNS meta-analysis only | 4 (Tasl, Parp14, Iigp1, Il1a) |

**Triple-validated candidates** = Tier 1 + Tier 2 genes. **This run: 0 triple-validated candidates** (Stream 2 GSEA returned 0 results due to insufficient gene set size after ortholog filtering; human dataset underpowered at N=14).

**Note on tier labels:** The master integration table uses automated tier assignment (Tier 4/5 above). The literature convergence table (Section 8 of nb08) uses manually curated tiers (Tier 3/4/4+) based on published evidence. Both are valid but serve different purposes.

---

## Reproducibility Notes

### R Package Reproducibility (renv)

This pipeline uses [renv](https://rstudio.github.io/renv/) to lock exact R package
versions. The lockfile (`renv.lock`) pins all 35 R packages to their
Bioconductor 3.18 / CRAN versions, ensuring identical results across machines and
over time.

**One-time setup (run before any analysis script):**

```bash
# 1. Activate the conda environment
conda env create -f environment.yml
conda activate salmon-rnaseq

# 2. Restore the exact R package library from renv.lock
Rscript renv_setup.R
```

**What renv does:**
- Installs all R packages into a project-local library (`renv/library/`)
- Does NOT modify your system R library
- Takes 20–60 minutes on first run (downloads ~35 packages + dependencies)
- Subsequent runs are instant (packages already installed)

**Kaggle Notebook users:**
Paste the following into a code cell before running any analysis script:
```r
source("/kaggle/working/renv_setup.R")
```

**After setup, activate renv at the start of each R session:**
```r
renv::activate()
```
Or open the project in RStudio — renv activates automatically via `.Rprofile`.

**Key pinned versions (Bioconductor 3.18, R ≥ 4.3):**

| Package | Version | Role |
|---------|---------|------|
| DESeq2 | 1.42.1 | Differential expression |
| apeglm | 1.24.0 | LFC shrinkage |
| tximeta | 1.20.3 | Salmon import |
| metafor | 4.6-0 | Random-effects meta-analysis |
| fgsea | 1.28.0 | GSEA (100,000 permutations) |
| biomaRt | 2.58.2 | Ortholog mapping |
| clusterProfiler | 4.10.1 | ORA enrichment |
| ComplexHeatmap | 2.18.0 | Publication heatmaps |
| ggplot2 | 3.5.1 | All figures |

Full lockfile: `renv.lock`

---

### Sample Exclusion Documentation

Of the 867 samples quantified by Salmon, **125 were excluded** before analysis
on the basis of biological incompatibility, and **3 further samples failed
sequencing QC** (1 in GSE149805, 2 in GSE144738), leaving **739 analyzed
samples** (516 Stratum A + 173 Stratum B + 14 Stratum C + 36 Stratum D).

| Study | N excluded | Reason |
|-------|-----------|--------|
| GSE149805 | 10 | *noCre* negative controls (no cell-type specificity) |
| GSE198063 | 101 | RiboTag IP fraction (ribosome-associated mRNA; incompatible with total RNA meta-analysis) |
| GSE124932 | 8 | SARM1−/− knockout arms (altered axon degeneration pathways) |
| GSE202275 | 6 | VU846-treated arms (ISR pharmacological intervention) |

Full scientific justification, per-study narratives, and a ready-to-paste
manuscript Methods paragraph are in:
**`docs/docs_sample_exclusion_rationale.md`**

---

### General Reproducibility Notes

1. **R session info** is written to `*_session_info.txt` by Scripts 05, 07, and 08.
   Script 06 writes session info within `stream2_summary.txt`.

2. **Random seed** is set with `set.seed(42)` in Scripts 06 and 07 before any
   stochastic operations (GSEA permutations, power simulations).

3. **BioMart ortholog cache** (`resources/cache/GSE214373_orthologs_cache.rds`)
   is written on first run and reused thereafter. Delete this file to force a
   fresh BioMart query.

4. **tximeta transcript-level cache** (`resources/cache/`) stores Ensembl
   reference metadata. Delete to force re-download.

5. **Ensembl release:** All scripts use Ensembl Release 111 (GRCm39 / GRCh38.p14).
   The Salmon index (Script 02), tximeta import (Scripts 05–07), and BioMart
   ortholog mapping (Script 06) all reference r111 for consistency.

6. **NCBI API key:** Set the `NCBI_API_KEY` environment variable before running
   Script 04 to increase the Entrez rate limit from 3 to 10 requests/second:
   ```bash
   export NCBI_API_KEY="your_key_here"
   Rscript scripts/04_build_coldata.R
   ```

7. **Memory management:** Script 05 processes studies sequentially and calls
   `rm() + gc()` after each study to keep peak RAM below ~6 GB. Do not modify
   this to parallel processing without increasing available RAM to ≥ 32 GB.

8. **Salmon locale fix:** Scripts 02 and 03 set `LC_ALL=en_US.UTF-8` to prevent
   Salmon ≥ 1.10 index write failures in minimal Docker/conda environments.

---

## Known Limitations

1. **GSE214373 sample size (N=14):** Power to detect small LFC is limited.
   The power analysis in Script 06 (Section 9) quantifies this explicitly.
   **Minimum detectable LFC at 80% power: ≥ 1.50 at median dispersion.** All
   observed MCPS gene LFCs in human are < 0.20 — well below the detection
   threshold. The human dataset is underpowered to detect the effect sizes seen
   in mouse; this is not evidence of non-replication.

2. **Age confound in GSE214373:** sCJD patients are significantly younger than
   controls (mean 71.6 vs. 83.3 years; p = 0.024). Age is included as a covariate
   in the primary model (VIF = 2.06). The sensitivity analysis (Section 5)
   quantifies the impact: 93% of sensitivity-model DEGs are age-confounded.

3. **GSE214373 GSEA returned 0 results:** After ortholog mapping, fewer than 3
   MCPS genes had valid rank statistics in the human DEG list, falling below
   fgsea's minimum gene set size. The GSEA leading-edge criterion for evidence
   tiers cannot be applied. Pathway-level human validation (GSE124571, N=10+10)
   is the primary human evidence layer.

3. **Terminal timepoint wpi unknown:** GSE210128 GEO metadata lists terminal
   samples as "term wpi" without a numeric value. Terminal samples are excluded
   from the continuous LRT (Arm A) and used only in the terminal Wald test (Arm B).
   `TERMINAL_WPI = 22L` is used as a filter key only, not as a regression value.

4. **Blood cohort:** GSE210128 whole blood has only the main cohort (no validation
   cohort), giving 3 replicates per condition per timepoint. LRT power is low.

5. **Replication cascade Stage 2:** Reports mean per-study replication count
   (not union across Stratum D studies). This is a conservative estimate.
   See `report_integration.md` Section 6.4.6 for full methodology.

6. **Pathway convergence in Script 08:** Uses gene-overlap fraction as a proxy
   for pathway enrichment in Streams 2 and 3 (ORA is not re-run). This is an
   approximation; full ORA on Stream 2/3 DEGs is available in Scripts 06/07.

7. **Cross-tissue multiple testing:** Tissues in Stream 3 are analysed
   independently. No cross-tissue FDR correction is applied (each tissue is
   a separate biological question).

---

## Citation

If you use this pipeline or its results, please cite the associated manuscript
(citation to be added upon publication) and this repository:

`https://github.com/[USERNAME]/prion_rnaseq_meta_analysis` — *[REPOSITORY LINK — TO BE INSERTED BY AUTHORS]*

Archived intermediates (Salmon quantifications, QC reports): Zenodo DOI
*[TO BE INSERTED AFTER ARCHIVAL]* — see [Large data files](#large-data-files).

## License

MIT — see [LICENSE](LICENSE).
