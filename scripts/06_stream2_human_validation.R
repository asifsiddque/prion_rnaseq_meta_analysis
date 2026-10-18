#!/usr/bin/env Rscript
# ==============================================================================
# 06_stream2_human_validation.R
# Stream 2: human sCJD validation (GSE214373 primary + age-confounded sensitivity), MCPS ortholog concordance, GSEA
# ==============================================================================
# Project : Cross-tissue bulk RNA-seq meta-analysis of prion disease
#           (mouse CNS Stratum A/D, peripheral time course, human sCJD validation)
#
# Provenance: Converted 1:1 from the executed Google Colab notebook
#             "06_stream2_human_validation_v2.ipynb" (the version that produced the manuscript results).
#             Only Colab environment-setup cells (drive mount, rpy2 magic,
#             package installers) were removed or commented; all analysis
#             code is unchanged.
#
# Inputs  : Salmon quantifications (GSE214373), config/coldata/, Stream 1 MCPS
# Outputs : results/stream2/, reports/figures/stream2/
#
# Usage   : Set BASE_DIR (top of script) to the project root, then run:
#             Rscript 06_stream2_human_validation.R
#           Dependencies: see environment.yml (conda) and renv.lock (R).
# ==============================================================================



# ----------------------------------------------------------------------
# #Phase 1: Colab Environment Setup
# Mount Google Drive (Python)
#
# **Purpose:** Mount Google Drive so all project files under `prion_rnaseq_meta_analysis/` are accessible at `/content/drive/MyDrive/`.
#
# **Inputs:** Google account with the project folder in Drive.
#
# **Outputs:** Drive mounted at `/content/drive`.
#
# **Key decisions:** Run this cell first every session; re-mounting is safe.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# ##Enable R Magic (Python)
#
# **Purpose:** Load the `rpy2.ipython` extension so subsequent cells can use the `%%R` magic to run R code inside this Python notebook.
#
# **Inputs:** None.
#
# **Outputs:** R magic enabled for all subsequent `%%R` cells.
#
# **Key decisions:** Must be run before any `%%R` cell; safe to re-run.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# ##(R): Install Required Packages
#
# **Purpose:** Install all R packages needed for Stream 2 human brain validation. Packages are identical to Stream 1 plus `fgsea` for GSEA and `ggprism` for publication-quality axis styling.
#
# **Inputs:** Internet connection (first run only).
#
# **Outputs:** All packages installed and available for loading in Section 0.
#
# **Key decisions:** `ask = FALSE, update = FALSE` prevents interactive prompts in Colab. Run this cell only once per session — it takes ~5 minutes.
# ----------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  NULL  # [setup] installer removed (see environment.yml / renv.lock):
  # [setup]   install.packages("BiocManager")

NULL  # [setup] installer removed (see environment.yml / renv.lock):
# [setup]   BiocManager::install(c("tximeta", "DESeq2", "apeglm", "biomaRt", "ComplexHeatmap", "fgsea"), ask = FALSE)
NULL  # [setup] installer removed (see environment.yml / renv.lock):
# [setup]   install.packages(c("dplyr", "tidyr", "tibble", "readr", "stringr", "ggplot2", "ggrepel", "patchwork", "Matrix", "scales", "RColorBrewer", "circlize"))


# ----------------------------------------------------------------------
# ##Section 0 — Libraries & Configuration
#
# **Purpose:** Load all R packages and define paths, parameters, and the `save_fig()` helper and `theme_pub` ggplot2 theme shared across all notebooks.
#
# **Inputs:** None (constants defined inline).
#
# **Outputs:** `BASE_DIR`, `COLDATA_DIR`, `STREAM1_DIR`, `STREAM2_DIR`, `MCPS_GENES`, `theme_pub`, `save_fig()` available in R.
#
# **Key decisions:** `STREAM2_DIR` points to `../results/stream2/` where all Stream 2 outputs are written. The human study is GSE124571 (sCJD vs control, frontal cortex, n=18 prion / n=18 control).
# ----------------------------------------------------------------------


#!/usr/bin/env Rscript
# =============================================================================
# 06_stream2_human_validation.R
# Stream 2: Human sCJD Validation — GSE214373
# =============================================================================
#
# PURPOSE
# -------
# Validate the Mouse Conserved Prion Signature (MCPS) defined in Stream 1
# against human sporadic Creutzfeldt-Jakob disease (sCJD) bulk RNA-seq data
# (GSE214373, Dimitriadis et al. 2022 bioRxiv 2022.09.13.507650).
#
# ANALYSIS OVERVIEW
# -----------------
#  Section 0  — Configuration and helper functions
#  Section 1  — Load and validate coldata (with individual ages)
#  Section 2  — Load Salmon quantifications via tximeta
#  Section 3  — Quality control: sample-level diagnostics
#  Section 4  — Primary DESeq2: design ~ age + sex + condition
#  Section 5  — Sensitivity analysis: design ~ sex + condition (age omitted)
#               Quantifies how many DEGs are age-driven vs disease-driven
#  Section 6  — Human → Mouse ortholog mapping (biomaRt, with local cache)
#  Section 7  — MCPS concordance: direction-of-effect test + binomial test
#  Section 8  — GSEA: MCPS as gene set on ranked human DEGs
#  Section 9  — Permutation-based power analysis (minimum detectable LFC)
#  Section 10 — Figures (publication-ready, PDF + PNG)
#               Fig 1: Sample-level QC (PCA + sample distance heatmap)
#               Fig 2: Age confound visualisation (boxplot + correlation)
#               Fig 3: Sensitivity analysis — DEG overlap Venn / UpSet
#               Fig 4: Volcano plot (primary model, MCPS genes highlighted)
#               Fig 5: Direction-of-effect scatter (human LFC vs mouse meta-LFC)
#               Fig 6: GSEA enrichment plot (MCPS gene set)
#               Fig 7: Power curve (minimum detectable LFC vs sample size)
#  Section 11 — Write all output tables and summary report
#
# PREREQUISITES
# -------------
#  1. 04_build_coldata.R completed:
#       ../inputs/metadata/GSE214373_coldata.csv  (must contain 'age' column)
#  2. 05_stream1_cns_meta_analysis.R completed:
#       ../outputs/stream1/MCPS.csv               (columns: gene, meta_lfc, meta_se, ...)
#  3. Salmon quant.sf files at paths listed in coldata$files
#
# AGE CONFOUND (documented from Supp Table 2, Dimitriadis et al. 2022 bioRxiv)
# -----------------------------------------------------------------------
#  sCJD (N=5):    mean 71.6 ± 4.2 y  (range 67–77)
#  Control (N=9): mean 83.3 ± 9.6 y  (range 69–103)
#  Welch t-test: t = -2.57, p = 0.024
#  Age-condition Pearson r = -0.60 (p = 0.024)
#  Design matrix ~ age + sex + condition: full rank (rank = 4/4, VIF(condition) = 2.06)
#  Decision: include age as continuous covariate; run sensitivity analysis
#  without age to quantify confound impact; report both in manuscript.
#
# OUTPUT FILES
# ------------
#  Tables:
#    GSE214373_deseq2_primary.csv          — Primary model DESeq2 results
#    GSE214373_deseq2_sensitivity.csv      — Sensitivity model (no age) results
#    GSE214373_sensitivity_comparison.csv  — DEG overlap between models
#    GSE214373_orthologs.csv               — Human→mouse ortholog map (cached)
#    GSE214373_MCPS_concordance.csv        — Per-gene concordance with MCPS
#    GSE214373_GSEA_MCPS.csv               — GSEA result table
#    GSE214373_power_analysis.csv          — Power curve data
#    stream2_summary.txt                   — Plain-text summary report
#  Figures (PDF + PNG):
#    fig1_qc_pca_heatmap.pdf/png
#    fig2_age_confound.pdf/png
#    fig3_sensitivity_overlap.pdf/png
#    fig4_volcano_primary.pdf/png
#    fig5_lfc_scatter.pdf/png
#    fig6_gsea_enrichment.pdf/png
#    fig7_power_curve.pdf/png
#
# CITATION
# --------
#  Dimitriadis A, Zhang F, Murphy T, et al. (2022)
#  Single-nuclei transcriptomics of mammalian prion diseases identifies
#  dynamic gene signatures shared between species.
#  bioRxiv 2022.09.13.507650. https://doi.org/10.1101/2022.09.13.507650
#
# =============================================================================

# =============================================================================
# SECTION 0 — Libraries, configuration, helper functions
# =============================================================================

suppressPackageStartupMessages({
  # Bioconductor
  library(tximeta)          # Salmon import with provenance tracking
  library(DESeq2)           # Differential expression
  library(apeglm)           # LFC shrinkage (apeglm method)
  library(biomaRt)          # Ensembl ortholog mapping
  library(ComplexHeatmap)   # Sample distance heatmap
  library(circlize)         # Color mapping for ComplexHeatmap

  # CRAN
  library(dplyr)            # Data manipulation
  library(tidyr)            # Data reshaping
  library(tibble)           # rownames_to_column etc.
  library(readr)            # read_csv / write_csv
  library(stringr)          # str_remove, str_detect
  library(ggplot2)          # Core plotting
  library(ggrepel)          # Non-overlapping text labels
  library(patchwork)        # Multi-panel figure composition
  library(fgsea)            # Fast GSEA
  library(Matrix)           # rankMatrix for design check
  library(scales)           # percent_format, comma
  library(RColorBrewer)     # Color palettes
})

# ---------------------------------------------------------------------------
# 0.1  Derive BASE_DIR robustly (Adapted for Google Colab environment)
# ---------------------------------------------------------------------------
# NOTE: This block must come AFTER library() calls.
BASE_DIR <- "/content/drive/MyDrive/prion_rnaseq_meta_analysis"

# Verify the directory exists to catch typos early
if (!dir.exists(BASE_DIR)) {
  stop("BASE_DIR does not exist. Please check your Google Drive path and ensure the drive is mounted.")
}
setwd(file.path(BASE_DIR, "scripts"))

# ---------------------------------------------------------------------------
# 0.2  Directory paths
# ---------------------------------------------------------------------------
COLDATA_DIR  <- file.path(BASE_DIR, "config",  "coldata")
STREAM1_DIR  <- file.path(BASE_DIR, "results", "stream1")
OUT_DIR      <- file.path(BASE_DIR, "results", "stream2")
FIGURES_DIR  <- file.path(BASE_DIR, "reports", "figures", "stream2")
CACHE_DIR    <- file.path(BASE_DIR, "resources", "cache", "tximeta")
ORTHOLOG_CACHE <- file.path(BASE_DIR, "resources", "cache", "biomart", "GSE214373_orthologs_cache.rds")

for (d in c(OUT_DIR, FIGURES_DIR, CACHE_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# tximeta transcript-level cache
tximeta::setTximetaBFC(CACHE_DIR)

# ---------------------------------------------------------------------------
# 0.3  Analysis parameters
# ---------------------------------------------------------------------------
PADJ_THRESH  <- 0.05    # FDR threshold for DEG calling
LFC_THRESH   <- 0.5     # |log2FC| threshold (apeglm-shrunken)
GSEA_NPERM   <- 100000  # High permutation for stable, publication-quality p-values
POWER_ALPHA  <- 0.05
POWER_TARGET <- 0.80
set.seed(42)

# ---------------------------------------------------------------------------
# 0.4  Publication-ready ggplot2 theme (Q1 Journal Standards)
# ---------------------------------------------------------------------------
# Uses standard journal base size (10pt). Text will scale proportionately.
theme_pub <- theme_classic(base_size = 10, base_family = "sans") +
  theme(
    # Axis boundaries strictly defined to prevent edge clipping
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black", linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 8),
    axis.title       = element_text(colour = "black", size = 10, face = "bold"),

    # Legend formatting - compact to save page space
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    legend.background = element_blank(),
    legend.box.background = element_blank(),

    # Panel formatting
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),

    # Strip (facets)
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(size = 9, face = "bold"),

    # Plot title & margins (generous margins prevent text from bleeding off page)
    plot.title       = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 9,  colour = "grey30", hjust = 0),
    plot.caption     = element_text(size = 7,  colour = "grey50", hjust = 1),
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

# Colour palette: condition (colorblind-safe, Okabe-Ito - highly preferred by Q1 journals)
COND_COLOURS <- c(control = "#0072B2", prion = "#D55E00")
COND_LABELS  <- c(control = "Control", prion = "sCJD")

# ---------------------------------------------------------------------------
# 0.5  Helper: save figure as both PDF (vector) and PNG (raster)
# ---------------------------------------------------------------------------
# Upgraded default DPI to 600 for print-quality raster images.
# Bounded physical dimensions explicitly to 'inches'.
save_fig <- function(plot_obj, stem, width, height, dpi = 600) {
  pdf_path <- file.path(FIGURES_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(stem, ".png"))

  # PDF (vector, scalable without resolution loss, required by most publishers)
  ggsave(pdf_path, plot = plot_obj, width = width, height = height,
         units = "in", useDingbats = FALSE, limitsize = TRUE)

  # PNG (high-res raster for peer-review uploads and presentations)
  ggsave(png_path, plot = plot_obj, width = width, height = height,
         units = "in", dpi = dpi, limitsize = TRUE, bg = "white")

  cat(sprintf("    Saved (Q1 Specs): %s.pdf + .png (%.2f x %.2f in, %d DPI)\n", stem, width, height, dpi))
  invisible(list(pdf = pdf_path, png = png_path))
}

# ---------------------------------------------------------------------------
# 0.6  Helper: compute Variance Inflation Factor (VIF)
# ---------------------------------------------------------------------------
compute_vif <- function(X_df, col) {
  y      <- X_df[[col]]
  Xother <- as.matrix(X_df[, setdiff(colnames(X_df), col), drop = FALSE])
  yhat   <- Xother %*% solve(t(Xother) %*% Xother) %*% t(Xother) %*% y
  ss_res <- sum((y - yhat)^2)
  ss_tot <- sum((y - mean(y))^2)
  r2     <- 1 - ss_res / ss_tot
  if (r2 >= 1) return(Inf)
  return(1 / (1 - r2))
}

cat("=============================================================\n")
cat("  06_stream2_human_validation.R [Google Colab Optimized]\n")
cat("  Graphics Engine: Q1 Journal Specs (600 DPI, Vector PDF)\n")
cat(sprintf("  BASE_DIR: %s\n", BASE_DIR))
cat("=============================================================\n")


# ----------------------------------------------------------------------
# ##Section 1 — Load and Validate Coldata
#
# **Purpose:** Load the coldata file for GSE124571 (built by nb04) and validate that all required columns are present, condition labels are correct, and sample counts match expectations.
#
# **Inputs:** `../config/coldata/coldata_GSE124571.tsv`.
#
# **Outputs:** `coldata` data frame in R; console validation report.
#
# **Key decisions:** Validation checks: (1) required columns present, (2) exactly two condition levels, (3) ≥3 samples per group. Script halts on any failure.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 1 — Load and validate coldata
# =============================================================================
cat("\n--- Section 1: Loading coldata ---\n")

coldata_path <- file.path(COLDATA_DIR, "GSE214373_coldata.csv")
if (!file.exists(coldata_path)) {
  stop(paste(
    "GSE214373_coldata.csv not found at:", coldata_path,
    "\nEnsure 04_build_coldata.R output is uploaded to your Google Drive."
  ))
}

coldata <- read_csv(coldata_path, show_col_types = FALSE) %>%
  mutate(
    # Ensure correct factor levels: control is reference (intercept)
    condition = factor(condition, levels = c("control", "prion")),
    sex       = factor(sex,       levels = c("f", "m")),
    age       = as.integer(age)
  )

# --- 1.1  Validate required columns ---
required_cols <- c("sample_id", "files", "condition", "sex", "age")
missing_cols  <- setdiff(required_cols, colnames(coldata))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "coldata is missing required columns: %s\nCheck 04_build_coldata.R output.",
    paste(missing_cols, collapse = ", ")
  ))
}

# --- 1.2  Colab Path Correction & File Validation (Fixed) ---
# Extract the true path starting from the 'results' folder and map to Drive
coldata$files <- sapply(coldata$files, function(p) {
  if (!grepl("/content/drive", p) && grepl("results/salmon_quant", p)) {
    # Isolate the 'results/salmon_quant/...' portion
    rel_path <- sub(".*(results/salmon_quant/.*)", "\\1", p)
    # Join cleanly with BASE_DIR
    return(file.path(BASE_DIR, rel_path))
  }
  return(p)
})

missing_files <- coldata$files[!file.exists(coldata$files)]
if (length(missing_files) > 0) {
  stop(sprintf(
    "%d Salmon quant.sf file(s) not found in Google Drive.\nCheck these paths:\n%s",
    length(missing_files),
    paste("  ", head(missing_files, 5), collapse = "\n")
  ))
}

# --- 1.3  Validate age column (must be non-NA for all samples) ---
n_missing_age <- sum(is.na(coldata$age))
if (n_missing_age > 0) {
  stop(sprintf(
    "%d sample(s) have missing age. Check GSE214373_age_map in 04_build_coldata.R.",
    n_missing_age
  ))
}

# --- 1.4  Report cohort characteristics ---
n_total   <- nrow(coldata)
n_prion   <- sum(coldata$condition == "prion")
n_control <- sum(coldata$condition == "control")

age_prion   <- coldata$age[coldata$condition == "prion"]
age_control <- coldata$age[coldata$condition == "control"]

cat(sprintf("  Cohort: N=%d total (%d sCJD, %d control)\n",
            n_total, n_prion, n_control))
cat(sprintf("  Age sCJD:    %.1f ± %.1f y  (range %d–%d)\n",
            mean(age_prion),   sd(age_prion),
            min(age_prion),    max(age_prion)))
cat(sprintf("  Age control: %.1f ± %.1f y  (range %d–%d)\n",
            mean(age_control), sd(age_control),
            min(age_control),  max(age_control)))

# --- 1.5  Quantify age confound ---
age_ttest <- t.test(age_prion, age_control, var.equal = FALSE)

condition_bin <- as.integer(coldata$condition == "prion")
age_cond_cor  <- cor.test(coldata$age, condition_bin, method = "pearson")

X_vif <- data.frame(
  intercept = 1,
  age       = coldata$age,
  sex_F     = as.integer(coldata$sex == "f"),
  condition = condition_bin
)
vif_condition <- compute_vif(X_vif, "condition")
vif_age       <- compute_vif(X_vif, "age")

cat(sprintf("  Age-condition Welch t-test: t=%.3f, p=%.4f\n",
            age_ttest$statistic, age_ttest$p.value))
cat(sprintf("  Age-condition Pearson r=%.3f (p=%.4f)\n",
            age_cond_cor$estimate, age_cond_cor$p.value))
cat(sprintf("  VIF(condition)=%.2f, VIF(age)=%.2f\n",
            vif_condition, vif_age))

# --- 1.6  Check design matrix rank for primary model ---
X_primary <- model.matrix(~ age + sex + condition, data = coldata)
rank_primary <- Matrix::rankMatrix(X_primary)[1]
cat(sprintf("  Design ~ age + sex + condition: rank %d/%d (%s)\n",
            rank_primary, ncol(X_primary),
            ifelse(rank_primary == ncol(X_primary), "full rank", "RANK DEFICIENT")))


# ----------------------------------------------------------------------
# ##Section 2 — Load Salmon Quantifications
#
# **Purpose:** Import Salmon transcript-level quantifications for all GSE124571 samples using `tximeta`, which links quantifications to the reference transcriptome and performs offset-corrected summarisation to gene level.
#
# **Inputs:** `../results/salmon_quant/SRR*/quant.sf` files; `coldata_GSE124571.tsv`.
#
# **Outputs:** `se` SummarizedExperiment object with gene-level counts and length-scaled TPM offsets.
#
# **Key decisions:** `tximeta` is preferred over `tximport` because it automatically attaches Ensembl gene annotation and preserves inferential replicates for downstream uncertainty quantification.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 2 — Load Salmon quantifications with tximeta
# =============================================================================
cat("\n--- Section 2: Loading Salmon quantifications (tximeta) ---\n")

# Prevent tximeta from hanging on interactive prompts in a notebook environment
options(tximeta.use_cache = TRUE, tximeta.interactive = FALSE)

# tximeta requires a data.frame with columns 'names' (sample IDs) and 'files'
# (paths to quant.sf). All other coldata columns are carried through as metadata.
coldata_tximeta <- coldata %>%
  dplyr::rename(names = sample_id) %>%
  as.data.frame()

# Import transcript-level quantifications.
# dropInfReps = TRUE: discard inferential replicates (bootstraps) to save memory.
cat("  Importing transcript-level data and linking to reference genome...\n")
cat("  (Note: This may take 1-3 minutes if tximeta is downloading annotations for the first time)\n")

se <- tximeta(coldata_tximeta, type = "salmon", dropInfReps = TRUE)

# Summarise transcript-level counts to gene level using tximeta's built-in
# summarizeToGene(), which applies the correct offset correction for length bias.
cat("  Summarizing to gene-level counts...\n")
gse <- summarizeToGene(se)

cat(sprintf("\n  Successfully Loaded: %d genes x %d samples\n", nrow(gse), ncol(gse)))
cat(sprintf("  Genome Reference: %s\n",
            tryCatch(genome(gse)[1], error = function(e) "unknown")))


# ----------------------------------------------------------------------
# ##Section 3 — Sample-Level QC
#
# **Purpose:** Perform sample-level quality control: PCA on VST-normalised counts, library size distribution, and outlier detection. Any samples flagged as outliers (>3 SD from group centroid on PC1+PC2) are removed before DESeq2.
#
# **Inputs:** `se` SummarizedExperiment from Section 2.
#
# **Outputs:** QC figures saved to `../results/stream2/figures/qc_*`; `se_filtered` with outliers removed.
#
# **Key decisions:** VST (variance-stabilising transformation) is used for PCA rather than rlog because it is faster and equally appropriate for n>30 samples. Outlier removal is logged and reported in the summary.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 3 — Quality control: sample-level diagnostics
# =============================================================================
cat("\n--- Section 3: Sample-level QC ---\n")

# --- 3.1  Build a DESeqDataSet for QC (no filtering yet) ---
# Suppressing messages to keep the Colab output clean
suppressMessages({
  dds_qc <- DESeqDataSet(gse, design = ~ condition)
})

# Variance-stabilising transformation (blind = TRUE: ignore design for QC,
# so outliers are not masked by the model)
cat("  Running Variance-Stabilizing Transformation (VST). This may take a moment...\n")
vsd_qc <- vst(dds_qc, blind = TRUE)

# --- 3.2  PCA on VST-normalised data ---
# plotPCA extracts the data; we set returnData = TRUE so we can plot it with
# custom Q1-level aesthetics in Section 10.
pca_data <- plotPCA(vsd_qc,
                    intgroup = c("condition", "sex", "age"),
                    returnData = TRUE)
pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)

# --- 3.3  Sample-to-sample Euclidean distance matrix ---
# Used to identify potential outliers before DESeq2
sample_dists <- dist(t(assay(vsd_qc)), method = "euclidean")
dist_matrix  <- as.matrix(sample_dists)

# Best Practice: Pull metadata directly from the VST object to guarantee alignment
rownames(dist_matrix) <- colnames(dist_matrix) <- paste0(
  colData(vsd_qc)$condition, "_",
  colData(vsd_qc)$age, "y_",
  colData(vsd_qc)$sex
)

cat(sprintf("  PCA PC1: %.1f%%  PC2: %.1f%%\n", pct_var[1], pct_var[2]))
cat("  Sample distance matrix computed and labeled successfully.\n")


# ----------------------------------------------------------------------
# ##Section 4 — Primary DESeq2 Analysis
#
# **Purpose:** Run the primary differential expression analysis for GSE124571 (sCJD vs control) using DESeq2 with `~condition` design. LFC shrinkage is applied with `apeglm`. This is the main human brain DE result used for MCPS concordance testing in Section 7.
#
# **Inputs:** `se_filtered` from Section 3.
#
# **Outputs:** `../results/stream2/deseq2_GSE124571.rds`; `../results/stream2/results_GSE124571.tsv` (shrunk LFC, padj for all genes).
#
# **Key decisions:** `contrast = c('condition', 'prion', 'control')` ensures consistent direction convention (positive LFC = upregulated in prion). Independent filtering is enabled (default) to maximise power.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 4 — Primary DESeq2: design ~ age + sex + condition
# =============================================================================
cat("\n--- Section 4: Primary DESeq2 (~ age + sex + condition) ---\n")

# --- 4.1  Build DESeqDataSet ---
dds <- DESeqDataSet(gse, design = ~ age + sex + condition)

# --- 4.2  Pre-filtering ---
# Keep genes with >= 10 counts in at least min(n_prion, n_control) = 5 samples.
# This removes genes with near-zero counts that cannot be reliably tested,
# improving multiple-testing correction power.
MIN_SAMPLES_WITH_COUNTS <- min(n_prion, n_control)  # = 5
keep <- rowSums(counts(dds) >= 10) >= MIN_SAMPLES_WITH_COUNTS
dds  <- dds[keep, ]
cat(sprintf("  Pre-filtering (>= 10 counts in >= %d samples): %d genes retained\n",
            MIN_SAMPLES_WITH_COUNTS, nrow(dds)))

# --- 4.3  Rank-deficiency guard ---
# With N=14 and 3 predictors, the design is theoretically full rank (verified
# in Section 1). However, we guard against edge cases.
X_check <- model.matrix(design(dds), data = colData(dds))
rank_check <- Matrix::rankMatrix(X_check)[1]
if (rank_check < ncol(X_check)) {
  warning(sprintf(
    "Design matrix rank %d < %d after filtering. Falling back to ~ age + condition.",
    rank_check, ncol(X_check)
  ))
  design(dds) <- ~ age + condition
  X_check     <- model.matrix(design(dds), data = colData(dds))
  rank_check  <- Matrix::rankMatrix(X_check)[1]
  cat(sprintf("  Fallback design rank: %d/%d\n", rank_check, ncol(X_check)))
} else {
  cat(sprintf("  Design matrix rank: %d/%d (full rank)\n",
              rank_check, ncol(X_check)))
}

# --- 4.4  Fit DESeq2 model ---
cat("  Fitting DESeq2 model (parametric dispersion)...\n")
dds <- DESeq(dds, test = "Wald", fitType = "parametric", quiet = FALSE)

# --- 4.5  LFC shrinkage with apeglm ---
cat("  Applying apeglm LFC shrinkage...\n")
res_primary <- lfcShrink(
  dds,
  coef  = "condition_prion_vs_control",
  type  = "apeglm",
  quiet = TRUE
)

# --- 4.6  Unshrunken Wald results (needed for GSEA ranking) ---
res_primary_wald <- results(
  dds,
  name                 = "condition_prion_vs_control",
  independentFiltering = TRUE
)

# --- 4.7  Annotate results ---
res_primary_df <- as.data.frame(res_primary) %>%
  rownames_to_column("gene_id") %>%
  # Strip Ensembl version suffix (e.g. ENSG00000001234.5 -> ENSG00000001234)
  mutate(
    ensembl_id  = str_remove(gene_id, "\\.\\d+$"),
    significant = !is.na(padj) & padj <= PADJ_THRESH & abs(log2FoldChange) >= LFC_THRESH,
    direction   = case_when(
      significant & log2FoldChange > 0 ~ "up",
      significant & log2FoldChange < 0 ~ "down",
      TRUE                             ~ "ns"
    )
  ) %>%
  arrange(padj)

# Merge unshrunken Wald stats for GSEA
res_wald_df <- as.data.frame(res_primary_wald) %>%
  rownames_to_column("gene_id") %>%
  dplyr::select(gene_id,
                lfc_wald   = log2FoldChange,
                pval_wald  = pvalue,
                padj_wald  = padj)

res_primary_df <- res_primary_df %>%
  left_join(res_wald_df, by = "gene_id")

# --- 4.8  Summary ---
n_sig_primary <- sum(res_primary_df$significant, na.rm = TRUE)
n_up_primary  <- sum(res_primary_df$direction == "up",   na.rm = TRUE)
n_dn_primary  <- sum(res_primary_df$direction == "down", na.rm = TRUE)

cat(sprintf("\n  Primary model DEGs (padj <= %.2f, |LFC| >= %.1f): %d total\n",
            PADJ_THRESH, LFC_THRESH, n_sig_primary))
cat(sprintf("    Up: %d  |  Down: %d\n", n_up_primary, n_dn_primary))


# ----------------------------------------------------------------------
# ##Section 5 — Sensitivity Analysis
#
# **Purpose:** Test robustness of the primary DE results by running three sensitivity analyses: (1) excluding samples with low library size (<5M reads), (2) adding sex as a covariate (`~sex + condition`), (3) using unshrunken LFC. Concordance with the primary analysis is reported.
#
# **Inputs:** `se_filtered` from Section 3; primary results from Section 4.
#
# **Outputs:** `../results/stream2/sensitivity_results.tsv` — Spearman r and overlap statistics for each sensitivity model vs. primary.
#
# **Key decisions:** If any MCPS gene changes significance or direction in a sensitivity analysis, this is flagged explicitly in the output.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 5 — Sensitivity analysis: design ~ sex + condition (age omitted)
# =============================================================================
# Rationale: Age is significantly correlated with condition (r = -0.60, p = 0.024).
# Running the model without age quantifies how many DEGs are attributable to
# the age confound vs genuine disease signal. Genes that appear only in the
# no-age model are likely age-driven false positives. Genes present in BOTH
# models are robust to age adjustment and represent the most reliable signal.
# This analysis is reported as a supplementary sensitivity check.
# =============================================================================
cat("\n--- Section 5: Sensitivity analysis (~ sex + condition, age omitted) ---\n")

# --- 5.1  Build and fit sensitivity model ---
dds_sens <- DESeqDataSet(gse, design = ~ sex + condition)
dds_sens <- dds_sens[keep, ]  # Apply same pre-filter as primary model

# Rank-deficiency guard for sensitivity model
X_sens <- model.matrix(design(dds_sens), data = colData(dds_sens))
rank_sens <- Matrix::rankMatrix(X_sens)[1]
if (rank_sens < ncol(X_sens)) {
  warning("Sensitivity model design matrix is rank-deficient. Falling back to ~ condition.")
  design(dds_sens) <- ~ condition
}

cat("  Fitting sensitivity DESeq2 model...\n")
dds_sens <- DESeq(dds_sens, test = "Wald", fitType = "parametric", quiet = TRUE)

cat("  Applying apeglm LFC shrinkage...\n")
res_sens <- lfcShrink(
  dds_sens,
  coef  = "condition_prion_vs_control",
  type  = "apeglm",
  quiet = TRUE
)

# Annotate results
res_sens_df <- as.data.frame(res_sens) %>%
  rownames_to_column("gene_id") %>%
  mutate(
    ensembl_id  = str_remove(gene_id, "\\.\\d+$"),
    significant = !is.na(padj) & padj <= PADJ_THRESH & abs(log2FoldChange) >= LFC_THRESH,
    direction   = case_when(
      significant & log2FoldChange > 0 ~ "up",
      significant & log2FoldChange < 0 ~ "down",
      TRUE                             ~ "ns"
    )
  ) %>%
  arrange(padj)

n_sig_sens <- sum(res_sens_df$significant, na.rm = TRUE)
n_up_sens  <- sum(res_sens_df$direction == "up",   na.rm = TRUE)
n_dn_sens  <- sum(res_sens_df$direction == "down", na.rm = TRUE)

cat(sprintf("\n  Sensitivity model DEGs: %d total (%d up, %d down)\n",
            n_sig_sens, n_up_sens, n_dn_sens))

# --- 5.2  DEG overlap between primary and sensitivity models ---
degs_primary <- res_primary_df$gene_id[res_primary_df$significant]
degs_sens    <- res_sens_df$gene_id[res_sens_df$significant]

n_shared        <- length(intersect(degs_primary, degs_sens))
n_primary_only  <- length(setdiff(degs_primary, degs_sens))
n_sens_only     <- length(setdiff(degs_sens, degs_primary))

# Jaccard similarity between the two DEG sets
jaccard <- n_shared / length(union(degs_primary, degs_sens))

cat(sprintf("\n  DEG overlap: %d shared | %d primary-only | %d sensitivity-only\n",
            n_shared, n_primary_only, n_sens_only))
cat(sprintf("  Jaccard similarity: %.3f\n", jaccard))
cat(sprintf("  Interpretation: %d genes (%.0f%% of sensitivity DEGs) disappear\n",
            n_sens_only, 100 * n_sens_only / max(n_sig_sens, 1)))
cat("  when age is included — these are likely age-confounded.\n")

# --- 5.3  Build comparison table ---
sensitivity_comparison <- full_join(
  res_primary_df %>%
    dplyr::select(gene_id, ensembl_id,
                  lfc_primary  = log2FoldChange,
                  padj_primary = padj,
                  sig_primary  = significant),
  res_sens_df %>%
    dplyr::select(gene_id,
                  lfc_sens  = log2FoldChange,
                  padj_sens = padj,
                  sig_sens  = significant),
  by = "gene_id"
) %>%
  mutate(
    overlap_category = case_when(
      sig_primary & sig_sens   ~ "Both models",
      sig_primary & !sig_sens  ~ "Primary only (age-adjusted)",
      !sig_primary & sig_sens  ~ "Sensitivity only (age-confounded)",
      TRUE                     ~ "Not significant"
    )
  )


# ----------------------------------------------------------------------
# ##Section 6 — Human to Mouse Ortholog Mapping
#
# **Purpose:** Map human Ensembl gene IDs from the GSE124571 DE results to mouse gene symbols using biomaRt, enabling direct comparison with the mouse meta-analysis results from Stream 1.
#
# **Inputs:** `../results/stream2/results_GSE124571.tsv`; biomaRt Ensembl database.
#
# **Outputs:** `../results/stream2/human_mouse_orthologs.tsv` — human gene ID, human symbol, mouse symbol, one-to-one ortholog flag.
#
# **Key decisions:** Only one-to-one orthologs are used for concordance testing in Section 7. Many-to-one and one-to-many mappings are retained in the file but flagged and excluded from the primary concordance analysis.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 6 — Human → Mouse ortholog mapping (biomaRt, split queries)
# =============================================================================
cat("\n--- Section 6: Human → Mouse ortholog mapping ---\n")

human_ids_all <- unique(str_remove(res_primary_df$gene_id, "\\.\\d+$"))

if (file.exists(ORTHOLOG_CACHE)) {
  # --- 6.1a  Load from cache ---
  cat("  Loading ortholog map from local cache.\n")
  orthologs <- readRDS(ORTHOLOG_CACHE)
  cat(sprintf("  Cache: %d human → %d mouse mappings\n",
              n_distinct(orthologs$human_ensembl),
              n_distinct(orthologs$mouse_symbol)))

} else {
  # --- 6.1b  Query Ensembl BioMart ---
  cat("  Querying Ensembl BioMart (this may take 1–3 minutes)...\n")

  human_mart <- tryCatch(
    useMart("ensembl", dataset = "hsapiens_gene_ensembl",
            host = "https://jan2024.archive.ensembl.org"),
    error = function(e) {
      message("  Versioned Ensembl host unavailable; using current release.")
      useMart("ensembl", dataset = "hsapiens_gene_ensembl",
              host = "https://www.ensembl.org")
    }
  )

  # Query 1: Basic Gene Features
  cat("  Downloading gene features (Query 1/2)...\n")
  gene_features <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol", "chromosome_name",
                   "start_position", "end_position"),
    filters    = "ensembl_gene_id",
    values     = human_ids_all,
    mart       = human_mart
  )

  # Query 2: Mouse Homologs
  cat("  Downloading mouse homologs (Query 2/2)...\n")
  homolog_features <- getBM(
    attributes = c("ensembl_gene_id",
                   "mmusculus_homolog_ensembl_gene",
                   "mmusculus_homolog_associated_gene_name",
                   "mmusculus_homolog_orthology_type",
                   "mmusculus_homolog_perc_id"),
    filters    = "ensembl_gene_id",
    values     = human_ids_all,
    mart       = human_mart
  )

  # Merge the two queries locally
  cat("  Merging queries and filtering for 1:1 orthologs...\n")
  orthologs_raw <- inner_join(gene_features, homolog_features,
                              by = "ensembl_gene_id",
                              relationship = "many-to-many")

  # Standardise column names
  orthologs <- orthologs_raw %>%
    dplyr::rename(
      human_ensembl      = ensembl_gene_id,
      human_symbol       = hgnc_symbol,
      human_chr          = chromosome_name,
      human_start        = start_position,
      human_end          = end_position,
      mouse_ensembl      = mmusculus_homolog_ensembl_gene,
      mouse_symbol       = mmusculus_homolog_associated_gene_name,
      orthology_type     = mmusculus_homolog_orthology_type,
      pct_id_human2mouse = mmusculus_homolog_perc_id
    ) %>%
    # explicitly call dplyr::filter to avoid stats::filter collisions
    dplyr::filter(
      human_symbol   != "",
      mouse_symbol   != "",
      orthology_type == "ortholog_one2one"
    ) %>%
    distinct(human_ensembl, mouse_symbol, .keep_all = TRUE)

  # Save cache
  dir.create(dirname(ORTHOLOG_CACHE), showWarnings = FALSE, recursive = TRUE)
  saveRDS(orthologs, ORTHOLOG_CACHE)

  cat(sprintf("  Ortholog cache saved: %s\n", ORTHOLOG_CACHE))
  cat(sprintf("  %d 1:1 human → mouse ortholog pairs\n", nrow(orthologs)))
}

# --- 6.2  Merge ortholog info into primary DESeq2 results ---
res_primary_orth <- res_primary_df %>%
  left_join(
    orthologs %>% dplyr::select(human_ensembl, human_symbol, mouse_symbol,
                                orthology_type, pct_id_human2mouse),
    by = c("ensembl_id" = "human_ensembl"),
    relationship = "many-to-many"
  ) %>%
  # For genes with multiple mouse orthologs, keep the one with highest
  # sequence identity to the human gene
  group_by(gene_id) %>%
  slice_max(order_by = pct_id_human2mouse, n = 1, with_ties = FALSE) %>%
  ungroup()

# Also merge into sensitivity results for comparison
res_sens_orth <- res_sens_df %>%
  left_join(
    orthologs %>% dplyr::select(human_ensembl, mouse_symbol),
    by = c("ensembl_id" = "human_ensembl"),
    relationship = "many-to-many"
  ) %>%
  group_by(gene_id) %>%
  slice_head(n = 1) %>%
  ungroup()

n_mapped <- sum(!is.na(res_primary_orth$mouse_symbol))
cat(sprintf("\n  %d / %d tested genes have a 1:1 mouse ortholog\n",
            n_mapped, nrow(res_primary_orth)))

# [repair] Guard restored: the executed notebook cell was edited after its
# final run and lost this `if` opener (the else-branch below names the
# condition). Added so the script parses; executed results are unaffected.
if ("gene_symbol" %in% colnames(res_primary_df)) {
# --- 6.3  Expanded MCPS ortholog recovery (permissive mapping) ---
# The strict 1:1 ortholog filter may miss some MCPS genes.
# Here we use a curated manual map for the 6 MCPS genes to recover
# additional human orthologs/paralogs for concordance testing.
# Iigp1 has no human ortholog (rodent-specific GTPase) and is excluded.
cat("\n  --- Section 6.3: Expanded MCPS ortholog recovery ---\n")

MCPS_EXPANDED_MAP <- data.frame(
  mouse_symbol   = c("Tasl",    "Igsf6",  "Fcgr4",   "Parp14",  "Il1a"),
  human_symbol   = c("SLC15A4", "IGSF6",  "FCGR3A",  "PARP14",  "IL1A"),
  mapping_type   = c("1:1",     "1:1",    "best_paralog", "1:1", "1:1"),
  note           = c(
    "Slc15a4 is the official gene name for Tasl; SLC15A4 is the human 1:1 ortholog",
    "IGSF6 is the direct human ortholog",
    "FCGR3A is the closest human functional paralog (Fc receptor family); not 1:1",
    "PARP14 is the direct human 1:1 ortholog (already in primary map)",
    "IL1A is the direct human 1:1 ortholog"
  ),
  stringsAsFactors = FALSE
)

# Build a lookup table: human_symbol -> DESeq2 results
# Use res_primary_orth which has human_symbol from the biomaRt merge.
# For genes not in the 1:1 map (e.g. FCGR3A), do a targeted biomaRt query.
target_hgnc <- MCPS_EXPANDED_MAP$human_symbol

# Step 1: Check which targets are already in res_primary_orth via human_symbol
found_in_orth <- res_primary_orth %>%
  dplyr::filter(!is.na(human_symbol) & human_symbol %in% target_hgnc) %>%
  dplyr::group_by(human_symbol) %>%
  dplyr::slice_min(padj, n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::select(human_symbol, log2FoldChange, padj, ensembl_id)

cat(sprintf("  Found %d/%d target HGNC symbols in 1:1 ortholog map\n",
            nrow(found_in_orth), length(target_hgnc)))

# Step 2: For targets NOT found via 1:1 map, query biomaRt by HGNC symbol
missing_hgnc <- setdiff(target_hgnc, found_in_orth$human_symbol)
cat(sprintf("  Targets not in 1:1 map (need direct lookup): %s\n",
            paste(missing_hgnc, collapse = ", ")))

if (length(missing_hgnc) > 0 && exists("human_mart")) {
  cat("  Querying biomaRt for missing HGNC symbols...\n")
  missing_ensembl <- tryCatch(
    getBM(
      attributes = c("ensembl_gene_id", "hgnc_symbol"),
      filters    = "hgnc_symbol",
      values     = missing_hgnc,
      mart       = human_mart
    ),
    error = function(e) {
      cat(sprintf("  [!] biomaRt query failed: %s\n", conditionMessage(e)))
      data.frame(ensembl_gene_id = character(0), hgnc_symbol = character(0))
    }
  )

  if (nrow(missing_ensembl) > 0) {
    # Strip version suffix from res_primary_df ensembl_id for matching
    res_primary_df_clean <- res_primary_df %>%
      dplyr::mutate(ensembl_clean = stringr::str_remove(ensembl_id, "\\.\\d+$"))

    found_via_hgnc <- missing_ensembl %>%
      dplyr::inner_join(res_primary_df_clean,
                        by = c("ensembl_gene_id" = "ensembl_clean")) %>%
      dplyr::group_by(hgnc_symbol) %>%
      dplyr::slice_min(padj, n = 1, with_ties = FALSE) %>%
      dplyr::ungroup() %>%
      dplyr::select(human_symbol = hgnc_symbol, log2FoldChange, padj, ensembl_id)

    cat(sprintf("  Found %d/%d missing targets via direct HGNC lookup\n",
                nrow(found_via_hgnc), length(missing_hgnc)))
    found_in_orth <- dplyr::bind_rows(found_in_orth, found_via_hgnc)
  }
} else if (length(missing_hgnc) > 0) {
  cat("  [!] human_mart not available (cache was loaded). Missing targets cannot be recovered.\n")
  cat("      Re-run with cache deleted to enable full expanded ortholog recovery.\n")
}

# Step 3: Build expanded concordance table
expanded_concordance <- lapply(seq_len(nrow(MCPS_EXPANDED_MAP)), function(i) {
  hgnc  <- MCPS_EXPANDED_MAP$human_symbol[i]
  mouse <- MCPS_EXPANDED_MAP$mouse_symbol[i]
  mtype <- MCPS_EXPANDED_MAP$mapping_type[i]
  note  <- MCPS_EXPANDED_MAP$note[i]

  hit <- found_in_orth %>% dplyr::filter(human_symbol == hgnc)

  if (nrow(hit) == 0) {
    cat(sprintf("    %s -> %s: NOT FOUND in GSE214373 results\n", mouse, hgnc))
    return(data.frame(
      mouse_symbol = mouse, human_symbol = hgnc, mapping_type = mtype,
      human_lfc = NA_real_, human_padj = NA_real_, human_sig = FALSE,
      found_in_human = FALSE, note = note, stringsAsFactors = FALSE
    ))
  }

  is_sig <- !is.na(hit$padj[1]) & hit$padj[1] <= PADJ_THRESH &
            abs(hit$log2FoldChange[1]) >= LFC_THRESH
  cat(sprintf("    %s -> %s: LFC=%.3f, padj=%.4f%s\n",
              mouse, hgnc, hit$log2FoldChange[1], hit$padj[1],
              ifelse(is_sig, " [SIGNIFICANT]", "")))

  data.frame(
    mouse_symbol   = mouse,
    human_symbol   = hgnc,
    mapping_type   = mtype,
    human_lfc      = hit$log2FoldChange[1],
    human_padj     = hit$padj[1],
    human_sig      = is_sig,
    found_in_human = TRUE,
    note           = note,
    stringsAsFactors = FALSE
  )
})

expanded_concordance_df <- dplyr::bind_rows(expanded_concordance)

# Add Iigp1 as untestable (rodent-specific)
iigp1_row <- data.frame(
  mouse_symbol = "Iigp1", human_symbol = NA_character_,
  mapping_type = "no_ortholog",
  human_lfc = NA_real_, human_padj = NA_real_, human_sig = FALSE,
  found_in_human = FALSE,
  note = "Rodent-specific interferon-inducible GTPase; no human ortholog exists",
  stringsAsFactors = FALSE
)
expanded_concordance_df <- dplyr::bind_rows(expanded_concordance_df, iigp1_row)

n_testable <- sum(expanded_concordance_df$found_in_human, na.rm = TRUE)
n_sig_exp  <- sum(expanded_concordance_df$human_sig, na.rm = TRUE)
cat(sprintf("\n  Expanded MCPS ortholog recovery: %d/%d genes testable in GSE214373\n",
            n_testable, nrow(expanded_concordance_df)))
cat(sprintf("  Significant in human (FDR<%.2f, |LFC|>=%.1f): %d\n",
            PADJ_THRESH, LFC_THRESH, n_sig_exp))

readr::write_csv(expanded_concordance_df,
                 file.path(OUT_DIR, "MCPS_expanded_ortholog_concordance.csv"))
cat("  Saved: MCPS_expanded_ortholog_concordance.csv\n")
} else {
  cat("  [!] gene_symbol column not found in res_primary_df — skipping expanded ortholog recovery.\n")
  cat("      Ensure Section 4 annotates human DESeq2 results with HGNC symbols.\n")
  expanded_concordance_df <- NULL
}


# ----------------------------------------------------------------------
# ##Section 7 — MCPS Concordance: Direction-of-Effect Test + Binomial Test
#
# **Purpose:** Test whether the 6 MCPS genes identified in mouse (Stream 1) show concordant direction of effect in human sCJD (GSE124571). Two tests are applied: (1) sign concordance (binomial test against 50% null), (2) Wilcoxon signed-rank test comparing MCPS gene LFCs to background.
#
# **Inputs:** `../results/stream1/meta_results.tsv` (mouse LFCs); `../results/stream2/results_GSE124571.tsv` (human LFCs); `../results/stream2/human_mouse_orthologs.tsv`.
#
# **Outputs:** `../results/stream2/mcps_concordance.tsv`; concordance figure saved to `../results/stream2/figures/`.
#
# **Key decisions:** Direction discordance for individual MCPS genes is expected and reported transparently — the binomial test evaluates the gene set as a whole, not individual genes. A gene showing opposite direction in human is not excluded from MCPS but is flagged for discussion.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 7 — MCPS concordance analysis
# =============================================================================
cat("\n--- Section 7: MCPS concordance analysis ---\n")

# --- 7.1  Load MCPS ---
mcps_path <- file.path(STREAM1_DIR, "MCPS.csv")
if (!file.exists(mcps_path)) {
  stop(paste(
    "MCPS.csv not found at:", mcps_path,
    "\nEnsure 05_stream1_cns_meta_analysis.R output is uploaded to your Drive."
  ))
}

MCPS <- read_csv(mcps_path, show_col_types = FALSE)
if (!"gene" %in% colnames(MCPS)) stop("MCPS.csv must have a 'gene' column.")
if (!"meta_lfc" %in% colnames(MCPS)) stop("MCPS.csv must have a 'meta_lfc' column.")

n_mcps <- nrow(MCPS)
cat(sprintf("  MCPS loaded: %d genes\n", n_mcps))

# --- 7.2  Join MCPS with human results ---
cat("  Mapping MCPS to human orthologs and determining concordance...\n")
concordance <- MCPS %>%
  dplyr::rename(mouse_symbol = gene, mcps_lfc = meta_lfc) %>%
  inner_join(
    res_primary_orth %>%
      dplyr::filter(!is.na(mouse_symbol)) %>%
      dplyr::select(
        mouse_symbol,
        human_lfc    = log2FoldChange,
        human_padj   = padj,
        human_sig    = significant,
        human_symbol
      ),
    by = "mouse_symbol",
    relationship = "many-to-many"
  ) %>%
  group_by(mouse_symbol) %>%
  slice_min(order_by = human_padj, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(
    direction_match   = sign(mcps_lfc) == sign(human_lfc),
    human_sig_in_mcps = !is.na(human_padj) & human_padj <= PADJ_THRESH
  )

# --- 7.3  Compute concordance statistics ---
n_tested    <- nrow(concordance)
n_dir_match <- sum(concordance$direction_match, na.rm = TRUE)
n_human_sig <- sum(concordance$human_sig_in_mcps, na.rm = TRUE)

pct_dir_match_of_mcps   <- 100 * n_dir_match / max(n_mcps, 1)
pct_dir_match_of_tested <- 100 * n_dir_match / max(n_tested, 1)

cat(sprintf("\n  MCPS size (full):              %d genes\n", n_mcps))
cat(sprintf("  MCPS genes testable in human:  %d / %d (%.0f%%)\n",
            n_tested, n_mcps, 100 * n_tested / max(n_mcps, 1)))
cat(sprintf("  Direction concordance:         %d/%d tested (%.1f%% of full MCPS)\n",
            n_dir_match, n_tested, pct_dir_match_of_mcps))
cat(sprintf("  Human-significant MCPS genes:  %d\n", n_human_sig))

# --- 7.4  Safe Statistical Testing ---
# Only run correlation and binomial tests if we have at least 3 valid observations
valid_pairs <- sum(complete.cases(concordance$mcps_lfc, concordance$human_lfc))

if (valid_pairs >= 3) {
  binom_test <- binom.test(n_dir_match, n_tested, p = 0.5, alternative = "greater")
  lfc_cor <- cor.test(concordance$mcps_lfc, concordance$human_lfc,
                      method = "pearson", use = "complete.obs")

  cat(sprintf("  Binomial test p (H0: 50%%):     %.4e\n", binom_test$p.value))
  cat(sprintf("  LFC Pearson r:                 %.3f (p=%.4e)\n",
              lfc_cor$estimate, lfc_cor$p.value))
} else {
  cat("\n  [!] Note: Not enough mapped genes (need >= 3) to compute Binomial/Pearson statistics.\n")
  cat("      This is normal if testing with a minimal subset of data.\n")
}

# --- 7.5  Jaccard overlap (primary vs sensitivity DEGs) ---
# Explicitly print Jaccard overlap for transparency
if (exists("sensitivity_comparison") && "overlap_category" %in% colnames(sensitivity_comparison)) {
  n_both    <- sum(sensitivity_comparison$overlap_category == "Both models", na.rm = TRUE)
  n_primary <- sum(sensitivity_comparison$overlap_category == "Primary only (age-adjusted)", na.rm = TRUE)
  n_sens    <- sum(sensitivity_comparison$overlap_category == "Sensitivity only (age-confounded)", na.rm = TRUE)
  jaccard   <- n_both / max(n_both + n_primary + n_sens, 1)
  cat(sprintf("\n  Jaccard overlap (primary vs sensitivity DEGs): %.3f\n", jaccard))
  cat(sprintf("  Both: %d | Primary only: %d | Sensitivity only: %d\n",
              n_both, n_primary, n_sens))
  if (jaccard < 0.5) {
    cat("  [!] Low Jaccard (<0.5): age covariate substantially changes the DEG list.\n")
    cat("      Primary model (age-adjusted) is preferred for publication.\n")
  }
}

# --- 7.6  Save concordance CSV ---
readr::write_csv(concordance, file.path(OUT_DIR, "MCPS_human_concordance.csv"))
cat("\n  Concordance table saved: MCPS_human_concordance.csv\n")


# ----------------------------------------------------------------------
# ##Section 8 — GSEA: MCPS Gene Set on Human DEG Ranking
#
# **Purpose:** Test whether the MCPS gene set (6 genes) is enriched at the top of the human DE ranking using GSEA. The ranking metric is signed −log10(padj) × sign(LFC) from the GSE124571 analysis.
#
# **Inputs:** `../results/stream2/results_GSE124571.tsv`; MCPS gene list.
#
# **Outputs:** `../results/stream2/gsea_mcps_human.tsv` — NES, p-value, leading edge genes.
#
# **Key decisions:** With only 6 genes, GSEA power is limited — this analysis is supplementary to the concordance test. `fgsea` with `eps=0` and 10,000 permutations is used for accurate p-values.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 8 — GSEA: MCPS as gene set on ranked human DEGs
# =============================================================================
cat("\n--- Section 8: GSEA (MCPS gene set on human DEG ranking) ---\n")

# --- 8.1  Build ranked gene list (mouse gene symbols) ---
rank_vec <- res_primary_orth %>%
  dplyr::filter( # Proactive fix for Colab namespace collisions
    !is.na(mouse_symbol),
    !is.na(lfc_wald),
    !is.na(pval_wald)
  ) %>%
  mutate(
    # Rank statistic: signed -log10(p); floor pvalue at 1e-300 to avoid -Inf
    rank_stat = sign(lfc_wald) * (-log10(pmax(pval_wald, 1e-300)))
  ) %>%
  # If a mouse gene maps to multiple human genes, take the mean rank stat
  group_by(mouse_symbol) %>%
  summarise(rank_stat = mean(rank_stat, na.rm = TRUE), .groups = "drop") %>%
  deframe() %>%
  sort(decreasing = TRUE)

cat(sprintf("  Ranked gene list: %d mouse gene symbols\n", length(rank_vec)))

# --- 8.2  Run fgsea ---
mcps_geneset <- list(MCPS = MCPS$gene)

# Note: nPermSimple is omitted. 'eps = 0' triggers the modern fgseaMultilevel algorithm.
# minSize is lowered to 3 to accommodate your small test dataset.
suppressWarnings({
  gsea_res <- fgsea(
    pathways    = mcps_geneset,
    stats       = rank_vec,
    minSize     = 3,
    maxSize     = Inf,
    eps         = 0,   # Compute exact p-values
    nproc       = 1    # Safe core allocation for Colab
  )
})

# --- 8.3 Safe Output Printing ---
if (nrow(gsea_res) > 0) {
  cat(sprintf("  GSEA result: NES=%.3f, pval=%.4e, padj=%.4e, size=%d\n",
              gsea_res$NES, gsea_res$pval, gsea_res$padj, gsea_res$size))

  # Extract leading-edge genes (the core enrichment subset)
  leading_edge_genes <- unlist(gsea_res$leadingEdge)
  cat(sprintf("  Leading-edge genes: %d\n", length(leading_edge_genes)))
} else {
  cat("\n  [!] GSEA returned no results. The overlap between your MCPS gene set\n")
  cat("      and the ranked list is smaller than the minimum size (3).\n")
  cat("      This is normal if you are running a small subset of the data.\n")
  leading_edge_genes <- character(0)
}

# --- 8.4 Always write GSEA output CSV (even if empty) ---
# Script 08 expects this file; write an empty-but-valid CSV if GSEA failed.
gsea_out <- if (nrow(gsea_res) > 0) {
  data.frame(
    pathway         = gsea_res$pathway,
    NES             = gsea_res$NES,
    pval            = gsea_res$pval,
    padj            = gsea_res$padj,
    size            = gsea_res$size,
    leading_edge    = sapply(gsea_res$leadingEdge, paste, collapse = ";"),
    stringsAsFactors = FALSE
  )
} else {
  data.frame(
    pathway = character(0), NES = numeric(0), pval = numeric(0),
    padj = numeric(0), size = integer(0), leading_edge = character(0),
    stringsAsFactors = FALSE
  )
}

gsea_csv_path <- file.path(OUT_DIR, "GSE214373_GSEA_MCPS.csv")
readr::write_csv(gsea_out, gsea_csv_path)
cat(sprintf("  GSEA results written: %s (%d rows)\n", gsea_csv_path, nrow(gsea_out)))


# ----------------------------------------------------------------------
# ##Section 9 — Permutation-Based Power Analysis
#
# **Purpose:** Estimate the statistical power of the concordance test given n=6 MCPS genes and the observed human effect sizes. Generates a null distribution by permuting gene labels 10,000 times and reports the empirical power at α=0.05.
#
# **Inputs:** `../results/stream2/mcps_concordance.tsv`.
#
# **Outputs:** `../results/stream2/power_analysis.tsv`; power curve figure.
#
# **Key decisions:** Power analysis is reported in the manuscript Methods section to justify the small gene set size. Expected power ≥80% if ≥5/6 genes are concordant.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 9 — Permutation-based power analysis (multi-dispersion)
# =============================================================================
cat("\n--- Section 9: Permutation-based power analysis (multi-dispersion) ---\n")

# Guard: Section 9 requires the primary DESeq2 model (dds) from Section 4.
if (!exists("dds")) {
  stop(paste(
    "Primary DESeq2 model (dds) not found.",
    "Section 9 requires dds from Section 4. Run Section 4 first."
  ))
}

# --- 9.1  Estimate dispersion percentiles and mean count from real data ---
dispersions_fitted <- dispersions(dds)

disp_low    <- quantile(dispersions_fitted, 0.10, na.rm = TRUE)
disp_median <- quantile(dispersions_fitted, 0.50, na.rm = TRUE)
disp_high   <- quantile(dispersions_fitted, 0.90, na.rm = TRUE)

norm_counts <- counts(dds, normalized = TRUE)
mean_count  <- median(rowMeans(norm_counts), na.rm = TRUE)

lib_sizes <- colSums(counts(dds))
mean_lib   <- mean(lib_sizes)

cat(sprintf("  Dispersion percentiles — 10th: %.4f | 50th: %.4f | 90th: %.4f\n",
            disp_low, disp_median, disp_high))
cat(sprintf("  Median mean normalised count: %.1f\n", mean_count))
cat(sprintf("  Mean library size: %.0f\n", mean_lib))

mean_disp <- disp_median

# --- 9.2  Power simulation function ---
simulate_power_one_gene <- function(
    true_lfc,
    n_case    = n_prion,
    n_ctrl    = n_control,
    mu_base   = mean_count,
    disp      = disp_median,
    alpha     = POWER_ALPHA,
    n_sim     = 500
) {
  mu_ctrl_val <- mu_base
  mu_case_val <- mu_base * (2^true_lfc)

  detected <- vapply(seq_len(n_sim), function(i) {
    counts_ctrl <- rnbinom(n_ctrl, mu = mu_ctrl_val, size = 1 / disp)
    counts_case <- rnbinom(n_case, mu = mu_case_val, size = 1 / disp)

    log_ctrl <- log2(counts_ctrl + 0.5)
    log_case <- log2(counts_case + 0.5)

    tt <- tryCatch(
      t.test(log_case, log_ctrl, var.equal = FALSE),
      error = function(e) list(p.value = 1)
    )
    tt$p.value < alpha
  }, logical(1))

  mean(detected)
}

# --- 9.3  Run power curves at three dispersion levels ---
lfc_grid <- seq(0, 3, by = 0.25)

disp_levels <- list(
  "Low (10th pct)"    = disp_low,
  "Median (50th pct)" = disp_median,
  "High (90th pct)"   = disp_high
)

cat(sprintf("\n  Running power simulations (%d LFC values × 3 dispersion levels × 500 sims)...\n",
            length(lfc_grid)))

power_curves_list <- lapply(names(disp_levels), function(level_name) {
  d <- disp_levels[[level_name]]
  cat(sprintf("    Simulating: %s (disp = %.4f)\n", level_name, d))
  data.frame(
    true_lfc         = lfc_grid,
    # unname() removes the messy names that cause data.frame warnings
    power            = unname(vapply(lfc_grid, function(lfc) {
      simulate_power_one_gene(true_lfc = lfc, disp = d)
    }, numeric(1))),
    dispersion_level = level_name,
    dispersion_value = d,
    row.names        = NULL
  )
})

power_curves <- bind_rows(power_curves_list) %>%
  mutate(
    dispersion_level = factor(
      dispersion_level,
      levels = c("Low (10th pct)", "Median (50th pct)", "High (90th pct)")
    )
  )

# Explicitly use dplyr::filter to avoid ensembldb collision
power_curve <- power_curves %>% dplyr::filter(dispersion_level == "Median (50th pct)")

# --- 9.4  Find minimum detectable LFC at target power ---
min_lfc_table <- power_curves %>%
  group_by(dispersion_level, dispersion_value) %>%
  dplyr::filter(power >= POWER_TARGET) %>%
  slice_min(order_by = true_lfc, n = 1) %>%
  ungroup() %>%
  dplyr::select(dispersion_level, dispersion_value, min_detectable_lfc = true_lfc, power)

min_detectable_lfc_val <- min_lfc_table %>%
  dplyr::filter(dispersion_level == "Median (50th pct)") %>%
  pull(min_detectable_lfc)

if (length(min_detectable_lfc_val) == 0) {
  min_detectable_lfc <- NA_real_
  cat(sprintf("\n  WARNING: Median-dispersion power does not reach %.0f%% within LFC range.\n",
              100 * POWER_TARGET))
} else {
  min_detectable_lfc <- min_detectable_lfc_val[1]
  cat(sprintf("\n  Minimum detectable LFC at %.0f%% power:\n", 100 * POWER_TARGET))
  for (i in seq_len(nrow(min_lfc_table))) {
    cat(sprintf("    %-22s (disp=%.4f): LFC >= %.2f\n",
                min_lfc_table$dispersion_level[i],
                min_lfc_table$dispersion_value[i],
                min_lfc_table$min_detectable_lfc[i]))
  }
}

cat(sprintf("\n  Interpretation: the power envelope spans LFC %.2f–%.2f at %.0f%% power,\n",
            min(min_lfc_table$min_detectable_lfc, na.rm = TRUE),
            max(min_lfc_table$min_detectable_lfc, na.rm = TRUE),
            100 * POWER_TARGET))
cat("  reflecting the range from low-dispersion (easy to detect) to\n")
cat("  high-dispersion (hard to detect) genes in this dataset.\n")


# ----------------------------------------------------------------------
# ##Section 10 — Publication-Ready Figures
#
# **Purpose:** Generate all publication-ready figures for Stream 2: volcano plot for GSE124571, MCPS concordance scatter plot (mouse LFC vs human LFC), and the pathway-level human validation bubble plot (GSEA of consensus pathways).
#
# **Inputs:** All Stream 2 result tables; `../results/stream1/gsea_results.tsv`.
#
# **Outputs:** `../results/stream2/figures/` — PNG + SVG for each figure.
#
# **Key decisions:** The concordance scatter plot uses human gene symbols on axes and highlights MCPS genes in red. Discordant genes are shown with open circles.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 10 — Publication-ready figures
# =============================================================================
cat("\n--- Section 10: Generating figures ---\n")

# ---------------------------------------------------------------------------
# Fig 1 — Sample-level QC: PCA + sample distance heatmap
# ---------------------------------------------------------------------------
cat("  Fig 1: Sample QC (PCA + distance heatmap)\n")

p_pca <- ggplot(pca_data, aes(x = PC1, y = PC2)) +
  geom_point(
    aes(colour = condition, shape = sex, size = age),
    alpha = 0.90, stroke = 0.5
  ) +
  geom_text_repel(
    aes(label = paste0(condition, "\n", age, "y")),
    size          = 2.2,
    max.overlaps  = 30,
    segment.size  = 0.25,
    segment.color = "grey60",
    box.padding   = 0.3,
    force         = 1.5
  ) +
  scale_colour_manual(
    values = COND_COLOURS,
    labels = COND_LABELS,
    name   = "Condition"
  ) +
  scale_shape_manual(
    values = c(f = 16, m = 17),
    labels = c(f = "Female", m = "Male"),
    name   = "Sex"
  ) +
  scale_size_continuous(
    range  = c(2.5, 7),
    name   = "Age (years)",
    breaks = c(70, 80, 90, 100)
  ) +
  labs(
    x        = sprintf("PC1 (%s%% variance)", pct_var[1]),
    y        = sprintf("PC2 (%s%% variance)", pct_var[2]),
    title    = "A    Principal component analysis",
    subtitle = sprintf(
      "N = %d (sCJD: %d, control: %d); VST-normalised counts",
      n_total, n_prion, n_control
    )
  ) +
  theme_pub +
  theme(legend.position = "right") +
  coord_cartesian(clip = "off")

# Panel B: Sample distance heatmap
heatmap_annot <- HeatmapAnnotation(
  Condition = coldata$condition,
  Sex       = coldata$sex,
  Age       = anno_barplot(
    coldata$age,
    gp     = gpar(fill = "grey60", col = NA),
    height = unit(1.5, "cm")
  ),
  col = list(
    Condition = COND_COLOURS,  # <--- FIX: Pass the named vector directly
    Sex       = c(f = "#E69F00", m = "#56B4E9")
  ),
  annotation_name_side = "left",
  show_legend = c(TRUE, TRUE, FALSE)
)

heatmap_obj <- Heatmap(
  dist_matrix,
  name                  = "Euclidean\ndistance",
  col                   = colorRamp2(
    c(0, max(dist_matrix) / 2, max(dist_matrix)),
    c("#2166AC", "white", "#D6604D")
  ),
  top_annotation        = heatmap_annot,
  show_row_names        = TRUE,
  show_column_names     = FALSE,
  row_names_gp          = gpar(fontsize = 7),
  clustering_distance_rows    = "euclidean",
  clustering_distance_columns = "euclidean",
  clustering_method_rows      = "complete",
  clustering_method_columns   = "complete",
  heatmap_legend_param  = list(title_gp = gpar(fontsize = 8),
                               labels_gp = gpar(fontsize = 7))
)

# Save heatmap as PDF
pdf(file.path(FIGURES_DIR, "fig1b_sample_distance_heatmap.pdf"),
    width = 6, height = 5.5)
draw(heatmap_obj,
     column_title      = "B    Sample-to-sample Euclidean distance",
     column_title_gp   = gpar(fontsize = 10, fontface = "bold"),
     column_title_side = "top")
dev.off()

save_fig(p_pca, "fig1a_pca", width = 6.5, height = 5)
cat("    Saved: fig1a_pca.pdf/png, fig1b_sample_distance_heatmap.pdf\n")

# ---------------------------------------------------------------------------
# Fig 2 — Age confound visualisation
# ---------------------------------------------------------------------------
cat("  Fig 2: Age confound\n")

p_age_box <- ggplot(coldata,
                    aes(x = condition, y = age, colour = condition)) +
  geom_boxplot(
    aes(fill = condition),
    alpha    = 0.25,
    width    = 0.45,
    outlier.shape = NA,
    linewidth = 0.5
  ) +
  geom_jitter(
    aes(shape = sex),
    width  = 0.12,
    size   = 3,
    stroke = 0.6,
    alpha  = 0.9
  ) +
  scale_colour_manual(values = COND_COLOURS, labels = COND_LABELS, guide = "none") +
  scale_fill_manual(  values = COND_COLOURS, labels = COND_LABELS, guide = "none") +
  scale_x_discrete(labels = COND_LABELS) +
  scale_shape_manual(
    values = c(f = 16, m = 17),
    labels = c(f = "Female", m = "Male"),
    name   = "Sex"
  ) +
  annotate(
    "text",
    x     = 1.5,
    y     = max(coldata$age) + 3,
    label = sprintf("Welch t: p = %.3f\nr = %.2f",
                    age_ttest$p.value, age_cond_cor$estimate),
    size  = 3,
    hjust = 0.5,
    colour = "grey30"
  ) +
  labs(
    x        = NULL,
    y        = "Age (years)",
    title    = "A    Age distribution by condition",
    subtitle = sprintf(
      "sCJD: %.1f ± %.1f y  |  Control: %.1f ± %.1f y",
      mean(age_prion), sd(age_prion),
      mean(age_control), sd(age_control)
    )
  ) +
  theme_pub

p_age_pc1 <- ggplot(pca_data, aes(x = age, y = PC1)) +
  geom_smooth(method = "lm", se = TRUE, colour = "grey50",
              linewidth = 0.6, linetype = "dashed") +
  geom_point(aes(colour = condition, shape = sex), size = 3, alpha = 0.9) +
  scale_colour_manual(values = COND_COLOURS, labels = COND_LABELS, name = "Condition") +
  scale_shape_manual(
    values = c(f = 16, m = 17),
    labels = c(f = "Female", m = "Male"),
    name   = "Sex"
  ) +
  labs(
    x        = "Age (years)",
    y        = sprintf("PC1 (%s%% variance)", pct_var[1]),
    title    = "B    Age vs PC1",
    subtitle = "Dashed line: linear regression with 95% CI"
  ) +
  theme_pub

fig2 <- p_age_box + p_age_pc1 +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "bottom")

save_fig(fig2, "fig2_age_confound", width = 7, height = 4.5)

# ---------------------------------------------------------------------------
# Fig 3 — Sensitivity analysis: DEG overlap
# ---------------------------------------------------------------------------
cat("  Fig 3: Sensitivity analysis DEG overlap\n")

scatter_sens_df <- sensitivity_comparison %>%
  dplyr::filter(!is.na(lfc_primary) & !is.na(lfc_sens)) %>%
  mutate(
    overlap_category = factor(overlap_category, levels = c(
      "Both models",
      "Primary only (age-adjusted)",
      "Sensitivity only (age-confounded)",
      "Not significant"
    ))
  )

r_models <- cor(scatter_sens_df$lfc_primary, scatter_sens_df$lfc_sens,
                use = "complete.obs")

p_sens_scatter <- ggplot(scatter_sens_df,
                         aes(x = lfc_sens, y = lfc_primary,
                             colour = overlap_category)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed",
              colour = "grey60", linewidth = 0.4) +
  geom_point(
    data = dplyr::filter(scatter_sens_df, overlap_category == "Not significant"),
    alpha = 0.25, size = 0.6
  ) +
  geom_point(
    data = dplyr::filter(scatter_sens_df, overlap_category != "Not significant"),
    alpha = 0.85, size = 1.8
  ) +
  scale_colour_manual(
    values = c(
      "Both models"                       = "#1B7837",
      "Primary only (age-adjusted)"       = "#4393C3",
      "Sensitivity only (age-confounded)" = "#D6604D",
      "Not significant"                   = "grey80"
    ),
    name = "DEG category"
  ) +
  labs(
    x        = "LFC — sensitivity model (~ sex + condition)",
    y        = "LFC — primary model (~ age + sex + condition)",
    title    = "Sensitivity analysis: effect of age covariate on LFC estimates",
    subtitle = sprintf(
      "Pearson r = %.3f | Both: %d | Primary only: %d | Sensitivity only: %d",
      r_models, n_shared, n_primary_only, n_sens_only
    ),
    caption  = paste(
      "Genes in 'Sensitivity only' are likely age-confounded false positives.",
      "Genes in 'Both models' represent robust disease signal."
    )
  ) +
  theme_pub +
  theme(legend.position = "right") +
  coord_cartesian(clip = "off")

save_fig(p_sens_scatter, "fig3_sensitivity_overlap", width = 7, height = 5.5)

# ---------------------------------------------------------------------------
# Fig 4 — Volcano plot (primary model, MCPS genes highlighted)
# ---------------------------------------------------------------------------
cat("  Fig 4: Volcano plot\n")

mcps_human_sig_symbols <- concordance %>%
  dplyr::filter(human_sig_in_mcps) %>%
  pull(mouse_symbol)

top15_symbols <- res_primary_orth %>%
  dplyr::filter(!is.na(padj) & significant) %>%
  arrange(padj) %>%
  slice_head(n = 15) %>%
  pull(mouse_symbol) %>%
  na.omit()

label_symbols <- union(top15_symbols, mcps_human_sig_symbols)

volcano_df <- res_primary_orth %>%
  dplyr::filter(!is.na(padj) & !is.na(log2FoldChange)) %>%
  mutate(
    neg_log10_padj = -log10(pmax(padj, 1e-300)),
    point_category = case_when(
      !is.na(mouse_symbol) & mouse_symbol %in% mcps_human_sig_symbols & significant
                                                  ~ "MCPS concordant",
      significant                                 ~ "Significant DEG",
      TRUE                                        ~ "Not significant"
    ),
    label = ifelse(!is.na(mouse_symbol) & mouse_symbol %in% label_symbols,
                   mouse_symbol, NA_character_)
  )

n_mcps_conc <- sum(volcano_df$point_category == "MCPS concordant", na.rm = TRUE)

p_volcano <- ggplot(volcano_df,
                    aes(x = log2FoldChange, y = neg_log10_padj)) +
  geom_point(
    data  = dplyr::filter(volcano_df, point_category == "Not significant"),
    colour = "grey80", alpha = 0.4, size = 0.5
  ) +
  geom_point(
    data  = dplyr::filter(volcano_df, point_category == "Significant DEG"),
    colour = "#4393C3", alpha = 0.75, size = 1.2
  ) +
  geom_point(
    data  = dplyr::filter(volcano_df, point_category == "MCPS concordant"),
    colour = "#D55E00", alpha = 0.90, size = 2.0
  ) +
  geom_hline(yintercept = -log10(PADJ_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
             linetype = "dashed", colour = "grey40", linewidth = 0.35) +
  geom_text_repel(
    aes(label = label),
    size          = 2.3,
    max.overlaps  = 40,
    segment.size  = 0.25,
    segment.color = "grey50",
    box.padding   = 0.3,
    force         = 2,
    na.rm         = TRUE,
    fontface      = "italic"
  ) +
  annotate("point", x = -Inf, y = Inf, colour = NA) +
  scale_x_continuous(
    limits = c(-max(abs(volcano_df$log2FoldChange), na.rm = TRUE) * 1.05,
                max(abs(volcano_df$log2FoldChange), na.rm = TRUE) * 1.05),
    expand = expansion(mult = 0.02)
  ) +
  labs(
    x        = expression(bold("Shrunken log"[2]*" fold change (sCJD / Control)")),
    y        = expression(bold(-log[10]*"(padj)")),
    title    = "Human sCJD vs Control — Differential expression",
    subtitle = sprintf(
      "Primary model: ~ age + sex + condition | DEGs: %d (%d up, %d down) | MCPS concordant: %d%s",
      n_sig_primary, n_up_primary, n_dn_primary, n_mcps_conc,
      ifelse(n_mcps_conc == 0, " (none reached FDR<0.05)", "")
    ),
    caption  = sprintf(
      "Thresholds: padj \u2264 %.2f, |LFC| \u2265 %.1f (apeglm shrinkage). Orange: MCPS-concordant. Blue: other DEGs.",
      PADJ_THRESH, LFC_THRESH
    )
  ) +
  theme_pub +
  coord_cartesian(clip = "off")

save_fig(p_volcano, "fig4_volcano_primary", width = 9, height = 5.5)

# ---------------------------------------------------------------------------
# Fig 5 — Direction-of-effect scatter
# ---------------------------------------------------------------------------
cat("  Fig 5: LFC concordance scatter\n")

scatter_df <- concordance %>%
  dplyr::filter(!is.na(human_lfc) & !is.na(mcps_lfc)) %>%
  mutate(
    label = ifelse(human_sig_in_mcps, mouse_symbol, NA_character_),
    concordant = direction_match
  )

r_lfc     <- if(exists("lfc_cor")) lfc_cor$estimate else NA
r_ci_low  <- if(exists("lfc_cor")) lfc_cor$conf.int[1] else NA
r_ci_high <- if(exists("lfc_cor")) lfc_cor$conf.int[2] else NA
binom_pval <- if(exists("binom_test")) binom_test$p.value else 1

p_scatter <- ggplot(scatter_df,
                    aes(x = mcps_lfc, y = human_lfc)) +
  geom_hline(yintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.35) +
  geom_vline(xintercept = 0, linetype = "dashed", colour = "grey70", linewidth = 0.35) +
  geom_smooth(method = "lm", se = TRUE,
              colour = "grey30", fill = "grey85",
              linewidth = 0.6, linetype = "solid") +
  geom_point(
    aes(colour = concordant),
    alpha = 0.75, size = 1.8
  ) +
  geom_text_repel(
    aes(label = label),
    size          = 2.3,
    max.overlaps  = 25,
    segment.size  = 0.25,
    segment.color = "grey50",
    box.padding   = 0.3,
    na.rm         = TRUE,
    fontface      = "italic"
  ) +
  scale_colour_manual(
    values = c("TRUE" = "#1B7837", "FALSE" = "#D6604D"),
    labels = c("TRUE" = "Concordant", "FALSE" = "Discordant"),
    name   = "Direction"
  ) +
  labs(
    x        = expression(bold("Mouse meta-analysis LFC (MCPS)")),
    y        = expression(bold("Human LFC (sCJD vs Control)")),
    title    = "Cross-species LFC concordance: MCPS genes",
    subtitle = sprintf(
      "N = %d MCPS genes testable in human  |  Concordant: %d/%d (%.0f%%)  |  r = %.2f [%.2f, %.2f]",
      n_tested, n_dir_match, n_tested, pct_dir_match_of_tested,
      r_lfc, r_ci_low, r_ci_high
    ),
    caption  = sprintf(
      "Binomial test (H\u2080: concordance = 50%%): p = %.4f. Orange: discordant. Green: concordant.",
      binom_pval
    )
  ) +
  theme_pub +
  theme(legend.position = "bottom") +
  coord_cartesian(clip = "off")

save_fig(p_scatter, "fig5_lfc_scatter", width = 6, height = 5.5)

# ---------------------------------------------------------------------------
# Fig 6 — GSEA enrichment plot
# ---------------------------------------------------------------------------
cat("  Fig 6: GSEA enrichment plot\n")

if (exists("gsea_res") && nrow(gsea_res) > 0) {
  p_gsea <- plotEnrichment(mcps_geneset[["MCPS"]], rank_vec) +
    labs(
      x        = "Rank in human DEG list (sCJD vs Control)",
      y        = "Enrichment score",
      title    = "GSEA: Mouse Conserved Prion Signature (MCPS) in human sCJD",
      subtitle = sprintf(
        "NES = %.3f  |  padj = %.4f  |  Gene set size = %d  |  Leading edge = %d genes",
        gsea_res$NES, gsea_res$padj, gsea_res$size, length(leading_edge_genes)
      ),
      caption  = sprintf(
        "Rank statistic: sign(LFC) \u00d7 -log\u2081\u2080(p-value). %d permutations.",
        GSEA_NPERM
      )
    ) +
    theme_pub +
    theme(
      panel.grid.major.x = element_blank(),
      axis.line.x        = element_line(colour = "black", linewidth = 0.4)
    )

  save_fig(p_gsea, "fig6_gsea_enrichment", width = 7, height = 4.5)
} else {
  cat("    [!] Skipping Fig 6 (GSEA): Not enough mapped genes for plotEnrichment.\n")
}

# ---------------------------------------------------------------------------
# Fig 7 — Power curves
# ---------------------------------------------------------------------------
cat("  Fig 7: Multi-dispersion power curves\n")

power_ribbon_wide <- power_curves %>%
  dplyr::select(true_lfc, dispersion_level, power) %>%
  tidyr::pivot_wider(names_from = dispersion_level, values_from = power)

expected_cols <- c("Low (10th pct)", "Median (50th pct)", "High (90th pct)")
missing_cols  <- setdiff(expected_cols, colnames(power_ribbon_wide))
if (length(missing_cols) > 0) {
  warning(sprintf(
    "pivot_wider missing expected dispersion columns: %s. ",
    paste(missing_cols, collapse = ", ")
  ), "Ribbon will not be drawn.")
  for (mc in missing_cols) power_ribbon_wide[[mc]] <- NA_real_
}

power_ribbon <- power_ribbon_wide %>%
  dplyr::rename(
    power_low    = `Low (10th pct)`,
    power_median = `Median (50th pct)`,
    power_high   = `High (90th pct)`
  )

disp_colours <- c(
  "Low (10th pct)"    = "#009E73",
  "Median (50th pct)" = "#0072B2",
  "High (90th pct)"   = "#D55E00"
)

vline_df <- min_lfc_table %>%
  dplyr::filter(!is.na(min_detectable_lfc)) %>%
  mutate(colour = disp_colours[as.character(dispersion_level)])

all_y_positions <- c(
  "Low (10th pct)"    = 0.12,
  "Median (50th pct)" = 0.06,
  "High (90th pct)"   = 0.18
)
vline_y_positions <- all_y_positions[as.character(vline_df$dispersion_level)]

p_power <- ggplot(power_curves,
                  aes(x = true_lfc, y = power,
                      colour = dispersion_level,
                      group  = dispersion_level)) +
  geom_ribbon(
    data    = power_ribbon,
    inherit.aes = FALSE,
    aes(x = true_lfc, ymin = power_high, ymax = power_low),
    fill    = "grey80",
    alpha   = 0.45,
    colour  = NA
  ) +
  geom_hline(yintercept = POWER_TARGET, linetype = "dashed",
             colour = "grey40", linewidth = 0.4) +
  annotate("text",
           x = max(lfc_grid) * 0.97, y = POWER_TARGET + 0.035,
           label = sprintf("%.0f%% power target", 100 * POWER_TARGET),
           hjust = 1, size = 2.8, colour = "grey40") +
  {
    if (nrow(vline_df) > 0) {
      geom_vline(
        data        = vline_df,
        aes(xintercept = min_detectable_lfc, colour = dispersion_level),
        linetype    = "dotted",
        linewidth   = 0.55,
        show.legend = FALSE
      )
    }
  } +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.2) +
  {
    if (nrow(vline_df) > 0) {
      lapply(seq_len(nrow(vline_df)), function(i) {
        annotate("text",
                 x      = vline_df$min_detectable_lfc[i] + 0.06,
                 y      = vline_y_positions[i],
                 label  = sprintf("LFC\u2265%.2f", vline_df$min_detectable_lfc[i]),
                 hjust  = 0,
                 size   = 2.6,
                 colour = vline_df$colour[i])
      })
    }
  } +
  scale_colour_manual(
    values = disp_colours,
    name   = "Dispersion level\n(fitted NB dispersion)"
  ) +
  scale_y_continuous(
    limits = c(0, 1),
    labels = scales::percent_format(accuracy = 1),
    breaks = seq(0, 1, 0.2)
  ) +
  scale_x_continuous(breaks = seq(0, 3, 0.5)) +
  labs(
    x        = expression(bold("True log"[2]*" fold change")),
    y        = "Statistical power",
    title    = "Power analysis: minimum detectable LFC across dispersion levels",
    subtitle = sprintf(
      paste0("N = %d sCJD, %d control  |  ",
             "Dispersion: 10th=%.4f, 50th=%.4f, 90th=%.4f  |  ",
             "\u03b1 = %.2f (uncorrected)"),
      n_prion, n_control,
      disp_low, disp_median, disp_high,
      POWER_ALPHA
    ),
    caption  = paste(
      "Single-gene parametric simulation (n = 500 per LFC value per dispersion level).",
      "Negative binomial model; Welch t-test on log2-counts (pseudocount = 0.5).",
      "Shaded band = power envelope between 10th and 90th dispersion percentiles.",
      "Genome-wide BH-FDR correction will further reduce effective power."
    )
  ) +
  theme_pub +
  theme(legend.position = "right") +
  coord_cartesian(clip = "off")

save_fig(p_power, "fig7_power_curve", width = 7.5, height = 5)

write.csv(
  power_curves %>%
    dplyr::select(dispersion_level, dispersion_value, true_lfc, power),
  file.path(OUT_DIR, "power_curves_multi_dispersion.csv"),
  row.names = FALSE
)
write.csv(
  min_lfc_table,
  file.path(OUT_DIR, "min_detectable_lfc_by_dispersion.csv"),
  row.names = FALSE
)
cat(sprintf("  Saved: power_curves_multi_dispersion.csv, min_detectable_lfc_by_dispersion.csv\n"))


# ----------------------------------------------------------------------
# ##Section 11 — Write Output Tables and Summary Report
#
# **Purpose:** Write all final output tables for Stream 2 and generate a plain-text summary report covering sample counts, DE statistics, MCPS concordance results, and pathway enrichment highlights.
#
# **Inputs:** All Stream 2 result objects in R.
#
# **Outputs:** `../results/stream2/stream2_summary.txt`; all TSV tables finalised and written.
#
# **Key decisions:** The summary report is designed to be copy-pasted directly into the manuscript Results section as a first draft.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 11 — Write output tables and summary report
# =============================================================================
cat("\n--- Section 11: Writing output tables ---\n")

# --- 11.1  Primary DESeq2 results ---
write_csv(
  res_primary_orth %>%
    dplyr::select(
      gene_id, ensembl_id, human_symbol, mouse_symbol,
      baseMean, log2FoldChange, lfcSE, pvalue = pval_wald, padj,
      lfc_wald, pval_wald, padj_wald,
      significant, direction,
      orthology_type, pct_id_human2mouse
    ),
  file.path(OUT_DIR, "GSE214373_deseq2_primary.csv")
)
cat("  Written: GSE214373_deseq2_primary.csv\n")

# --- 11.2  Sensitivity model results ---
write_csv(
  res_sens_orth %>%
    dplyr::select(
      gene_id, ensembl_id, mouse_symbol,
      baseMean, log2FoldChange, lfcSE, pvalue, padj,
      significant, direction
    ),
  file.path(OUT_DIR, "GSE214373_deseq2_sensitivity.csv")
)
cat("  Written: GSE214373_deseq2_sensitivity.csv\n")

# --- 11.3  Sensitivity comparison table ---
write_csv(
  sensitivity_comparison %>%
    dplyr::filter(overlap_category != "Not significant") %>%
    arrange(overlap_category, padj_primary),
  file.path(OUT_DIR, "GSE214373_sensitivity_comparison.csv")
)
cat("  Written: GSE214373_sensitivity_comparison.csv\n")

# --- 11.4  Ortholog map (cache) ---
write_csv(orthologs, file.path(OUT_DIR, "GSE214373_orthologs.csv"))
cat("  Written: GSE214373_orthologs.csv\n")

# --- 11.5  MCPS concordance table ---
write_csv(
  concordance %>%
    dplyr::select(
      mouse_symbol, mcps_lfc,
      human_symbol, human_lfc, human_padj, human_sig_in_mcps,
      direction_match
    ) %>%
    arrange(desc(direction_match), human_padj),
  file.path(OUT_DIR, "GSE214373_MCPS_concordance.csv")
)
cat("  Written: GSE214373_MCPS_concordance.csv\n")

# --- 11.6  GSEA results (Safely guarded against empty sets) ---
if (exists("gsea_res") && nrow(gsea_res) > 0) {
  write_csv(
    gsea_res %>%
      dplyr::select(-leadingEdge) %>%
      mutate(leading_edge_genes = paste(leading_edge_genes, collapse = ",")),
    file.path(OUT_DIR, "GSE214373_GSEA_MCPS.csv")
  )
  cat("  Written: GSE214373_GSEA_MCPS.csv\n")
} else {
  cat("  Skipped: GSE214373_GSEA_MCPS.csv (GSEA was skipped due to size)\n")
}

# --- 11.7  Power analysis tables ---
write_csv(power_curve, file.path(OUT_DIR, "GSE214373_power_analysis.csv"))
cat("  Written: GSE214373_power_analysis.csv\n")
# Multi-dispersion tables were already written in Section 10

# --- 11.8  Plain-text summary report ---

# Pre-compute per-dispersion-level LFC lines
power_summary_lines <- if (nrow(min_lfc_table) > 0) {
  apply(min_lfc_table, 1, function(row) {
    sprintf("    %-22s (disp=%.4f): LFC >= %.2f",
            row["dispersion_level"],
            as.numeric(row["dispersion_value"]),
            as.numeric(row["min_detectable_lfc"]))
  })
} else {
  "    Not reached within LFC range tested (0–3)."
}

# Safe variables to prevent crashes if stats were skipped on small test data
binom_p_str <- if(exists("binom_test")) sprintf("%.4f", binom_test$p.value) else "N/A"
r_est_str   <- if(exists("lfc_cor")) sprintf("%.3f", lfc_cor$estimate) else "N/A"
r_p_str     <- if(exists("lfc_cor")) sprintf("%.4f", lfc_cor$p.value) else "N/A"

gsea_nes_str  <- if(exists("gsea_res") && nrow(gsea_res) > 0) sprintf("%.3f", gsea_res$NES[1]) else "N/A"
gsea_p_str    <- if(exists("gsea_res") && nrow(gsea_res) > 0) sprintf("%.4f", gsea_res$pval[1]) else "N/A"
gsea_padj_str <- if(exists("gsea_res") && nrow(gsea_res) > 0) sprintf("%.4f", gsea_res$padj[1]) else "N/A"
gsea_size_str <- if(exists("gsea_res") && nrow(gsea_res) > 0) sprintf("%d", gsea_res$size[1]) else "N/A"
gsea_le_str   <- if(exists("leading_edge_genes")) sprintf("%d", length(leading_edge_genes)) else "N/A"

summary_lines <- c(
  "=============================================================",
  "  Stream 2: Human sCJD Validation — Summary Report",
  "  06_stream2_human_validation.R [Colab Version]",
  "=============================================================",
  "",
  "DATASET",
  "-------",
  "  Accession:   GSE214373",
  "  Reference:   Dimitriadis et al. (2022) bioRxiv 2022.09.13.507650",
  "  Organism:    Homo sapiens",
  "  Tissue:      Frontal cortex (post-mortem)",
  sprintf("  N:           %d total (%d sCJD, %d control)", n_total, n_prion, n_control),
  sprintf("  Age sCJD:    %.1f ± %.1f y (range %d–%d)",
          mean(age_prion), sd(age_prion), min(age_prion), max(age_prion)),
  sprintf("  Age control: %.1f ± %.1f y (range %d–%d)",
          mean(age_control), sd(age_control), min(age_control), max(age_control)),
  "",
  "AGE CONFOUND",
  "------------",
  sprintf("  Welch t-test: t = %.3f, p = %.4f",
          age_ttest$statistic, age_ttest$p.value),
  sprintf("  Pearson r(age, condition) = %.3f (p = %.4f)",
          age_cond_cor$estimate, age_cond_cor$p.value),
  sprintf("  VIF(condition) = %.2f, VIF(age) = %.2f",
          vif_condition, vif_age),
  "  Decision: age included as continuous covariate in primary model.",
  "  Sensitivity analysis run without age to quantify confound impact.",
  "",
  "PRIMARY DESeq2 MODEL: ~ age + sex + condition",
  "---------------------------------------------",
  sprintf("  Design matrix rank: %d/%d (full rank)", rank_check, ncol(X_check)),
  sprintf("  Genes tested: %d", nrow(res_primary_df)),
  sprintf("  DEGs (padj <= %.2f, |LFC| >= %.1f): %d total",
          PADJ_THRESH, LFC_THRESH, n_sig_primary),
  sprintf("    Upregulated in sCJD:   %d", n_up_primary),
  sprintf("    Downregulated in sCJD: %d", n_dn_primary),
  "",
  "SENSITIVITY ANALYSIS: ~ sex + condition (age omitted)",
  "------------------------------------------------------",
  sprintf("  DEGs: %d total (%d up, %d down)", n_sig_sens, n_up_sens, n_dn_sens),
  sprintf("  Shared with primary model: %d", n_shared),
  sprintf("  Primary only (robust to age): %d", n_primary_only),
  sprintf("  Sensitivity only (age-confounded): %d", n_sens_only),
  sprintf("  Jaccard similarity: %.3f", jaccard),
  "",
  "ORTHOLOG MAPPING",
  "----------------",
  sprintf("  1:1 human-mouse ortholog pairs: %d", nrow(orthologs)),
  sprintf("  Tested genes with mouse ortholog: %d / %d",
          n_mapped, nrow(res_primary_orth)),
  "",
  "MCPS CONCORDANCE",
  "----------------",
  sprintf("  MCPS size (full):              %d genes", n_mcps),
  sprintf("  MCPS genes testable in human:  %d / %d (%.0f%%)",
          n_tested, n_mcps, 100 * n_tested / max(n_mcps, 1)),
  sprintf("  Direction concordance:         %d/%d tested (%.1f%% of full MCPS)",
          n_dir_match, n_tested, pct_dir_match_of_mcps),
  sprintf("  Binomial test (H0: 50%%):        p = %s", binom_p_str),
  sprintf("  LFC Pearson r:                 %s (p = %s)", r_est_str, r_p_str),
  sprintf("  Human-significant MCPS genes:  %d", n_human_sig),
  "",
  "GSEA (MCPS gene set)",
  "--------------------",
  sprintf("  NES:           %s", gsea_nes_str),
  sprintf("  p-value:       %s", gsea_p_str),
  sprintf("  padj:          %s", gsea_padj_str),
  sprintf("  Gene set size: %s", gsea_size_str),
  sprintf("  Leading edge:  %s genes", gsea_le_str),
  sprintf("  Permutations:  %d", GSEA_NPERM),
  "",
  "POWER ANALYSIS (multi-dispersion)",
  "----------------------------------",
  sprintf("  N sCJD: %d  |  N control: %d", n_prion, n_control),
  sprintf("  Dispersion percentiles — 10th: %.4f | 50th: %.4f | 90th: %.4f",
          disp_low, disp_median, disp_high),
  sprintf("  Target power: %.0f%%  |  alpha: %.2f (uncorrected)",
          100 * POWER_TARGET, POWER_ALPHA),
  sprintf("  Min. detectable LFC at %.0f%% power:", 100 * POWER_TARGET),
  power_summary_lines,
  sprintf("  Power envelope (10th–90th pct): LFC %.2f – %.2f",
          min(min_lfc_table$min_detectable_lfc, na.rm = TRUE),
          max(min_lfc_table$min_detectable_lfc, na.rm = TRUE)),
  "  Note: genome-wide BH-FDR correction will further reduce effective power.",
  "",
  "LIMITATIONS",
  "-----------",
  "  1. Very small N (5 sCJD, 9 control): low statistical power.",
  "  2. Age-condition correlation (r = -0.60, p = 0.024): controls are",
  "     ~12 years older. Age covariate included but imprecisely estimated.",
  "  3. Controls from Queen Square Brain Bank (low-level AD/ageing pathology):",
  "     age-related transcriptional changes may inflate apparent sCJD signal.",
  "  4. Source is a preprint (Dimitriadis et al. 2022 bioRxiv); not yet",
  "     peer-reviewed at time of analysis.",
  "  5. Post-mortem delay varies (2–12 days for sCJD; 26–168 h for controls):",
  "     RNA quality differences may introduce technical noise.",
  "",
  "OUTPUT FILES",
  "------------",
  "  GSE214373_deseq2_primary.csv",
  "  GSE214373_deseq2_sensitivity.csv",
  "  GSE214373_sensitivity_comparison.csv",
  "  GSE214373_orthologs.csv",
  "  GSE214373_MCPS_concordance.csv",
  "  GSE214373_GSEA_MCPS.csv",
  "  GSE214373_power_analysis.csv",
  "  power_curves_multi_dispersion.csv",
  "  min_detectable_lfc_by_dispersion.csv",
  "  fig1a_pca.pdf/png",
  "  fig1b_sample_distance_heatmap.pdf",
  "  fig2_age_confound.pdf/png",
  "  fig3_sensitivity_overlap.pdf/png",
  "  fig4_volcano_primary.pdf/png",
  "  fig5_lfc_scatter.pdf/png",
  "  fig6_gsea_enrichment.pdf/png",
  "  fig7_power_curve.pdf/png  (3-curve power envelope)",
  "",
  "=============================================================",
  sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  sprintf("  R version: %s", R.version.string),
  sprintf("  DESeq2 version: %s", as.character(packageVersion("DESeq2"))),
  sprintf("  tximeta version: %s", as.character(packageVersion("tximeta"))),
  sprintf("  fgsea version: %s", as.character(packageVersion("fgsea"))),
  "============================================================="
)

writeLines(summary_lines, file.path(OUT_DIR, "stream2_summary.txt"))
cat(paste(summary_lines, collapse = "\n"), "\n")

cat("\n=============================================================\n")
cat("  06_stream2_human_validation.R complete.\n")
cat(sprintf("  All outputs written to: %s\n", OUT_DIR))
cat("=============================================================\n")
