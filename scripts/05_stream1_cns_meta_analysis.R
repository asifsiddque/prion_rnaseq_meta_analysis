#!/usr/bin/env Rscript
# ==============================================================================
# 05_stream1_cns_meta_analysis.R
# Stream 1: mouse CNS per-study DESeq2 + DerSimonian-Laird random-effects meta-analysis; MCPS derivation; LOO sensitivity; replication; ORA/GSEA
# ==============================================================================
# Project : Cross-tissue bulk RNA-seq meta-analysis of prion disease
#           (mouse CNS Stratum A/D, peripheral time course, human sCJD validation)
#
# Provenance: Converted 1:1 from the executed Google Colab notebook
#             "05_stream1_cns_meta_analysis_v2.ipynb" (the version that produced the manuscript results).
#             Only Colab environment-setup cells (drive mount, rpy2 magic,
#             package installers) were removed or commented; all analysis
#             code is unchanged.
#
# Inputs  : Salmon quantifications (results/salmon_quant/), config/coldata/
# Outputs : results/stream1/, results/meta_analysis/, reports/figures/stream1/
#
# Usage   : Set BASE_DIR (top of script) to the project root, then run:
#             Rscript 05_stream1_cns_meta_analysis.R
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
# ##Section 0 — Package Installation & Loading
#
# **Purpose:** Install and load all R packages required for the CNS meta-analysis: tximeta/DESeq2 for differential expression, apeglm for LFC shrinkage, biomaRt for gene annotation, ComplexHeatmap for visualisation, and metafor for random-effects meta-analysis.
#
# **Inputs:** Internet connection (first run only).
#
# **Outputs:** All packages loaded into the R session.
#
# **Key decisions:** BiocParallel is configured to use `MulticoreParam` on Linux (Colab) and `SerialParam` on Windows — detected automatically.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 0: PACKAGE INSTALLATION AND LOADING
# =============================================================================

# Increase timeouts for Colab downloads (crucial for large Bioconductor packages)
options(timeout = 600)

# Install Bioconductor packages if not already installed
if (!requireNamespace("BiocManager", quietly = TRUE))
  NULL  # [setup] installer removed (see environment.yml / renv.lock):
  # [setup]   install.packages("BiocManager")

bioc_pkgs <- c(
  "tximeta",          # Salmon import with transcript-length offsets
  "tximport",         # Fallback import if tximeta fails
  "DESeq2",           # Differential expression
  "apeglm",           # Log2FC shrinkage (apeglm method)
  "BiocParallel",     # Parallel backend control
  "SummarizedExperiment",
  "clusterProfiler",  # ORA pathway enrichment
  "fgsea",            # GSEA
  "org.Mm.eg.db",     # Mouse gene annotation (Entrez IDs)
  "AnnotationDbi",    # Gene ID mapping
  "enrichplot",       # Enrichment visualizations
  "msigdbr",          # MSigDB gene sets (Hallmark, KEGG, Reactome)
  "biomaRt",          # Added biomaRt based on your Colab fallback requirements
  "ComplexHeatmap"    # Added for Q1 quality heatmaps
)

cran_pkgs <- c(
  "metafor",          # Random-effects meta-analysis (DerSimonian-Laird)
  "ggplot2",          # Plotting
  "ggrepel",          # Non-overlapping labels on plots
  "dplyr",            # Data manipulation
  "tidyr",            # Data reshaping
  "readr",            # Fast CSV reading/writing
  "purrr",            # Functional programming utilities
  "stringr",          # String manipulation
  "scales",           # Plot axis formatting
  "RColorBrewer",     # Color palettes
  "circlize",         # Color mapping for ComplexHeatmap
  "Matrix"            # Rank check for GSE201249 design matrix
)

# Install packages without asking for updates to prevent Colab hanging
NULL  # [setup] installer removed (see environment.yml / renv.lock):
# [setup]   BiocManager::install(bioc_pkgs, ask = FALSE, update = FALSE)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    NULL  # [setup] installer removed (see environment.yml / renv.lock):
    # [setup]   install.packages(pkg)
}

# Load all packages
suppressPackageStartupMessages({
  library(tximeta)
  library(tximport)
  library(DESeq2)
  library(apeglm)
  library(BiocParallel)
  library(SummarizedExperiment)
  library(clusterProfiler)
  library(fgsea)
  library(org.Mm.eg.db)
  library(AnnotationDbi)
  library(enrichplot)
  library(msigdbr)
  library(metafor)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(purrr)
  library(stringr)
  library(scales)
  library(RColorBrewer)
  library(ComplexHeatmap)
  library(circlize)
  library(Matrix)
  library(biomaRt)
})

# Use single-core BiocParallel to avoid memory duplication from parallel workers
register(SerialParam())

cat("All packages loaded successfully.\n")


# ----------------------------------------------------------------------
# ##Section 1 — Global Configuration
#
# **Purpose:** Define all paths, analysis parameters, and the `save_fig()` helper used throughout this notebook. Also defines `theme_pub` — the shared ggplot2 theme applied to all publication figures.
#
# **Inputs:** None (constants defined inline).
#
# **Outputs:** `BASE_DIR`, `COLDATA_DIR`, `STREAM1_DIR`, `OUT_DIR`, `MCPS_GENES`, `theme_pub`, `save_fig()` available in R.
#
# **Key decisions:** `MCPS_GENES` lists the 6 multi-criteria prioritised signature genes (Tasl/Slc15a4, Igsf6, Fcgr4, Parp14, Iigp1, Il1a) that are the primary outcome of this entire pipeline. `FDR_THRESH = 0.05`, `LFC_THRESH = 0.5`.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 1: GLOBAL CONFIGURATION
# =============================================================================

# --- Directory paths ---
# Define the root of your project directory in Google Drive.
# IMPORTANT: Verify that "prion_rnaseq_meta_analysis" perfectly matches your folder name in Drive.
BASE_DIR <- "/content/drive/MyDrive/prion_rnaseq_meta_analysis"

# Set working directory to scripts to mimic local execution
setwd(file.path(BASE_DIR, "scripts"))

# Use absolute paths linked to the Colab Drive mount
SALMON_DIR   <- file.path(BASE_DIR, "results/salmon_quant") # <GSE>/<SRR>/quant.sf
METADATA_DIR <- file.path(BASE_DIR, "config/coldata")       # <GSE>_coldata.csv
DESEQ2_DIR   <- file.path(BASE_DIR, "results/deseq2")
META_DIR     <- file.path(BASE_DIR, "results/meta_analysis")
REPLIC_DIR   <- file.path(BASE_DIR, "results/stream1/replication")
ENRICH_DIR   <- file.path(BASE_DIR, "results/stream1/enrichment")
FIGURES_DIR  <- file.path(BASE_DIR, "reports/figures/stream1")

for (d in c(DESEQ2_DIR, META_DIR, REPLIC_DIR, ENRICH_DIR, FIGURES_DIR)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# --- Significance thresholds ---
PADJ_THRESH  <- 0.05   # Benjamini-Hochberg FDR threshold (inclusive)
LFC_THRESH   <- 0.5    # |log2FoldChange| threshold (inclusive)
I2_THRESH    <- 75     # I² threshold: genes with I² > 75% excluded from MCPS
MIN_STUDIES  <- 3      # Minimum studies a gene must appear in for meta-analysis

# --- Study registries ---

# Stratum A: CNS Discovery (5 studies — primary meta-analysis)
STRATUM_A <- c("GSE149805", "GSE198063", "GSE144738", "GSE189527", "GSE201249")

# Stratum D: Held-Out Replication (4 studies)
STRATUM_D <- c("GSE90977", "GSE124932", "GSE207251", "GSE202275")

# --- DESeq2 design formulas per study ---
# NOTE: GSE149805 is handled specially (stratified by cell_type)
STUDY_DESIGNS <- list(
  GSE149805 = "~ time_point + condition",   # applied per cell_type stratum
  GSE198063 = "~ condition",                # check sex covariate in coldata
  GSE144738 = "~ time_point + layout + condition",
  GSE189527 = "~ time_point + condition",
  GSE201249 = "~ brain_region + cell_type + condition",
  # Stratum D
  GSE90977  = "~ time_point + condition",
  GSE124932 = "~ condition",                # WT arms only
  GSE207251 = "~ condition",
  GSE202275 = "~ condition"                 # prion+vehicle vs control+vehicle only
)

# --- LFC shrinkage coefficient names (DESeq2 resultsNames convention) ---
# These are the contrast coefficient names produced by DESeq2 for each study.
# Adjust if your coldata uses different factor level names.
# Convention: condition_<treatment>_vs_<reference>
# Reference level should always be the control/mock/uninfected group.
STUDY_COEF <- list(
  GSE149805 = "condition_prion_vs_control",
  GSE198063 = "condition_prion_vs_control",
  GSE144738 = "condition_prion_vs_control",
  GSE189527 = "condition_prion_vs_control",
  GSE201249 = "condition_prion_vs_control",
  GSE90977  = "condition_prion_vs_control",
  GSE124932 = "condition_prion_vs_control",
  GSE207251 = "condition_prion_vs_control",
  GSE202275 = "condition_prion_vs_control"
)

cat("Configuration loaded.\n")
cat(sprintf("  Stratum A studies : %s\n", paste(STRATUM_A, collapse = ", ")))
cat(sprintf("  Stratum D studies : %s\n", paste(STRATUM_D, collapse = ", ")))
cat(sprintf("  PADJ threshold    : %.2f\n", PADJ_THRESH))
cat(sprintf("  |LFC| threshold   : %.1f\n", LFC_THRESH))
cat(sprintf("  I2 exclusion      : > %d%%\n", I2_THRESH))


# ----------------------------------------------------------------------
# ##Section 2 — Helper Functions
#
# **Purpose:** Define reusable functions for the meta-analysis: `run_deseq2_study()` runs a full DESeq2 pipeline for one study and returns shrunk LFC results; `extract_meta_input()` formats per-study results for metafor; `run_re_meta()` fits a DerSimonian-Laird random-effects model.
#
# **Inputs:** None (function definitions only).
#
# **Outputs:** Three functions available in the R environment.
#
# **Key decisions:** LFC shrinkage uses `apeglm` (adaptive t-prior), which is more accurate than the original `normal` shrinkage for genes with low counts. Studies with <3 samples per group are skipped with a warning.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 2: HELPER FUNCTIONS (Colab & Q1 Publication Optimized)
# =============================================================================

#' Load coldata CSV for a given study.
#' Expected columns: sample_id, files (path to quant.sf), condition,
#' and any covariates (time_point, cell_type, brain_region, sex, etc.)
#' condition must be a factor with the CONTROL level set as reference.
load_coldata <- function(gse, salmon_dir = SALMON_DIR,
                         metadata_dir = METADATA_DIR) {
  meta_path <- file.path(metadata_dir, paste0(gse, "_coldata.csv"))
  if (!file.exists(meta_path)) {
    stop(sprintf(
      "Coldata file not found: %s\n  Create a CSV with columns: sample_id, files, condition [+ covariates]",
      meta_path
    ))
  }
  cd <- read_csv(meta_path, show_col_types = FALSE)

  # Validate required columns
  required_cols <- c("sample_id", "files", "condition")
  missing <- setdiff(required_cols, colnames(cd))
  if (length(missing) > 0) {
    stop(sprintf("Coldata for %s is missing columns: %s", gse, paste(missing, collapse = ", ")))
  }

  # --- COLAB FIX: Dynamically rewrite file paths to point to mounted Google Drive ---
  cd$files <- file.path(salmon_dir, gse, cd$sample_id, "quant.sf")

  # Verify all quant.sf files exist
  missing_files <- cd$files[!file.exists(cd$files)]
  if (length(missing_files) > 0) {
    stop(sprintf(
      "%d quant.sf file(s) not found for %s:\n  %s",
      length(missing_files), gse,
      paste(head(missing_files, 5), collapse = "\n  ")
    ))
  }

  # Ensure condition is a factor; set control as reference level
  cd$condition <- as.factor(cd$condition)
  ctrl_patterns <- c("^control$", "^mock$", "^uninfected$", "^nbh$", "^vehicle$",
                     "control", "mock", "uninfected")
  ctrl_level <- NA
  for (pat in ctrl_patterns) {
    matches <- grep(pat, levels(cd$condition), ignore.case = TRUE, value = TRUE)
    if (length(matches) > 0) { ctrl_level <- matches[1]; break }
  }
  if (!is.na(ctrl_level)) {
    cd$condition <- relevel(cd$condition, ref = ctrl_level)
    cat(sprintf("  [%s] Reference level set to: '%s'\n", gse, ctrl_level))
  } else {
    warning(sprintf(
      "[%s] Could not auto-detect control level. Current reference: '%s'. Check manually.",
      gse, levels(cd$condition)[1]
    ))
  }

  # Convert covariate columns to factors if present
  covariate_cols <- c("time_point", "cell_type", "brain_region", "sex", "age_group", "sample_type", "layout")
  for (col in intersect(covariate_cols, colnames(cd))) {
    cd[[col]] <- as.factor(cd[[col]])
  }

  return(cd)
}


#' Import Salmon quantification with tximeta.
#' Returns a gene-level SummarizedExperiment with length offsets.
#' Includes tximport fallback for Colab API timeouts.
import_with_tximeta <- function(coldata, gse) {
  cat(sprintf("  [%s] Importing with tximeta...\n", gse))

  coldata_tximeta <- as.data.frame(coldata)
  coldata_tximeta$names <- as.character(coldata$sample_id)
  coldata_tximeta$files <- as.character(coldata$files)

  cache_dir <- file.path(BASE_DIR, "resources", "cache", "tximeta")
  dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
  tximeta::setTximetaBFC(cache_dir)

  # --- COLAB FIX: Try tximeta first, fallback to tximport+biomaRt if API fails ---
  gse_gene <- tryCatch({
    se_tx <- tximeta::tximeta(coldata_tximeta, type = "salmon", dropInfReps = TRUE)
    tximeta::summarizeToGene(se_tx)
  }, error = function(e) {
    cat(sprintf("  [!] tximeta failed: %s\n  [!] Falling back to tximport + biomaRt...\n", e$message))
    mart <- tryCatch(
      useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://nov2022.archive.ensembl.org"),
      error = function(err) useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")
    )
    tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)
    txi <- tximport(coldata_tximeta$files, type = "salmon", tx2gene = tx2gene, ignoreTxVersion = TRUE, dropInfReps = TRUE)
    counts_mat <- round(txi$counts)
    storage.mode(counts_mat) <- "integer"
    SummarizedExperiment(
      assays = list(counts = counts_mat, abundance = txi$abundance, length = txi$length),
      colData = DataFrame(coldata_tximeta)
    )
  })

  cat(sprintf("  [%s] Imported: %d genes x %d samples\n",
              gse, nrow(gse_gene), ncol(gse_gene)))
  return(gse_gene)
}


#' Quality filter: remove low-count genes.
#' Keeps genes with >= MIN_COUNTS counts in >= MIN_FRAC fraction of samples.
filter_low_counts <- function(gse_obj, gse_id,
                              min_counts = 10, min_frac = 0.20) {
  n_before <- nrow(gse_obj)
  keep <- rowSums(assay(gse_obj, "counts") >= min_counts) >=
    ceiling(min_frac * ncol(gse_obj))
  gse_obj <- gse_obj[keep, ]
  n_after <- nrow(gse_obj)
  cat(sprintf("  [%s] Low-count filter: %d -> %d genes (removed %d)\n",
              gse_id, n_before, n_after, n_before - n_after))
  return(gse_obj)
}


#' PCA-based outlier detection on VST-normalized counts.
#' Flags samples > N_SD standard deviations from the centroid on PC1+PC2.
#' Returns a list: $gse_obj (filtered), $outliers (character vector of sample IDs)
detect_pca_outliers <- function(gse_obj, gse_id, n_sd = 3) {
  vst_design <- if ("condition" %in% colnames(colData(gse_obj))) ~ condition else ~ 1
  dds_for_pca <- DESeqDataSet(gse_obj, design = vst_design)

  # --- COLAB FIX: Prevent VST crash on highly purified datasets (<1000 genes) ---
  vst_obj <- tryCatch({
    if (nrow(dds_for_pca) < 1000) {
      varianceStabilizingTransformation(dds_for_pca, blind = TRUE)
    } else {
      vst(dds_for_pca, blind = TRUE)
    }
  }, error = function(e) {
    cat(sprintf("  [%s] VST fallback triggered. Using normTransform.\n", gse_id))
    normTransform(dds_for_pca)
  })

  pca_data <- plotPCA(vst_obj, intgroup = "condition", returnData = TRUE)

  pc_mat <- as.matrix(pca_data[, c("PC1", "PC2")])
  centroid <- colMeans(pc_mat)
  dists <- sqrt(rowSums(sweep(pc_mat, 2, centroid)^2))
  threshold <- mean(dists) + n_sd * sd(dists)
  outlier_names <- rownames(pca_data)[dists > threshold]

  if (length(outlier_names) > 0) {
    cat(sprintf("  [%s] PCA outliers detected (>%d SD): %s\n",
                gse_id, n_sd, paste(outlier_names, collapse = ", ")))
    gse_obj <- gse_obj[, !colnames(gse_obj) %in% outlier_names]
    cat(sprintf("  [%s] Samples after outlier removal: %d\n",
                gse_id, ncol(gse_obj)))
  } else {
    cat(sprintf("  [%s] No PCA outliers detected.\n", gse_id))
  }

  return(list(gse_obj = gse_obj, outliers = outlier_names))
}


#' Run DESeq2 + apeglm shrinkage for a single study.
#' Returns a data.frame with gene_id, log2FoldChange, lfcSE, pvalue, padj.
#' Saves results CSV to disk and removes dds from memory.
run_deseq2 <- function(gse_obj, gse_id, design_formula, coef_name,
                       out_dir = DESEQ2_DIR) {
  cat(sprintf("  [%s] Running DESeq2 with design: %s\n", gse_id, design_formula))

  dds_fit <- DESeq(
    DESeqDataSet(gse_obj, design = as.formula(design_formula)),
    fitType = "parametric",
    quiet   = TRUE
  )

  available_coefs <- resultsNames(dds_fit)
  cat(sprintf("  [%s] Available coefficients: %s\n",
              gse_id, paste(available_coefs, collapse = ", ")))

  if (!coef_name %in% available_coefs) {
    stop(sprintf(
      "[%s] Coefficient '%s' not found in resultsNames.\n  Available: %s\n  Check condition factor levels in coldata.",
      gse_id, coef_name, paste(available_coefs, collapse = ", ")
    ))
  }

  res_shrunk <- lfcShrink(dds_fit, coef = coef_name, type = "apeglm", quiet = TRUE)

  res_df <- as.data.frame(res_shrunk)
  res_df$gene_id <- rownames(res_df)
  res_df <- res_df[, c("gene_id", "baseMean", "log2FoldChange", "lfcSE",
                        "pvalue", "padj")]

  # Save results
  out_path <- file.path(out_dir, paste0(gse_id, "_deseq2_results.csv"))
  write_csv(res_df, out_path)
  cat(sprintf("  [%s] Results saved: %s\n", gse_id, out_path))

  # --- COLAB FIX: Apply the same VST fallback for the output RDS file ---
  vst_obj <- tryCatch({
    if (nrow(dds_fit) < 1000) {
      varianceStabilizingTransformation(dds_fit, blind = FALSE)
    } else {
      vst(dds_fit, blind = FALSE)
    }
  }, error = function(e) {
    normTransform(dds_fit)
  })

  vst_path <- file.path(out_dir, paste0(gse_id, "_vst.rds"))
  saveRDS(vst_obj, vst_path)

  # CRITICAL: remove large objects from memory before returning
  rm(dds_fit, res_shrunk, vst_obj)
  gc()

  return(res_df)
}


#' Map Ensembl gene IDs to gene symbols and Entrez IDs (mouse).
annotate_mouse_genes <- function(res_df) {
  ensembl_ids <- sub("\\..*$", "", res_df$gene_id)
  symbols <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "SYMBOL", keytype = "ENSEMBL", multiVals = "first")
  entrez <- mapIds(org.Mm.eg.db, keys = ensembl_ids, column = "ENTREZID", keytype = "ENSEMBL", multiVals = "first")
  res_df$gene_symbol <- symbols[ensembl_ids]
  res_df$entrez_id   <- entrez[ensembl_ids]
  res_df$ensembl_id  <- ensembl_ids
  return(res_df)
}


#' Save a ggplot to both PDF and PNG.
#' Q1 PUBLICATION UPGRADE: Vector PDFs & 600 DPI PNGs
save_plot <- function(p, filename_base, width = 8, height = 6) {
  pdf_path <- file.path(FIGURES_DIR, paste0(filename_base, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(filename_base, ".png"))

  # useDingbats = FALSE is crucial for opening PDFs in Adobe Illustrator without font errors
  ggsave(pdf_path, plot = p, width = width, height = height, useDingbats = FALSE)

  # 600 DPI ensures crisp, publication-ready raster images
  ggsave(png_path, plot = p, width = width, height = height, dpi = 600)

  cat(sprintf("  Saved: %s (Q1 Format: 600 DPI & Vector PDF)\n", filename_base))
}

cat("Helper functions defined.\n")


# ----------------------------------------------------------------------
# ##Section 3 — Per-Study DESeq2 (Stratum A Mouse Brain)
#
# **Purpose:** Run DESeq2 independently for each of the five Stratum A mouse brain studies (GSE149805, GSE198063, GSE144738, GSE189527, GSE201249). Each study is analysed with `~condition` as the design formula. Results are saved as per-study RDS files for downstream meta-analysis.
#
# **Inputs:** `../config/coldata/coldata_GSE*.tsv`; `../results/salmon_quant/` Salmon output directories.
#
# **Outputs:** `../results/stream1/deseq2_GSE*.rds` (one per study); `../results/stream1/results_GSE*.tsv` (shrunk LFC tables).
#
# **Key decisions:** `tximeta` is used for import (preserves transcript-level uncertainty and links to reference annotation). Genes with <10 total counts across all samples are pre-filtered before DESeq2 to improve power.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 3: PER-STUDY DESeq2 — STRATUM A (CNS DISCOVERY)
# =============================================================================

all_stratum_a_results <- list()  # stores per-study results data.frames

# -----------------------------------------------------------------------------
# 3.1  GSE149805 — RML scrapie, N=186, 4 cell types, HiSeq 2500
# -----------------------------------------------------------------------------
cat("\n=== GSE149805 (stratified by cell_type) ===\n")

gse149805_coldata <- load_coldata("GSE149805")

# Verify cell_type column exists
if (!"cell_type" %in% colnames(gse149805_coldata)) {
  stop("GSE149805 coldata must contain a 'cell_type' column with values: CamKIIa, Cx3Cr1, GFAP, PV")
}

cell_types_149805 <- levels(as.factor(gse149805_coldata$cell_type))
cat(sprintf("  Cell types found: %s\n", paste(cell_types_149805, collapse = ", ")))

gse149805_results_list <- list()

for (ct in cell_types_149805) {
  cat(sprintf("\n  --- GSE149805 | cell_type = %s ---\n", ct))

  # Subset coldata to this cell type
  cd_sub <- gse149805_coldata[gse149805_coldata$cell_type == ct, ]
  cat(sprintf("  Samples in this stratum: %d\n", nrow(cd_sub)))

  # Import with tximeta
  gse_sub <- import_with_tximeta(cd_sub, paste0("GSE149805_", ct))

  # Low-count filter
  gse_sub <- filter_low_counts(gse_sub, paste0("GSE149805_", ct))

  # PCA outlier detection
  outlier_res <- detect_pca_outliers(gse_sub, paste0("GSE149805_", ct))
  gse_sub <- outlier_res$gse_obj

  # DESeq2
  res_ct <- run_deseq2(
    gse_obj        = gse_sub,
    gse_id         = paste0("GSE149805_", ct),
    design_formula = STUDY_DESIGNS[["GSE149805"]],
    coef_name      = STUDY_COEF[["GSE149805"]],
    out_dir        = DESEQ2_DIR
  )

  res_ct <- annotate_mouse_genes(res_ct)
  res_ct$cell_type <- ct
  gse149805_results_list[[ct]] <- res_ct

  rm(gse_sub, cd_sub, outlier_res)
  gc()
}

# --- BUG FIX: Explicitly namespace dplyr to avoid AnnotationDbi conflicts ---
gse149805_combined <- dplyr::bind_rows(gse149805_results_list) %>%
  dplyr::filter(!is.na(log2FoldChange) & !is.na(lfcSE) & lfcSE > 0) %>%
  dplyr::group_by(gene_id, gene_symbol, entrez_id, ensembl_id) %>%
  dplyr::summarise(
    log2FoldChange = sum(log2FoldChange / lfcSE^2) / sum(1 / lfcSE^2),
    lfcSE          = sqrt(1 / sum(1 / lfcSE^2)),
    baseMean       = mean(baseMean, na.rm = TRUE),
    n_cell_types   = dplyr::n(),
    .groups = "drop"
  ) %>%
  dplyr::mutate(
    z_stat = log2FoldChange / lfcSE,
    pvalue = 2 * pnorm(-abs(z_stat)),
    padj   = p.adjust(pvalue, method = "BH")
  ) %>%
  dplyr::select(-z_stat)

# Re-save combined results (already includes gene_symbol from annotate_mouse_genes)
write_csv(gse149805_combined, file.path(DESEQ2_DIR, "GSE149805_combined_results.csv"))
all_stratum_a_results[["GSE149805"]] <- gse149805_combined
cat("\n  GSE149805 combined across cell types: saved.\n")
rm(gse149805_results_list, gse149805_combined)
gc()


# -----------------------------------------------------------------------------
# 3.2  GSE198063 — FFI D178N knock-in, N=159, pre-clinical, NovaSeq 6000
# -----------------------------------------------------------------------------
cat("\n=== GSE198063 ===\n")

gse198063_coldata <- load_coldata("GSE198063")

# Auto-detect sex covariate
design_198063 <- if ("sex" %in% colnames(gse198063_coldata)) {
  cat("  Sex covariate detected — using design: ~ sex + condition\n")
  "~ sex + condition"
} else {
  STUDY_DESIGNS[["GSE198063"]]
}

gse198063 <- import_with_tximeta(gse198063_coldata, "GSE198063")
gse198063 <- filter_low_counts(gse198063, "GSE198063")
outlier_res <- detect_pca_outliers(gse198063, "GSE198063")
gse198063 <- outlier_res$gse_obj

res_198063 <- run_deseq2(
  gse_obj        = gse198063,
  gse_id         = "GSE198063",
  design_formula = design_198063,
  coef_name      = STUDY_COEF[["GSE198063"]],
  out_dir        = DESEQ2_DIR
)
res_198063 <- annotate_mouse_genes(res_198063)
readr::write_csv(res_198063, file.path(DESEQ2_DIR, "GSE198063_deseq2_results.csv"))
all_stratum_a_results[["GSE198063"]] <- res_198063

rm(gse198063, outlier_res)
gc()


# -----------------------------------------------------------------------------
# 3.3  GSE144738 — RML scrapie, N=137, pre-clinical, hippocampus, HiSeq 4000
# -----------------------------------------------------------------------------
cat("\n=== GSE144738 ===\n")

gse144738_coldata <- load_coldata("GSE144738")

if (!"time_point" %in% colnames(gse144738_coldata)) {
  warning("GSE144738 coldata missing 'time_point' column. Falling back to ~ condition.")
  design_144738 <- "~ condition"
} else if (!"layout" %in% colnames(gse144738_coldata)) {
  warning("GSE144738 coldata missing 'layout' column. Re-run 04_build_coldata.R to add it. Falling back to ~ time_point + condition.")
  design_144738 <- "~ time_point + condition"
} else {
  layout_condition_tab <- table(gse144738_coldata$layout, gse144738_coldata$condition)
  cat("  [GSE144738] Layout x condition table:\n")
  print(layout_condition_tab)

  full_formula <- as.formula(STUDY_DESIGNS[["GSE144738"]])
  test_mm <- tryCatch(
    model.matrix(full_formula, data = gse144738_coldata),
    error = function(e) NULL
  )

  if (is.null(test_mm)) {
    warning("GSE144738: Could not build model matrix for ~ time_point + layout + condition. Falling back to ~ time_point + condition.")
    design_144738 <- "~ time_point + condition"
  } else if (Matrix::rankMatrix(test_mm)[1] < ncol(test_mm)) {
    warning(sprintf(
      "GSE144738: Design matrix rank-deficient (rank %d < %d cols). layout is collinear with condition — all single-end samples may share one condition. Falling back to ~ time_point + condition.",
      Matrix::rankMatrix(test_mm)[1], ncol(test_mm)
    ))
    design_144738 <- "~ time_point + condition"
  } else {
    cat(sprintf("  [GSE144738] Design matrix full rank (%d x %d). Using ~ time_point + layout + condition.\n",
                nrow(test_mm), ncol(test_mm)))
    design_144738 <- STUDY_DESIGNS[["GSE144738"]]
  }
  rm(full_formula, test_mm)
}

gse144738 <- import_with_tximeta(gse144738_coldata, "GSE144738")
gse144738 <- filter_low_counts(gse144738, "GSE144738")
outlier_res <- detect_pca_outliers(gse144738, "GSE144738")
gse144738 <- outlier_res$gse_obj

res_144738 <- run_deseq2(
  gse_obj        = gse144738,
  gse_id         = "GSE144738",
  design_formula = design_144738,
  coef_name      = STUDY_COEF[["GSE144738"]],
  out_dir        = DESEQ2_DIR
)
res_144738 <- annotate_mouse_genes(res_144738)
readr::write_csv(res_144738, file.path(DESEQ2_DIR, "GSE144738_deseq2_results.csv"))
all_stratum_a_results[["GSE144738"]] <- res_144738

rm(gse144738, outlier_res)
gc()


# -----------------------------------------------------------------------------
# 3.4  GSE189527 — RML scrapie, N=98, pre-clinical to terminal, HiSeq 2500
# -----------------------------------------------------------------------------
cat("\n=== GSE189527 ===\n")

gse189527_coldata <- load_coldata("GSE189527")

if (!"time_point" %in% colnames(gse189527_coldata)) {
  warning("GSE189527 coldata missing 'time_point' column. Falling back to ~ condition.")
  design_189527 <- "~ condition"
} else {
  design_189527 <- STUDY_DESIGNS[["GSE189527"]]
}

gse189527 <- import_with_tximeta(gse189527_coldata, "GSE189527")
gse189527 <- filter_low_counts(gse189527, "GSE189527")
outlier_res <- detect_pca_outliers(gse189527, "GSE189527")
gse189527 <- outlier_res$gse_obj

res_189527 <- run_deseq2(
  gse_obj        = gse189527,
  gse_id         = "GSE189527",
  design_formula = design_189527,
  coef_name      = STUDY_COEF[["GSE189527"]],
  out_dir        = DESEQ2_DIR
)
res_189527 <- annotate_mouse_genes(res_189527)
readr::write_csv(res_189527, file.path(DESEQ2_DIR, "GSE189527_deseq2_results.csv"))
all_stratum_a_results[["GSE189527"]] <- res_189527

rm(gse189527, outlier_res)
gc()


# -----------------------------------------------------------------------------
# 3.5  GSE201249 — RML scrapie, N=50, terminal, CA1 + Thalamus, NextSeq 550
# -----------------------------------------------------------------------------
cat("\n=== GSE201249 ===\n")

gse201249_coldata <- load_coldata("GSE201249")

missing_covs <- setdiff(c("brain_region", "cell_type"), colnames(gse201249_coldata))
if (length(missing_covs) > 0) {
  warning(sprintf(
    "GSE201249 coldata missing columns: %s. Falling back to ~ condition.",
    paste(missing_covs, collapse = ", ")
  ))
  design_201249 <- "~ condition"
} else {
  test_formula <- as.formula(STUDY_DESIGNS[["GSE201249"]])
  test_mm <- tryCatch(
    model.matrix(test_formula, data = gse201249_coldata),
    error = function(e) NULL
  )
  if (is.null(test_mm)) {
    warning("GSE201249: Could not build model matrix. Falling back to ~ brain_region + condition.")
    design_201249 <- "~ brain_region + condition"
  } else if (Matrix::rankMatrix(test_mm)[1] < ncol(test_mm)) {
    warning(sprintf(
      "GSE201249: Design matrix is rank-deficient. Falling back to ~ brain_region + condition.",
      Matrix::rankMatrix(test_mm)[1], ncol(test_mm)
    ))
    design_201249 <- "~ brain_region + condition"
  } else {
    cat(sprintf("  GSE201249: Design matrix is full rank (%d x %d). Using full design.\n",
                nrow(test_mm), ncol(test_mm)))
    design_201249 <- STUDY_DESIGNS[["GSE201249"]]
  }
  rm(test_formula, test_mm)
}

gse201249 <- import_with_tximeta(gse201249_coldata, "GSE201249")
gse201249 <- filter_low_counts(gse201249, "GSE201249")
outlier_res <- detect_pca_outliers(gse201249, "GSE201249")
gse201249 <- outlier_res$gse_obj

res_201249 <- run_deseq2(
  gse_obj        = gse201249,
  gse_id         = "GSE201249",
  design_formula = design_201249,
  coef_name      = STUDY_COEF[["GSE201249"]],
  out_dir        = DESEQ2_DIR
)
res_201249 <- annotate_mouse_genes(res_201249)
readr::write_csv(res_201249, file.path(DESEQ2_DIR, "GSE201249_deseq2_results.csv"))
all_stratum_a_results[["GSE201249"]] <- res_201249

rm(gse201249, outlier_res)
gc()

cat("\n=== All Stratum A DESeq2 runs complete ===\n")
cat(sprintf("Studies processed: %s\n", paste(names(all_stratum_a_results), collapse = ", ")))


# ----------------------------------------------------------------------
# ##Section 4 — Random-Effects Meta-Analysis (DerSimonian-Laird)
#
# **Purpose:** Combine per-study LFC estimates across all five Stratum A studies using a DerSimonian-Laird random-effects model (metafor). This accounts for between-study heterogeneity (τ²) and produces a pooled effect size with heterogeneity statistics (I², Q-test).
#
# **Inputs:** Per-study shrunk LFC + SE from Section 3.
#
# **Outputs:** `../results/stream1/meta_results.tsv` — pooled LFC, SE, z, p, padj (BH), I², τ² for all genes tested in ≥2 studies.
#
# **Key decisions:** Genes must appear in ≥2 studies to enter meta-analysis. Multiple testing correction uses Benjamini-Hochberg across all tested genes. High heterogeneity (I²>75%) is flagged but not excluded — reported transparently.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 4: RANDOM-EFFECTS META-ANALYSIS (DerSimonian-Laird)
# =============================================================================
# Operates on the 5 per-study results data.frames (small, in-memory).
# For each gene present in >= MIN_STUDIES (3) studies:
#   rma(yi = log2FC_vector, sei = lfcSE_vector, method = "DL")
# Apply Benjamini-Hochberg FDR correction across all genes.
# MCPS = genes with meta-padj <= 0.05, |meta_log2FC| >= 0.5, I2 <= 75%
# =============================================================================

cat("\n=== SECTION 4: Random-Effects Meta-Analysis ===\n")

# --- 4.1 Build unified gene universe ---
all_gene_ids <- Reduce(union, lapply(all_stratum_a_results, function(r) r$gene_id))
cat(sprintf("Total unique genes across all studies: %d\n", length(all_gene_ids)))

# Count how many studies each gene appears in with valid LFC + SE
gene_study_counts <- sapply(all_gene_ids, function(g) {
  sum(sapply(all_stratum_a_results, function(r) {
    idx <- which(r$gene_id == g)
    if (length(idx) == 0) return(FALSE)
    !is.na(r$log2FoldChange[idx]) & !is.na(r$lfcSE[idx]) & r$lfcSE[idx] > 0
  }))
})

meta_genes <- names(gene_study_counts)[gene_study_counts >= MIN_STUDIES]
cat(sprintf("Genes present in >= %d studies: %d\n", MIN_STUDIES, length(meta_genes)))

# --- 4.2 Run meta-analysis per gene ---
cat("Running DerSimonian-Laird random-effects meta-analysis...\n")

meta_results_list <- lapply(meta_genes, function(g) {
  yi_vec <- sei_vec <- numeric(0)
  study_vec <- character(0)

  for (study in names(all_stratum_a_results)) {
    r   <- all_stratum_a_results[[study]]
    idx <- which(r$gene_id == g)
    if (length(idx) == 0) next
    lfc <- r$log2FoldChange[idx]
    se  <- r$lfcSE[idx]
    if (is.na(lfc) || is.na(se) || se <= 0) next
    yi_vec    <- c(yi_vec, lfc)
    sei_vec   <- c(sei_vec, se)
    study_vec <- c(study_vec, study)
  }

  if (length(yi_vec) < 2) return(NULL)

  fit <- tryCatch(
    metafor::rma(yi = yi_vec, sei = sei_vec, method = "DL", verbose = FALSE),
    error = function(e) NULL
  )
  if (is.null(fit)) return(NULL)

  # Get gene annotation from first available study
  gene_info <- NULL
  for (study in names(all_stratum_a_results)) {
    r   <- all_stratum_a_results[[study]]
    idx <- which(r$gene_id == g)
    if (length(idx) > 0) {
      gene_info <- r[idx, intersect(c("gene_id","gene_symbol","entrez_id","ensembl_id"),
                                    colnames(r))]
      break
    }
  }

  data.frame(
    gene_id      = g,
    gene_symbol  = if (!is.null(gene_info) && "gene_symbol" %in% colnames(gene_info))
                     gene_info$gene_symbol else NA_character_,
    entrez_id    = if (!is.null(gene_info) && "entrez_id"   %in% colnames(gene_info))
                     gene_info$entrez_id   else NA_character_,
    ensembl_id   = if (!is.null(gene_info) && "ensembl_id"  %in% colnames(gene_info))
                     gene_info$ensembl_id  else NA_character_,
    meta_log2FC  = as.numeric(fit$beta),
    meta_SE      = as.numeric(fit$se),
    meta_zval    = as.numeric(fit$zval),
    meta_pvalue  = as.numeric(fit$pval),
    Q_stat       = as.numeric(fit$QE),
    Q_pvalue     = as.numeric(fit$QEp),
    I2           = as.numeric(fit$I2),
    tau2         = as.numeric(fit$tau2),
    n_studies    = length(yi_vec),
    studies      = paste(study_vec, collapse = ";"),
    stringsAsFactors = FALSE
  )
})

# COLAB FIX: Explicitly namespace dplyr here as well to avoid conflicts
meta_results <- dplyr::bind_rows(Filter(Negate(is.null), meta_results_list))
rm(meta_results_list); gc()
cat(sprintf("Meta-analysis complete: %d genes tested\n", nrow(meta_results)))

# --- 4.3 Multiple testing correction ---
meta_results$meta_padj <- p.adjust(meta_results$meta_pvalue, method = "BH")

# --- 4.4 Save full meta-analysis table ---
readr::write_csv(meta_results, file.path(META_DIR, "meta_analysis_full_results.csv"))
cat(sprintf("Full meta-analysis results saved: %d genes\n", nrow(meta_results)))

# --- 4.5 Define MCPS ---
MCPS <- meta_results %>%
  dplyr::filter(
    meta_padj        <= PADJ_THRESH,
    abs(meta_log2FC) >= LFC_THRESH,
    I2               <= I2_THRESH
  ) %>%
  dplyr::arrange(meta_padj)

MCPS_up <- MCPS %>% dplyr::filter(meta_log2FC >= LFC_THRESH)
MCPS_dn <- MCPS %>% dplyr::filter(meta_log2FC <= -LFC_THRESH)

cat(sprintf("\nMouse Conserved Prion Signature (MCPS):\n"))
cat(sprintf("  Total MCPS genes : %d\n",  nrow(MCPS)))
cat(sprintf("  Upregulated      : %d\n",  nrow(MCPS_up)))
cat(sprintf("  Downregulated    : %d\n",  nrow(MCPS_dn)))
cat(sprintf("  Excluded (I2>75%%): %d\n",
            sum(meta_results$meta_padj <= PADJ_THRESH &
                  abs(meta_results$meta_log2FC) >= LFC_THRESH &
                  meta_results$I2 > I2_THRESH, na.rm = TRUE)))

readr::write_csv(MCPS, file.path(META_DIR, "MCPS_genes.csv"))
cat(sprintf("MCPS gene list saved.\n"))

# --- 4.6 Per-study LFC matrix for MCPS genes ---
mcps_lfc_matrix <- lapply(names(all_stratum_a_results), function(study) {
  r     <- all_stratum_a_results[[study]]
  r_sub <- r[r$gene_id %in% MCPS$gene_id,
             c("gene_id", "log2FoldChange", "lfcSE", "padj")]
  colnames(r_sub)[2:4] <- paste0(study, c("_log2FC", "_lfcSE", "_padj"))
  r_sub
})

# COLAB FIX: Explicitly namespace join functions
mcps_lfc_wide <- Reduce(function(a, b) dplyr::full_join(a, b, by = "gene_id"),
                         mcps_lfc_matrix)
mcps_lfc_wide <- dplyr::left_join(
  mcps_lfc_wide,
  MCPS[, c("gene_id","gene_symbol","meta_log2FC","meta_padj","I2","n_studies")],
  by = "gene_id"
)
readr::write_csv(mcps_lfc_wide, file.path(META_DIR, "MCPS_per_study_lfc_matrix.csv"))
cat("MCPS per-study LFC matrix saved.\n")


# ----------------------------------------------------------------------
# ##Section 5 — Held-Out Replication (Stratum D)
#
# **Purpose:** Test whether the meta-analysis top hits replicate in four independent held-out datasets (Stratum D: GSE90977, GSE124932, GSE207251, GSE202275) that were not used in the discovery meta-analysis. This is the primary replication validation for the MCPS gene set.
#
# **Inputs:** `../results/stream1/meta_results.tsv`; Salmon quants for Stratum D studies.
#
# **Outputs:** `../results/stream1/replication_results.tsv` — per-study and pooled replication statistics; `replication_summary.tsv`.
#
# **Key decisions:** Replication criterion: same direction of effect AND nominal p<0.05 in ≥3/4 held-out studies. This is a pre-specified criterion that defines the MCPS gene set.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 5: HELD-OUT REPLICATION (STRATUM D)
# =============================================================================

cat("\n=== SECTION 5: Held-Out Replication (Stratum D) ===\n")

stratum_d_results <- list()

# --- 5.1 GSE90977 ---
cat("\n--- GSE90977 ---\n")
gse90977_coldata <- load_coldata("GSE90977")
design_90977 <- if ("time_point" %in% colnames(gse90977_coldata))
  STUDY_DESIGNS[["GSE90977"]] else "~ condition"
gse90977 <- import_with_tximeta(gse90977_coldata, "GSE90977")
gse90977 <- filter_low_counts(gse90977, "GSE90977")
res_90977 <- run_deseq2(gse90977, "GSE90977", design_90977,
                         STUDY_COEF[["GSE90977"]], REPLIC_DIR)
res_90977 <- annotate_mouse_genes(res_90977)
stratum_d_results[["GSE90977"]] <- res_90977
rm(gse90977); gc()

# --- 5.2 GSE124932 (WT arms only) ---
cat("\n--- GSE124932 (WT arms only) ---\n")
gse124932_coldata <- load_coldata("GSE124932")
if ("genotype" %in% colnames(gse124932_coldata)) {
  # COLAB FIX: Explicit dplyr namespace
  gse124932_coldata <- gse124932_coldata %>%
    dplyr::filter(grepl("^WT$|^wild.type$|^wildtype$", genotype, ignore.case = TRUE))
  cat(sprintf("  WT samples retained: %d\n", nrow(gse124932_coldata)))
} else {
  warning("GSE124932: 'genotype' column missing — verify WT-only samples manually.")
}
gse124932 <- import_with_tximeta(gse124932_coldata, "GSE124932")
gse124932 <- filter_low_counts(gse124932, "GSE124932")
res_124932 <- run_deseq2(gse124932, "GSE124932", STUDY_DESIGNS[["GSE124932"]],
                          STUDY_COEF[["GSE124932"]], REPLIC_DIR)
res_124932 <- annotate_mouse_genes(res_124932)
stratum_d_results[["GSE124932"]] <- res_124932
rm(gse124932); gc()

# --- 5.3 GSE207251 ---
cat("\n--- GSE207251 ---\n")
gse207251_coldata <- load_coldata("GSE207251")
gse207251 <- import_with_tximeta(gse207251_coldata, "GSE207251")
gse207251 <- filter_low_counts(gse207251, "GSE207251")
res_207251 <- run_deseq2(gse207251, "GSE207251", STUDY_DESIGNS[["GSE207251"]],
                          STUDY_COEF[["GSE207251"]], REPLIC_DIR)
res_207251 <- annotate_mouse_genes(res_207251)
stratum_d_results[["GSE207251"]] <- res_207251
rm(gse207251); gc()

# --- 5.4 GSE202275 (vehicle arms only) ---
cat("\n--- GSE202275 (vehicle arms only) ---\n")
gse202275_coldata <- load_coldata("GSE202275")
if ("treatment" %in% colnames(gse202275_coldata)) {
  # COLAB FIX: Explicit dplyr namespace
  gse202275_coldata <- gse202275_coldata %>%
    dplyr::filter(grepl("vehicle|veh", treatment, ignore.case = TRUE))
  cat(sprintf("  Vehicle-arm samples retained: %d\n", nrow(gse202275_coldata)))
} else {
  warning("GSE202275: 'treatment' column missing — verify vehicle arms manually.")
}
gse202275 <- import_with_tximeta(gse202275_coldata, "GSE202275")
gse202275 <- filter_low_counts(gse202275, "GSE202275")
res_202275 <- run_deseq2(gse202275, "GSE202275", STUDY_DESIGNS[["GSE202275"]],
                          STUDY_COEF[["GSE202275"]], REPLIC_DIR)
res_202275 <- annotate_mouse_genes(res_202275)
stratum_d_results[["GSE202275"]] <- res_202275
rm(gse202275); gc()

# --- 5.5 Replication statistics ---
cat("\n--- Replication statistics ---\n")

replication_summary <- lapply(names(stratum_d_results), function(study) {
  r <- stratum_d_results[[study]]

  # COLAB FIX: Explicit dplyr namespace
  mcps_check <- MCPS %>%
    dplyr::left_join(r[, c("gene_id","log2FoldChange","padj")],
              by = "gene_id", suffix = c("_meta","_replic"))

  n_tested   <- sum(!is.na(mcps_check$log2FoldChange))
  n_sig      <- sum(!is.na(mcps_check$padj) &
                      mcps_check$padj <= PADJ_THRESH, na.rm = TRUE)
  n_same_dir <- sum(!is.na(mcps_check$padj) &
                      mcps_check$padj <= PADJ_THRESH &
                      sign(mcps_check$log2FoldChange) == sign(mcps_check$meta_log2FC),
                    na.rm = TRUE)
  # Denominator is nrow(MCPS) (all MCPS genes), NOT n_tested.
  # Genes absent from the replication study (low counts / not expressed) are
  # genuine non-replications and must not be excluded from the denominator.
  replic_rate <- round(n_same_dir / nrow(MCPS) * 100, 1)

  cat(sprintf("  %s: %d/%d MCPS genes tested, %d sig & same-direction (%.1f%% replication)\n",
              study, n_tested, nrow(MCPS), n_same_dir, replic_rate))

  data.frame(study, mcps_total = nrow(MCPS), n_tested, n_significant = n_sig,
             n_same_direction = n_same_dir, replication_rate_pct = replic_rate,
             stringsAsFactors = FALSE)
})

# COLAB FIX: Explicit namespaces
replication_summary_df <- dplyr::bind_rows(replication_summary)
readr::write_csv(replication_summary_df, file.path(REPLIC_DIR, "replication_summary.csv"))
cat("Replication summary saved.\n")

# --- 5.6 Leave-One-Out (LOO) sensitivity analysis ---
cat("\n--- Leave-One-Out sensitivity analysis ---\n")

loo_results <- lapply(STRATUM_A, function(leave_out) {
  cat(sprintf("  LOO: excluding %s\n", leave_out))
  loo_studies <- setdiff(names(all_stratum_a_results), leave_out)

  loo_meta_list <- lapply(meta_genes, function(g) {
    yi_vec <- sei_vec <- numeric(0)
    for (study in loo_studies) {
      r   <- all_stratum_a_results[[study]]
      idx <- which(r$gene_id == g)
      if (length(idx) == 0) next
      lfc <- r$log2FoldChange[idx]; se <- r$lfcSE[idx]
      if (is.na(lfc) || is.na(se) || se <= 0) next
      yi_vec <- c(yi_vec, lfc); sei_vec <- c(sei_vec, se)
    }
    if (length(yi_vec) < 2) return(NULL)

    # COLAB FIX: Explicit metafor namespace
    fit <- tryCatch(metafor::rma(yi = yi_vec, sei = sei_vec, method = "DL", verbose = FALSE),
                    error = function(e) NULL)

    if (is.null(fit)) return(NULL)
    data.frame(gene_id = g, meta_log2FC = as.numeric(fit$beta),
               meta_pvalue = as.numeric(fit$pval), I2 = as.numeric(fit$I2),
               stringsAsFactors = FALSE)
  })

  # COLAB FIX: Explicit namespaces
  loo_df <- dplyr::bind_rows(Filter(Negate(is.null), loo_meta_list))
  loo_df$meta_padj <- p.adjust(loo_df$meta_pvalue, method = "BH")
  loo_mcps <- loo_df %>%
    dplyr::filter(meta_padj <= PADJ_THRESH, abs(meta_log2FC) >= LFC_THRESH, I2 <= I2_THRESH)

  overlap     <- sum(loo_mcps$gene_id %in% MCPS$gene_id)
  overlap_pct <- round(overlap / nrow(MCPS) * 100, 1)
  cat(sprintf("    LOO-%s: %d MCPS genes, %.1f%% overlap with full MCPS\n",
              leave_out, nrow(loo_mcps), overlap_pct))

  data.frame(excluded_study = leave_out, loo_mcps_size = nrow(loo_mcps),
             full_mcps_size = nrow(MCPS), overlap_n = overlap,
             overlap_pct = overlap_pct, stringsAsFactors = FALSE)
})

# COLAB FIX: Explicit namespaces
loo_summary <- dplyr::bind_rows(loo_results)
readr::write_csv(loo_summary, file.path(REPLIC_DIR, "LOO_sensitivity_summary.csv"))
cat("LOO sensitivity analysis saved.\n")

# --- LOO dominance warning ---
dominant_studies <- loo_summary %>%
  dplyr::filter(overlap_pct < 20) %>%
  dplyr::arrange(overlap_pct)

if (nrow(dominant_studies) > 0) {
  cat("\n*** LOO DOMINANCE WARNING ***\n")
  for (i in seq_len(nrow(dominant_studies))) {
    cat(sprintf(
      "  Excluding %s reduces MCPS overlap to %.1f%% (%d/%d genes).\n",
      dominant_studies$excluded_study[i],
      dominant_studies$overlap_pct[i],
      dominant_studies$overlap_n[i],
      dominant_studies$full_mcps_size[i]
    ))
  }
  cat("  ACTION REQUIRED: These studies are dominant drivers of the MCPS.\n")
  cat("  Manuscript must disclose this sensitivity finding.\n")
  cat("*****************************\n")
} else {
  cat("  LOO check passed: all studies show >=20%% MCPS overlap when excluded.\n")
}


# ----------------------------------------------------------------------
# ##Section 6 — Pathway Enrichment: ORA (per study) + GSEA (meta-analysis)
#
# **Purpose:** Characterise the biological processes enriched among prion-regulated genes. ORA (over-representation analysis) is run per study using significant DEGs; GSEA is run on the meta-analysis ranking (signed −log10 p × sign(LFC)) using MSigDB gene sets (Hallmark, GOBP, KEGG).
#
# **Inputs:** `../results/stream1/meta_results.tsv`; MSigDB GMT files.
#
# **Outputs:** `../results/stream1/gsea_results.tsv`; `../results/stream1/ora_results.tsv`.
#
# **Key decisions:** GSEA uses `fgsea` with 1,000 permutations and `eps=0` for accurate p-values at the extreme tail. Minimum gene set size = 15, maximum = 500.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 6: PATHWAY ENRICHMENT — ORA (per study) + GSEA (meta-analysis)
# =============================================================================

cat("\n=== SECTION 6: Pathway Enrichment ===\n")

# --- 6.1 Load MSigDB gene sets for mouse (Mus musculus) ---
cat("Loading MSigDB gene sets (mouse)...\n")

# COLAB FIX 1: Use 'collection' instead of 'category' (msigdbr >= 10.0.0)
# COLAB FIX 2: Strict dplyr namespaces to prevent AnnotationDbi crashes
# COLAB FIX 3: 'entrez_gene' was renamed to 'ncbi_gene' in msigdbr 10.0.0.
#              We select and rename it back to 'entrez_gene' for compatibility.
msig_hallmark  <- msigdbr(species = "Mus musculus", collection = "H") %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

# COLAB FIX 4: Fetch all C2 curated sets, then manually filter for KEGG and Reactome
# to bypass the "Unknown subcollection" error caused by MSigDB renaming CP:KEGG.
msig_c2 <- msigdbr(species = "Mus musculus", collection = "C2")

msig_kegg <- msig_c2 %>%
  dplyr::filter(grepl("^KEGG", gs_name)) %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

msig_reactome <- msig_c2 %>%
  dplyr::filter(grepl("^REACTOME", gs_name)) %>%
  dplyr::select(gs_name, entrez_gene = ncbi_gene) %>%
  dplyr::mutate(entrez_gene = as.character(entrez_gene))

# Convert to named list format for fgsea
make_gene_set_list <- function(msig_df) {
  split(msig_df$entrez_gene, msig_df$gs_name)
}

hallmark_list  <- make_gene_set_list(msig_hallmark)
kegg_list      <- make_gene_set_list(msig_kegg)
reactome_list  <- make_gene_set_list(msig_reactome)

cat(sprintf("  Hallmark gene sets : %d\n",  length(hallmark_list)))
cat(sprintf("  KEGG gene sets     : %d\n",  length(kegg_list)))
cat(sprintf("  Reactome gene sets : %d\n",  length(reactome_list)))

# --- 6.2 ORA helper function ---
run_ora <- function(gene_entrez_ids, background_entrez_ids, gse_id, direction) {
  cat(sprintf("  ORA [%s | %s]: %d genes\n", gse_id, direction, length(gene_entrez_ids)))

  if (length(gene_entrez_ids) < 5) {
    cat(sprintf("  [SKIP] Too few genes for ORA (%d < 5)\n", length(gene_entrez_ids)))
    return(NULL)
  }

  results_list <- list()

  # GO Biological Process
  go_res <- tryCatch(
    enrichGO(gene          = gene_entrez_ids,
             universe      = background_entrez_ids,
             OrgDb         = org.Mm.eg.db,
             ont           = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05,
             qvalueCutoff  = 0.2,
             readable      = TRUE),
    error = function(e) { cat(sprintf("  GO-BP error: %s\n", e$message)); NULL }
  )
  if (!is.null(go_res) && nrow(as.data.frame(go_res)) > 0) {
    df <- as.data.frame(go_res)
    df$database  <- "GO_BP"
    df$direction <- direction
    df$study     <- gse_id
    results_list[["GO_BP"]] <- df
  }

  # KEGG
  kegg_res <- tryCatch(
    enrichKEGG(gene          = gene_entrez_ids,
               universe      = background_entrez_ids,
               organism      = "mmu",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05),
    error = function(e) { cat(sprintf("  KEGG error: %s\n", e$message)); NULL }
  )
  if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    df <- as.data.frame(kegg_res)
    df$database  <- "KEGG"
    df$direction <- direction
    df$study     <- gse_id
    results_list[["KEGG"]] <- df
  }

  # MSigDB Hallmark
  hallmark_res <- tryCatch(
    enricher(gene          = gene_entrez_ids,
             universe      = background_entrez_ids,
             TERM2GENE     = msig_hallmark,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05),
    error = function(e) { cat(sprintf("  Hallmark error: %s\n", e$message)); NULL }
  )
  if (!is.null(hallmark_res) && nrow(as.data.frame(hallmark_res)) > 0) {
    df <- as.data.frame(hallmark_res)
    df$database  <- "Hallmark"
    df$direction <- direction
    df$study     <- gse_id
    results_list[["Hallmark"]] <- df
  }

  # Reactome
  reactome_res <- tryCatch(
    enricher(gene          = gene_entrez_ids,
             universe      = background_entrez_ids,
             TERM2GENE     = msig_reactome,
             pAdjustMethod = "BH",
             pvalueCutoff  = 0.05),
    error = function(e) { cat(sprintf("  Reactome error: %s\n", e$message)); NULL }
  )
  if (!is.null(reactome_res) && nrow(as.data.frame(reactome_res)) > 0) {
    df <- as.data.frame(reactome_res)
    df$database  <- "Reactome"
    df$direction <- direction
    df$study     <- gse_id
    results_list[["Reactome"]] <- df
  }

  if (length(results_list) == 0) return(NULL)
  dplyr::bind_rows(results_list)
}

# --- 6.3 Run ORA for each Stratum A study ---
all_ora_results <- list()

for (study in names(all_stratum_a_results)) {
  cat(sprintf("\n  ORA for %s\n", study))
  r <- all_stratum_a_results[[study]]

  background_ids <- r$entrez_id[!is.na(r$padj) & !is.na(r$entrez_id)]
  background_ids <- unique(as.character(background_ids))

  up_ids <- r$entrez_id[!is.na(r$padj) & r$padj <= PADJ_THRESH &
                           r$log2FoldChange >= LFC_THRESH & !is.na(r$entrez_id)]
  up_ids <- unique(as.character(up_ids))

  dn_ids <- r$entrez_id[!is.na(r$padj) & r$padj <= PADJ_THRESH &
                           r$log2FoldChange <= -LFC_THRESH & !is.na(r$entrez_id)]
  dn_ids <- unique(as.character(dn_ids))

  ora_up <- run_ora(up_ids, background_ids, study, "up")
  ora_dn <- run_ora(dn_ids, background_ids, study, "down")

  ora_combined <- dplyr::bind_rows(ora_up, ora_dn)
  if (!is.null(ora_combined) && nrow(ora_combined) > 0) {
    readr::write_csv(ora_combined, file.path(ENRICH_DIR, paste0(study, "_ORA_results.csv")))
    all_ora_results[[study]] <- ora_combined
    cat(sprintf("  %s ORA: %d significant pathways saved\n", study, nrow(ora_combined)))
  }
}

# --- 6.4 Consensus pathways: enriched in >= 3 of 5 Discovery studies ---
cat("\nIdentifying consensus pathways (enriched in >= 3/5 studies)...\n")

if (length(all_ora_results) >= 3) {
  all_ora_df <- dplyr::bind_rows(all_ora_results)

  pathway_counts <- all_ora_df %>%
    dplyr::filter(p.adjust <= 0.05) %>%
    dplyr::group_by(ID, Description, database) %>%
    dplyr::summarise(
      n_studies_significant = dplyr::n_distinct(study),
      studies               = paste(sort(unique(study)), collapse = ";"),
      mean_pvalue           = mean(pvalue, na.rm = TRUE),
      mean_padj             = mean(p.adjust, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::filter(n_studies_significant >= 3) %>%
    dplyr::arrange(dplyr::desc(n_studies_significant), mean_padj)

  readr::write_csv(pathway_counts, file.path(ENRICH_DIR, "consensus_pathways_3of5.csv"))
  cat(sprintf("Consensus pathways (>= 3/5 studies): %d\n", nrow(pathway_counts)))
} else {
  cat("Not enough ORA results for consensus pathway analysis.\n")
}

# --- 6.5 GSEA on meta-analysis ranking ---
cat("\nRunning GSEA on meta-analysis ranking...\n")

gsea_input <- meta_results %>%
  dplyr::filter(!is.na(meta_pvalue) & !is.na(meta_log2FC) & !is.na(entrez_id)) %>%
  dplyr::mutate(
    rank_metric = sign(meta_log2FC) * (-log10(pmax(meta_pvalue, 1e-300))),
    entrez_id   = as.character(entrez_id)
  ) %>%
  dplyr::filter(!is.na(entrez_id) & entrez_id != "NA") %>%
  dplyr::group_by(entrez_id) %>%
  dplyr::slice_max(abs(rank_metric), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup()

ranked_list <- setNames(gsea_input$rank_metric, gsea_input$entrez_id)
ranked_list <- sort(ranked_list, decreasing = TRUE)

cat(sprintf("  Ranked gene list: %d genes\n", length(ranked_list)))

run_gsea_safe <- function(gene_sets, set_name, ranked_vec) {
  cat(sprintf("  GSEA [%s]...\n", set_name))
  res <- tryCatch(
    fgsea(pathways   = gene_sets,
          stats      = ranked_vec,
          minSize    = 15,
          maxSize    = 500,
          eps        = 0.0),
    error = function(e) { cat(sprintf("  GSEA error: %s\n", e$message)); NULL }
  )
  if (is.null(res)) return(NULL)
  res <- res %>%
    dplyr::filter(padj <= 0.05) %>%
    dplyr::arrange(padj) %>%
    dplyr::mutate(database = set_name)
  res$leadingEdge <- sapply(res$leadingEdge, paste, collapse = ";")
  res
}

gsea_hallmark  <- run_gsea_safe(hallmark_list,  "Hallmark",  ranked_list)
gsea_kegg      <- run_gsea_safe(kegg_list,      "KEGG",      ranked_list)
gsea_reactome  <- run_gsea_safe(reactome_list,  "Reactome",  ranked_list)

gsea_all <- dplyr::bind_rows(gsea_hallmark, gsea_kegg, gsea_reactome)
if (nrow(gsea_all) > 0) {
  readr::write_csv(gsea_all, file.path(ENRICH_DIR, "meta_GSEA_results.csv"))
  cat(sprintf("GSEA results: %d significant gene sets saved\n", nrow(gsea_all)))
} else {
  cat("No significant GSEA results at padj <= 0.05\n")
}


# ----------------------------------------------------------------------
# ##Section 7 — Visualisations
#
# **Purpose:** Generate all publication-ready figures for Stream 1: volcano plots per study, forest plot for MCPS genes, heterogeneity summary, pathway bubble plot, and MCPS gene expression heatmap.
#
# **Inputs:** `../results/stream1/meta_results.tsv`; `../results/stream1/gsea_results.tsv`; per-study DESeq2 results.
#
# **Outputs:** `../results/stream1/figures/` — PNG + SVG for each figure.
#
# **Key decisions:** All figures use `theme_pub` for consistent styling. MCPS genes are highlighted in red on volcano plots. Forest plot shows per-study LFC ± 95% CI with pooled estimate.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 7: VISUALIZATIONS
# =============================================================================

cat("\n=== SECTION 7: Visualizations ===\n")

# --- Color palette ---
PRION_COLORS <- c(
  up       = "#D73027",   # red — upregulated
  down     = "#4575B4",   # blue — downregulated
  ns       = "#BDBDBD",   # grey — not significant
  low_i2   = "#1A9850",   # green — low heterogeneity
  high_i2  = "#FD8D3C"    # orange — high heterogeneity
)

# ---- Figure 1A: Meta-analysis volcano plot ----
cat("  Plotting Fig 1A: Meta-analysis volcano plot\n")

volcano_df <- meta_results %>%
  dplyr::mutate(
    neg_log10_padj = -log10(pmax(meta_padj, 1e-300)),
    sig_class = dplyr::case_when(
      meta_padj <= PADJ_THRESH & meta_log2FC >= LFC_THRESH  & I2 <= I2_THRESH ~ "MCPS_up",
      meta_padj <= PADJ_THRESH & meta_log2FC <= -LFC_THRESH & I2 <= I2_THRESH ~ "MCPS_down",
      meta_padj <= PADJ_THRESH & abs(meta_log2FC) >= LFC_THRESH & I2 > I2_THRESH ~ "High_I2",
      TRUE ~ "NS"
    )
  )

top_labels <- volcano_df %>%
  dplyr::filter(sig_class %in% c("MCPS_up","MCPS_down")) %>%
  dplyr::slice_min(meta_padj, n = 15) %>%
  dplyr::filter(!is.na(gene_symbol))

p_volcano <- ggplot(volcano_df, aes(x = meta_log2FC, y = neg_log10_padj)) +
  geom_point(data = dplyr::filter(volcano_df, sig_class == "NS"),
             color = PRION_COLORS["ns"], alpha = 0.3, size = 0.8) +
  geom_point(data = dplyr::filter(volcano_df, sig_class == "High_I2"),
             color = PRION_COLORS["high_i2"], alpha = 0.6, size = 1.2) +
  geom_point(data = dplyr::filter(volcano_df, sig_class == "MCPS_up"),
             color = PRION_COLORS["up"], alpha = 0.8, size = 1.5) +
  geom_point(data = dplyr::filter(volcano_df, sig_class == "MCPS_down"),
             color = PRION_COLORS["down"], alpha = 0.8, size = 1.5) +
  geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_hline(yintercept = -log10(PADJ_THRESH),
             linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_text_repel(data = top_labels,
                  aes(label = gene_symbol),
                  size = 2.8, max.overlaps = 20,
                  segment.color = "grey50", segment.size = 0.3) +
  labs(
    title    = "Prion Disease Meta-Analysis: Mouse Conserved Prion Signature",
    subtitle = sprintf("N=5 CNS Discovery studies | MCPS: %d genes (padj<=%.2f, |LFC|>=%.1f, I²<=75%%)",
                       nrow(MCPS), PADJ_THRESH, LFC_THRESH),
    x        = "Meta-analysis log2 Fold Change (prion vs control)",
    y        = expression(-log[10](adjusted~p-value)),
    caption  = paste0("Grey: NS | Red: MCPS upregulated | Blue: MCPS downregulated | ",
                      "Orange: Significant but I²>75% (excluded from MCPS)")
  ) +
  theme_bw(base_size = 11) +
  theme(plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9),
        plot.caption  = element_text(size = 7, color = "grey50"))

save_plot(p_volcano, "Fig1A_meta_volcano", width = 9, height = 7)


# ---- Figure 1B: Forest plots for top 20 MCPS genes ----
cat("  Plotting Fig 1B: Forest plots (top 20 MCPS genes)\n")

top20_mcps <- MCPS %>%
  dplyr::slice_min(meta_padj, n = 20) %>%
  dplyr::arrange(meta_log2FC)

forest_df <- lapply(top20_mcps$gene_id, function(g) {
  sym <- top20_mcps$gene_symbol[top20_mcps$gene_id == g]
  sym <- if (is.na(sym) || sym == "") g else sym

  per_study <- lapply(names(all_stratum_a_results), function(study) {
    r   <- all_stratum_a_results[[study]]
    idx <- which(r$gene_id == g)
    if (length(idx) == 0) return(NULL)
    data.frame(gene_symbol = sym, study = study,
               log2FC = r$log2FoldChange[idx], se = r$lfcSE[idx],
               type = "study", stringsAsFactors = FALSE)
  })

  meta_row <- data.frame(
    gene_symbol = sym,
    study       = "META",
    log2FC      = top20_mcps$meta_log2FC[top20_mcps$gene_id == g],
    se          = top20_mcps$meta_SE[top20_mcps$gene_id == g],
    type        = "meta",
    stringsAsFactors = FALSE
  )

  dplyr::bind_rows(c(Filter(Negate(is.null), per_study), list(meta_row)))
})

forest_df <- dplyr::bind_rows(forest_df) %>%
  dplyr::mutate(
    ci_lo = log2FC - 1.96 * se,
    ci_hi = log2FC + 1.96 * se,
    study = factor(study, levels = c(STRATUM_A, "META"))
  )

p_forest <- ggplot(forest_df, aes(x = log2FC, y = study,
                                   color = type, size = type)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.4) +
  # COLAB FIX: Updated geom_errorbar to modern ggplot2 standards (removed geom_errorbarh warning)
  geom_errorbar(aes(xmin = ci_lo, xmax = ci_hi), width = 0.3, linewidth = 0.5) +
  geom_point() +
  scale_color_manual(values = c(study = "steelblue", meta = "firebrick"),
                     labels = c(study = "Per-study estimate", meta = "Meta-estimate")) +
  scale_size_manual(values  = c(study = 2, meta = 3.5)) +
  facet_wrap(~ gene_symbol, ncol = 4, scales = "free_x") +
  labs(
    title  = "Forest Plots: Top 20 MCPS Genes",
    x      = "log2 Fold Change (95% CI)",
    y      = NULL,
    color  = NULL, size = NULL
  ) +
  theme_bw(base_size = 9) +
  theme(legend.position = "bottom",
        strip.text = element_text(face = "bold", size = 8),
        axis.text.y = element_text(size = 7))

save_plot(p_forest, "Fig1B_forest_plots_top20", width = 14, height = 10)


# ---- Figure 1C: I2 distribution histogram ----
cat("  Plotting Fig 1C: I2 distribution\n")

i2_df <- meta_results %>%
  dplyr::filter(!is.na(I2)) %>%
  dplyr::mutate(sig = meta_padj <= PADJ_THRESH & abs(meta_log2FC) >= LFC_THRESH)

p_i2 <- ggplot(i2_df, aes(x = I2, fill = sig)) +
  geom_histogram(binwidth = 5, color = "white", linewidth = 0.2) +
  geom_vline(xintercept = I2_THRESH, linetype = "dashed",
             color = "firebrick", linewidth = 0.8) +
  annotate("text", x = I2_THRESH + 2, y = Inf,
           label = sprintf("I²=%d%% exclusion threshold", I2_THRESH),
           hjust = 0, vjust = 1.5, size = 3, color = "firebrick") +
  scale_fill_manual(values = c("FALSE" = "grey70", "TRUE" = "steelblue"),
                    labels = c("FALSE" = "Not significant", "TRUE" = "Significant DEG")) +
  labs(
    title = "Heterogeneity Distribution (I²) Across Meta-Analysis Genes",
    x     = "I² (%)",
    y     = "Number of genes",
    fill  = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(legend.position = "top")

save_plot(p_i2, "Fig1C_I2_distribution", width = 7, height = 5)


# ---- Figure 2: Per-study DEG count barplot ----
cat("  Plotting Fig 2: Per-study DEG counts\n")

deg_counts <- lapply(names(all_stratum_a_results), function(study) {
  r <- all_stratum_a_results[[study]]
  data.frame(
    study = study,
    direction = c("Upregulated", "Downregulated"),
    n_degs = c(
      sum(!is.na(r$padj) & r$padj <= PADJ_THRESH & r$log2FoldChange >= LFC_THRESH),
      -sum(!is.na(r$padj) & r$padj <= PADJ_THRESH & r$log2FoldChange <= -LFC_THRESH)
    ),
    stringsAsFactors = FALSE
  )
})
deg_counts_df <- dplyr::bind_rows(deg_counts)

p_deg_counts <- ggplot(deg_counts_df,
                        aes(x = study, y = n_degs, fill = direction)) +
  geom_col(width = 0.6) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_fill_manual(values = c(Upregulated = PRION_COLORS["up"],
                                Downregulated = PRION_COLORS["down"])) +
  scale_y_continuous(labels = function(x) abs(x)) +
  labs(
    title = "DEG Counts per CNS Discovery Study",
    subtitle = sprintf("padj <= %.2f, |log2FC| >= %.1f", PADJ_THRESH, LFC_THRESH),
    x = NULL, y = "Number of DEGs", fill = NULL
  ) +
  theme_bw(base_size = 11) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1),
        legend.position = "top")

save_plot(p_deg_counts, "Fig2_per_study_DEG_counts", width = 7, height = 5)


# ---- Figure 3: MCPS heatmap (top 50 genes by |meta_log2FC|) ----
if (exists("mcps_lfc_wide") && nrow(mcps_lfc_wide) > 0) {

cat("  Plotting Fig 3: MCPS heatmap\n")

top50_mcps <- MCPS %>%
  dplyr::slice_max(abs(meta_log2FC), n = 50) %>%
  dplyr::arrange(dplyr::desc(meta_log2FC))

# COLAB FIX: Removed redundant left_join to avoid .x and .y suffix crash
lfc_cols <- paste0(STRATUM_A, "_log2FC")
heatmap_mat <- mcps_lfc_wide %>%
  dplyr::filter(gene_id %in% top50_mcps$gene_id) %>%
  dplyr::arrange(dplyr::desc(meta_log2FC))

row_labels <- ifelse(is.na(heatmap_mat$gene_symbol) | heatmap_mat$gene_symbol == "",
                     heatmap_mat$gene_id, heatmap_mat$gene_symbol)

mat <- as.matrix(heatmap_mat[, lfc_cols])
rownames(mat) <- row_labels
colnames(mat) <- STRATUM_A

mat[is.na(mat)] <- 0

max_lfc <- max(abs(mat), na.rm = TRUE)
col_fun <- colorRamp2(
  c(-max_lfc, -LFC_THRESH, 0, LFC_THRESH, max_lfc),
  c("#4575B4", "#ABD9E9", "white", "#FDAE61", "#D73027")
)

meta_lfc_vec <- heatmap_mat$meta_log2FC
ha_row <- rowAnnotation(
  `Meta log2FC` = anno_barplot(
    meta_lfc_vec,
    gp = gpar(fill = ifelse(meta_lfc_vec > 0, PRION_COLORS["up"], PRION_COLORS["down"]),
              col  = NA),
    width = unit(2, "cm")
  )
)

pdf(file.path(FIGURES_DIR, "Fig3_MCPS_heatmap.pdf"), width = 10, height = 12, useDingbats = FALSE)
draw(Heatmap(
  mat,
  name            = "log2FC",
  col             = col_fun,
  cluster_rows    = FALSE,
  cluster_columns = FALSE,
  show_row_names  = TRUE,
  row_names_gp    = gpar(fontsize = 8),
  column_names_gp = gpar(fontsize = 9, fontface = "bold"),
  column_title    = "Top 50 MCPS Genes by |meta-log2FC|",
  column_title_gp = gpar(fontsize = 11, fontface = "bold"),
  right_annotation = ha_row,
  heatmap_legend_param = list(title = "log2FC", direction = "horizontal")
))
dev.off()
cat(sprintf("  Saved: %s\n", file.path(FIGURES_DIR, "Fig3_MCPS_heatmap.pdf")))

} else {
  cat("  Skipping Fig 3: mcps_lfc_wide not available.\n")
}


# ---- Figure 4: Consensus pathway dotplot ----
cat("  Plotting Fig 4: Consensus pathway dotplot\n")

if (exists("pathway_counts") && nrow(pathway_counts) > 0) {
  top_pathways <- pathway_counts %>%
    dplyr::slice_min(mean_padj, n = 30) %>%
    dplyr::mutate(
      Description_short = stringr::str_trunc(Description, 55),
      Description_short = factor(Description_short,
                                  levels = rev(unique(Description_short)))
    )

  p_pathways <- ggplot(top_pathways,
                        aes(x = n_studies_significant,
                            y = Description_short,
                            color = mean_padj,
                            size  = n_studies_significant)) +
    geom_point() +
    scale_color_gradient(low = "firebrick", high = "steelblue",
                         name = "Mean adj. p-value") +
    scale_size_continuous(range = c(3, 8), name = "Studies significant") +
    scale_x_continuous(breaks = 1:5, limits = c(0.5, 5.5)) +
    geom_vline(xintercept = 3, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    facet_wrap(~ database, scales = "free_y", ncol = 1) +
    labs(
      title    = "Consensus Pathways Enriched in >= 3/5 CNS Discovery Studies",
      subtitle = "ORA: padj <= 0.05 per study",
      x        = "Number of studies with significant enrichment",
      y        = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(axis.text.y = element_text(size = 8),
          strip.text  = element_text(face = "bold"))

  save_plot(p_pathways, "Fig4_consensus_pathways", width = 10, height = 14)
} else {
  cat("  Skipping Fig 4: no consensus pathway data available.\n")
}


# ---- Figure 5: GSEA enrichment barplot (top Hallmark gene sets) ----
cat("  Plotting Fig 5: GSEA Hallmark enrichment\n")

if (exists("gsea_hallmark") && !is.null(gsea_hallmark) && nrow(gsea_hallmark) > 0) {
  gsea_plot_df <- gsea_hallmark %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 20) %>%
    dplyr::mutate(
      pathway_short = stringr::str_remove(pathway, "^HALLMARK_"),
      pathway_short = stringr::str_replace_all(pathway_short, "_", " "),
      pathway_short = factor(pathway_short, levels = rev(unique(pathway_short))),
      direction     = ifelse(NES > 0, "Upregulated", "Downregulated")
    )

  p_gsea <- ggplot(gsea_plot_df,
                    aes(x = NES, y = pathway_short, fill = direction)) +
    geom_col(width = 0.7) +
    geom_text(aes(label = sprintf("padj=%.3f", padj)),
              hjust = ifelse(gsea_plot_df$NES > 0, -0.1, 1.1),
              size = 2.8) +
    scale_fill_manual(values = c(Upregulated   = PRION_COLORS["up"],
                                  Downregulated = PRION_COLORS["down"])) +
    geom_vline(xintercept = 0, linewidth = 0.5) +
    labs(
      title    = "GSEA: Hallmark Gene Sets (Meta-Analysis Ranking)",
      subtitle = "Ranked by sign(meta_log2FC) × -log10(meta_pvalue)",
      x        = "Normalized Enrichment Score (NES)",
      y        = NULL,
      fill     = NULL
    ) +
    theme_bw(base_size = 10) +
    theme(legend.position = "top",
          axis.text.y = element_text(size = 9))

  save_plot(p_gsea, "Fig5_GSEA_Hallmark", width = 9, height = 8)
} else {
  cat("  Skipping Fig 5: no significant GSEA Hallmark results.\n")
}


# ---- Figure 6: Replication barplot (Stratum D) ----
cat("  Plotting Fig 6: Replication barplot\n")

if (exists("replication_summary_df") && nrow(replication_summary_df) > 0) {
  p_replic <- ggplot(replication_summary_df,
                      aes(x = study, y = replication_rate_pct, fill = study)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d/%d)",
                                   replication_rate_pct,
                                   n_same_direction, n_tested)),
              vjust = -0.3, size = 3.2) +
    geom_hline(yintercept = 50, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    annotate("text", x = 0.5, y = 52, label = "50% (chance)", hjust = 0,
             size = 3, color = "grey40") +
    scale_fill_brewer(palette = "Set2") +
    scale_y_continuous(limits = c(0, 110), labels = function(x) paste0(x, "%")) +
    labs(
      title    = "MCPS Replication in Held-Out Studies (Stratum D)",
      subtitle = sprintf("MCPS genes: %d | Threshold: padj<=%.2f, same direction",
                         nrow(MCPS), PADJ_THRESH),
      x        = NULL,
      y        = "Replication rate (%)"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  save_plot(p_replic, "Fig6_replication_barplot", width = 7, height = 5)
}


# ---- Figure 7: LOO sensitivity plot ----
cat("  Plotting Fig 7: LOO sensitivity\n")

if (exists("loo_summary") && nrow(loo_summary) > 0) {
  p_loo <- ggplot(loo_summary,
                   aes(x = excluded_study, y = overlap_pct, fill = excluded_study)) +
    geom_col(width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%\n(%d genes)", overlap_pct, loo_mcps_size)),
              vjust = -0.3, size = 3.2) +
    geom_hline(yintercept = 80, linetype = "dashed", color = "grey40", linewidth = 0.5) +
    scale_fill_brewer(palette = "Pastel1") +
    scale_y_continuous(limits = c(0, 115), labels = function(x) paste0(x, "%")) +
    labs(
      title    = "Leave-One-Out Sensitivity: MCPS Stability",
      subtitle = "% overlap between LOO-MCPS and full MCPS",
      x        = "Excluded study",
      y        = "MCPS overlap (%)"
    ) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 20, hjust = 1))

  save_plot(p_loo, "Fig7_LOO_sensitivity", width = 7, height = 5)
}

cat("\nAll figures saved to:", FIGURES_DIR, "\n")


# ----------------------------------------------------------------------
# ##Section 8 — Summary & Gene Symbol Aliases
#
# **Purpose:** Write a plain-text summary of Stream 1 results (sample counts, DEG counts, MCPS gene statistics, pathway top hits) and resolve all gene symbol aliases to ensure consistent naming across streams (e.g. Slc15a4 ↔ Tasl, Iigp1 ↔ Ifi47).
#
# **Inputs:** All Stream 1 result tables.
#
# **Outputs:** `../results/stream1/stream1_summary.txt`; `../results/stream1/gene_alias_map.tsv`.
#
# **Key decisions:** Alias resolution uses biomaRt `external_synonym` field. The alias map is used by nb08 to harmonise gene names across streams.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 8: FINAL SUMMARY REPORT
# =============================================================================

cat("\n=== SECTION 8: Final Summary ===\n")

summary_lines <- c(
  "=============================================================",
  "  PRION RNA-SEQ META-ANALYSIS — STREAM 1 SUMMARY",
  "=============================================================",
  "",
  sprintf("  Stratum A studies processed : %d", length(STRATUM_A)),
  sprintf("  Stratum D studies processed : %d", length(STRATUM_D)),
  "",
  "  --- Per-Study DEG Counts (padj<=0.05, |LFC|>=0.5) ---"
)

for (study in names(all_stratum_a_results)) {
  r  <- all_stratum_a_results[[study]]
  nu <- sum(!is.na(r$padj) & r$padj <= PADJ_THRESH & r$log2FoldChange >= LFC_THRESH)
  nd <- sum(!is.na(r$padj) & r$padj <= PADJ_THRESH & r$log2FoldChange <= -LFC_THRESH)
  summary_lines <- c(summary_lines,
    sprintf("    %-12s : %d up, %d down (%d total)", study, nu, nd, nu + nd))
}

summary_lines <- c(summary_lines,
  "",
  "  --- Meta-Analysis ---",
  sprintf("    Genes tested (>= %d studies) : %d", MIN_STUDIES, nrow(meta_results)),
  sprintf("    MCPS total                   : %d", nrow(MCPS)),
  sprintf("    MCPS upregulated             : %d", nrow(MCPS_up)),
  sprintf("    MCPS downregulated           : %d", nrow(MCPS_dn)),
  sprintf("    Excluded (I2 > %d%%)          : %d",
          I2_THRESH,
          sum(meta_results$meta_padj <= PADJ_THRESH &
                abs(meta_results$meta_log2FC) >= LFC_THRESH &
                meta_results$I2 > I2_THRESH, na.rm = TRUE)),
  "",
  "  --- Replication (Stratum D) ---"
)

if (exists("replication_summary_df")) {
  for (i in seq_len(nrow(replication_summary_df))) {
    r <- replication_summary_df[i, ]
    summary_lines <- c(summary_lines,
      sprintf("    %-12s : %.1f%% (%d/%d same direction, padj<=0.05)",
              r$study, r$replication_rate_pct, r$n_same_direction, r$n_tested))
  }
}

summary_lines <- c(summary_lines,
  "",
  "  --- LOO Sensitivity ---"
)

if (exists("loo_summary")) {
  for (i in seq_len(nrow(loo_summary))) {
    l <- loo_summary[i, ]
    summary_lines <- c(summary_lines,
      sprintf("    LOO (excl. %-12s): %.1f%% MCPS overlap (%d genes)",
              l$excluded_study, l$overlap_pct, l$loo_mcps_size))
  }
}

summary_lines <- c(summary_lines,
  "",
  "  --- Pathway Enrichment ---"
)

if (exists("pathway_counts")) {
  summary_lines <- c(summary_lines,
    sprintf("    Consensus pathways (>= 3/5 studies) : %d", nrow(pathway_counts)))
}
if (exists("gsea_all") && nrow(gsea_all) > 0) {
  summary_lines <- c(summary_lines,
    sprintf("    GSEA significant gene sets          : %d", nrow(gsea_all)))
}

summary_lines <- c(summary_lines,
  "",
  "  --- Output Files ---",
  sprintf("    DESeq2 results  : %s/", DESEQ2_DIR),
  sprintf("    Meta-analysis   : %s/", META_DIR),
  sprintf("    Replication     : %s/", REPLIC_DIR),
  sprintf("    Enrichment      : %s/", ENRICH_DIR),
  sprintf("    Figures         : %s/", FIGURES_DIR),
  "",
  "=============================================================",
  "  Stream 1 complete. Proceed to Stream 2 (Human Validation)",
  "  and Stream 3 (Peripheral Analysis).",
  "============================================================="
)

cat(paste(summary_lines, collapse = "\n"), "\n")

# Save summary to text file
writeLines(summary_lines, file.path(META_DIR, "stream1_summary.txt"))
cat(sprintf("\nSummary saved: %s\n", file.path(META_DIR, "stream1_summary.txt")))

# Save full session info for reproducibility
sink(file.path(META_DIR, "session_info.txt"))
print(sessionInfo())
sink()
cat("Session info saved.\n")

# =============================================================================
# STREAM 1 ALIAS EXPORT
# Write a canonical directory consumed by Stream 2 & 3.
# =============================================================================
cat("\n--- Writing stream1 alias outputs ---\n")

# COLAB FIX: Use absolute BASE_DIR path instead of relative "../results/stream1"
STREAM1_ALIAS_DIR <- file.path(BASE_DIR, "results/stream1")
dir.create(STREAM1_ALIAS_DIR, showWarnings = FALSE, recursive = TRUE)

# --- MCPS_genes.csv: full MCPS table with canonical Stream 1 column names ---
readr::write_csv(MCPS, file.path(STREAM1_ALIAS_DIR, "MCPS_genes.csv"))
cat(sprintf("  Written: stream1/MCPS_genes.csv  (%d genes)\n", nrow(MCPS)))

# --- MCPS.csv: column-renamed alias for Stream 2 ---
MCPS_alias <- MCPS %>%
  dplyr::rename(gene     = gene_symbol,
                meta_lfc = meta_log2FC) %>%
  dplyr::select(gene, meta_lfc, dplyr::everything())

readr::write_csv(MCPS_alias, file.path(STREAM1_ALIAS_DIR, "MCPS.csv"))
cat(sprintf("  Written: stream1/MCPS.csv  (%d genes, columns: gene, meta_lfc, ...)\n",
            nrow(MCPS_alias)))

# --- meta_analysis_full_results.csv: full meta-analysis table ---
meta_full_path <- file.path(META_DIR, "meta_analysis_full_results.csv")
if (file.exists(meta_full_path)) {
  file.copy(meta_full_path,
            file.path(STREAM1_ALIAS_DIR, "meta_analysis_full_results.csv"),
            overwrite = TRUE)
  cat("  Copied: stream1/meta_analysis_full_results.csv\n")
}

# --- MCPS_per_study_lfc_matrix.csv ---
mcps_matrix_path <- file.path(META_DIR, "MCPS_per_study_lfc_matrix.csv")
if (file.exists(mcps_matrix_path)) {
  file.copy(mcps_matrix_path,
            file.path(STREAM1_ALIAS_DIR, "MCPS_per_study_lfc_matrix.csv"),
            overwrite = TRUE)
  cat("  Copied: stream1/MCPS_per_study_lfc_matrix.csv\n")
}

cat("  Stream 1 alias export complete.\n")

cat("\n=== Stream 1 complete ===\n")
