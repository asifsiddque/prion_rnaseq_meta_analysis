#!/usr/bin/env Rscript
# ==============================================================================
# 07_stream3_peripheral_analysis.R
# Stream 3: peripheral time-course (GSE210128 per-timepoint Wald + LRT), terminal-stage replication, cross-tissue overlap, evidence tiers
# ==============================================================================
# Project : Cross-tissue bulk RNA-seq meta-analysis of prion disease
#           (mouse CNS Stratum A/D, peripheral time course, human sCJD validation)
#
# Provenance: Converted 1:1 from the executed Google Colab notebook
#             "07_stream3_peripheral_analysis_v2.ipynb" (the version that produced the manuscript results).
#             Only Colab environment-setup cells (drive mount, rpy2 magic,
#             package installers) were removed or commented; all analysis
#             code is unchanged.
#
# Inputs  : Salmon quantifications (GSE210128 + Stratum D), config/coldata/, Stream 1 MCPS, Stream 2 outputs
# Outputs : results/stream3/, reports/figures/stream3/
#
# Usage   : Set BASE_DIR (top of script) to the project root, then run:
#             Rscript 07_stream3_peripheral_analysis.R
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
# ##Install Required Packages (R)
#
# **Purpose:** Install all R packages needed for Stream 3 peripheral analysis. Identical package set to Streams 1–2 plus `survival` for time-to-event modelling in the first-detection analysis (Section 8).
#
# **Inputs:** Internet connection (first run only).
#
# **Outputs:** All packages installed and available for loading in Section 0.
#
# **Key decisions:** Run this cell only once per session — it takes ~5 minutes. `ask = FALSE, update = FALSE` prevents interactive prompts in Colab.
# ----------------------------------------------------------------------

if (!requireNamespace("BiocManager", quietly = TRUE))
  NULL  # [setup] installer removed (see environment.yml / renv.lock):
  # [setup]   install.packages("BiocManager")

NULL  # [setup] installer removed (see environment.yml / renv.lock):
# [setup]   BiocManager::install(c("tximeta", "DESeq2", "apeglm", "BiocParallel",
# [setup]   "SummarizedExperiment", "biomaRt", "AnnotationDbi",
# [setup]   "ComplexHeatmap"), ask = FALSE, update = FALSE)

NULL  # [setup] installer removed (see environment.yml / renv.lock):
# [setup]   install.packages(c("ggplot2", "ggrepel", "patchwork", "dplyr", "tidyr",
# [setup]   "tibble", "readr", "stringr", "scales", "RColorBrewer",
# [setup]   "viridis", "Matrix", "circlize"), repos = "https://cloud.r-project.org")


# ----------------------------------------------------------------------
# ##Phase 2: The Script
# Section 0 — Libraries, Configuration & Helper Functions
#
# **Purpose:** Load all R packages, define paths and parameters, configure BiocParallel, and define the `save_fig()` helper and `theme_pub` ggplot2 theme. Also defines the three analysis arms: Arm A (longitudinal LRT), Arm B (terminal Wald), Arm C (concordance + first-detection).
#
# **Inputs:** None (constants defined inline).
#
# **Outputs:** `BASE_DIR`, `COLDATA_DIR`, `STREAM1_DIR`, `STREAM3_DIR`, `MCPS_GENES`, `theme_pub`, `save_fig()` available in R.
#
# **Key decisions:** Stream 3 uses GSE149805 (longitudinal mouse blood, 4 time points: 60/90/120/150 dpi). The three-arm design separates trajectory (LRT), endpoint (Wald), and cross-tissue concordance analyses.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 0 — Libraries, configuration, helper functions
# =============================================================================

cat("=============================================================\n")
cat("  07_stream3_peripheral_analysis.R\n")
cat("  Stream 3: Peripheral Tissue Longitudinal Analysis\n")
cat("=============================================================\n")

suppressPackageStartupMessages({
  library(tximeta)
  library(DESeq2)
  library(apeglm)
  library(BiocParallel)
  library(SummarizedExperiment)
  library(biomaRt)          # Gene annotation (Section 7)
  library(ComplexHeatmap)
  library(circlize)
  library(dplyr)
  library(tidyr)
  library(tibble)
  library(readr)
  library(stringr)
  library(ggplot2)
  library(ggrepel)
  library(patchwork)
  library(Matrix)
  library(scales)
  library(RColorBrewer)
  library(viridis)
})

# ---------------------------------------------------------------------------
# 0.2  Parallel backend
# ---------------------------------------------------------------------------
n_cores <- min(4L, parallel::detectCores(logical = FALSE))
if (.Platform$OS.type == "windows" || n_cores <= 1L) {
  register(SerialParam())
  cat("  BiocParallel: SerialParam (Windows or single core)\n")
} else {
  register(MulticoreParam(workers = n_cores))
  cat(sprintf("  BiocParallel: MulticoreParam with %d workers\n", n_cores))
}

# ---------------------------------------------------------------------------
# 0.3  BASE_DIR resolution (Adapted for Colab)
# ---------------------------------------------------------------------------
BASE_DIR <- "/content/drive/MyDrive/prion_rnaseq_meta_analysis"

COLDATA_DIR <- file.path(BASE_DIR, "config",  "coldata")
STREAM1_DIR <- file.path(BASE_DIR, "results", "stream1")
OUT_DIR     <- file.path(BASE_DIR, "results", "stream3")
FIGURES_DIR <- file.path(BASE_DIR, "reports", "figures", "stream3")
CACHE_DIR   <- file.path(BASE_DIR, "resources", "cache", "tximeta")

for (d in c(OUT_DIR, FIGURES_DIR, CACHE_DIR)) dir.create(d, showWarnings = FALSE, recursive = TRUE)
tximeta::setTximetaBFC(CACHE_DIR)

# ---------------------------------------------------------------------------
# 0.4  Analysis parameters
# ---------------------------------------------------------------------------
PADJ_THRESH   <- 0.05    # Adjusted p-value threshold (inclusive)
LFC_THRESH    <- 0.5     # |log2FoldChange| threshold (inclusive)
TERMINAL_WPI  <- 22L     # Numeric encoding for "term wpi" (filter key only)
TOP_LRT_GENES <- 50L     # Genes shown in trajectory heatmap
MIN_REPS      <- 3L      # Minimum replicates per group per timepoint

TISSUES <- c("whole_blood", "spleen", "skeletal_muscle")
TISSUE_LABELS <- c(
  whole_blood     = "Whole blood",
  spleen          = "Spleen",
  skeletal_muscle = "Skeletal muscle"
)

# Timepoint levels and display labels (22 = terminal)
WPI_LEVELS <- c(4L, 8L, 12L, 14L, 16L, 18L, 20L, 22L)
WPI_LABELS <- c("4", "8", "12", "14", "16", "18", "20", "Term")

# Reproducibility seed (used before all DESeq2 calls and simulations)
set.seed(42)

# ---------------------------------------------------------------------------
# 0.5  Publication theme (Q1 Journal Standards)
# ---------------------------------------------------------------------------
theme_pub <- theme_classic(base_size = 10, base_family = "sans") +
  theme(
    axis.line        = element_line(colour = "black", linewidth = 0.4),
    axis.ticks       = element_line(colour = "black", linewidth = 0.3),
    axis.text        = element_text(colour = "black", size = 9),
    axis.title       = element_text(colour = "black", size = 10, face = "bold"),
    legend.title     = element_text(size = 9, face = "bold"),
    legend.text      = element_text(size = 8),
    legend.key.size  = unit(0.4, "cm"),
    legend.background     = element_blank(),
    legend.box.background = element_blank(),
    panel.grid.major = element_line(colour = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    panel.border     = element_blank(),
    strip.background = element_rect(fill = "grey95", colour = NA),
    strip.text       = element_text(size = 9, face = "bold"),
    plot.title       = element_text(size = 11, face = "bold", hjust = 0),
    plot.subtitle    = element_text(size = 9, colour = "grey30", hjust = 0),
    plot.caption     = element_text(size = 7, colour = "grey50", hjust = 1),
    # Strict padding to ensure no text/labels extend beyond the plotting area
    plot.margin      = margin(t = 10, r = 10, b = 10, l = 10, unit = "pt")
  )

# ---------------------------------------------------------------------------
# 0.6  Color palettes (colorblind-safe, Wong 2011)
# ---------------------------------------------------------------------------
COND_COLOURS <- c(control = "#0072B2", prion = "#D55E00")
COND_LABELS  <- c(control = "NBH (control)", prion = "RML6 (prion)")

TISSUE_COLOURS <- c(
  whole_blood     = "#E69F00",   # amber
  spleen          = "#56B4E9",   # sky blue
  skeletal_muscle = "#009E73"    # bluish green
)

# ---------------------------------------------------------------------------
# 0.7  Helper: save ggplot figure as PDF + PNG (Publication Quality)
# ---------------------------------------------------------------------------
# Default DPI increased to 600 for Q1 journal line-art standards.
# cairo_pdf used to ensure proper font rendering and transparency in final PDFs.
save_fig <- function(plot_obj, stem, width, height, dpi = 600) {
  pdf_path <- file.path(FIGURES_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(stem, ".png"))
  tryCatch({
    # Use cairo_pdf if cairo graphics are available (preferred: embeds fonts, no Dingbats).
    # Fall back to standard pdf(useDingbats=FALSE) on systems without cairo (e.g. some Colab instances).
    pdf_device <- if (capabilities("cairo")) {
      cairo_pdf
    } else {
      function(filename, ...) pdf(filename, ..., useDingbats = FALSE)
    }
    ggsave(pdf_path, plot = plot_obj, width = width, height = height,
           device = pdf_device, fallback_resolution = dpi, scale = 1)
    ggsave(png_path, plot = plot_obj, width = width, height = height,
           dpi = dpi, scale = 1)
    cat(sprintf("    Saved: %s.pdf + .png\n", stem))
  }, error = function(e) {
    warning(sprintf("save_fig failed for '%s': %s", stem, conditionMessage(e)))
  })
  invisible(list(pdf = pdf_path, png = png_path))
}

# ---------------------------------------------------------------------------
# 0.8  Helper: save ComplexHeatmap as PDF + PNG (Publication Quality)
# ---------------------------------------------------------------------------
save_heatmap <- function(ht_obj, stem, width, height, dpi = 600,
                         column_title = NULL, column_title_gp = NULL) {
  pdf_path <- file.path(FIGURES_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(stem, ".png"))

  draw_args <- list(ht_obj)
  if (!is.null(column_title))    draw_args$column_title    <- column_title
  if (!is.null(column_title_gp)) draw_args$column_title_gp <- column_title_gp

  tryCatch({
    # Using cairo_pdf equivalent for base graphics if available, else standard pdf
    pdf(pdf_path, width = width, height = height, useDingbats = FALSE)
    do.call(draw, draw_args)
    dev.off()
    cat(sprintf("    Saved: %s.pdf\n", stem))
  }, error = function(e) {
    warning(sprintf("save_heatmap PDF failed for '%s': %s", stem, conditionMessage(e)))
    if (dev.cur() > 1) dev.off()
  })

  tryCatch({
    # High-res PNG rendering
    png(png_path, width = width * dpi, height = height * dpi,
        res = dpi, units = "px", type = "cairo")
    do.call(draw, draw_args)
    dev.off()
    cat(sprintf("    Saved: %s.png\n", stem))
  }, error = function(e) {
    warning(sprintf("save_heatmap PNG failed for '%s': %s", stem, conditionMessage(e)))
    if (dev.cur() > 1) dev.off()
  })

  invisible(list(pdf = pdf_path, png = png_path))
}

# ---------------------------------------------------------------------------
# 0.9  Helper: run per-timepoint Wald test
# ---------------------------------------------------------------------------
run_wald_at_wpi <- function(dds_full_tissue, wpi_val, min_reps = MIN_REPS) {

  idx <- which(dds_full_tissue$time_point_wpi == wpi_val)
  if (length(idx) == 0) return(NULL)

  sub <- dds_full_tissue[, idx]
  n_prion   <- sum(sub$condition == "prion")
  n_control <- sum(sub$condition == "control")

  if (n_prion < min_reps || n_control < min_reps) {
    cat(sprintf("      wpi=%2d: skipped (prion=%d, control=%d; need >=%d each)\n",
                wpi_val, n_prion, n_control, min_reps))
    return(NULL)
  }

  has_both_cohorts <- length(unique(sub$cohort)) > 1
  des <- if (has_both_cohorts) ~ cohort + condition else ~ condition

  X <- tryCatch(
    model.matrix(des, data = as.data.frame(colData(sub))),
    error = function(e) NULL
  )
  if (is.null(X) || Matrix::rankMatrix(X)[1] < ncol(X)) {
    des <- ~ condition
    cat(sprintf("      wpi=%2d: cohort rank-deficient; using ~ condition\n", wpi_val))
  }
  design(sub) <- des

  sub <- tryCatch({
    dispersionFunction(sub) <- dispersionFunction(dds_full_tissue)
    sub
  }, error = function(e) {
    cat(sprintf("      wpi=%2d: dispersionFunction transfer failed; using independent estimation\n",
                wpi_val))
    NULL
  })

  if (is.null(sub)) {
    sub <- dds_full_tissue[, idx]
    design(sub) <- des
    sub <- tryCatch(
      estimateDispersions(sub, fitType = "parametric", quiet = TRUE),
      error = function(e) {
        tryCatch(
          estimateDispersions(sub, fitType = "local", quiet = TRUE),
          error = function(e2) NULL
        )
      }
    )
  } else {
    sub <- tryCatch(
      estimateDispersionsGeneEst(sub, quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(sub)) return(NULL)

    sub <- tryCatch(
      estimateDispersionsFit(sub, fitType = "parametric", quiet = TRUE),
      error = function(e) {
        tryCatch(
          estimateDispersionsFit(sub, fitType = "mean", quiet = TRUE),
          error = function(e2) NULL
        )
      }
    )
    if (is.null(sub)) return(NULL)

    sub <- tryCatch(
      estimateDispersionsMAP(sub, quiet = TRUE),
      error = function(e) NULL
    )
    if (is.null(sub)) return(NULL)
  }

  if (is.null(sub)) return(NULL)

  sub <- tryCatch(
    nbinomWaldTest(sub, quiet = TRUE),
    error = function(e) NULL
  )
  if (is.null(sub)) return(NULL)

  res <- tryCatch(
    results(sub, name = "condition_prion_vs_control",
            independentFiltering = TRUE),
    error = function(e) NULL
  )
  if (is.null(res)) return(NULL)

  cat(sprintf("      wpi=%2d: %d genes tested (%d prion, %d control)\n",
              wpi_val, sum(!is.na(res$padj)), n_prion, n_control))

  as.data.frame(res) %>%
    rownames_to_column("gene_id") %>%
    mutate(wpi = wpi_val) %>%
    dplyr::select(gene_id, log2FoldChange, lfcSE, pvalue, padj, wpi)
}

cat(sprintf("  BASE_DIR : %s\n", BASE_DIR))
cat(sprintf("  OUT_DIR  : %s\n", OUT_DIR))


# ----------------------------------------------------------------------
# ##Section 1 — Coldata & CONSORT
#
# **Purpose:** Load the GSE149805 coldata, validate sample assignments, and generate a CONSORT-style flow diagram showing sample inclusion/exclusion at each time point.
#
# **Inputs:** `../config/coldata/coldata_GSE149805.tsv`.
#
# **Outputs:** `coldata` data frame; CONSORT figure saved to `../results/stream3/figures/consort_flow.png`.
#
# **Key decisions:** Samples with missing time-point annotation are excluded. The CONSORT diagram is required for the manuscript Methods figure.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 1 — Load and validate coldata; CONSORT sample flow
# =============================================================================
cat("\n--- Section 1: Loading and validating coldata ---\n")

coldata_path <- file.path(COLDATA_DIR, "GSE210128_coldata.csv")
if (!file.exists(coldata_path)) {
  stop(paste(
    "GSE210128_coldata.csv not found at:", coldata_path,
    "\nRun 04_build_coldata.R first."
  ))
}

# --- 1.1  Load raw coldata ---
coldata_raw <- readr::read_csv(coldata_path, show_col_types = FALSE)
n_raw <- nrow(coldata_raw)
cat(sprintf("  Raw coldata: %d rows\n", n_raw))

# --- 1.2  Validate required columns ---
required_cols <- c("sample_id", "files", "condition", "tissue",
                   "time_point_wpi", "cohort")
missing_cols <- setdiff(required_cols, colnames(coldata_raw))
if (length(missing_cols) > 0) {
  stop(sprintf(
    "coldata is missing required columns: %s\nCheck 04_build_coldata.R output.",
    paste(missing_cols, collapse = ", ")
  ))
}

# --- 1.3  Type coercion, factor encoding, and COLAB PATH FIX ---
coldata <- coldata_raw %>%
  dplyr::mutate(
    condition      = factor(condition,      levels = c("control", "prion")),
    tissue         = factor(tissue,         levels = TISSUES),
    cohort         = factor(cohort,         levels = c("main", "validation")),
    time_point_wpi = as.integer(time_point_wpi),
    # FIXED FOR COLAB: Reconstruct exact absolute paths to Google Drive
    files          = file.path(BASE_DIR, "results", "salmon_quant", "GSE210128", sample_id, "quant.sf")
  )

# --- 1.4  CONSORT-style sample flow ---
# Track how many samples are excluded at each QC step.
n_after_coerce <- nrow(coldata)

# Remove samples with missing tissue, timepoint, or condition
coldata <- coldata %>%
  dplyr::filter(!is.na(tissue), !is.na(time_point_wpi), !is.na(condition))
n_after_na_filter <- nrow(coldata)

# Remove samples with unrecognised tissue
coldata <- coldata %>% dplyr::filter(tissue %in% TISSUES)
n_after_tissue_filter <- nrow(coldata)

# Remove samples with unrecognised timepoint
coldata <- coldata %>% dplyr::filter(time_point_wpi %in% WPI_LEVELS)
n_after_wpi_filter <- nrow(coldata)

cat(sprintf("  CONSORT sample flow:\n"))
cat(sprintf("    Raw:                    %d\n", n_raw))
cat(sprintf("    After NA filter:        %d  (excluded: %d)\n",
            n_after_na_filter, n_raw - n_after_na_filter))
cat(sprintf("    After tissue filter:    %d  (excluded: %d)\n",
            n_after_tissue_filter, n_after_na_filter - n_after_tissue_filter))
cat(sprintf("    After timepoint filter: %d  (excluded: %d)\n",
            n_after_wpi_filter, n_after_tissue_filter - n_after_wpi_filter))
cat(sprintf("    Final N:                %d\n", n_after_wpi_filter))

# --- 1.5  Validate quant.sf paths ---
missing_files <- coldata$files[!file.exists(coldata$files)]
if (length(missing_files) > 0) {
  warning(sprintf(
    "%d quant.sf file(s) not found. DESeq2 will fail if files are absent.\n%s%s",
    length(missing_files),
    paste(" ", head(missing_files, 5), collapse = "\n"),
    if (length(missing_files) > 5)
      sprintf("\n  ... and %d more", length(missing_files) - 5) else ""
  ))
} else {
  cat("  [QC PASS] All quant.sf paths verified successfully in Google Drive.\n")
}

# --- 1.6  Add display factor for timepoint ---
coldata <- coldata %>%
  dplyr::mutate(
    wpi_factor = factor(time_point_wpi, levels = WPI_LEVELS, labels = WPI_LABELS)
  )

# --- 1.7  Report cohort characteristics per tissue ---
cat(sprintf("\n  Sample breakdown (N = %d total):\n", nrow(coldata)))
for (tis in TISSUES) {
  sub <- coldata[coldata$tissue == tis, ]
  cat(sprintf("  %-20s N=%3d  (%d prion, %d control)  cohorts: %s  timepoints: %s\n",
              TISSUE_LABELS[tis],
              nrow(sub),
              sum(sub$condition == "prion"),
              sum(sub$condition == "control"),
              paste(sort(unique(as.character(sub$cohort))), collapse = "+"),
              paste(sort(unique(sub$time_point_wpi)), collapse = ",")))
}


# ----------------------------------------------------------------------
# ##Section 2 — Robust Salmon Import
#
# **Purpose:** Import Salmon quantifications for all GSE149805 samples using `tximeta` with error handling for missing quant files. Summarises to gene level and applies pre-filtering (≥10 counts in ≥3 samples).
#
# **Inputs:** `../results/salmon_quant/SRR*/quant.sf` for GSE149805 samples.
#
# **Outputs:** `se` SummarizedExperiment object; import log with any missing files.
#
# **Key decisions:** Missing quant files are logged and the corresponding samples are dropped rather than halting — allows partial runs during development.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 2 — Load Salmon quantifications with tximeta; summarizeToGene
# =============================================================================
cat("\n--- Section 2: Loading Salmon quantifications ---\n")

# Extend timeout for large transcriptome downloads in Colab
options(timeout = 600)

# Clean up Colab cache locks that often cause tximeta to freeze
lock_file <- file.path(CACHE_DIR, "BiocFileCache.sqlite.LOCK")
if (file.exists(lock_file)) {
  cat("  [QC] Removed lingering BiocFileCache lock file.\n")
  file.remove(lock_file)
}

# tximeta requires a data.frame with columns 'names' (sample IDs) and 'files'
# (The 'files' paths were already fixed to absolute Google Drive paths in Section 1)
coldata_tximeta <- coldata %>%
  dplyr::rename(names = sample_id) %>%
  as.data.frame()

# --- Robust Import Strategy ---
gse <- tryCatch({
  cat("  Attempting import with tximeta (strict provenance tracking)...\n")
  se <- tximeta::tximeta(coldata_tximeta, type = "salmon", dropInfReps = TRUE)

  cat("  Attempting to summarize to gene level...\n")
  tximeta::summarizeToGene(se)

}, error = function(e) {
  cat(sprintf("\n  [!] tximeta import failed. Error: %s\n", e$message))
  cat("  [!] Falling back to robust tximport + biomaRt (bypassing hash checks)...\n")

  library(tximport)
  library(biomaRt)

  cat("      Fetching tx2gene mapping from Ensembl...\n")
  mart <- tryCatch(
    useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://nov2022.archive.ensembl.org"),
    error = function(err) useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = "https://www.ensembl.org")
  )
  tx2gene <- getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id"), mart = mart)

  cat("      Running tximport...\n")
  txi <- tximport(coldata_tximeta$files, type = "salmon", tx2gene = tx2gene,
                  ignoreTxVersion = TRUE, dropInfReps = TRUE)

  cat("      Constructing final SummarizedExperiment...\n")
  counts_mat <- round(txi$counts)
  storage.mode(counts_mat) <- "integer"

  SummarizedExperiment(
    assays = list(counts = counts_mat, abundance = txi$abundance, length = txi$length),
    colData = DataFrame(coldata_tximeta)
  )
})

cat(sprintf("\n  Loaded: %d genes x %d samples\n", nrow(gse), ncol(gse)))
cat(sprintf("  Genome: %s\n", tryCatch(genome(gse)[1], error = function(e) "unknown")))

# Propagate coldata factors to the SummarizedExperiment colData
# (tximeta may not preserve factor levels from the input data.frame)
colData(gse)$condition      <- factor(colData(gse)$condition, levels = c("control", "prion"))
colData(gse)$tissue         <- factor(colData(gse)$tissue, levels = TISSUES)
colData(gse)$cohort         <- factor(colData(gse)$cohort, levels = c("main", "validation"))
colData(gse)$time_point_wpi <- as.integer(colData(gse)$time_point_wpi)
colData(gse)$wpi_factor     <- factor(colData(gse)$time_point_wpi, levels = WPI_LEVELS, labels = WPI_LABELS)


# ----------------------------------------------------------------------
# ##Sections 3 & 4 — Sample QC & MCPS Gene Loading
#
# **Purpose:** Section 3: PCA-based sample QC with outlier detection. Section 4: Load the MCPS gene list from Stream 1 results and confirm all 6 MCPS genes are present in the GSE149805 count matrix.
#
# **Inputs:** `se` from Section 2; `../results/stream1/meta_results.tsv`.
#
# **Outputs:** `se_filtered` (QC-passed samples); `mcps_present` logical vector; QC figures in `../results/stream3/figures/`.
#
# **Key decisions:** MCPS gene presence check uses both primary symbols and aliases (from `gene_alias_map.tsv`) to handle cross-study naming differences.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 3 — Sample-level QC: library size + VST PCA per tissue
# =============================================================================
cat("\n--- Section 3: Sample-level QC ---\n")

# --- 3.1  Library size summary ---
lib_sizes <- colSums(counts(DESeqDataSet(gse, design = ~ condition)))
lib_size_df <- data.frame(
  sample_id      = names(lib_sizes),
  library_size   = lib_sizes,
  condition      = colData(gse)$condition,
  tissue         = colData(gse)$tissue,
  time_point_wpi = colData(gse)$time_point_wpi,
  cohort         = colData(gse)$cohort
)

cat(sprintf("  Library size: median = %.1f M  |  range = %.1f – %.1f M\n",
            median(lib_sizes) / 1e6,
            min(lib_sizes) / 1e6,
            max(lib_sizes) / 1e6))

# Flag low-library-size samples (< 1 M reads) for reporting
low_lib <- lib_size_df %>% dplyr::filter(library_size < 1e6)
if (nrow(low_lib) > 0) {
  warning(sprintf("%d sample(s) have library size < 1 M reads:\n%s",
                  nrow(low_lib),
                  paste(low_lib$sample_id, collapse = ", ")))
}

# --- 3.2  VST PCA per tissue (blind = TRUE) ---
# Build a single DESeqDataSet for VST across all samples.
# Pre-filter: >= 10 counts in >= MIN_REPS samples (global filter for VST).
dds_all  <- DESeqDataSet(gse, design = ~ condition)
keep_all <- rowSums(counts(dds_all) >= 10) >= MIN_REPS
dds_all  <- dds_all[keep_all, ]
vsd_all  <- vst(dds_all, blind = TRUE)

cat(sprintf("  Genes retained after global pre-filter (>= 10 counts in >= %d samples): %d\n",
            MIN_REPS, nrow(dds_all)))

# Per-tissue PCA (stored for Fig 10)
pca_list <- lapply(TISSUES, function(tis) {
  idx      <- which(colData(vsd_all)$tissue == tis)
  vsd_sub  <- vsd_all[, idx]
  pca_data <- plotPCA(vsd_sub,
                      intgroup = c("condition", "time_point_wpi", "cohort"),
                      returnData = TRUE)
  pct_var  <- round(100 * attr(pca_data, "percentVar"), 1)
  cat(sprintf("  %s: PC1 = %.1f%%, PC2 = %.1f%%\n",
              TISSUE_LABELS[tis], pct_var[1], pct_var[2]))
  list(pca_data = pca_data, pct_var = pct_var, tissue = tis)
})
names(pca_list) <- TISSUES

# =============================================================================
# SECTION 4 — Load and validate MCPS gene list from Stream 1
# =============================================================================
cat("\n--- Section 4: Loading MCPS gene list ---\n")

# Stream 1 (05_stream1_cns_meta_analysis.R) writes MCPS_genes.csv with columns:
#   gene_id, gene_symbol, meta_log2FC, meta_padj, I2, n_studies
# gene_symbol = mouse MGI symbol (used for matching against GSE210128 results)
mcps_path <- file.path(STREAM1_DIR, "MCPS_genes.csv")
if (!file.exists(mcps_path)) {
  stop(paste(
    "MCPS_genes.csv not found at:", mcps_path,
    "\nRun 05_stream1_cns_meta_analysis.R first."
  ))
}

MCPS <- readr::read_csv(mcps_path, show_col_types = FALSE)

# --- 4.1  Validate required columns ---
if (!"gene_symbol" %in% colnames(MCPS))
  stop("MCPS_genes.csv must have a 'gene_symbol' column. Check Stream 1 output.")
if (!"meta_log2FC" %in% colnames(MCPS))
  stop("MCPS_genes.csv must have a 'meta_log2FC' column. Check Stream 1 output.")

# --- 4.2  Remove duplicates (keep highest |meta_log2FC| per symbol) ---
n_mcps_raw <- nrow(MCPS)
MCPS <- MCPS %>%
  dplyr::group_by(gene_symbol) %>%
  dplyr::slice_max(order_by = abs(meta_log2FC), n = 1, with_ties = FALSE) %>%
  dplyr::ungroup() %>%
  dplyr::filter(!is.na(gene_symbol), gene_symbol != "")

n_mcps    <- nrow(MCPS)
n_mcps_up <- sum(MCPS$meta_log2FC > 0, na.rm = TRUE)
n_mcps_dn <- sum(MCPS$meta_log2FC < 0, na.rm = TRUE)

if (n_mcps_raw != n_mcps)
  cat(sprintf("  MCPS: removed %d duplicate gene symbols\n", n_mcps_raw - n_mcps))

cat(sprintf("  MCPS loaded: %d genes (%d up, %d down in CNS)\n",
            n_mcps, n_mcps_up, n_mcps_dn))


# ----------------------------------------------------------------------
# ##Section 5 — Arm A: Longitudinal LRT
#
# **Purpose:** Test for genes with significant expression trajectories over prion disease progression using a likelihood ratio test (LRT). Design: `~time_point + condition + time_point:condition`; reduced model: `~time_point + condition`. The interaction term captures genes whose trajectory differs between prion and control.
#
# **Inputs:** `se_filtered` from Section 3.
#
# **Outputs:** `../results/stream3/lrt_results.tsv` — LRT statistic, p, padj for all genes; `../results/stream3/lrt_trajectories.tsv` — per-time-point normalised counts for top LRT genes.
#
# **Key decisions:** LRT is preferred over Wald test for trajectory analysis because it tests the full interaction term jointly across all time points. FDR threshold = 0.05 for LRT significance.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 5 — Arm A: Longitudinal LRT per tissue
# =============================================================================
# Tests whether gene expression trajectories DIVERGE between RML6 and NBH
# over time (condition × time interaction).
#
# Primary LRT (continuous time):
#   Full:    ~ cohort + condition + time_point_wpi + condition:time_point_wpi
#   Reduced: ~ cohort + condition + time_point_wpi
#   Rationale: Both models include condition main effect and shared time trend.
#   Only the interaction term (diverging trajectories) is dropped in the
#   reduced model. This is the correct symmetric formulation.
#
# Sensitivity LRT (factor time):
#   Full:    ~ cohort + condition + time_point_fac + condition:time_point_fac
#   Reduced: ~ cohort + condition + time_point_fac
#   Rationale: Same interaction-only test but allows non-monotonic patterns.
#   Run on a SEPARATE dds object (dds_fac) to preserve the primary dds.
#
# IMPORTANT: Terminal timepoint (time_point_wpi == TERMINAL_WPI) is EXCLUDED
# from both LRTs. Its exact numeric wpi is not available in GEO metadata.
# Including an approximate value in a continuous interaction term would distort
# the slope estimate. Terminal samples are used only in Arm B (Section 6).
# =============================================================================
cat("\n--- Section 5: Longitudinal LRT (Arm A) ---\n")

lrt_results_list  <- list()   # per-tissue LRT results data.frames
dds_tissue_list   <- list()   # per-tissue fitted DESeqDataSet (size factors for Arm C)
lrt_sig_counts    <- integer(length(TISSUES))
names(lrt_sig_counts) <- TISSUES
lrt_model_used    <- character(length(TISSUES))
names(lrt_model_used) <- TISSUES

for (tis in TISSUES) {
  cat(sprintf("\n  === %s ===\n", TISSUE_LABELS[tis]))

  # --- 5.1  Subset to this tissue, EXCLUDING terminal timepoint ---
  idx <- which(colData(gse)$tissue == tis &
                 colData(gse)$time_point_wpi != TERMINAL_WPI)
  if (length(idx) == 0) {
    cat(sprintf("    No non-terminal samples for %s; skipping LRT.\n",
                TISSUE_LABELS[tis]))
    lrt_results_list[[tis]] <- NULL
    next
  }
  gse_sub <- gse[, idx]
  cat(sprintf("    Samples for LRT (terminal excluded): %d\n", length(idx)))

  # --- 5.2  Pre-filter: >= 10 counts in >= MIN_REPS samples ---
  dds <- DESeqDataSet(gse_sub, design = ~ condition)
  keep <- rowSums(counts(dds) >= 10) >= MIN_REPS
  dds  <- dds[keep, ]
  cat(sprintf("    Genes after pre-filter: %d\n", nrow(dds)))

  # --- 5.3  Add covariates to colData ---
  colData(dds)$time_point_wpi <- as.numeric(colData(dds)$time_point_wpi)
  colData(dds)$time_point_fac <- factor(colData(dds)$time_point_wpi)
  colData(dds)$cohort         <- factor(colData(dds)$cohort,
                                        levels = c("main", "validation"))

  # Determine whether cohort is informative (> 1 level present)
  has_cohort <- length(unique(colData(dds)$cohort)) > 1

  # --- 5.4  Define full and reduced models ---
  if (has_cohort) {
    design_full    <- ~ cohort + condition + time_point_wpi + condition:time_point_wpi
    design_reduced <- ~ cohort + condition + time_point_wpi
    model_desc     <- "~ cohort + condition + wpi + condition:wpi"
  } else {
    design_full    <- ~ condition + time_point_wpi + condition:time_point_wpi
    design_reduced <- ~ condition + time_point_wpi
    model_desc     <- "~ condition + wpi + condition:wpi"
  }
  lrt_model_used[tis] <- model_desc

  # --- 5.5  Rank guard for full model ---
  X_full <- tryCatch(
    model.matrix(design_full, data = as.data.frame(colData(dds))),
    error = function(e) NULL
  )
  if (is.null(X_full) || Matrix::rankMatrix(X_full)[1] < ncol(X_full)) {
    warning(sprintf("  %s: full LRT design is rank-deficient; skipping LRT.\n",
                    TISSUE_LABELS[tis]))
    lrt_results_list[[tis]] <- NULL
    next
  }
  cat(sprintf("    LRT design rank: %d/%d (full rank)\n",
              Matrix::rankMatrix(X_full)[1], ncol(X_full)))
  cat(sprintf("    Model: %s\n", model_desc))

  # --- 5.6  Fit primary LRT (continuous time) ---
  design(dds) <- design_full
  set.seed(42)
  dds <- tryCatch(
    DESeq(dds, test = "LRT", reduced = design_reduced,
          fitType = "parametric", quiet = FALSE),
    error = function(e) {
      message(sprintf("  DESeq2 LRT failed for %s: %s", tis, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dds)) {
    lrt_results_list[[tis]] <- NULL
    next
  }

  # Store fitted dds for downstream use (size factors needed for Arm C)
  dds_tissue_list[[tis]] <- dds

  # --- 5.7  Extract primary LRT results ---
  res_lrt <- results(dds, independentFiltering = TRUE)
  res_lrt_df <- as.data.frame(res_lrt) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      ensembl_id  = stringr::str_remove(gene_id, "\\.\\d+$"),
      significant = !is.na(padj) & padj <= PADJ_THRESH,
      tissue      = tis
    ) %>%
    dplyr::arrange(padj)

  n_sig <- sum(res_lrt_df$significant, na.rm = TRUE)
  lrt_sig_counts[tis] <- n_sig
  cat(sprintf("    LRT significant genes (padj <= %.2f): %d\n",
              PADJ_THRESH, n_sig))

  # --- 5.8  Factor-based sensitivity LRT (SEPARATE dds object) ---
  # Run on dds_fac (a fresh copy) to preserve the primary dds above.
  # This avoids overwriting dds_tissue_list[[tis]] with the factor model.
  if (has_cohort) {
    design_full_fac    <- ~ cohort + condition + time_point_fac + condition:time_point_fac
    design_reduced_fac <- ~ cohort + condition + time_point_fac
  } else {
    design_full_fac    <- ~ condition + time_point_fac + condition:time_point_fac
    design_reduced_fac <- ~ condition + time_point_fac
  }

  X_fac <- tryCatch(
    model.matrix(design_full_fac, data = as.data.frame(colData(dds))),
    error = function(e) NULL
  )

  if (!is.null(X_fac) && Matrix::rankMatrix(X_fac)[1] == ncol(X_fac)) {
    dds_fac <- dds                       # copy preserves size factors + gene universe
    design(dds_fac) <- design_full_fac
    set.seed(42)
    dds_fac <- tryCatch(
      DESeq(dds_fac, test = "LRT", reduced = design_reduced_fac,
            fitType = "parametric", quiet = TRUE),
      error = function(e) NULL
    )
    if (!is.null(dds_fac)) {
      res_fac_df <- as.data.frame(results(dds_fac)) %>%
        tibble::rownames_to_column("gene_id") %>%
        dplyr::select(gene_id,
                      lrt_fac_pvalue = pvalue,
                      lrt_fac_padj   = padj)

      res_lrt_df <- res_lrt_df %>%
        dplyr::left_join(res_fac_df, by = "gene_id")

      n_sig_fac <- sum(!is.na(res_lrt_df$lrt_fac_padj) &
                         res_lrt_df$lrt_fac_padj <= PADJ_THRESH, na.rm = TRUE)
      n_overlap_fac <- sum(res_lrt_df$significant &
                             !is.na(res_lrt_df$lrt_fac_padj) &
                             res_lrt_df$lrt_fac_padj <= PADJ_THRESH, na.rm = TRUE)

      cat(sprintf("    Factor LRT significant genes: %d\n", n_sig_fac))
      cat(sprintf("    Overlap (continuous AND factor LRT): %d (%.0f%% of continuous)\n",
                  n_overlap_fac, 100 * n_overlap_fac / max(n_sig, 1)))

      # Jaccard similarity between continuous and factor LRT gene sets
      set_cont <- res_lrt_df$gene_id[res_lrt_df$significant]
      set_fac  <- res_lrt_df$gene_id[!is.na(res_lrt_df$lrt_fac_padj) &
                                       res_lrt_df$lrt_fac_padj <= PADJ_THRESH]
      jaccard_lrt <- if (length(union(set_cont, set_fac)) > 0)
        length(intersect(set_cont, set_fac)) / length(union(set_cont, set_fac))
      else NA_real_
      cat(sprintf("    Jaccard (continuous vs factor LRT): %.3f\n", jaccard_lrt))

      rm(dds_fac)
      gc(verbose = FALSE)
    }
  } else {
    cat("    Factor LRT: design rank-deficient; skipping.\n")
  }

  lrt_results_list[[tis]] <- res_lrt_df

  rm(gse_sub, X_full)
  gc(verbose = FALSE)
}

cat("\n  LRT summary:\n")
for (tis in TISSUES) {
  cat(sprintf("    %-20s %d significant genes  |  model: %s\n",
              TISSUE_LABELS[tis], lrt_sig_counts[tis],
              lrt_model_used[tis]))
}


# ----------------------------------------------------------------------
# ##Section 6 — Arm B: Terminal Wald Test
#
# **Purpose:** Identify genes differentially expressed at the terminal time point (150 dpi) using a standard Wald test with `~condition` design on the 150 dpi subset only. This captures endpoint-specific changes that may not show a significant trajectory in the LRT.
#
# **Inputs:** `se_filtered` subset to 150 dpi samples.
#
# **Outputs:** `../results/stream3/wald_terminal_results.tsv`.
#
# **Key decisions:** Terminal Wald results are used in Section 8 for cross-tissue concordance with Stream 1 (which also uses endpoint comparisons). LFC shrinkage with `apeglm` is applied.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 6 — Arm B: Terminal-timepoint Wald test per tissue
# =============================================================================
# Standard DESeq2 Wald test at terminal timepoint only.
# Design: ~ cohort + condition (cohort included if both levels present).
# apeglm LFC shrinkage applied for publication-quality effect size estimates.
#
# CRITICAL DESIGN DECISION:
#   Size factors are estimated from the FULL tissue dataset (all timepoints)
#   and transferred to the terminal subset. This ensures consistent
#   normalisation with the LRT analysis and avoids re-estimating size factors
#   on a small subset (which can be unstable with N=6-12 samples).
#
# Gene universe: same pre-filter as LRT (>= 10 counts in >= MIN_REPS samples
# across the FULL tissue dataset), ensuring consistent gene sets.
# =============================================================================
cat("\n--- Section 6: Terminal Wald test (Arm B) ---\n")

terminal_results_list <- list()
terminal_sig_counts   <- integer(length(TISSUES))
names(terminal_sig_counts) <- TISSUES

for (tis in TISSUES) {
  cat(sprintf("\n  === %s ===\n", TISSUE_LABELS[tis]))

  # --- 6.1  Build full-tissue DESeqDataSet (all timepoints) for size factors ---
  idx_full <- which(colData(gse)$tissue == tis)
  gse_full <- gse[, idx_full]

  dds_full <- DESeqDataSet(gse_full, design = ~ condition)
  keep_full <- rowSums(counts(dds_full) >= 10) >= MIN_REPS
  dds_full  <- dds_full[keep_full, ]

  set.seed(42)
  dds_full <- tryCatch(
    estimateSizeFactors(dds_full),
    error = function(e) {
      message(sprintf("  estimateSizeFactors failed for %s: %s", tis, conditionMessage(e)))
      NULL
    }
  )
  if (is.null(dds_full)) {
    cat(sprintf("    Size factor estimation failed for %s; skipping.\n",
                TISSUE_LABELS[tis]))
    terminal_results_list[[tis]] <- NULL
    next
  }

  # ALWAYS overwrite dds_tissue_list with the full-tissue dds (all timepoints,
  # including terminal). The LRT dds from Section 5 EXCLUDES terminal samples
  # (it was built with time_point_wpi != TERMINAL_WPI). If we kept the LRT dds
  # here, run_wald_at_wpi() in Section 9 would find zero terminal samples and
  # silently return NULL for the terminal timepoint. The full-tissue dds has
  # size factors estimated from all timepoints and includes terminal samples.
  dds_tissue_list[[tis]] <- dds_full

  # --- 6.2  Subset to terminal timepoint ---
  idx_term <- which(colData(dds_full)$time_point_wpi == TERMINAL_WPI)
  if (length(idx_term) < (2 * MIN_REPS)) {
    cat(sprintf("    Insufficient samples at terminal (N=%d; need >= %d); skipping.\n",
                length(idx_term), 2 * MIN_REPS))
    terminal_results_list[[tis]] <- NULL
    next
  }

  dds_term <- dds_full[, idx_term]
  n_prion_term   <- sum(dds_term$condition == "prion")
  n_control_term <- sum(dds_term$condition == "control")
  cat(sprintf("    Terminal samples: %d (%d prion, %d control)\n",
              length(idx_term), n_prion_term, n_control_term))
  cat(sprintf("    Genes in universe: %d\n", nrow(dds_term)))

  # --- 6.3  Choose design ---
  colData(dds_term)$cohort <- factor(colData(dds_term)$cohort,
                                     levels = c("main", "validation"))
  has_cohort_term <- length(unique(colData(dds_term)$cohort)) > 1
  des_term <- if (has_cohort_term) ~ cohort + condition else ~ condition

  X_term <- tryCatch(
    model.matrix(des_term, data = as.data.frame(colData(dds_term))),
    error = function(e) NULL
  )
  if (is.null(X_term) || Matrix::rankMatrix(X_term)[1] < ncol(X_term)) {
    des_term <- ~ condition
    cat("    Cohort term rank-deficient; falling back to ~ condition.\n")
  }
  design(dds_term) <- des_term
  cat(sprintf("    Terminal model: %s\n", deparse(des_term, width.cutoff = 60)))

  # --- 6.4  Estimate dispersions and run Wald test ---
  # Size factors already set (transferred from full-tissue dds).
  set.seed(42)
  dds_term <- tryCatch({
    dds_term <- estimateDispersions(dds_term, fitType = "parametric", quiet = TRUE)
    nbinomWaldTest(dds_term, quiet = TRUE)
  }, error = function(e) {
    message(sprintf("  Terminal Wald failed for %s: %s", tis, conditionMessage(e)))
    NULL
  })
  if (is.null(dds_term)) {
    terminal_results_list[[tis]] <- NULL
    next
  }

  # --- 6.5  apeglm LFC shrinkage ---
  res_term_shrunk <- tryCatch(
    lfcShrink(dds_term, coef = "condition_prion_vs_control",
              type = "apeglm", quiet = TRUE),
    error = function(e) {
      message(sprintf("  apeglm failed for %s terminal; using unshrunken.", tis))
      results(dds_term, name = "condition_prion_vs_control")
    }
  )

  # Unshrunken Wald stats (for volcano plot y-axis: use unshrunken pvalue)
  res_term_wald <- results(dds_term, name = "condition_prion_vs_control",
                           independentFiltering = TRUE)

  # --- 6.6  Build results data.frame ---
  res_term_df <- as.data.frame(res_term_shrunk) %>%
    tibble::rownames_to_column("gene_id") %>%
    dplyr::mutate(
      ensembl_id  = stringr::str_remove(gene_id, "\\.\\d+$"),
      significant = !is.na(padj) & padj <= PADJ_THRESH &
                    abs(log2FoldChange) >= LFC_THRESH,
      direction   = dplyr::case_when(
        significant & log2FoldChange > 0 ~ "up",
        significant & log2FoldChange < 0 ~ "down",
        TRUE                             ~ "ns"
      ),
      tissue = tis
    ) %>%
    dplyr::left_join(
      as.data.frame(res_term_wald) %>%
        tibble::rownames_to_column("gene_id") %>%
        dplyr::select(gene_id,
                      pval_wald = pvalue,
                      lfc_wald  = log2FoldChange),
      by = "gene_id"
    ) %>%
    dplyr::arrange(padj)

  n_sig_term <- sum(res_term_df$significant, na.rm = TRUE)
  n_up_term  <- sum(res_term_df$direction == "up",   na.rm = TRUE)
  n_dn_term  <- sum(res_term_df$direction == "down", na.rm = TRUE)
  terminal_sig_counts[tis] <- n_sig_term

  cat(sprintf("    Terminal DEGs (padj<=%.2f, |LFC|>=%.1f): %d (%d up, %d down)\n",
              PADJ_THRESH, LFC_THRESH, n_sig_term, n_up_term, n_dn_term))

  terminal_results_list[[tis]] <- res_term_df

  rm(gse_full, dds_full, dds_term, X_term)
  gc(verbose = FALSE)
}

cat("\n  Terminal Wald summary:\n")
for (tis in TISSUES) {
  res <- terminal_results_list[[tis]]
  if (is.null(res)) {
    cat(sprintf("    %-20s SKIPPED\n", TISSUE_LABELS[tis]))
  } else {
    n_up <- sum(res$direction == "up",   na.rm = TRUE)
    n_dn <- sum(res$direction == "down", na.rm = TRUE)
    cat(sprintf("    %-20s %d DEGs (%d up, %d down)\n",
                TISSUE_LABELS[tis], terminal_sig_counts[tis], n_up, n_dn))
  }
}

# --- 6.7  Cross-tissue terminal DEG overlap (Jaccard similarity) ---
# Quantifies how tissue-specific the terminal transcriptional response is.
terminal_deg_sets <- lapply(terminal_results_list, function(r) {
  if (is.null(r)) return(character(0))
  r$gene_id[r$significant]
})

cross_tissue_overlap <- data.frame()
tissue_pairs <- combn(TISSUES, 2, simplify = FALSE)
for (pair in tissue_pairs) {
  t1 <- pair[1]; t2 <- pair[2]
  s1 <- terminal_deg_sets[[t1]]; s2 <- terminal_deg_sets[[t2]]
  n_shared <- length(intersect(s1, s2))
  n_union  <- length(union(s1, s2))
  jaccard  <- if (n_union > 0) n_shared / n_union else NA_real_
  cross_tissue_overlap <- dplyr::bind_rows(cross_tissue_overlap, data.frame(
    tissue_1  = TISSUE_LABELS[t1],
    tissue_2  = TISSUE_LABELS[t2],
    n_degs_t1 = length(s1),
    n_degs_t2 = length(s2),
    n_shared  = n_shared,
    jaccard   = round(jaccard, 3),
    stringsAsFactors = FALSE
  ))
}

cat("\n  Cross-tissue terminal DEG overlap (Jaccard similarity):\n")
print(cross_tissue_overlap)


# ----------------------------------------------------------------------
# ##Section 7 — Gene Annotation
#
# **Purpose:** Annotate all DE results with gene symbols, Entrez IDs, and gene descriptions using biomaRt. Also flags MCPS genes and resolves aliases for cross-stream harmonisation.
#
# **Inputs:** `../results/stream3/lrt_results.tsv`; `../results/stream3/wald_terminal_results.tsv`; biomaRt Ensembl database.
#
# **Outputs:** Annotated versions of both result tables; `../results/stream3/annotation_map.tsv`.
#
# **Key decisions:** Genes without biomaRt annotation are retained with Ensembl ID as symbol — they are not dropped, as some may be novel lncRNAs.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 7 — Gene symbol annotation (biomaRt, with local RDS cache)
# =============================================================================
# Maps Ensembl gene IDs to mouse MGI gene symbols using Ensembl BioMart.
# Results are cached locally (ANNOT_CACHE) to avoid repeated network calls
# and to ensure reproducibility if the Ensembl server is unavailable.
# A versioned Ensembl archive host is tried first (Nov 2022 release, matching
# the GEO submission date); falls back to the current release if unavailable.
#
# Annotation includes: gene_symbol (MGI), gene_biotype, chromosome, coordinates.
# Biotype is used to annotate concordance tables (protein_coding vs lncRNA etc.)
# =============================================================================
cat("\n--- Section 7: Gene symbol annotation (biomaRt) ---\n")

# Collect all unique Ensembl IDs across all tissues and analyses
all_ensembl_ids <- unique(c(
  unlist(lapply(lrt_results_list,      function(r) if (!is.null(r)) r$ensembl_id)),
  unlist(lapply(terminal_results_list, function(r) if (!is.null(r)) r$ensembl_id))
))
all_ensembl_ids <- all_ensembl_ids[!is.na(all_ensembl_ids) & all_ensembl_ids != ""]
cat(sprintf("  Unique Ensembl IDs to annotate: %d\n", length(all_ensembl_ids)))

ANNOT_CACHE <- file.path(BASE_DIR, "resources", "cache", "biomart", "GSE210128_gene_annotation_cache.rds")

# Safely create the biomart cache directory if it doesn't exist
annot_cache_dir <- dirname(ANNOT_CACHE)
if (!dir.exists(annot_cache_dir)) {
  dir.create(annot_cache_dir, recursive = TRUE, showWarnings = FALSE)
}

if (file.exists(ANNOT_CACHE)) {
  cat("  Loading gene annotation from local cache.\n")
  gene_annot <- readRDS(ANNOT_CACHE)
  cat(sprintf("  Cache: %d genes annotated\n", nrow(gene_annot)))
} else {
  cat("  Querying Ensembl BioMart for mouse gene symbols...\n")
  library(biomaRt)

  # Try versioned archive first (Nov 2022, matching GEO submission date)
  mouse_mart <- tryCatch(
    useMart("ensembl", dataset = "mmusculus_gene_ensembl",
            host = "https://nov2022.archive.ensembl.org"),
    error = function(e) {
      message("  Versioned host unavailable; trying current release.")
      tryCatch(
        useMart("ensembl", dataset = "mmusculus_gene_ensembl",
                host = "https://www.ensembl.org"),
        error = function(e2) {
          stop(paste("BioMart unavailable:", conditionMessage(e2)))
        }
      )
    }
  )

  gene_annot <- getBM(
    attributes = c("ensembl_gene_id", "mgi_symbol", "gene_biotype",
                   "chromosome_name", "start_position", "end_position"),
    filters    = "ensembl_gene_id",
    values     = all_ensembl_ids,
    mart       = mouse_mart
  ) %>%
    dplyr::rename(
      ensembl_id  = ensembl_gene_id,
      gene_symbol = mgi_symbol
    ) %>%
    # Remove entries with empty gene symbols
    dplyr::filter(gene_symbol != "") %>%
    # Keep one row per Ensembl ID (prefer protein_coding if duplicated)
    dplyr::group_by(ensembl_id) %>%
    dplyr::arrange(dplyr::desc(gene_biotype == "protein_coding")) %>%
    dplyr::slice_head(n = 1) %>%
    dplyr::ungroup()

  saveRDS(gene_annot, ANNOT_CACHE)
  cat(sprintf("  Annotation cache saved: %d genes\n", nrow(gene_annot)))
}

# --- 7.1  Validate symbol format consistency with MCPS ---
# MCPS uses MGI symbols; gene_annot should also use MGI symbols.
# Check for format mismatches (e.g., all-caps human symbols vs mixed-case mouse).
mcps_symbols_upper <- sum(MCPS$gene_symbol == toupper(MCPS$gene_symbol), na.rm = TRUE)
annot_symbols_upper <- sum(gene_annot$gene_symbol == toupper(gene_annot$gene_symbol),
                           na.rm = TRUE)
if (mcps_symbols_upper > nrow(MCPS) * 0.5 &&
    annot_symbols_upper < nrow(gene_annot) * 0.5) {
  warning(paste(
    "MCPS gene symbols appear to be in ALL-CAPS (human format) but",
    "gene_annot uses mixed-case (mouse MGI format).",
    "MCPS concordance matching may fail.",
    "Check that Stream 1 outputs mouse MGI symbols, not human HGNC symbols."
  ))
}

# --- 7.2  Annotate LRT and terminal results ---
lrt_results_list <- lapply(lrt_results_list, function(r) {
  if (is.null(r)) return(NULL)
  r %>% dplyr::left_join(
    gene_annot %>% dplyr::select(ensembl_id, gene_symbol, gene_biotype),
    by = "ensembl_id"
  )
})

terminal_results_list <- lapply(terminal_results_list, function(r) {
  if (is.null(r)) return(NULL)
  r %>% dplyr::left_join(
    gene_annot %>% dplyr::select(ensembl_id, gene_symbol, gene_biotype),
    by = "ensembl_id"
  )
})

# Report annotation coverage
for (tis in TISSUES) {
  r <- terminal_results_list[[tis]]
  if (is.null(r)) next
  n_annotated <- sum(!is.na(r$gene_symbol) & r$gene_symbol != "", na.rm = TRUE)
  cat(sprintf("  %s: %d / %d genes annotated (%.0f%%)\n",
              TISSUE_LABELS[tis], n_annotated, nrow(r),
              100 * n_annotated / nrow(r)))
}


# ----------------------------------------------------------------------
# ##Sections 8 & 9 — Arm C: Cross-Tissue Concordance & First-Detection
#
# **Purpose:** Section 8: Test whether Stream 3 peripheral DE genes are concordant with Stream 1 CNS meta-analysis results (direction-of-effect test). Section 9: For MCPS genes, determine the earliest time point at which each gene is significantly dysregulated in peripheral blood (first-detection analysis).
#
# **Inputs:** `../results/stream3/wald_terminal_results.tsv`; `../results/stream1/meta_results.tsv`.
#
# **Outputs:** `../results/stream3/concordance_results.tsv`; `../results/stream3/first_detection.tsv`.
#
# **Key decisions:** First-detection uses per-time-point Wald tests with Bonferroni correction across time points (not BH) to control FWER for the earliest-detection claim.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 8 — Arm C: MCPS concordance at terminal timepoint (per tissue)
# =============================================================================
# For each tissue, tests whether MCPS genes are concordantly dysregulated
# at the terminal timepoint (same direction as CNS meta-analysis).
#
# Metrics:
#   (a) Direction concordance: proportion of testable MCPS genes with
#       sign(peripheral LFC) == sign(MCPS meta-LFC).
#   (b) Binomial test: H0 = concordance = 50% (random), one-sided (greater).
#   (c) LFC Pearson correlation: quantifies linear agreement between
#       peripheral terminal LFC and MCPS meta-LFC.
#
# Denominator: nrow(MCPS) (all MCPS genes), consistent with Stream 1 and 2.
# Deduplication: if a gene maps to multiple Ensembl IDs, keep the one with
# the lowest terminal padj (most significant result).
# =============================================================================
cat("\n--- Section 8: MCPS concordance at terminal (Arm C) ---\n")

concordance_list    <- list()
concordance_summary <- data.frame()

for (tis in TISSUES) {
  cat(sprintf("\n  === %s ===\n", TISSUE_LABELS[tis]))

  res_term <- terminal_results_list[[tis]]
  if (is.null(res_term)) {
    cat("    No terminal results; skipping MCPS concordance.\n")
    next
  }

  # --- 8.1  Join terminal results with MCPS by gene_symbol ---
  concordance <- MCPS %>%
    dplyr::rename(mcps_lfc = meta_log2FC) %>%
    dplyr::select(gene_symbol, mcps_lfc,
                  dplyr::any_of(c("meta_padj", "I2", "n_studies"))) %>%
    dplyr::inner_join(
      res_term %>%
        dplyr::filter(!is.na(gene_symbol), gene_symbol != "") %>%
        dplyr::select(gene_id, gene_symbol, ensembl_id, gene_biotype,
                      terminal_lfc  = log2FoldChange,
                      terminal_padj = padj,
                      terminal_sig  = significant,
                      direction),
      by = "gene_symbol"
    ) %>%
    # Deduplicate: if multiple Ensembl IDs per symbol, keep lowest padj
    dplyr::group_by(gene_symbol) %>%
    dplyr::slice_min(order_by = terminal_padj, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      direction_match = !is.na(mcps_lfc) & !is.na(terminal_lfc) &
                        sign(mcps_lfc) == sign(terminal_lfc),
      tissue = tis
    )

  n_tested      <- nrow(concordance)
  n_dir_match   <- sum(concordance$direction_match, na.rm = TRUE)
  n_term_sig    <- sum(concordance$terminal_sig, na.rm = TRUE)
  pct_of_mcps   <- 100 * n_dir_match / n_mcps
  pct_of_tested <- 100 * n_dir_match / max(n_tested, 1)

  # Binomial test: H0 = concordance = 50%
  binom_test <- binom.test(n_dir_match, n_tested, p = 0.5,
                           alternative = "greater")

  # LFC Pearson correlation
  lfc_cor <- tryCatch(
    cor.test(concordance$mcps_lfc, concordance$terminal_lfc,
             method = "pearson", use = "complete.obs"),
    error = function(e) list(estimate = NA_real_, p.value = NA_real_,
                             conf.int = c(NA_real_, NA_real_))
  )

  cat(sprintf("    MCPS genes testable:     %d / %d (%.0f%%)\n",
              n_tested, n_mcps, 100 * n_tested / n_mcps))
  cat(sprintf("    Direction concordance:   %d/%d tested (%.1f%% of full MCPS)\n",
              n_dir_match, n_tested, pct_of_mcps))
  cat(sprintf("    Binomial p (H0=50%%):     %.4f\n", binom_test$p.value))
  cat(sprintf("    LFC Pearson r:           %.3f (p=%.4f)\n",
              lfc_cor$estimate, lfc_cor$p.value))
  cat(sprintf("    Terminal-sig MCPS genes: %d\n", n_term_sig))

  concordance_list[[tis]] <- concordance

  concordance_summary <- dplyr::bind_rows(concordance_summary, data.frame(
    tissue            = tis,
    tissue_label      = TISSUE_LABELS[tis],
    n_mcps_total      = n_mcps,
    n_tested          = n_tested,
    n_dir_match       = n_dir_match,
    pct_of_mcps       = round(pct_of_mcps, 1),
    pct_of_tested     = round(pct_of_tested, 1),
    binom_pvalue      = binom_test$p.value,
    binom_ci_lower    = binom_test$conf.int[1],
    lfc_pearson_r     = round(lfc_cor$estimate, 3),
    lfc_pearson_p     = lfc_cor$p.value,
    n_terminal_sig    = n_term_sig,
    stringsAsFactors  = FALSE
  ))
}

cat("\n  MCPS concordance summary:\n")
if (nrow(concordance_summary) > 0) {
  print(concordance_summary[, c("tissue_label", "n_tested", "n_dir_match",
                                 "pct_of_mcps", "binom_pvalue", "lfc_pearson_r")])
}

# =============================================================================
# SECTION 9 — MCPS first-detection: per-timepoint Wald tests
# =============================================================================
# For each tissue, runs per-timepoint Wald tests (RML6 vs NBH) and identifies
# the earliest wpi at which >= 1 MCPS gene reaches padj <= PADJ_THRESH.
# This provides a timeline of when peripheral prion signal first emerges.
#
# CRITICAL DESIGN DECISION:
#   Size factors are estimated from the FULL tissue dataset (all timepoints)
#   and transferred to each per-timepoint subset via run_wald_at_wpi().
#   This ensures consistent normalisation across timepoints.
#   Dispersions are re-estimated within each subset (unavoidable with small N).
#
# The per-timepoint Wald test uses the same gene universe as the full-tissue
# analysis (pre-filtered genes from dds_tissue_list[[tis]]).
# =============================================================================
cat("\n--- Section 9: MCPS first-detection analysis ---\n")

first_detection_list  <- list()
per_timepoint_wald_all <- data.frame()   # all per-timepoint results (supplementary)

for (tis in TISSUES) {
  cat(sprintf("\n  === %s ===\n", TISSUE_LABELS[tis]))

  # Use the full-tissue dds (with size factors already estimated)
  dds_tis <- dds_tissue_list[[tis]]
  if (is.null(dds_tis)) {
    cat(sprintf("    No full-tissue dds for %s; skipping first-detection.\n",
                TISSUE_LABELS[tis]))
    first_detection_list[[tis]] <- list(
      detection_df  = data.frame(),
      first_det_wpi = NA_integer_,
      per_wpi_res   = list()
    )
    next
  }

  # Add cohort factor if not already present
  if (!"cohort" %in% colnames(colData(dds_tis))) {
    colData(dds_tis)$cohort <- factor(colData(dds_tis)$cohort,
                                      levels = c("main", "validation"))
  }

  wpi_vals <- sort(unique(colData(dds_tis)$time_point_wpi))
  cat(sprintf("    Timepoints to test: %s\n",
              paste(wpi_vals, collapse = ", ")))

  # --- 9.1  Run per-timepoint Wald tests ---
  # run_wald_at_wpi() transfers size factors from dds_tis (full tissue model)
  per_wpi_results <- lapply(wpi_vals, function(w) {
    run_wald_at_wpi(dds_tis, w)
  })
  names(per_wpi_results) <- as.character(wpi_vals)

  # --- 9.2  Annotate per-timepoint results with gene symbols ---
  per_wpi_annotated <- lapply(names(per_wpi_results), function(w_chr) {
    res <- per_wpi_results[[w_chr]]
    if (is.null(res)) return(NULL)
    res %>%
      dplyr::mutate(ensembl_id = stringr::str_remove(gene_id, "\\.\\d+$")) %>%
      dplyr::left_join(gene_annot %>% dplyr::select(ensembl_id, gene_symbol, gene_biotype),
                by = "ensembl_id") %>%
      dplyr::mutate(tissue = tis)
  })
  names(per_wpi_annotated) <- names(per_wpi_results)

  # Collect all per-timepoint results for supplementary table
  per_timepoint_wald_all <- dplyr::bind_rows(
    per_timepoint_wald_all,
    dplyr::bind_rows(per_wpi_annotated)
  )

  # --- 9.3  Count MCPS genes significant at each timepoint ---
  mcps_symbols <- MCPS$gene_symbol

  detection_df <- lapply(wpi_vals, function(w) {
    res <- per_wpi_annotated[[as.character(w)]]
    if (is.null(res)) {
      return(data.frame(wpi = w, n_mcps_sig = 0L, n_mcps_tested = 0L,
                        tissue = tis, tissue_label = TISSUE_LABELS[tis]))
    }
    n_tested_mcps <- sum(res$gene_symbol %in% mcps_symbols, na.rm = TRUE)
    n_sig_mcps    <- sum(
      res$gene_symbol %in% mcps_symbols &
        !is.na(res$padj) & res$padj <= PADJ_THRESH,
      na.rm = TRUE
    )
    data.frame(wpi = w, n_mcps_sig = n_sig_mcps,
               n_mcps_tested = n_tested_mcps,
               tissue = tis, tissue_label = TISSUE_LABELS[tis])
  }) %>% dplyr::bind_rows()

  # --- 9.4  First detection = earliest wpi with >= 1 MCPS gene significant ---
  first_det_wpi <- detection_df %>%
    dplyr::filter(n_mcps_sig >= 1) %>%
    dplyr::slice_min(order_by = wpi, n = 1) %>%
    dplyr::pull(wpi)
  first_det_wpi <- if (length(first_det_wpi) == 0) NA_integer_ else first_det_wpi[1]

  if (is.na(first_det_wpi)) {
    cat(sprintf("    First MCPS detection: NOT DETECTED within tested timepoints\n"))
  } else {
    wpi_label <- WPI_LABELS[match(first_det_wpi, WPI_LEVELS)]
    cat(sprintf("    First MCPS detection: %s wpi\n",
                ifelse(is.na(wpi_label), as.character(first_det_wpi), wpi_label)))
  }

  detection_df$first_detection <- !is.na(first_det_wpi) &
                                    detection_df$wpi == first_det_wpi

  # --- 9.5  Gene-level first-detection table ---
  # Which specific MCPS genes are first detected, and at what timepoint?
  mcps_first_genes <- lapply(wpi_vals, function(w) {
    res <- per_wpi_annotated[[as.character(w)]]
    if (is.null(res)) return(NULL)
    res %>%
      dplyr::filter(gene_symbol %in% mcps_symbols,
             !is.na(padj), padj <= PADJ_THRESH) %>%
      dplyr::mutate(first_detected_wpi = w) %>%
      dplyr::select(gene_symbol, first_detected_wpi,
                    log2FoldChange, padj, tissue)
  }) %>%
    dplyr::bind_rows() %>%
    # Keep only the FIRST timepoint each gene is detected
    dplyr::group_by(gene_symbol) %>%
    dplyr::slice_min(order_by = first_detected_wpi, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::left_join(MCPS %>% dplyr::select(gene_symbol, mcps_lfc = meta_log2FC),
              by = "gene_symbol") %>%
    dplyr::mutate(
      direction_match = !is.na(mcps_lfc) & !is.na(log2FoldChange) &
                        sign(mcps_lfc) == sign(log2FoldChange),
      wpi_label = WPI_LABELS[match(first_detected_wpi, WPI_LEVELS)]
    ) %>%
    dplyr::arrange(first_detected_wpi, padj)

  first_detection_list[[tis]] <- list(
    detection_df      = detection_df,
    first_det_wpi     = first_det_wpi,
    per_wpi_res       = per_wpi_results,
    mcps_first_genes  = mcps_first_genes
  )
}

# --- 9.6  Combine detection data frames ---
first_detection_combined <- dplyr::bind_rows(
  lapply(first_detection_list, function(x) x$detection_df)
) %>%
  dplyr::mutate(
    wpi_label = factor(wpi, levels = WPI_LEVELS, labels = WPI_LABELS),
    tissue    = factor(tissue, levels = TISSUES)
  )

cat("\n  First MCPS detection summary:\n")
for (tis in TISSUES) {
  fd  <- first_detection_list[[tis]]$first_det_wpi
  lbl <- if (!is.na(fd)) {
    wl <- WPI_LABELS[match(fd, WPI_LEVELS)]
    sprintf("%s wpi", ifelse(is.na(wl), as.character(fd), wl))
  } else "not detected"
  cat(sprintf("    %-20s %s\n", TISSUE_LABELS[tis], lbl))
}

# =============================================================================
# SECTION 9.7 — Cross-stream integration: Stream 2 leading-edge cross-reference
# =============================================================================
# Identifies MCPS genes that satisfy ALL THREE of the following criteria:
#   (1) In the MCPS (CNS meta-analysis, Stream 1)
#   (2) First detected in at least one peripheral mouse tissue (Stream 3, above)
#   (3) In the GSEA leading edge in human sCJD (Stream 2, GSE214373)
#
# These triple-validated genes are the highest-priority cross-species peripheral
# biomarker candidates in the entire pipeline. A gene that is dysregulated in
# mouse CNS, detectable in mouse peripheral tissues early in disease, AND
# enriched in the human sCJD transcriptome represents a convergent signal
# across species, tissues, and disease stages.
#
# GRACEFUL FALLBACK: If Stream 2 outputs are not yet available (i.e., the
# pipeline is being run for the first time and Stream 2 has not completed),
# this section is skipped with a warning. The cross-reference can be re-run
# after Stream 2 completes by re-running this script.
#
# INPUT:  ../outputs/stream2/GSE214373_GSEA_MCPS.csv
#         Column 'leading_edge_genes': comma-separated mouse MGI symbols
# OUTPUT: GSE210128_triple_validated_candidates.csv
# =============================================================================
cat("\n--- Section 9.7: Stream 2 leading-edge cross-reference ---\n")

# ---------------------------------------------------------------------------
# 9.7.0  Build mcps_first_genes_all (needed here AND in Section 11.6)
# ---------------------------------------------------------------------------
# Aggregate per-tissue first-detection gene tables into a single data frame.
# This object is constructed here (before Section 9.7 uses it) and reused
# in Section 11.6 when writing the output CSV.
mcps_first_genes_all <- dplyr::bind_rows(
  lapply(TISSUES, function(tis) {
    fg <- first_detection_list[[tis]]$mcps_first_genes
    if (is.null(fg) || nrow(fg) == 0) return(NULL)
    fg
  })
)
cat(sprintf("  mcps_first_genes_all: %d rows across %d tissues\n",
            nrow(mcps_first_genes_all), length(TISSUES)))

# Initialise output objects (populated below if Stream 2 is available)
triple_validated      <- data.frame()
stream2_leading_genes <- character(0)
stream2_available     <- FALSE

STREAM2_DIR   <- file.path(BASE_DIR, "results", "stream2")
gsea_csv_path <- file.path(STREAM2_DIR, "GSE214373_GSEA_MCPS.csv")

if (!file.exists(gsea_csv_path)) {
  cat(sprintf(
    "  Stream 2 GSEA output not found at: %s\n",
    gsea_csv_path
  ))
  cat("  Skipping cross-reference. Re-run after 06_stream2_human_validation.R completes.\n")
} else {
  stream2_available <- TRUE
  cat(sprintf("  Loading Stream 2 GSEA results from: %s\n", gsea_csv_path))

  gsea_stream2 <- tryCatch(
    readr::read_csv(gsea_csv_path, show_col_types = FALSE),
    error = function(e) {
      warning(sprintf("Failed to read Stream 2 GSEA CSV: %s", conditionMessage(e)))
      NULL
    }
  )

  if (!is.null(gsea_stream2) && "leadingEdge" %in% colnames(gsea_stream2)) {

    # --- 9.7.1  Parse leading-edge gene list ---
    # The leading_edge_genes column is a comma-separated string of mouse MGI symbols.
    # These are the MCPS genes that contribute most to the enrichment score in
    # human sCJD (i.e., they are both in the MCPS and significantly ranked in
    # the human DEG list).
    leading_edge_raw <- gsea_stream2$leadingEdge[1]
    # Script 06 writes leadingEdge as semicolon-separated string (fgsea default).
    stream2_leading_genes <- stringr::str_split(leading_edge_raw, ";")[[1]] %>%
      stringr::str_trim() %>%
      .[. != ""]

    cat(sprintf("  Stream 2 GSEA leading-edge genes: %d\n",
                length(stream2_leading_genes)))
    cat(sprintf("  Stream 2 GSEA NES: %.3f  |  padj: %.4f\n",
                gsea_stream2$NES[1], gsea_stream2$padj[1]))

    # --- 9.7.2  Cross-reference with peripheral first-detection genes ---
    # mcps_first_genes_all: MCPS genes first detected in peripheral tissues
    # (built in Section 9.5, one row per gene per tissue)
    if (nrow(mcps_first_genes_all) > 0) {

      triple_validated <- mcps_first_genes_all %>%
        dplyr::filter(gene_symbol %in% stream2_leading_genes) %>%
        dplyr::left_join(
          # Add Stream 2 human LFC if available
          tryCatch({
            s2_conc_path <- file.path(STREAM2_DIR, "GSE214373_MCPS_concordance.csv")
            if (file.exists(s2_conc_path)) {
              readr::read_csv(s2_conc_path, show_col_types = FALSE) %>%
                dplyr::select(gene_symbol = mouse_symbol,
                              human_lfc   = human_lfc,
                              human_padj  = human_padj) %>%
                dplyr::distinct(gene_symbol, .keep_all = TRUE)
            } else {
              data.frame(gene_symbol = character(0),
                         human_lfc   = numeric(0),
                         human_padj  = numeric(0))
            }
          }, error = function(e) {
            data.frame(gene_symbol = character(0),
                       human_lfc   = numeric(0),
                       human_padj  = numeric(0))
          }),
          by = "gene_symbol"
        ) %>%
        dplyr::mutate(
          # Summarise evidence tier:
          # Tier 1 = concordant in peripheral tissue AND in human leading edge
          # Tier 2 = detected in peripheral tissue AND in human leading edge
          #          (direction may differ from MCPS)
          evidence_tier = dplyr::case_when(
            direction_match ~ "Tier 1: concordant peripheral + human leading edge",
            TRUE            ~ "Tier 2: detected peripheral + human leading edge"
          )
        ) %>%
        dplyr::arrange(evidence_tier, first_detected_wpi, padj)

      n_triple <- length(unique(triple_validated$gene_symbol))
      n_tier1  <- sum(triple_validated$evidence_tier ==
                        "Tier 1: concordant peripheral + human leading edge" &
                        !duplicated(triple_validated$gene_symbol))

      cat(sprintf("  Triple-validated candidates: %d unique genes\n", n_triple))
      cat(sprintf("    Tier 1 (concordant peripheral + human leading edge): %d\n",
                  n_tier1))
      cat(sprintf("    Tier 2 (detected peripheral + human leading edge):   %d\n",
                  n_triple - n_tier1))

      if (n_triple > 0) {
        cat("  Top Tier 1 candidates:\n")
        tier1_top <- triple_validated %>%
          dplyr::filter(evidence_tier == "Tier 1: concordant peripheral + human leading edge") %>%
          dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
          dplyr::slice_head(n = 10)
        for (i in seq_len(nrow(tier1_top))) {
          cat(sprintf("    %-15s  tissue=%-15s  wpi=%s  padj=%.4f\n",
                      tier1_top$gene_symbol[i],
                      tier1_top$tissue[i],
                      tier1_top$wpi_label[i],
                      tier1_top$padj[i]))
        }
      }

    } else {
      cat("  No peripheral first-detection genes available; skipping cross-reference.\n")
    }

  } else {
    cat("  Stream 2 GSEA CSV missing 'leadingEdge' column; skipping.\n")
  }
}


# ----------------------------------------------------------------------
# ##Section 10 — Figures 1–10
#
# **Purpose:** Generate the main publication figures for Stream 3: trajectory plots for top LRT genes, volcano plot for terminal Wald, MCPS gene trajectory heatmap, concordance scatter plot, and first-detection timeline.
#
# **Inputs:** All Stream 3 result tables.
#
# **Outputs:** `../results/stream3/figures/fig_stream3_01.png` through `fig_stream3_10.png` (PNG + SVG).
#
# **Key decisions:** Trajectory plots show mean ± SEM per time point with individual sample points overlaid. MCPS genes are highlighted in red.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 10 — Publication-ready figures
# =============================================================================
# All figures use theme_pub (defined in Section 0).
# Saved as both PDF (vector, for journal submission) and PNG (for inspection).
# Figure dimensions follow typical journal guidelines:
#   Single column: 3.5 in wide
#   Double column: 7.0 in wide
#   Max height:    9.0 in
# ComplexHeatmap figures use save_heatmap() (base graphics, not ggsave).
# =============================================================================
cat("\n--- Section 10: Generating figures ---\n")

# ---------------------------------------------------------------------------
# Fig 1: Sample overview — N per tissue × timepoint × condition
# ---------------------------------------------------------------------------
cat("  Fig 1: Sample overview\n")

sample_counts <- coldata %>%
  dplyr::count(tissue, time_point_wpi, condition) %>%
  dplyr::mutate(
    tissue_label = factor(TISSUE_LABELS[as.character(tissue)],
                          levels = TISSUE_LABELS),
    wpi_label    = factor(time_point_wpi, levels = WPI_LEVELS, labels = WPI_LABELS),
    condition    = factor(condition, levels = c("control", "prion"))
  )

p_overview <- ggplot(sample_counts,
                     aes(x = wpi_label, y = n, fill = condition)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6,
           colour = "white", linewidth = 0.2) +
  geom_text(aes(label = n),
            position = position_dodge(width = 0.7),
            vjust = -0.4, size = 2.4, colour = "grey30") +
  facet_wrap(~ tissue_label, ncol = 1) +
  scale_fill_manual(values = COND_COLOURS, labels = COND_LABELS,
                    name = "Condition") +
  scale_y_continuous(expand = expansion(mult = c(0, 0.18)),
                     breaks = seq(0, 12, 2)) +
  labs(
    x        = "Timepoint (wpi)",
    y        = "Number of samples",
    title    = "GSE210128 — Sample overview",
    subtitle = sprintf(
      "N = %d total  |  3 tissues  |  8 timepoints  |  RML6 prion model (C57BL/6)",
      nrow(coldata)
    ),
    caption  = paste(
      "Term = terminal (~22 wpi). Main + validation cohorts combined.",
      "Bars show N per condition per timepoint."
    )
  ) +
  theme_pub +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())

save_fig(p_overview, "fig1_sample_overview", width = 7, height = 7.5)

# ---------------------------------------------------------------------------
# Fig 2: Library size QC + VST PCA per tissue (2-panel)
# ---------------------------------------------------------------------------
cat("  Fig 2: Library size QC + VST PCA\n")

# Panel A: Library size violin
lib_size_plot_df <- lib_size_df %>%
  dplyr::mutate(
    tissue_label = factor(TISSUE_LABELS[as.character(tissue)],
                          levels = TISSUE_LABELS),
    condition    = factor(condition, levels = c("control", "prion"))
  )

p_libsize <- ggplot(lib_size_plot_df,
                    aes(x = tissue_label, y = library_size / 1e6,
                        fill = condition)) +
  geom_violin(alpha = 0.5, colour = NA, position = position_dodge(0.8)) +
  geom_jitter(aes(colour = condition),
              position = position_jitterdodge(jitter.width = 0.15,
                                              dodge.width  = 0.8),
              size = 1.5, alpha = 0.8) +
  scale_fill_manual(values  = COND_COLOURS, labels = COND_LABELS, name = NULL) +
  scale_colour_manual(values = COND_COLOURS, labels = COND_LABELS, name = NULL) +
  scale_y_continuous(labels = label_number(suffix = " M")) +
  labs(
    x        = NULL,
    y        = "Library size (million reads)",
    title    = "A   Library size distribution",
    subtitle = sprintf("N = %d samples  |  Median = %.1f M reads",
                       nrow(lib_size_df), median(lib_size_df$library_size) / 1e6)
  ) +
  theme_pub +
  theme(legend.position = "bottom",
        panel.grid.major.x = element_blank())

# Panel B: VST PCA (all tissues combined, coloured by tissue + condition)
pca_global <- plotPCA(vsd_all,
                      intgroup = c("condition", "tissue"),
                      returnData = TRUE)
pct_var_global <- round(100 * attr(pca_global, "percentVar"), 1)

p_pca_global <- ggplot(pca_global, aes(x = PC1, y = PC2)) +
  geom_point(aes(colour = tissue, shape = condition),
             size = 2.5, alpha = 0.85, stroke = 0.4) +
  scale_colour_manual(values = TISSUE_COLOURS, labels = TISSUE_LABELS,
                      name = "Tissue") +
  scale_shape_manual(values = c(control = 16, prion = 17),
                     labels = COND_LABELS, name = "Condition") +
  labs(
    x        = sprintf("PC1 (%s%% variance)", pct_var_global[1]),
    y        = sprintf("PC2 (%s%% variance)", pct_var_global[2]),
    title    = "B   Global VST PCA (all tissues)",
    subtitle = "Blind VST normalisation. Shape = condition; colour = tissue."
  ) +
  theme_pub +
  theme(legend.position = "right")

fig2 <- p_libsize + p_pca_global +
  plot_layout(ncol = 2, guides = "keep") &
  theme(legend.position = "bottom")

save_fig(fig2, "fig2_qc_libsize_pca", width = 10, height = 5)

# ---------------------------------------------------------------------------
# Fig 3: Trajectory heatmaps — top LRT genes per tissue (VST z-score)
# ---------------------------------------------------------------------------
cat("  Fig 3: Trajectory heatmaps (top LRT genes)\n")

for (tis in TISSUES) {
  res_lrt <- lrt_results_list[[tis]]
  if (is.null(res_lrt) || sum(res_lrt$significant, na.rm = TRUE) == 0) {
    cat(sprintf("    %s: no LRT significant genes; skipping heatmap.\n",
                TISSUE_LABELS[tis]))
    next
  }

  top_genes <- res_lrt %>%
    dplyr::filter(significant) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = TOP_LRT_GENES) %>%
    dplyr::pull(gene_id)

  idx_tis <- which(colData(vsd_all)$tissue == tis)
  vsd_tis <- vsd_all[, idx_tis]
  mat     <- assay(vsd_tis)[top_genes[top_genes %in% rownames(vsd_tis)], ,
                            drop = FALSE]
  if (nrow(mat) == 0) {
    cat(sprintf("    %s: top genes not in VST matrix; skipping.\n",
                TISSUE_LABELS[tis]))
    next
  }

  mat_z <- t(scale(t(mat)))
  mat_z <- pmax(pmin(mat_z, 3), -3)

  row_label_df <- res_lrt %>%
    dplyr::filter(gene_id %in% rownames(mat_z)) %>%
    dplyr::select(gene_id, gene_symbol) %>%
    dplyr::distinct(gene_id, .keep_all = TRUE) %>%
    dplyr::mutate(label = ifelse(!is.na(gene_symbol) & gene_symbol != "",
                                 gene_symbol, gene_id))
  row_labels <- row_label_df$label[match(rownames(mat_z), row_label_df$gene_id)]
  row_labels[is.na(row_labels)] <- rownames(mat_z)[is.na(row_labels)]

  cd_tis <- as.data.frame(colData(vsd_tis))
  col_annot <- HeatmapAnnotation(
    Condition = cd_tis$condition,
    Cohort    = cd_tis$cohort,
    Timepoint = anno_barplot(
      cd_tis$time_point_wpi,
      gp     = gpar(fill = "grey55", col = NA),
      height = unit(1.0, "cm"),
      axis_param = list(gp = gpar(fontsize = 6))
    ),
    col = list(
      Condition = c(control = COND_COLOURS["control"],
                    prion   = COND_COLOURS["prion"]),
      Cohort    = c(main       = "#999999",
                    validation = "#333333")
    ),
    annotation_name_side = "left",
    annotation_name_gp   = gpar(fontsize = 8),
    show_legend = c(TRUE, TRUE, FALSE)
  )

  col_order <- order(cd_tis$condition, cd_tis$time_point_wpi)

  ht <- Heatmap(
    mat_z[, col_order],
    name                      = "VST\n(z-score)",
    col                       = colorRamp2(c(-3, 0, 3),
                                           c("#2166AC", "white", "#D6604D")),
    top_annotation            = col_annot[col_order],
    show_column_names         = FALSE,
    row_labels                = row_labels,
    row_names_gp              = gpar(fontsize = 6, fontface = "italic"),
    row_names_side            = "right",
    clustering_distance_rows  = "euclidean",
    clustering_method_rows    = "complete",
    cluster_columns           = FALSE,
    column_split              = cd_tis$condition[col_order],
    column_title_gp           = gpar(fontsize = 9, fontface = "bold"),
    heatmap_legend_param      = list(
      title_gp  = gpar(fontsize = 8, fontface = "bold"),
      labels_gp = gpar(fontsize = 7),
      legend_height = unit(3, "cm")
    )
  )

  stem <- sprintf("fig3_%s_trajectory_heatmap", tis)
  save_heatmap(
    ht, stem, width = 8, height = 9,
    column_title = sprintf(
      "%s — Top %d LRT genes (padj <= %.2f, continuous time model)",
      TISSUE_LABELS[tis], min(nrow(mat_z), TOP_LRT_GENES), PADJ_THRESH
    ),
    column_title_gp = gpar(fontsize = 10, fontface = "bold")
  )
}

# ---------------------------------------------------------------------------
# Fig 4: MCPS trajectory line plots (VST mean ± SE per timepoint)
# ---------------------------------------------------------------------------
cat("  Fig 4: MCPS trajectory line plots\n")

mcps_traj_plots <- list()

for (tis in TISSUES) {
  conc <- concordance_list[[tis]]
  if (is.null(conc) || nrow(conc) == 0) next

  mcps_ensembl <- gene_annot %>%
    dplyr::filter(gene_symbol %in% MCPS$gene_symbol) %>%
    dplyr::pull(ensembl_id)

  idx_tis <- which(colData(vsd_all)$tissue == tis)
  vsd_tis <- vsd_all[, idx_tis]

  vst_ensembl <- stringr::str_remove(rownames(vsd_tis), "\\.\\d+$")
  mcps_rows   <- which(vst_ensembl %in% mcps_ensembl)

  if (length(mcps_rows) == 0) {
    cat(sprintf("    %s: no MCPS genes in VST matrix; skipping trajectory plot.\n",
                TISSUE_LABELS[tis]))
    next
  }

  cd_tis <- as.data.frame(colData(vsd_tis))
  vst_mcps <- assay(vsd_tis)[mcps_rows, , drop = FALSE]

  traj_df <- lapply(WPI_LEVELS, function(w) {
    idx_w <- which(cd_tis$time_point_wpi == w)
    if (length(idx_w) == 0) return(NULL)
    lapply(c("control", "prion"), function(cond) {
      idx_c <- idx_w[cd_tis$condition[idx_w] == cond]
      if (length(idx_c) == 0) return(NULL)
      vals <- colMeans(vst_mcps[, idx_c, drop = FALSE])
      data.frame(
        wpi       = w,
        condition = cond,
        mean_vst  = mean(vals),
        se_vst    = sd(vals) / sqrt(length(idx_c)),
        n         = length(idx_c)
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows() %>%
    dplyr::mutate(
      wpi_label = factor(wpi, levels = WPI_LEVELS, labels = WPI_LABELS),
      condition = factor(condition, levels = c("control", "prion"))
    )

  if (nrow(traj_df) == 0) next

  p_traj <- ggplot(traj_df, aes(x = wpi_label, y = mean_vst,
                                colour = condition, group = condition)) +
    geom_ribbon(aes(ymin = mean_vst - se_vst,
                    ymax = mean_vst + se_vst,
                    fill = condition),
                alpha = 0.15, colour = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.5) +
    scale_colour_manual(values = COND_COLOURS, labels = COND_LABELS,
                        name = "Condition") +
    scale_fill_manual(values = COND_COLOURS, labels = COND_LABELS,
                      name = "Condition") +
    labs(
      x        = "Timepoint (wpi)",
      y        = "Mean VST expression (MCPS genes)",
      title    = TISSUE_LABELS[tis],
      subtitle = sprintf(
        "%d MCPS genes  |  Mean ± SE across genes and replicates",
        length(mcps_rows)
      )
    ) +
    theme_pub +
    theme(legend.position = "bottom")

  mcps_traj_plots[[tis]] <- p_traj
}

if (length(mcps_traj_plots) > 0) {
  combined_traj <- wrap_plots(mcps_traj_plots, ncol = 1) +
    plot_annotation(
      title    = "MCPS gene expression trajectories in peripheral tissues",
      subtitle = "Mean VST expression ± SE across MCPS genes. RML6 vs NBH.",
      caption  = sprintf(
        "MCPS = %d CNS-derived genes. Shaded band = ±1 SE.", n_mcps
      ),
      theme = theme_pub
    )
  save_fig(combined_traj, "fig4_mcps_trajectory_lines",
           width = 7, height = 4 * length(mcps_traj_plots))
}

# ---------------------------------------------------------------------------
# Fig 5: Terminal Wald volcano plots (3-panel combined)
# ---------------------------------------------------------------------------
cat("  Fig 5: Terminal volcano plots (3-panel)\n")

volcano_plots <- list()

for (tis in TISSUES) {
  res_term <- terminal_results_list[[tis]]
  if (is.null(res_term)) next

  mcps_term_sig <- if (!is.null(concordance_list[[tis]])) {
    concordance_list[[tis]] %>%
      dplyr::filter(terminal_sig) %>%
      dplyr::pull(gene_symbol)
  } else character(0)

  top15 <- res_term %>%
    dplyr::filter(significant) %>%
    dplyr::arrange(padj) %>%
    dplyr::slice_head(n = 15) %>%
    dplyr::pull(gene_symbol) %>%
    stats::na.omit()

  label_genes <- union(top15, mcps_term_sig)

  vol_df <- res_term %>%
    dplyr::filter(!is.na(padj), !is.na(log2FoldChange)) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(pmax(padj, 1e-300)),
      point_cat = dplyr::case_when(
        !is.na(gene_symbol) & gene_symbol %in% mcps_term_sig & significant
                                                      ~ "MCPS concordant",
        significant                                   ~ "Significant DEG",
        TRUE                                          ~ "Not significant"
      ),
      label = ifelse(!is.na(gene_symbol) & gene_symbol %in% label_genes,
                     gene_symbol, NA_character_)
    )

  n_sig_vol  <- sum(vol_df$point_cat != "Not significant")
  n_mcps_vol <- sum(vol_df$point_cat == "MCPS concordant")
  lfc_max    <- max(abs(vol_df$log2FoldChange), na.rm = TRUE)

  p_vol <- ggplot(vol_df, aes(x = log2FoldChange, y = neg_log10_padj)) +
    geom_point(data = dplyr::filter(vol_df, point_cat == "Not significant"),
               colour = "grey82", alpha = 0.30, size = 0.45) +
    geom_point(data = dplyr::filter(vol_df, point_cat == "Significant DEG"),
               colour = "#4393C3", alpha = 0.80, size = 1.3) +
    geom_point(data = dplyr::filter(vol_df, point_cat == "MCPS concordant"),
               colour = "#D55E00", alpha = 0.95, size = 2.2) +
    geom_hline(yintercept = -log10(PADJ_THRESH),
               linetype = "dashed", colour = "grey40", linewidth = 0.35) +
    geom_vline(xintercept = c(-LFC_THRESH, LFC_THRESH),
               linetype = "dashed", colour = "grey40", linewidth = 0.35) +
    geom_text_repel(
      aes(label = label),
      size          = 2.2,
      max.overlaps  = 40,
      segment.size  = 0.25,
      segment.color = "grey55",
      box.padding   = 0.3,
      force         = 2.5,
      na.rm         = TRUE,
      fontface      = "italic"
    ) +
    scale_x_continuous(
      limits = c(-lfc_max * 1.08, lfc_max * 1.08),
      expand = expansion(mult = 0.02)
    ) +
    labs(
      x        = expression(bold("Shrunken log"[2]*" fold change (RML6 / NBH)")),
      y        = expression(bold(-log[10]*"(padj)")),
      title    = TISSUE_LABELS[tis],
      subtitle = sprintf("DEGs: %d  |  MCPS concordant: %d",
                         n_sig_vol, n_mcps_vol),
      caption  = sprintf(
        "padj <= %.2f, |LFC| >= %.1f (apeglm). Orange = MCPS-concordant.",
        PADJ_THRESH, LFC_THRESH
      )
    ) +
    theme_pub

  volcano_plots[[tis]] <- p_vol
}

for (tis in names(volcano_plots)) {
  save_fig(volcano_plots[[tis]],
           sprintf("fig5_%s_terminal_volcano", tis),
           width = 7, height = 5)
}

if (length(volcano_plots) >= 2) {
  combined_vol <- wrap_plots(volcano_plots, ncol = length(volcano_plots)) +
    plot_annotation(
      title    = "Terminal-timepoint DEGs across peripheral tissues",
      subtitle = "RML6 prion vs NBH control at terminal (~22 wpi). apeglm LFC shrinkage.",
      theme    = theme_pub
    )
  save_fig(combined_vol, "fig5_terminal_volcanos",
           width = 5 * length(volcano_plots), height = 5.5)
}

# ---------------------------------------------------------------------------
# Fig 6: LFC scatter — peripheral terminal LFC vs MCPS meta-LFC
# ---------------------------------------------------------------------------
cat("  Fig 6: LFC scatter (peripheral vs MCPS)\n")

lfc_scatter_plots <- list()

for (tis in TISSUES) {
  conc <- concordance_list[[tis]]
  if (is.null(conc) || nrow(conc) == 0) next

  scatter_df <- conc %>%
    dplyr::filter(!is.na(mcps_lfc), !is.na(terminal_lfc)) %>%
    dplyr::mutate(
      concordance_cat = dplyr::case_when(
        terminal_sig & direction_match  ~ "Concordant & significant",
        terminal_sig & !direction_match ~ "Discordant & significant",
        direction_match                 ~ "Concordant (NS)",
        TRUE                            ~ "Discordant (NS)"
      ),
      concordance_cat = factor(concordance_cat, levels = c(
        "Concordant & significant",
        "Concordant (NS)",
        "Discordant & significant",
        "Discordant (NS)"
      ))
    )

  r_val <- cor(scatter_df$terminal_lfc, scatter_df$mcps_lfc,
               use = "complete.obs")
  binom_p <- concordance_summary$binom_pvalue[
    concordance_summary$tissue == tis]
  binom_p <- if (length(binom_p) == 0) NA_real_ else binom_p[1]

  p_scatter <- ggplot(scatter_df,
                      aes(x = terminal_lfc, y = mcps_lfc)) +
    geom_hline(yintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    geom_vline(xintercept = 0, linetype = "dashed",
               colour = "grey60", linewidth = 0.3) +
    geom_abline(slope = 1, intercept = 0, linetype = "dotted",
                colour = "grey50", linewidth = 0.4) +
    geom_smooth(method = "lm", se = TRUE, colour = "grey40",
                linewidth = 0.6, linetype = "dashed", alpha = 0.15) +
    geom_point(
      data = dplyr::filter(scatter_df, concordance_cat %in%
                      c("Concordant (NS)", "Discordant (NS)")),
      aes(colour = concordance_cat), alpha = 0.45, size = 1.2
    ) +
    geom_point(
      data = dplyr::filter(scatter_df, concordance_cat %in%
                      c("Concordant & significant", "Discordant & significant")),
      aes(colour = concordance_cat), alpha = 0.90, size = 2.0
    ) +
    geom_text_repel(
      data = dplyr::filter(scatter_df, terminal_sig),
      aes(label = gene_symbol),
      size = 2.2, max.overlaps = 25,
      segment.size = 0.25, segment.color = "grey55",
      box.padding = 0.3, force = 2, na.rm = TRUE,
      fontface = "italic"
    ) +
    scale_colour_manual(
      values = c(
        "Concordant & significant" = "#1B7837",
        "Concordant (NS)"          = "#74C476",
        "Discordant & significant" = "#D6604D",
        "Discordant (NS)"          = "#FDAE6B"
      ),
      name = "MCPS gene status"
    ) +
    labs(
      x        = expression(bold("Terminal LFC (RML6 / NBH, apeglm)")),
      y        = expression(bold("MCPS meta-LFC (CNS meta-analysis)")),
      title    = TISSUE_LABELS[tis],
      subtitle = sprintf(
        "Pearson r = %.3f  |  Binomial p (concordance > 50%%) = %.4f  |  N = %d genes",
        r_val, binom_p, nrow(scatter_df)
      ),
      caption  = "Dotted diagonal: identity line (perfect concordance)."
    ) +
    theme_pub +
    theme(legend.position = "right")

  lfc_scatter_plots[[tis]] <- p_scatter
}

for (tis in names(lfc_scatter_plots)) {
  save_fig(lfc_scatter_plots[[tis]],
           sprintf("fig6_%s_lfc_scatter", tis),
           width = 6.5, height = 5.5)
}

if (length(lfc_scatter_plots) >= 2) {
  combined_scatter <- wrap_plots(lfc_scatter_plots, ncol = 1) +
    plot_annotation(
      title    = "Peripheral terminal LFC vs MCPS meta-LFC",
      subtitle = "Each point = one MCPS gene. Concordance = same direction as CNS.",
      theme    = theme_pub
    )
  save_fig(combined_scatter, "fig6_lfc_scatter_peripheral_vs_mcps",
           width = 6.5, height = 5.5 * length(lfc_scatter_plots))
}

# ---------------------------------------------------------------------------
# Fig 7: MCPS concordance bar chart (per tissue)
# ---------------------------------------------------------------------------
cat("  Fig 7: MCPS concordance bar chart\n")

if (nrow(concordance_summary) > 0) {
  conc_long <- concordance_summary %>%
    dplyr::select(tissue_label, n_tested, n_dir_match) %>%
    dplyr::mutate(
      n_discordant = n_tested - n_dir_match,
      n_untested   = n_mcps - n_tested
    ) %>%
    tidyr::pivot_longer(cols = c(n_dir_match, n_discordant, n_untested),
                        names_to = "category", values_to = "count") %>%
    dplyr::mutate(
      category = factor(category,
                        levels = c("n_dir_match", "n_discordant", "n_untested"),
                        labels = c("Concordant", "Discordant", "Not testable")),
      tissue_label = factor(tissue_label, levels = TISSUE_LABELS)
    )

  conc_annot <- concordance_summary %>%
    dplyr::mutate(
      tissue_label = factor(tissue_label, levels = TISSUE_LABELS),
      y_pos        = n_mcps * 1.04,
      label        = sprintf("p=%.3f\nr=%.2f", binom_pvalue, lfc_pearson_r)
    )

  p_conc <- ggplot(conc_long,
                   aes(x = tissue_label, y = count, fill = category)) +
    geom_col(width = 0.6, colour = "white", linewidth = 0.3) +
    geom_text(
      data        = conc_annot,
      aes(x = tissue_label, y = y_pos, label = label),
      inherit.aes = FALSE,
      size = 2.8, hjust = 0.5, colour = "grey30", lineheight = 1.2
    ) +
    scale_fill_manual(
      values = c(
        Concordant      = "#1B7837",
        Discordant      = "#D6604D",
        `Not testable`  = "grey85"
      ),
      name = "MCPS gene status"
    ) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18)),
      breaks = scales::pretty_breaks(n = 6)
    ) +
    labs(
      x        = NULL,
      y        = "Number of MCPS genes",
      title    = "MCPS concordance in peripheral tissues (terminal timepoint)",
      subtitle = sprintf(
        "MCPS = %d CNS-derived genes  |  Concordance = same direction as CNS meta-analysis",
        n_mcps
      ),
      caption  = paste(
        "p: one-sided binomial test (H0: concordance = 50%).",
        "r: Pearson correlation of terminal LFC vs MCPS meta-LFC."
      )
    ) +
    theme_pub +
    theme(legend.position = "bottom",
          panel.grid.major.x = element_blank())

  save_fig(p_conc, "fig7_mcps_concordance_bar", width = 6, height = 5.5)
}

# ---------------------------------------------------------------------------
# Fig 8: MCPS first-detection timeline
# ---------------------------------------------------------------------------
cat("  Fig 8: MCPS first-detection timeline\n")

if (nrow(first_detection_combined) > 0) {
  fd_plot_df <- first_detection_combined %>%
    dplyr::mutate(
      pt_size  = ifelse(first_detection, 4.0, 2.5),
      pt_shape = ifelse(first_detection, 17L, 16L),
      tissue_label = factor(TISSUE_LABELS[as.character(tissue)],
                            levels = TISSUE_LABELS)
    )

  p_detect <- ggplot(fd_plot_df,
                     aes(x = wpi_label, y = n_mcps_sig,
                         colour = tissue, group = tissue)) +
    geom_line(linewidth = 0.85, alpha = 0.9) +
    geom_point(aes(size = pt_size, shape = pt_shape), alpha = 0.9) +
    scale_colour_manual(values = TISSUE_COLOURS, labels = TISSUE_LABELS,
                        name = "Tissue") +
    scale_size_identity() +
    scale_shape_identity() +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18)),
      breaks = scales::pretty_breaks(n = 5),
      limits = c(0, NA)
    ) +
    labs(
      x        = "Timepoint (wpi)",
      y        = sprintf("MCPS genes significant (padj <= %.2f)", PADJ_THRESH),
      title    = "MCPS first detection in peripheral tissues",
      subtitle = paste(
        "Per-timepoint Wald test (RML6 vs NBH).",
        "Triangle = first detection timepoint."
      ),
      caption  = sprintf(
        "MCPS = %d CNS-derived genes. Size factors from full tissue model.",
        n_mcps
      )
    ) +
    theme_pub +
    theme(legend.position = "bottom")

  save_fig(p_detect, "fig8_mcps_first_detection", width = 7, height = 4.5)
}

# ---------------------------------------------------------------------------
# Fig 9: Per-timepoint DEG count (Arm B across wpi)
# ---------------------------------------------------------------------------
cat("  Fig 9: Per-timepoint DEG count\n")

if (nrow(per_timepoint_wald_all) > 0) {
  per_wpi_deg_counts <- per_timepoint_wald_all %>%
    dplyr::filter(!is.na(padj)) %>%
    dplyr::group_by(tissue, wpi) %>%
    dplyr::summarise(
      n_sig = sum(padj <= PADJ_THRESH, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      tissue_label = factor(TISSUE_LABELS[as.character(tissue)],
                            levels = TISSUE_LABELS),
      wpi_label    = factor(wpi, levels = WPI_LEVELS, labels = WPI_LABELS)
    )

  p_deg_count <- ggplot(per_wpi_deg_counts,
                        aes(x = wpi_label, y = n_sig,
                            fill = tissue_label)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.2) +
    geom_text(aes(label = ifelse(n_sig > 0, n_sig, "")),
              vjust = -0.4, size = 2.4, colour = "grey30") +
    facet_wrap(~ tissue_label, ncol = 1, scales = "free_y") +
    scale_fill_manual(values = TISSUE_COLOURS, guide = "none") +
    scale_y_continuous(expand = expansion(mult = c(0, 0.18))) +
    labs(
      x        = "Timepoint (wpi)",
      y        = sprintf("Significant genes (padj <= %.2f)", PADJ_THRESH),
      title    = "Per-timepoint DEG counts in peripheral tissues",
      subtitle = "Per-timepoint Wald test (RML6 vs NBH). Size factors from full tissue model.",
      caption  = "No LFC threshold applied (LRT-equivalent: any significant change)."
    ) +
    theme_pub +
    theme(panel.grid.major.x = element_blank())

  save_fig(p_deg_count, "fig9_per_timepoint_deg_count",
           width = 7, height = 7)
}

# ---------------------------------------------------------------------------
# Fig 10: Longitudinal PCA per tissue (coloured by timepoint)
# ---------------------------------------------------------------------------
cat("  Fig 10: Longitudinal PCA per tissue\n")

pca_plots <- lapply(TISSUES, function(tis) {
  pca_info <- pca_list[[tis]]
  pca_df   <- pca_info$pca_data %>%
    dplyr::mutate(
      wpi_label = factor(time_point_wpi, levels = WPI_LEVELS, labels = WPI_LABELS),
      condition = factor(condition, levels = c("control", "prion"))
    )

  ggplot(pca_df, aes(x = PC1, y = PC2)) +
    geom_point(aes(colour = wpi_label, shape = condition),
               size = 3, alpha = 0.90, stroke = 0.5) +
    scale_colour_viridis_d(option = "plasma", name = "Timepoint (wpi)",
                           direction = -1) +
    scale_shape_manual(values = c(control = 16, prion = 17),
                       labels = COND_LABELS, name = "Condition") +
    labs(
      x        = sprintf("PC1 (%s%% variance)", pca_info$pct_var[1]),
      y        = sprintf("PC2 (%s%% variance)", pca_info$pct_var[2]),
      title    = TISSUE_LABELS[tis],
      subtitle = "Blind VST normalisation. Shape = condition; colour = timepoint."
    ) +
    theme_pub +
    theme(legend.position = "right")
})
names(pca_plots) <- TISSUES

for (tis in TISSUES) {
  save_fig(pca_plots[[tis]],
           sprintf("fig10_%s_pca_longitudinal", tis),
           width = 6.5, height = 5)
}

combined_pca <- wrap_plots(pca_plots, ncol = 1) +
  plot_annotation(
    title    = "Longitudinal PCA — peripheral tissues",
    subtitle = "GSE210128: RML6 prion vs NBH control. Blind VST normalisation.",
    theme    = theme_pub
  )
save_fig(combined_pca, "fig10_all_pca_longitudinal",
         width = 6.5, height = 5 * length(TISSUES))

# =============================================================================
# FIGURE 11 — Cross-tissue MCPS LFC heatmap
# =============================================================================
cat("  Generating Fig 11: cross-tissue MCPS LFC heatmap...\n")

tryCatch({

  lfc_mat_list <- lapply(TISSUES, function(tis) {
    df <- concordance_list[[tis]]
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df %>%
      dplyr::select(gene_symbol, log2FoldChange) %>%
      dplyr::rename(!!tis := log2FoldChange)
  })
  lfc_mat_list <- Filter(Negate(is.null), lfc_mat_list)

  if (length(lfc_mat_list) == 0) {
    cat("  Skipping Fig 11: no concordance data available.\n")
  } else {

    lfc_wide <- Reduce(function(a, b) dplyr::full_join(a, b, by = "gene_symbol"),
                       lfc_mat_list)

    lfc_wide <- lfc_wide %>% dplyr::filter(!is.na(gene_symbol))

    mcps_dir <- MCPS %>%
      dplyr::select(gene_symbol, meta_log2FC) %>%
      dplyr::distinct(gene_symbol, .keep_all = TRUE) %>%
      dplyr::mutate(mcps_direction = ifelse(meta_log2FC > 0, "Up in CNS", "Down in CNS"))

    lfc_wide <- lfc_wide %>%
      dplyr::left_join(mcps_dir, by = "gene_symbol")

    mat_genes <- lfc_wide$gene_symbol
    mat_dir   <- lfc_wide$mcps_direction
    mat       <- as.matrix(lfc_wide[, TISSUES, drop = FALSE])
    rownames(mat) <- mat_genes
    colnames(mat) <- TISSUE_LABELS[TISSUES]

    LFC_CAP <- 3
    mat_capped <- pmin(pmax(mat, -LFC_CAP), LFC_CAP)

    mean_abs_lfc <- rowMeans(abs(mat_capped), na.rm = TRUE)
    row_order    <- order(mat_dir, -mean_abs_lfc)
    mat_capped   <- mat_capped[row_order, , drop = FALSE]
    mat_dir_ord  <- mat_dir[row_order]

    col_fun <- circlize::colorRamp2(
      c(-LFC_CAP, 0, LFC_CAP),
      c("#2166AC", "white", "#D6604D")
    )

    dir_colours <- c("Up in CNS" = "#D6604D", "Down in CNS" = "#2166AC")
    row_ann <- rowAnnotation(
      `CNS direction` = mat_dir_ord,
      col = list(`CNS direction` = dir_colours),
      annotation_legend_param = list(
        `CNS direction` = list(title = "MCPS direction")
      ),
      show_annotation_name = TRUE
    )

    tis_colours <- setNames(
      RColorBrewer::brewer.pal(max(3, length(TISSUES)), "Set2")[seq_along(TISSUES)],
      TISSUE_LABELS[TISSUES]
    )
    col_ann <- HeatmapAnnotation(
      Tissue = colnames(mat_capped),
      col    = list(Tissue = tis_colours),
      show_annotation_name = FALSE
    )

    n_genes_hm <- nrow(mat_capped)
    show_row_names_hm <- n_genes_hm <= 80
    row_names_gp_hm   <- if (show_row_names_hm) {
      gpar(fontsize = max(4, min(8, 200 / n_genes_hm)))
    } else {
      gpar(fontsize = 0)
    }

    hm_fig11 <- Heatmap(
      mat_capped,
      name            = "LFC (capped ±3)",
      col             = col_fun,
      na_col          = "grey90",
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      top_annotation  = col_ann,
      left_annotation = row_ann,
      show_row_names  = show_row_names_hm,
      row_names_gp    = row_names_gp_hm,
      show_column_names = TRUE,
      column_names_gp   = gpar(fontsize = 10),
      column_title    = "Peripheral tissue",
      row_title       = sprintf("MCPS genes (n=%d)", n_genes_hm),
      row_title_gp    = gpar(fontsize = 10),
      heatmap_legend_param = list(
        title          = "apeglm LFC\n(capped ±3)",
        legend_height  = unit(4, "cm"),
        at             = c(-3, -2, -1, 0, 1, 2, 3)
      ),
      border = TRUE
    )

    hm_height <- min(14, max(5, 0.12 * n_genes_hm + 2))
    save_heatmap(hm_fig11, "fig11_mcps_cross_tissue_heatmap",
                 width = 6.5, height = hm_height)
    cat(sprintf("  Fig 11 saved (%d genes × %d tissues).\n",
                n_genes_hm, length(TISSUES)))
  }

}, error = function(e) {
  warning(sprintf("Fig 11 (cross-tissue heatmap) failed: %s", conditionMessage(e)))
  cat("  Fig 11 skipped due to error (see warning above).\n")
})

cat("  All figures generated.\n")


# ----------------------------------------------------------------------
# ##Section 10 — Figure 11: Cross-Tissue Heatmap
#
# **Purpose:** Generate the cross-tissue heatmap showing MCPS gene LFCs across all three streams (CNS meta-analysis, human brain, peripheral blood) side by side. This is the key integration figure for the manuscript.
#
# **Inputs:** `../results/stream1/meta_results.tsv`; `../results/stream2/results_GSE124571.tsv`; `../results/stream3/wald_terminal_results.tsv`.
#
# **Outputs:** `../results/stream3/figures/fig_cross_tissue_heatmap.png` (PNG + SVG).
#
# **Key decisions:** LFCs are z-scored within each stream for visual comparability. Asterisks indicate padj<0.05 within each stream. Slc15a4 direction discordance in human is shown explicitly.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# ##Section 11 — Outputs & Summary
#
# **Purpose:** Write all final output tables for Stream 3 and generate a plain-text summary report covering sample counts, LRT/Wald statistics, MCPS gene trajectories, concordance results, and first-detection findings.
#
# **Inputs:** All Stream 3 result objects in R.
#
# **Outputs:** `../results/stream3/stream3_summary.txt`; all TSV tables finalised and written.
#
# **Key decisions:** The first-detection results for MCPS genes are highlighted in the summary as a key translational finding (potential early biomarker signal).
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 11 — Write output tables and summary report
# =============================================================================
# All CSV outputs are written with tryCatch to prevent a single write failure
# from aborting the entire script. The summary report is written last so it
# can include all computed statistics.
# =============================================================================
cat("\n--- Section 11: Writing output tables ---\n")

# ---------------------------------------------------------------------------
# 11.1  LRT results (per tissue)
# ---------------------------------------------------------------------------
for (tis in TISSUES) {
  res <- lrt_results_list[[tis]]
  if (is.null(res)) next
  out_path <- file.path(OUT_DIR, sprintf("GSE210128_%s_LRT_results.csv", tis))
  tryCatch({
    write_csv(
      res %>%
        dplyr::select(
          gene_id, ensembl_id, gene_symbol, gene_biotype,
          baseMean, stat, pvalue, padj,
          any_of(c("lrt_fac_pvalue", "lrt_fac_padj")),
          significant, tissue
        ),
      out_path
    )
    cat(sprintf("  Written: GSE210128_%s_LRT_results.csv\n", tis))
  }, error = function(e) {
    warning(sprintf("Failed to write LRT results for %s: %s", tis, conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.2  Terminal Wald results (per tissue)
# ---------------------------------------------------------------------------
for (tis in TISSUES) {
  res <- terminal_results_list[[tis]]
  if (is.null(res)) next
  out_path <- file.path(OUT_DIR,
                        sprintf("GSE210128_%s_terminal_results.csv", tis))
  tryCatch({
    write_csv(
      res %>%
        dplyr::select(
          gene_id, ensembl_id, gene_symbol, gene_biotype,
          baseMean, log2FoldChange, lfcSE,
          pval_wald, lfc_wald, padj,
          significant, direction, tissue
        ),
      out_path
    )
    cat(sprintf("  Written: GSE210128_%s_terminal_results.csv\n", tis))
  }, error = function(e) {
    warning(sprintf("Failed to write terminal results for %s: %s", tis, conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.3  MCPS concordance (per tissue)
# ---------------------------------------------------------------------------
for (tis in TISSUES) {
  conc <- concordance_list[[tis]]
  if (is.null(conc)) next
  out_path <- file.path(OUT_DIR,
                        sprintf("GSE210128_%s_MCPS_concordance.csv", tis))
  tryCatch({
    write_csv(
      conc %>%
        dplyr::select(
          gene_symbol, gene_biotype,
          any_of(c("meta_padj", "I2", "n_studies")),
          mcps_lfc,
          terminal_lfc, terminal_padj, terminal_sig,
          direction_match, tissue
        ) %>%
        arrange(desc(direction_match), terminal_padj),
      out_path
    )
    cat(sprintf("  Written: GSE210128_%s_MCPS_concordance.csv\n", tis))
  }, error = function(e) {
    warning(sprintf("Failed to write concordance for %s: %s", tis, conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.4  MCPS concordance summary (cross-tissue)
# ---------------------------------------------------------------------------
if (nrow(concordance_summary) > 0) {
  tryCatch({
    write_csv(concordance_summary,
              file.path(OUT_DIR, "GSE210128_MCPS_concordance_summary.csv"))
    cat("  Written: GSE210128_MCPS_concordance_summary.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write concordance summary: %s", conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.5  MCPS first-detection timeline
# ---------------------------------------------------------------------------
if (nrow(first_detection_combined) > 0) {
  tryCatch({
    write_csv(
      first_detection_combined %>%
        dplyr::select(tissue, tissue_label, wpi, wpi_label,
                      n_mcps_tested, n_mcps_sig, first_detection),
      file.path(OUT_DIR, "GSE210128_MCPS_first_detection.csv")
    )
    cat("  Written: GSE210128_MCPS_first_detection.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write first-detection table: %s", conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.6  MCPS gene-level first-detection table
# ---------------------------------------------------------------------------
# mcps_first_genes_all was already constructed in Section 9.7.0 (before the
# Stream 2 cross-reference). It is reused here for writing the output CSV.
# No re-construction needed.

if (nrow(mcps_first_genes_all) > 0) {
  tryCatch({
    write_csv(
      mcps_first_genes_all %>%
        dplyr::select(tissue, gene_symbol, first_detected_wpi, wpi_label,
                      log2FoldChange, padj, mcps_lfc, direction_match) %>%
        arrange(tissue, first_detected_wpi, padj),
      file.path(OUT_DIR, "GSE210128_MCPS_first_detection_genes.csv")
    )
    cat("  Written: GSE210128_MCPS_first_detection_genes.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write gene-level first-detection: %s", conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.6b  Triple-validated candidates (Section 9.7 cross-reference)
# ---------------------------------------------------------------------------
# Written only when Stream 2 output was available and cross-reference succeeded.
if (stream2_available && nrow(triple_validated) > 0) {
  tryCatch({
    write_csv(
      triple_validated %>%
        dplyr::select(
          gene_symbol, tissue, first_detected_wpi, wpi_label,
          log2FoldChange, padj, mcps_lfc, direction_match,
          any_of(c("human_lfc", "human_padj")),
          evidence_tier
        ) %>%
        arrange(evidence_tier, first_detected_wpi, padj),
      file.path(OUT_DIR, "GSE210128_triple_validated_candidates.csv")
    )
    cat("  Written: GSE210128_triple_validated_candidates.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write triple-validated candidates: %s", conditionMessage(e)))
  })
} else if (stream2_available && nrow(triple_validated) == 0) {
  # Stream 2 was available but no genes passed all three criteria.
  # Write an empty-but-valid CSV so Script 08 doesn't report this file as MISSING.
  empty_triple <- data.frame(
    gene_symbol = character(0), tissue = character(0),
    first_detected_wpi = integer(0), wpi_label = character(0),
    log2FoldChange = numeric(0), padj = numeric(0),
    mcps_lfc = numeric(0), direction_match = logical(0),
    evidence_tier = character(0),
    stringsAsFactors = FALSE
  )
  tryCatch({
    write_csv(empty_triple,
              file.path(OUT_DIR, "GSE210128_triple_validated_candidates.csv"))
    cat("  Written: GSE210128_triple_validated_candidates.csv (0 rows — no triple-validated genes found)\n")
  }, error = function(e) {
    warning(sprintf("Failed to write empty triple-validated CSV: %s", conditionMessage(e)))
  })
} else if (!stream2_available) {
  cat("  Skipped: GSE210128_triple_validated_candidates.csv (Stream 2 not available)\n")
}

# ---------------------------------------------------------------------------
# 11.7  Per-timepoint Wald results (supplementary table)
# ---------------------------------------------------------------------------
if (nrow(per_timepoint_wald_all) > 0) {
  tryCatch({
    write_csv(
      per_timepoint_wald_all %>%
        dplyr::select(tissue, wpi, gene_id,
                      any_of(c("ensembl_id", "gene_symbol", "gene_biotype")),
                      log2FoldChange, lfcSE, pvalue, padj) %>%
        arrange(tissue, wpi, padj),
      file.path(OUT_DIR, "GSE210128_per_timepoint_wald_all.csv")
    )
    cat("  Written: GSE210128_per_timepoint_wald_all.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write per-timepoint Wald table: %s", conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.8  Cross-tissue terminal DEG overlap
# ---------------------------------------------------------------------------
if (nrow(cross_tissue_overlap) > 0) {
  tryCatch({
    write_csv(cross_tissue_overlap,
              file.path(OUT_DIR, "GSE210128_cross_tissue_deg_overlap.csv"))
    cat("  Written: GSE210128_cross_tissue_deg_overlap.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write cross-tissue overlap: %s", conditionMessage(e)))
  })
}

# ---------------------------------------------------------------------------
# 11.9  Plain-text summary report
# ---------------------------------------------------------------------------
summary_lines <- c(
  "=============================================================",
  "  Stream 3: Peripheral Tissue Longitudinal Analysis",
  "  07_stream3_peripheral_analysis.R",
  "=============================================================",
  "",
  "DATASET",
  "-------",
  "  Accession:  GSE210128",
  "  Reference:  Aguzzi A, Frontzek K, Koderman M, Caredio D (2022)",
  "              University Hospital of Zurich",
  "              https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE210128",
  "  Organism:   Mus musculus (C57BL/6)",
  "  Prion:      RML6 (Rocky Mountain Laboratory strain, passage 6)",
  "  Control:    NBH (non-infectious brain homogenate)",
  sprintf("  N total:    %d samples (after QC filters)", nrow(coldata)),
  "  Tissues:    whole_blood, spleen, skeletal_muscle",
  "  Timepoints: 4, 8, 12, 14, 16, 18, 20 wpi + terminal (encoded as 22 wpi)",
  "",
  "CONSORT SAMPLE FLOW",
  "-------------------",
  sprintf("  Raw coldata:            %d", n_raw),
  sprintf("  After NA filter:        %d  (excluded: %d)", n_after_na_filter,
          n_raw - n_after_na_filter),
  sprintf("  After tissue filter:    %d  (excluded: %d)", n_after_tissue_filter,
          n_after_na_filter - n_after_tissue_filter),
  sprintf("  After timepoint filter: %d  (excluded: %d)", n_after_wpi_filter,
          n_after_tissue_filter - n_after_wpi_filter),
  "",
  "SAMPLE BREAKDOWN",
  "----------------"
)

for (tis in TISSUES) {
  sub <- coldata[coldata$tissue == tis, ]
  summary_lines <- c(summary_lines,
    sprintf("  %-20s N=%3d  (%d prion, %d control)  cohorts: %s",
            TISSUE_LABELS[tis], nrow(sub),
            sum(sub$condition == "prion"),
            sum(sub$condition == "control"),
            paste(sort(unique(as.character(sub$cohort))), collapse = "+")))
}

summary_lines <- c(summary_lines,
  "",
  "ARM A — LONGITUDINAL LRT (continuous time)",
  "------------------------------------------",
  "  Test: condition x time interaction only (diverging trajectories)",
  "  NOTE: Terminal timepoint excluded (exact wpi unknown; used only in Arm B)",
  sprintf("  Threshold: padj <= %.2f", PADJ_THRESH)
)

for (tis in TISSUES) {
  summary_lines <- c(summary_lines,
    sprintf("  %-20s %d significant genes  |  model: %s",
            TISSUE_LABELS[tis], lrt_sig_counts[tis],
            lrt_model_used[tis]))
}

summary_lines <- c(summary_lines,
  "",
  "ARM B — TERMINAL WALD TEST",
  "--------------------------",
  sprintf("  Terminal timepoint: %d wpi (encoded; exact wpi not in GEO metadata)",
          TERMINAL_WPI),
  "  LFC shrinkage: apeglm",
  "  Size factors: estimated from full tissue dataset (all timepoints)",
  sprintf("  Thresholds: padj <= %.2f, |LFC| >= %.1f", PADJ_THRESH, LFC_THRESH)
)

for (tis in TISSUES) {
  res <- terminal_results_list[[tis]]
  if (is.null(res)) {
    summary_lines <- c(summary_lines,
      sprintf("  %-20s SKIPPED (insufficient samples)", TISSUE_LABELS[tis]))
  } else {
    n_up <- sum(res$direction == "up",   na.rm = TRUE)
    n_dn <- sum(res$direction == "down", na.rm = TRUE)
    summary_lines <- c(summary_lines,
      sprintf("  %-20s %d DEGs (%d up, %d down)",
              TISSUE_LABELS[tis], terminal_sig_counts[tis], n_up, n_dn))
  }
}

summary_lines <- c(summary_lines,
  "",
  "CROSS-TISSUE TERMINAL DEG OVERLAP",
  "----------------------------------"
)
if (nrow(cross_tissue_overlap) > 0) {
  for (i in seq_len(nrow(cross_tissue_overlap))) {
    r <- cross_tissue_overlap[i, ]
    summary_lines <- c(summary_lines,
      sprintf("  %s vs %s: shared=%d  Jaccard=%.3f",
              r$tissue_1, r$tissue_2, r$n_shared, r$jaccard))
  }
}

summary_lines <- c(summary_lines,
  "",
  "ARM C — MCPS CONCORDANCE (TERMINAL)",
  "------------------------------------",
  sprintf("  MCPS size (full): %d genes (%d up, %d down in CNS)",
          n_mcps, n_mcps_up, n_mcps_dn),
  "  Concordance = same direction as CNS meta-analysis LFC",
  "  Denominator = nrow(MCPS) (all MCPS genes)"
)

if (nrow(concordance_summary) > 0) {
  for (i in seq_len(nrow(concordance_summary))) {
    r <- concordance_summary[i, ]
    summary_lines <- c(summary_lines,
      sprintf("  %-20s tested=%d/%d  concordant=%d (%.1f%%)  binom_p=%.4f  r=%.3f",
              r$tissue_label, r$n_tested, r$n_mcps_total,
              r$n_dir_match, r$pct_of_mcps,
              r$binom_pvalue, r$lfc_pearson_r))
  }
}

summary_lines <- c(summary_lines,
  "",
  "MCPS FIRST DETECTION",
  "--------------------"
)

for (tis in TISSUES) {
  fd  <- first_detection_list[[tis]]$first_det_wpi
  lbl <- if (!is.na(fd)) {
    wl <- WPI_LABELS[match(fd, WPI_LEVELS)]
    sprintf("%s wpi", ifelse(is.na(wl), as.character(fd), wl))
  } else "not detected"
  summary_lines <- c(summary_lines,
    sprintf("  %-20s %s", TISSUE_LABELS[tis], lbl))
}

summary_lines <- c(summary_lines,
  "",
  "CROSS-STREAM INTEGRATION (SECTION 9.7)",
  "---------------------------------------"
)

if (stream2_available && nrow(triple_validated) > 0) {
  n_triple_sum <- length(unique(triple_validated$gene_symbol))
  n_tier1_sum  <- triple_validated %>%
    filter(evidence_tier == "Tier 1: concordant peripheral + human leading edge") %>%
    pull(gene_symbol) %>% unique() %>% length()
  summary_lines <- c(summary_lines,
    sprintf("  Stream 2 GSEA leading-edge genes: %d",
            length(stream2_leading_genes)),
    sprintf("  Triple-validated candidates:      %d unique genes",
            n_triple_sum),
    sprintf("    Tier 1 (concordant peripheral + human leading edge): %d",
            n_tier1_sum),
    sprintf("    Tier 2 (detected peripheral + human leading edge):   %d",
            n_triple_sum - n_tier1_sum),
    "  Output: GSE210128_triple_validated_candidates.csv"
  )
} else if (!stream2_available) {
  summary_lines <- c(summary_lines,
    "  Stream 2 output not available; cross-reference skipped.",
    "  Re-run after 06_stream2_human_validation.R completes."
  )
} else {
  summary_lines <- c(summary_lines,
    "  No triple-validated candidates found."
  )
}

summary_lines <- c(summary_lines,
  "",
  "LIMITATIONS",
  "-----------",
  "  1. Pseudo-replication: main + validation cohorts pooled; cohort blocked",
  "     but animals not individually tracked across timepoints.",
  "  2. Blood has only 3 replicates/group/timepoint (no validation cohort);",
  "     LRT power is low for blood.",
  "  3. Terminal timepoint encoded as 22 wpi (approximate); exact days vary",
  "     between animals and are not available in GEO metadata.",
  "  4. LRT p-values are not directional; use trajectory heatmaps and",
  "     per-timepoint Wald tests for directionality.",
  "  5. MCPS was derived from CNS (brain) studies; peripheral concordance",
  "     is expected to be lower than CNS replication (Stream 1 Stratum D).",
  "  6. Per-timepoint Wald tests have low power (small N per timepoint);",
  "     first-detection timepoints should be interpreted as lower bounds.",
  "  7. No cross-tissue multiple testing correction applied.",
  "",
  "OUTPUT FILES",
  "------------",
  "  GSE210128_<tissue>_LRT_results.csv          (3 files)",
  "  GSE210128_<tissue>_terminal_results.csv      (3 files)",
  "  GSE210128_<tissue>_MCPS_concordance.csv      (3 files)",
  "  GSE210128_MCPS_concordance_summary.csv",
  "  GSE210128_MCPS_first_detection.csv",
  "  GSE210128_MCPS_first_detection_genes.csv",
  "  GSE210128_per_timepoint_wald_all.csv",
  "  GSE210128_cross_tissue_deg_overlap.csv",
  "  fig1_sample_overview.pdf/png",
  "  fig2_qc_libsize_pca.pdf/png",
  "  fig3_<tissue>_trajectory_heatmap.pdf/png     (3 files)",
  "  fig4_mcps_trajectory_lines.pdf/png",
  "  fig5_<tissue>_terminal_volcano.pdf/png       (3 files)",
  "  fig5_terminal_volcanos.pdf/png               (combined 3-panel)",
  "  fig6_<tissue>_lfc_scatter.pdf/png            (3 files)",
  "  fig6_lfc_scatter_peripheral_vs_mcps.pdf/png  (combined)",
  "  fig7_mcps_concordance_bar.pdf/png",
  "  fig8_mcps_first_detection.pdf/png",
  "  fig9_per_timepoint_deg_count.pdf/png",
  "  fig10_<tissue>_pca_longitudinal.pdf/png      (3 files)",
  "  fig10_all_pca_longitudinal.pdf/png           (combined)",
  "  fig11_mcps_cross_tissue_heatmap.pdf/png",
  "  GSE210128_triple_validated_candidates.csv    (if Stream 2 available)",
  "",
  "PACKAGE VERSIONS",
  "----------------",
  sprintf("  R version:       %s", R.version.string),
  sprintf("  DESeq2:          %s", as.character(packageVersion("DESeq2"))),
  sprintf("  tximeta:         %s", as.character(packageVersion("tximeta"))),
  sprintf("  apeglm:          %s", as.character(packageVersion("apeglm"))),
  sprintf("  ComplexHeatmap:  %s", as.character(packageVersion("ComplexHeatmap"))),
  sprintf("  ggplot2:         %s", as.character(packageVersion("ggplot2"))),
  sprintf("  patchwork:       %s", as.character(packageVersion("patchwork"))),
  sprintf("  dplyr:           %s", as.character(packageVersion("dplyr"))),
  "",
  "=============================================================",
  sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "============================================================="
)

tryCatch({
  writeLines(summary_lines,
             file.path(OUT_DIR, "GSE210128_stream3_summary.txt"))
  cat("  Written: GSE210128_stream3_summary.txt\n")
}, error = function(e) {
  warning(sprintf("Failed to write summary report: %s", conditionMessage(e)))
})

cat(paste(summary_lines, collapse = "\n"), "\n")

# ---------------------------------------------------------------------------
# 11.10  Session info (for reproducibility)
# ---------------------------------------------------------------------------
tryCatch({
  session_info_path <- file.path(OUT_DIR, "GSE210128_session_info.txt")
  writeLines(capture.output(sessionInfo()), session_info_path)
  cat("  Written: GSE210128_session_info.txt\n")
}, error = function(e) {
  warning(sprintf("Failed to write session info: %s", conditionMessage(e)))
})

cat("\n=============================================================\n")
cat("  07_stream3_peripheral_analysis.R complete.\n")
cat(sprintf("  All outputs written to: %s\n", OUT_DIR))
cat("=============================================================\n")
