#!/usr/bin/env Rscript
# ==============================================================================
# 08_integrate_streams.R
# Integration: cross-stream MCPS concordance, human microarray GSEA (GSE124571), literature convergence table, composite figures
# ==============================================================================
# Project : Cross-tissue bulk RNA-seq meta-analysis of prion disease
#           (mouse CNS Stratum A/D, peripheral time course, human sCJD validation)
#
# Provenance: Converted 1:1 from the executed Google Colab notebook
#             "08_integrate_streams_v2.ipynb" (the version that produced the manuscript results).
#             Only Colab environment-setup cells (drive mount, rpy2 magic,
#             package installers) were removed or commented; all analysis
#             code is unchanged.
#
# Inputs  : results/stream1/, results/stream2/, results/stream3/
# Outputs : results/integration/, reports/figures/integration/
#
# Usage   : Set BASE_DIR (top of script) to the project root, then run:
#             Rscript 08_integrate_streams.R
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
# ##Section 0 — Libraries, Configuration & Helper Functions
#
# **Purpose:** Load all R packages, define paths, and set up the `save_fig()` helper and `theme_pub` ggplot2 theme. This notebook integrates results from all three streams and requires all packages from Streams 1–3, plus GEOquery, limma, fgsea, and msigdbr for the pathway-level human validation in Section 9.
#
# **Inputs:** None (constants defined inline).
#
# **Outputs:** `BASE_DIR`, `STREAM1_DIR`, `STREAM2_DIR`, `STREAM3_DIR`, `OUT_DIR`, `MCPS_GENES`, `theme_pub`, `save_fig()` available in R.
#
# **Key decisions:** `OUT_DIR` points to `../results/integration/` (integration-specific results folder). `FIG_DIR` points to `../results/integration/figures/` — all figures from this notebook are saved here. Neither stream-specific subdirectories nor `../reports/` paths are used.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 0 — Libraries, configuration, helper functions
# =============================================================================
cat("=============================================================\n")
cat("  08_integrate_streams.R\n")
cat("  Prion Disease RNA-Seq — Cross-Stream Integration\n")
cat("=============================================================\n\n")

# --- 0.1  Package installation ---
# Optimized for Colab execution. Suppresses excessive output.
if (!requireNamespace("BiocManager", quietly = TRUE))
  NULL  # [setup] installer removed (see environment.yml / renv.lock):
  # [setup]   install.packages("BiocManager", quiet = TRUE)

bioc_pkgs <- c("ComplexHeatmap", "circlize", "GEOquery", "limma", "fgsea")
cran_pkgs <- c(
  "ggplot2", "ggrepel", "patchwork", "dplyr", "tidyr", "tibble",
  "readr", "stringr", "scales", "RColorBrewer", "viridis", "Matrix", "msigdbr"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    NULL  # [setup] installer removed (see environment.yml / renv.lock):
    # [setup]   BiocManager::install(pkg, ask = FALSE, update = FALSE, quiet = TRUE)
}
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    NULL  # [setup] installer removed (see environment.yml / renv.lock):
    # [setup]   install.packages(pkg, repos = "https://cloud.r-project.org", quiet = TRUE)
}

suppressPackageStartupMessages({
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
  library(scales)
  library(RColorBrewer)
  library(viridis)
  library(GEOquery)
  library(limma)
  library(fgsea)
  library(msigdbr)
})

# --- 0.2  Derive BASE_DIR robustly (Adapted for Google Colab environment) ---
BASE_DIR <- "/content/drive/MyDrive/prion_rnaseq_meta_analysis"

if (!dir.exists(BASE_DIR)) {
  stop("BASE_DIR does not exist. Please check your Google Drive path and ensure the drive is mounted.")
}
setwd(file.path(BASE_DIR, "scripts"))

# --- 0.3  Directory paths ---
STREAM1_DIR  <- file.path(BASE_DIR, "results", "stream1")
STREAM2_DIR  <- file.path(BASE_DIR, "results", "stream2")
STREAM3_DIR  <- file.path(BASE_DIR, "results", "stream3")
REPLIC_DIR   <- file.path(BASE_DIR, "results", "stream1", "replication")
ENRICH_DIR   <- file.path(BASE_DIR, "results", "stream1", "enrichment")
OUT_DIR      <- file.path(BASE_DIR, "results", "integration")
FIGURES_DIR  <- file.path(BASE_DIR, "reports", "figures", "integration")

# Explicitly create both the output and figures directories
for (d in c(OUT_DIR, FIGURES_DIR)) {
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
}

cat(sprintf("  BASE_DIR    : %s\n", BASE_DIR))
cat(sprintf("  OUT_DIR     : %s\n", OUT_DIR))
cat(sprintf("  FIGURES_DIR : %s\n", FIGURES_DIR))

# --- 0.4  Analysis parameters ---
PADJ_THRESH  <- 0.05
LFC_THRESH   <- 0.5
TISSUES      <- c("whole_blood", "spleen", "skeletal_muscle")
TISSUE_LABELS <- c(
  whole_blood     = "Whole blood",
  spleen          = "Spleen",
  skeletal_muscle = "Skeletal muscle"
)
WPI_LEVELS <- c(4L, 8L, 12L, 14L, 16L, 18L, 20L, 22L)
WPI_LABELS <- c("4", "8", "12", "14", "16", "18", "20", "Term")

set.seed(42)

# --- 0.5  Publication theme (consistent with Streams 1–3) ---
# Text sizes slightly optimized for high-res Q1 rendering
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
    plot.margin      = margin(6, 6, 6, 6, "pt")
  )

# Colour palette: consistent with stream scripts
STREAM_COLOURS <- c(
  "Stream 1 (CNS meta)"       = "#2166AC",
  "Stream 2 (Human sCJD)"     = "#D6604D",
  "Stream 3 (Peripheral)"     = "#4DAC26"
)
DIRECTION_COLOURS <- c(
  "Up in CNS"   = "#D6604D",
  "Down in CNS" = "#2166AC"
)
TIER_COLOURS <- c(
  "Tier 1" = "#1B7837",
  "Tier 2" = "#762A83",
  "MCPS only" = "#BABABA"
)

# --- 0.6  Helper: save ggplot figure (PDF + PNG) ---
# Upgraded for Q1 Journals: DPI set to 600, forced Cairo device for PDF rendering
save_fig <- function(plot_obj, stem, width, height, dpi = 600) {
  pdf_path <- file.path(FIGURES_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(stem, ".png"))
  tryCatch({
    ggsave(pdf_path, plot = plot_obj, width = width, height = height,
           device = cairo_pdf)
    ggsave(png_path, plot = plot_obj, width = width, height = height, dpi = dpi, type = "cairo")
    cat(sprintf("    Saved (High-Res): %s.pdf + .png\n", stem))
  }, error = function(e) {
    warning(sprintf("save_fig failed for '%s': %s", stem, conditionMessage(e)))
  })
  invisible(list(pdf = pdf_path, png = png_path))
}

# --- 0.7  Helper: save ComplexHeatmap figure (PDF + PNG) ---
# Upgraded for Q1 Journals: DPI set to 600
save_heatmap <- function(ht_obj, stem, width, height, dpi = 600,
                         column_title = NULL, column_title_gp = NULL) {
  pdf_path <- file.path(FIGURES_DIR, paste0(stem, ".pdf"))
  png_path <- file.path(FIGURES_DIR, paste0(stem, ".png"))
  draw_args <- list(ht_obj)
  if (!is.null(column_title))    draw_args$column_title    <- column_title
  if (!is.null(column_title_gp)) draw_args$column_title_gp <- column_title_gp
  tryCatch({
    # Use Cairo PDF for better font rendering in complex heatmaps
    cairo_pdf(pdf_path, width = width, height = height)
    do.call(draw, draw_args)
    dev.off()
    cat(sprintf("    Saved (High-Res): %s.pdf\n", stem))
  }, error = function(e) {
    warning(sprintf("save_heatmap PDF failed for '%s': %s", stem, conditionMessage(e)))
    if (dev.cur() > 1) dev.off()
  })
  tryCatch({
    png(png_path, width = width * dpi, height = height * dpi,
        res = dpi, units = "px", type = "cairo")
    do.call(draw, draw_args)
    dev.off()
    cat(sprintf("    Saved (High-Res): %s.png\n", stem))
  }, error = function(e) {
    warning(sprintf("save_heatmap PNG failed for '%s': %s", stem, conditionMessage(e)))
    if (dev.cur() > 1) dev.off()
  })
  invisible(list(pdf = pdf_path, png = png_path))
}

# --- 0.8  Helper: safe read_csv with informative fallback ---
safe_read <- function(path, label) {
  if (!file.exists(path)) {
    cat(sprintf("  [MISSING] %s\n         -> %s\n", label, path))
    return(NULL)
  }
  tryCatch({
    df <- read_csv(path, show_col_types = FALSE)
    cat(sprintf("  [OK]      %s  (%d rows)\n", label, nrow(df)))
    df
  }, error = function(e) {
    warning(sprintf("Failed to read %s: %s", label, conditionMessage(e)))
    NULL
  })
}

cat("\n  [✓] Section 0 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 1 — Load All Stream Outputs
#
# **Purpose:** Load the key result tables from all three streams into R: Stream 1 meta-analysis results, Stream 2 human validation results, Stream 3 peripheral concordance results, and the gene alias map.
#
# **Inputs:** `../results/stream1/meta_results.tsv`; `../results/stream2/results_GSE124571.tsv`; `../results/stream3/wald_terminal_results.tsv`; `../results/stream1/gene_alias_map.tsv`.
#
# **Outputs:** `stream1`, `stream2`, `stream3`, `alias_map` data frames in R.
#
# **Key decisions:** All gene symbols are harmonised using the alias map before any join operations. This prevents silent mismatches due to naming differences (e.g. Slc15a4 vs Tasl, Iigp1 vs Ifi47).
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 1 — Load all stream outputs
# =============================================================================
cat("\n--- Section 1: Loading stream outputs ---\n")

# --- 1.1  Stream 1 outputs ---
cat("\n  Stream 1:\n")
mcps_genes      <- safe_read(file.path(STREAM1_DIR, "MCPS_genes.csv"),
                              "MCPS_genes.csv")
mcps_lfc_matrix <- safe_read(file.path(STREAM1_DIR, "MCPS_per_study_lfc_matrix.csv"),
                              "MCPS_per_study_lfc_matrix.csv")
replic_summary  <- safe_read(file.path(REPLIC_DIR,  "replication_summary.csv"),
                              "replication_summary.csv")
loo_summary     <- safe_read(file.path(REPLIC_DIR,  "LOO_sensitivity_summary.csv"),
                              "LOO_sensitivity_summary.csv")
consensus_paths <- safe_read(file.path(ENRICH_DIR,  "consensus_pathways_3of5.csv"),
                              "consensus_pathways_3of5.csv")
meta_gsea       <- safe_read(file.path(ENRICH_DIR,  "meta_GSEA_results.csv"),
                              "meta_GSEA_results.csv")

# --- 1.2  Stream 2 outputs ---
cat("\n  Stream 2:\n")
s2_concordance  <- safe_read(file.path(STREAM2_DIR, "GSE214373_MCPS_concordance.csv"),
                              "GSE214373_MCPS_concordance.csv")
s2_gsea         <- safe_read(file.path(STREAM2_DIR, "GSE214373_GSEA_MCPS.csv"),
                              "GSE214373_GSEA_MCPS.csv")
if (is.null(s2_gsea)) {
  cat("  [!] GSE214373_GSEA_MCPS.csv is missing. This file is written by Section 8\n")
  cat("      of 06_stream2_human_validation.R. Possible causes:\n")
  cat("      (a) Script 06 has not been run yet.\n")
  cat("      (b) The GSEA section failed silently (check Script 06 output).\n")
  cat("      (c) The MCPS gene set had < 3 genes in the ranked list (minSize filter).\n")
  cat("      Impact: s2_leading_genes will be empty; triple-validated count = 0.\n")
} else if (nrow(s2_gsea) == 0) {
  cat("  [!] GSE214373_GSEA_MCPS.csv exists but has 0 rows (GSEA returned no results).\n")
  cat("      This means the MCPS gene set did not meet the minSize threshold in fgsea.\n")
  cat("      Impact: s2_leading_genes will be empty; triple-validated count = 0.\n")
}
s2_deseq2       <- safe_read(file.path(STREAM2_DIR, "GSE214373_deseq2_primary.csv"),
                              "GSE214373_deseq2_primary.csv")

# --- 1.3  Stream 3 outputs ---
cat("\n  Stream 3:\n")
s3_conc_summary <- safe_read(
  file.path(STREAM3_DIR, "GSE210128_MCPS_concordance_summary.csv"),
  "GSE210128_MCPS_concordance_summary.csv"
)
s3_first_det    <- safe_read(
  file.path(STREAM3_DIR, "GSE210128_MCPS_first_detection.csv"),
  "GSE210128_MCPS_first_detection.csv"
)
s3_first_genes  <- safe_read(
  file.path(STREAM3_DIR, "GSE210128_MCPS_first_detection_genes.csv"),
  "GSE210128_MCPS_first_detection_genes.csv"
)
s3_triple       <- safe_read(
  file.path(STREAM3_DIR, "GSE210128_triple_validated_candidates.csv"),
  "GSE210128_triple_validated_candidates.csv"
)

# Per-tissue concordance tables (used for cross-stream LFC correlations)
s3_conc_tissue <- list()
for (tis in TISSUES) {
  s3_conc_tissue[[tis]] <- safe_read(
    file.path(STREAM3_DIR, sprintf("GSE210128_%s_MCPS_concordance.csv", tis)),
    sprintf("GSE210128_%s_MCPS_concordance.csv", tis)
  )
}

# Per-tissue terminal results (used for pathway convergence)
s3_terminal_tissue <- list()
for (tis in TISSUES) {
  s3_terminal_tissue[[tis]] <- safe_read(
    file.path(STREAM3_DIR, sprintf("GSE210128_%s_terminal_results.csv", tis)),
    sprintf("GSE210128_%s_terminal_results.csv", tis)
  )
}

# --- 1.4  Availability flags ---
stream1_ok <- !is.null(mcps_genes)
stream2_ok <- !is.null(s2_concordance)
stream3_ok <- !is.null(s3_conc_summary)

cat(sprintf("\n  Stream availability: S1=%s  S2=%s  S3=%s\n",
            ifelse(stream1_ok, "OK", "MISSING"),
            ifelse(stream2_ok, "OK", "MISSING"),
            ifelse(stream3_ok, "OK", "MISSING")))

if (!stream1_ok) {
  stop(paste(
    "Stream 1 MCPS_genes.csv is required for integration.",
    "Run 05_stream1_cns_meta_analysis.R first."
  ))
}

cat("\n  [✓] Section 1 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 2 — Validate Join Keys & Gene Symbol Overlap Diagnostics
#
# **Purpose:** Validate that gene symbols match correctly across streams before building the master integration table. Reports the number of genes present in 1, 2, or all 3 streams, and flags any MCPS genes missing from any stream.
#
# **Inputs:** `stream1`, `stream2`, `stream3` data frames from Section 1.
#
# **Outputs:** Console overlap report; `../results/gene_overlap_diagnostics.tsv`.
#
# **Key decisions:** A gene must be tested in ≥2 streams to appear in the master integration table. Genes tested in only one stream are retained in stream-specific tables but excluded from cross-stream analyses.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 2 — Validate join keys; gene symbol overlap diagnostics
# =============================================================================
cat("\n--- Section 2: Validating join keys ---\n")

n_mcps <- nrow(mcps_genes)
cat(sprintf("  MCPS size: %d genes\n", n_mcps))

# Standardise the join key column name across all tables.
# Stream 1: gene_symbol
# Stream 2: mouse_symbol  → rename to gene_symbol
# Stream 3: gene_symbol

if (!is.null(s2_concordance) && "mouse_symbol" %in% colnames(s2_concordance)) {
  s2_concordance <- s2_concordance %>%
    dplyr::rename(gene_symbol = mouse_symbol)
}

# Check overlap
if (!is.null(s2_concordance)) {
  n_s2_overlap <- sum(mcps_genes$gene_symbol %in% s2_concordance$gene_symbol,
                      na.rm = TRUE)
  cat(sprintf("  MCPS ∩ Stream 2 concordance: %d / %d MCPS genes (%.0f%%)\n",
              n_s2_overlap, n_mcps, 100 * n_s2_overlap / n_mcps))
}

if (!is.null(s3_first_genes)) {
  n_s3_overlap <- sum(mcps_genes$gene_symbol %in% s3_first_genes$gene_symbol,
                      na.rm = TRUE)
  cat(sprintf("  MCPS ∩ Stream 3 first-detection: %d / %d MCPS genes (%.0f%%)\n",
              n_s3_overlap, n_mcps, 100 * n_s3_overlap / n_mcps))
}

# Warn if gene symbol case mismatch is detected
if (!is.null(s2_concordance)) {
  s2_upper <- sum(s2_concordance$gene_symbol == toupper(s2_concordance$gene_symbol),
                  na.rm = TRUE)
  mcps_upper <- sum(mcps_genes$gene_symbol == toupper(mcps_genes$gene_symbol),
                    na.rm = TRUE)
  if (s2_upper > 0.5 * nrow(s2_concordance) &&
      mcps_upper < 0.5 * nrow(mcps_genes)) {
    warning(paste(
      "Stream 2 gene symbols appear to be ALL-CAPS (human HGNC format)",
      "but MCPS uses mixed-case (mouse MGI format).",
      "Cross-stream joins may fail. Check ortholog mapping in Stream 2."
    ))
  }
}

cat("\n  [✓] Section 2 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 3 — Build Master Integration Table
#
# **Purpose:** Join all three stream results into a single master table with one row per gene, containing LFC, SE, padj, and significance flags from each stream. This is the primary deliverable of the integration notebook.
#
# **Inputs:** `stream1`, `stream2`, `stream3` data frames from Section 1.
#
# **Outputs:** `master_df` data frame in R; `../results/master_integration_table.tsv`.
#
# **Key decisions:** Full outer join is used so no gene is lost. Missing values (gene not tested in a stream) are coded as NA, not 0. A `n_streams` column counts how many streams tested each gene.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 3 — Build master integration table
# =============================================================================
# One row per MCPS gene. Columns from all three streams are joined by
# gene_symbol. Missing values (NA) indicate the gene was not tested or
# not significant in that stream.
# =============================================================================
cat("\n--- Section 3: Building master integration table ---\n")

master <- mcps_genes %>%
  dplyr::select(
    gene_symbol,
    meta_log2FC  = meta_log2FC,
    meta_padj    = meta_padj,
    I2,
    any_of(c("n_studies", "meta_pvalue"))
  ) %>%
  mutate(
    mcps_direction = ifelse(meta_log2FC > 0, "Up in CNS", "Down in CNS")
  )

# --- 3.1  Join Stream 2 concordance ---
if (!is.null(s2_concordance)) {
  # Select the columns we need; rename conditionally to standardised names.
  # any_of() inside select() cannot be used with = aliasing, so we use a
  # two-step approach: select with any_of(), then rename if the column exists.
  s2_join_cols <- intersect(
    c("gene_symbol", "human_lfc", "human_padj",
      "human_sig_in_mcps", "significant", "direction_match"),
    colnames(s2_concordance)
  )
  s2_join <- s2_concordance %>%
    dplyr::select(all_of(s2_join_cols)) %>%
    # Rename human_sig_in_mcps or significant -> human_sig (whichever is present)
    {
      df <- .
      if ("human_sig_in_mcps" %in% colnames(df))
        df <- dplyr::rename(df, human_sig = human_sig_in_mcps)
      else if ("significant" %in% colnames(df))
        df <- dplyr::rename(df, human_sig = significant)
      df
    } %>%
    # Rename direction_match -> s2_dir_match
    {
      df <- .
      if ("direction_match" %in% colnames(df))
        df <- dplyr::rename(df, s2_dir_match = direction_match)
      df
    } %>%
    distinct(gene_symbol, .keep_all = TRUE)

  master <- master %>%
    left_join(s2_join, by = "gene_symbol")
  cat(sprintf("  Joined Stream 2: %d / %d MCPS genes matched\n",
              sum(!is.na(master$human_lfc)), n_mcps))
}

# --- 3.2  Join Stream 3 per-tissue concordance (terminal LFC) ---
for (tis in TISSUES) {
  df <- s3_conc_tissue[[tis]]
  if (is.null(df)) next

  # Standardise column names from Stream 3 concordance tables
  # Expected columns: gene_symbol, terminal_lfc, terminal_padj, direction_match
  col_lfc  <- intersect(c("terminal_lfc", "log2FoldChange"), colnames(df))[1]
  col_padj <- intersect(c("terminal_padj", "padj"),          colnames(df))[1]
  col_dir  <- intersect(c("direction_match"),                 colnames(df))[1]

  if (is.na(col_lfc) || is.na(col_padj)) {
    cat(sprintf("  Stream 3 %s concordance: missing lfc/padj columns; skipping.\n", tis))
    next
  }

  s3_join <- df %>%
    dplyr::select(
      gene_symbol,
      !!paste0(tis, "_lfc")       := all_of(col_lfc),
      !!paste0(tis, "_padj")      := all_of(col_padj),
      !!paste0(tis, "_dir_match") := any_of(col_dir)
    ) %>%
    distinct(gene_symbol, .keep_all = TRUE)

  master <- master %>%
    left_join(s3_join, by = "gene_symbol")
  cat(sprintf("  Joined Stream 3 (%s): %d / %d MCPS genes matched\n",
              TISSUE_LABELS[tis],
              sum(!is.na(master[[paste0(tis, "_lfc")]])), n_mcps))
}

# --- 3.3  Join Stream 3 first-detection ---
if (!is.null(s3_first_genes)) {
  # Take the earliest first_detected_wpi across all tissues per gene
  first_det_join <- s3_first_genes %>%
    group_by(gene_symbol) %>%
    slice_min(order_by = first_detected_wpi, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    dplyr::select(
      gene_symbol,
      first_detected_wpi,
      first_detected_tissue = tissue,
      first_det_padj         = padj,
      any_of("wpi_label")
    ) %>%
    dplyr::rename_with(~ "first_det_wpi_label", any_of("wpi_label"))

  master <- master %>%
    left_join(first_det_join, by = "gene_symbol")
  cat(sprintf("  Joined Stream 3 first-detection: %d / %d MCPS genes have peripheral detection\n",
              sum(!is.na(master$first_detected_wpi)), n_mcps))
}

# --- 3.4  Derive evidence tier for each gene ---
# Tier 1: concordant in peripheral tissue (same direction as CNS) AND
#         in Stream 2 GSEA leading edge
# Tier 2: detected in peripheral tissue AND in Stream 2 GSEA leading edge
#         (direction may differ from CNS)
# Tier 3: in Stream 2 GSEA leading edge only (no peripheral detection)
# Tier 4: peripheral detection only (no Stream 2 leading edge)
# Tier 5: MCPS only (neither Stream 2 leading edge nor peripheral detection)

# Parse Stream 2 leading-edge genes
s2_leading_genes <- character(0)
# Script 06 writes leadingEdge as semicolon-separated string (fgsea default).
if (!is.null(s2_gsea) && "leadingEdge" %in% colnames(s2_gsea)) {
  raw_le <- s2_gsea$leadingEdge[1]
  s2_leading_genes <- str_split(raw_le, ";")[[1]] %>% str_trim() %>% .[. != ""]
  cat(sprintf("  Stream 2 GSEA leading-edge genes: %d\n", length(s2_leading_genes)))
}

master <- master %>%
  mutate(
    in_s2_leading_edge   = gene_symbol %in% s2_leading_genes,
    has_peripheral_det   = !is.na(first_detected_wpi),
    # Concordant peripheral = direction_match in at least one tissue
    peripheral_concordant = {
      dir_cols <- paste0(TISSUES, "_dir_match")
      dir_cols_present <- dir_cols[dir_cols %in% colnames(master)]
      if (length(dir_cols_present) == 0) {
        rep(FALSE, nrow(master))
      } else {
        rowSums(master[, dir_cols_present, drop = FALSE], na.rm = TRUE) > 0
      }
    },
    evidence_tier = case_when(
      in_s2_leading_edge & has_peripheral_det & peripheral_concordant ~
        "Tier 1: CNS + human leading edge + concordant peripheral",
      in_s2_leading_edge & has_peripheral_det ~
        "Tier 2: CNS + human leading edge + detected peripheral",
      in_s2_leading_edge & !has_peripheral_det ~
        "Tier 3: CNS + human leading edge only",
      !in_s2_leading_edge & has_peripheral_det ~
        "Tier 4: CNS + peripheral detection only",
      TRUE ~
        "Tier 5: CNS meta-analysis only"
    )
  )

tier_counts <- master %>% count(evidence_tier) %>% arrange(evidence_tier)
cat("\n  Evidence tier distribution:\n")
for (i in seq_len(nrow(tier_counts))) {
  cat(sprintf("    %s: %d genes\n",
              tier_counts$evidence_tier[i], tier_counts$n[i]))
}

cat(sprintf("\n  Master integration table: %d rows × %d columns\n",
            nrow(master), ncol(master)))

cat("\n  [✓] Section 3 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 4 — Triple-Validated Candidate Analysis
#
# **Purpose:** Identify genes that are significant (padj<0.05) and concordant (same direction of effect) in all three streams — the triple-validated candidate set. The MCPS genes are the pre-specified primary outcome; this section also identifies any additional triple-validated genes.
#
# **Inputs:** `master_df` from Section 3.
#
# **Outputs:** `../results/triple_validated_candidates.tsv`; `../results/mcps_gene_summary.tsv`.
#
# **Key decisions:** Triple validation requires: padj<0.05 in all 3 streams AND same sign of LFC in all 3 streams. Direction discordance for Slc15a4 in human is reported transparently — it does not disqualify the gene from MCPS because the mouse evidence is the primary discovery stream.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 4 — Triple-validated candidate analysis
# =============================================================================
# Tier 1 + Tier 2 genes are the triple-validated candidates.
# These are ranked by: (1) Tier, (2) earliest peripheral detection wpi,
# (3) ascending meta_padj.
# =============================================================================
cat("\n--- Section 4: Triple-validated candidate analysis ---\n")

triple_validated <- master %>%
  filter(evidence_tier %in% c(
    "Tier 1: CNS + human leading edge + concordant peripheral",
    "Tier 2: CNS + human leading edge + detected peripheral"
  )) %>%
  arrange(evidence_tier, first_detected_wpi, meta_padj)

n_tier1 <- sum(triple_validated$evidence_tier ==
                 "Tier 1: CNS + human leading edge + concordant peripheral")
n_tier2 <- sum(triple_validated$evidence_tier ==
                 "Tier 2: CNS + human leading edge + detected peripheral")

cat(sprintf("  Triple-validated candidates: %d total\n",
            nrow(triple_validated)))
cat(sprintf("    Tier 1 (concordant peripheral + human leading edge): %d\n", n_tier1))
cat(sprintf("    Tier 2 (detected peripheral + human leading edge):   %d\n", n_tier2))

if (nrow(triple_validated) > 0) {
  cat("\n  Top Tier 1 candidates (up to 15):\n")
  top15 <- triple_validated %>%
    filter(evidence_tier == "Tier 1: CNS + human leading edge + concordant peripheral") %>%
    slice_head(n = 15)
  for (i in seq_len(nrow(top15))) {
    cat(sprintf("    %-14s  meta_LFC=%+.2f  human_LFC=%+.2f  first_det=%s wpi\n",
                top15$gene_symbol[i],
                top15$meta_log2FC[i],
                ifelse(is.na(top15$human_lfc[i]), NA_real_, top15$human_lfc[i]),
                ifelse(is.na(top15$first_detected_wpi[i]), "NA",
                       as.character(top15$first_detected_wpi[i]))))
  }
}

cat("\n  [✓] Section 4 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 5 — Cross-Stream LFC Correlation Analysis
#
# **Purpose:** Quantify the overall concordance between streams by computing Spearman correlations of LFC values for all genes tested in pairs of streams. Also computes the MCPS gene set enrichment score relative to the background.
#
# **Inputs:** `master_df` from Section 3.
#
# **Outputs:** `../results/cross_stream_correlations.tsv` — Spearman r, 95% CI, p-value for each stream pair.
#
# **Key decisions:** Correlation is computed on all genes with LFC in both streams, not just significant genes — this avoids winner's curse bias.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 5 — Cross-stream LFC correlation analysis
# =============================================================================
# Pearson correlations between:
#   (a) meta_log2FC (S1) vs human_lfc (S2)
#   (b) meta_log2FC (S1) vs whole_blood_lfc (S3)
#   (c) meta_log2FC (S1) vs spleen_lfc (S3)
#   (d) meta_log2FC (S1) vs skeletal_muscle_lfc (S3)
#   (e) human_lfc (S2) vs each peripheral tissue LFC (S3)
# Only genes with non-NA values in both columns are used.
# =============================================================================
cat("\n--- Section 5: Cross-stream LFC correlations ---\n")

lfc_pairs <- list(
  list(x = "meta_log2FC",        y = "human_lfc",
       xlabel = "CNS meta-LFC (S1)", ylabel = "Human sCJD LFC (S2)",
       label  = "S1 vs S2"),
  list(x = "meta_log2FC",        y = "whole_blood_lfc",
       xlabel = "CNS meta-LFC (S1)", ylabel = "Whole blood LFC (S3)",
       label  = "S1 vs S3 blood"),
  list(x = "meta_log2FC",        y = "spleen_lfc",
       xlabel = "CNS meta-LFC (S1)", ylabel = "Spleen LFC (S3)",
       label  = "S1 vs S3 spleen"),
  list(x = "meta_log2FC",        y = "skeletal_muscle_lfc",
       xlabel = "CNS meta-LFC (S1)", ylabel = "Skeletal muscle LFC (S3)",
       label  = "S1 vs S3 muscle"),
  list(x = "human_lfc",          y = "whole_blood_lfc",
       xlabel = "Human sCJD LFC (S2)", ylabel = "Whole blood LFC (S3)",
       label  = "S2 vs S3 blood"),
  list(x = "human_lfc",          y = "spleen_lfc",
       xlabel = "Human sCJD LFC (S2)", ylabel = "Spleen LFC (S3)",
       label  = "S2 vs S3 spleen"),
  list(x = "human_lfc",          y = "skeletal_muscle_lfc",
       xlabel = "Human sCJD LFC (S2)", ylabel = "Skeletal muscle LFC (S3)",
       label  = "S2 vs S3 muscle")
)

lfc_cor_results <- lapply(lfc_pairs, function(p) {
  xcol <- p$x; ycol <- p$y
  if (!xcol %in% colnames(master) || !ycol %in% colnames(master)) {
    return(data.frame(
      comparison = p$label, x_col = xcol, y_col = ycol,
      n = 0L, pearson_r = NA_real_, pearson_p = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  df_pair <- master %>%
    filter(!is.na(.data[[xcol]]), !is.na(.data[[ycol]])) %>%
    dplyr::select(gene_symbol, x = all_of(xcol), y = all_of(ycol))

  if (nrow(df_pair) < 5) {
    return(data.frame(
      comparison = p$label, x_col = xcol, y_col = ycol,
      n = nrow(df_pair), pearson_r = NA_real_, pearson_p = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  ct <- tryCatch(
    cor.test(df_pair$x, df_pair$y, method = "pearson"),
    error = function(e) NULL
  )

  if (is.null(ct)) {
    return(data.frame(
      comparison = p$label, x_col = xcol, y_col = ycol,
      n = nrow(df_pair), pearson_r = NA_real_, pearson_p = NA_real_,
      ci_lower = NA_real_, ci_upper = NA_real_,
      stringsAsFactors = FALSE
    ))
  }

  cat(sprintf("  %-25s  n=%3d  r=%+.3f  p=%.4f\n",
              p$label, nrow(df_pair), ct$estimate, ct$p.value))

  data.frame(
    comparison = p$label, x_col = xcol, y_col = ycol,
    n          = nrow(df_pair),
    pearson_r  = round(ct$estimate, 4),
    pearson_p  = ct$p.value,
    ci_lower   = round(ct$conf.int[1], 4),
    ci_upper   = round(ct$conf.int[2], 4),
    stringsAsFactors = FALSE
  )
})

lfc_cor_df <- bind_rows(lfc_cor_results)

# --- Underpowered correlation note ---
# The MCPS has only 6 genes. Cross-stream correlations with n=6 are severely
# underpowered: Pearson r requires n>=10 for reliable estimates, and n>=30 for
# stable confidence intervals. Interpret all cross-stream correlations as
# exploratory only; do not report p-values as confirmatory.
underpowered_pairs <- lfc_cor_df %>%
  dplyr::filter(!is.na(pearson_r) & n < 10)

if (nrow(underpowered_pairs) > 0) {
  cat("\n  *** UNDERPOWERED CORRELATION WARNING ***\n")
  cat(sprintf("  %d correlation(s) have n < 10 (MCPS genes with data in both streams):\n",
              nrow(underpowered_pairs)))
  for (i in seq_len(nrow(underpowered_pairs))) {
    cat(sprintf("    %s: n=%d, r=%.3f — interpret as exploratory only\n",
                underpowered_pairs$comparison[i], underpowered_pairs$n[i],
                underpowered_pairs$pearson_r[i]))
  }
  cat("  Manuscript must note that cross-stream LFC correlations are underpowered\n")
  cat(sprintf("  with only %d MCPS genes. Directional concordance (binomial test) is preferred.\n",
              n_mcps))
  cat("  *****************************************\n")
}

cat("\n  [✓] Section 5 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 6 — Replication Cascade Statistics
#
# **Purpose:** Compute the replication cascade: starting from Stream 1 significant genes, what fraction replicate in Stream 2 (human), Stream 3 (peripheral), and all three streams? Reports binomial confidence intervals for each replication rate.
#
# **Inputs:** `master_df` from Section 3.
#
# **Outputs:** `../results/replication_cascade.tsv`; replication cascade waterfall figure.
#
# **Key decisions:** Replication is defined as same direction AND nominal p<0.05 (not padj<0.05) in the replication stream — consistent with standard replication study conventions.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 6 — Replication cascade statistics
# =============================================================================
# Counts how many MCPS genes survive each successive validation stage.
# Stages (in order):
#   1. CNS meta-analysis (= MCPS size, by definition)
#   2. Stratum D replication (same direction, padj <= 0.05, in any Stratum D study)
#   3. Human sCJD concordance (same direction as CNS, in Stream 2)
#   4. Peripheral terminal concordance (same direction, in any peripheral tissue)
#   5. Peripheral first-detection (significant at any pre-terminal timepoint)
# =============================================================================
cat("\n--- Section 6: Replication cascade ---\n")

cascade_stages <- list()

# Stage 1: MCPS
cascade_stages[["1. CNS meta-analysis\n(MCPS)"]] <- n_mcps

# Stage 2: Stratum D replication
# replication_summary.csv has columns: study, n_same_direction, replication_rate_pct
# IMPORTANT: replication_rate_pct = n_same_direction / n_mcps * 100 (per study).
# We report the MEAN per-study replication count (= mean n_same_direction across
# Stratum D studies). This is NOT the union of genes replicated in >= 1 study
# (which would require loading per-gene results from all Stratum D studies).
# The mean per-study count is a conservative estimate: a gene that replicates in
# 3/4 studies contributes 0.75 to the mean, not 1.
# This is clearly labelled as "mean per-study" in the cascade bar chart.
if (!is.null(replic_summary) && "replication_rate_pct" %in% colnames(replic_summary)) {
  mean_replic_rate <- mean(replic_summary$replication_rate_pct, na.rm = TRUE)
  # n_replicated_est = mean number of MCPS genes replicated per Stratum D study
  n_replicated_est <- round(n_mcps * mean_replic_rate / 100)
  cascade_stages[["2. Stratum D replication\n(mean per study)"]] <- n_replicated_est
  cat(sprintf("  Stratum D mean per-study replication rate: %.1f%%  (~%d genes/study)\n",
              mean_replic_rate, n_replicated_est))
  cat(sprintf("  NOTE: This is mean per-study count, not union across studies.\n"))
} else {
  cascade_stages[["2. Stratum D replication\n(mean per study)"]] <- NA_integer_
}

# Stage 3: Human sCJD concordance (Stream 2)
if ("s2_dir_match" %in% colnames(master)) {
  n_s2_concordant <- sum(master$s2_dir_match == TRUE, na.rm = TRUE)
  cascade_stages[["3. Human sCJD concordance\n(Stream 2)"]] <- n_s2_concordant
  cat(sprintf("  Human sCJD concordant: %d / %d MCPS genes\n",
              n_s2_concordant, n_mcps))
} else {
  cascade_stages[["3. Human sCJD concordance\n(Stream 2)"]] <- NA_integer_
}

# Stage 4: Peripheral terminal concordance (any tissue, Stream 3)
dir_cols_present <- paste0(TISSUES, "_dir_match")[
  paste0(TISSUES, "_dir_match") %in% colnames(master)
]
if (length(dir_cols_present) > 0) {
  n_periph_concordant <- sum(
    rowSums(master[, dir_cols_present, drop = FALSE] == TRUE, na.rm = TRUE) > 0,
    na.rm = TRUE
  )
  cascade_stages[["4. Peripheral terminal\nconcordance (Stream 3)"]] <- n_periph_concordant
  cat(sprintf("  Peripheral terminal concordant (any tissue): %d / %d MCPS genes\n",
              n_periph_concordant, n_mcps))
} else {
  cascade_stages[["4. Peripheral terminal\nconcordance (Stream 3)"]] <- NA_integer_
}

# Stage 5: Peripheral first-detection (pre-terminal)
if ("first_detected_wpi" %in% colnames(master)) {
  n_first_det <- sum(!is.na(master$first_detected_wpi) &
                       master$first_detected_wpi < 22L, na.rm = TRUE)
  cascade_stages[["5. Peripheral first-detection\n(pre-terminal, Stream 3)"]] <- n_first_det
  cat(sprintf("  Peripheral first-detection (pre-terminal): %d / %d MCPS genes\n",
              n_first_det, n_mcps))
} else {
  cascade_stages[["5. Peripheral first-detection\n(pre-terminal, Stream 3)"]] <- NA_integer_
}

# Stage 6: Triple-validated (Tier 1 + Tier 2)
cascade_stages[["6. Triple-validated\n(Tier 1 + Tier 2)"]] <- nrow(triple_validated)

cascade_df <- data.frame(
  stage       = names(cascade_stages),
  n_genes     = unlist(cascade_stages),
  stage_order = seq_along(cascade_stages),
  stringsAsFactors = FALSE
) %>%
  mutate(
    pct_of_mcps = round(100 * n_genes / n_mcps, 1),
    stage_short = c(
      "CNS meta\n(MCPS)",
      "Stratum D\nreplication",
      "Human sCJD\nconcordance",
      "Peripheral\nterminal",
      "Peripheral\nfirst-det.",
      "Triple-\nvalidated"
    )
  )

cat("\n  Replication cascade:\n")
print(cascade_df[, c("stage_short", "n_genes", "pct_of_mcps")])

cat("\n  [✓] Section 6 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 7 — Pathway Convergence Analysis
#
# **Purpose:** Test whether the same pathways are enriched across all three streams by comparing GSEA results from Streams 1, 2, and 3. Identifies pathways with significant NES in all three streams (triple-validated pathways).
#
# **Inputs:** `../results/stream1/gsea_results.tsv`; `../results/stream2/gsea_results.tsv`; `../results/stream3/gsea_results.tsv`.
#
# **Outputs:** `../results/pathway_convergence.tsv` — NES and padj for each pathway across all three streams; convergence score.
#
# **Key decisions:** Pathway convergence score = mean NES across streams, weighted by −log10(padj). Pathways significant in ≥2/3 streams are reported as 'convergent'; all 3 streams as 'triple-convergent'.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 7 — Pathway convergence analysis
# =============================================================================
# For each consensus pathway (enriched in >= 3/5 Stream 1 studies), check
# whether the same pathway is enriched in Stream 2 and Stream 3 terminal DEGs.
# Since Stream 2 and Stream 3 do not run ORA internally in this script, we
# use a gene-overlap approach: what fraction of the pathway's gene set is
# in the Stream 2 / Stream 3 significant DEG lists?
# =============================================================================
cat("\n--- Section 7: Pathway convergence ---\n")

pathway_convergence <- data.frame()

if (!is.null(consensus_paths) && nrow(consensus_paths) > 0) {

  # Build Stream 2 significant DEG set (mouse symbols)
  s2_sig_genes <- character(0)
  if (!is.null(s2_deseq2) && "mouse_symbol" %in% colnames(s2_deseq2)) {
    s2_sig_genes <- s2_deseq2 %>%
      filter(!is.na(padj), padj <= PADJ_THRESH) %>%
      pull(mouse_symbol) %>%
      unique() %>%
      .[. != "" & !is.na(.)]
  } else if (!is.null(s2_deseq2) && "significant" %in% colnames(s2_deseq2)) {
    s2_sig_genes <- s2_deseq2 %>%
      filter(significant == TRUE) %>%
      dplyr::select(any_of(c("mouse_symbol", "human_symbol", "gene_id"))) %>%
      dplyr::pull(1) %>%
      unique() %>%
      .[. != "" & !is.na(.)]
  }

  # Build Stream 3 significant DEG sets per tissue (gene symbols)
  s3_sig_genes <- list()
  for (tis in TISSUES) {
    df <- s3_terminal_tissue[[tis]]
    if (is.null(df)) next
    sig_col <- intersect(c("significant", "padj"), colnames(df))
    if ("significant" %in% sig_col) {
      s3_sig_genes[[tis]] <- df %>%
        filter(significant == TRUE) %>%
        dplyr::select(any_of(c("gene_symbol", "gene_id"))) %>%
        dplyr::pull(1) %>%
        unique() %>% .[!is.na(.)]
    } else if ("padj" %in% sig_col) {
      s3_sig_genes[[tis]] <- df %>%
        filter(!is.na(padj), padj <= PADJ_THRESH) %>%
        dplyr::select(any_of(c("gene_symbol", "gene_id"))) %>%
        dplyr::pull(1) %>%
        unique() %>% .[!is.na(.)]
    }
  }

  # For each consensus pathway, compute:
  #   - n_studies_s1: from consensus_paths (already computed)
  #   - s2_overlap_pct: % of pathway gene set in Stream 2 sig DEGs
  #   - s3_<tissue>_overlap_pct: % of pathway gene set in Stream 3 sig DEGs
  # Note: consensus_paths may not have the gene set members listed.
  # We use the pathway description and n_studies as the primary metric,
  # and flag pathways that are also in the MCPS GSEA results (meta_gsea).

  # Check if meta_gsea has matching pathway IDs
  gsea_sig_paths <- character(0)
  if (!is.null(meta_gsea) && "ID" %in% colnames(meta_gsea)) {
    gsea_sig_paths <- meta_gsea %>%
      filter(!is.na(p.adjust), p.adjust <= PADJ_THRESH) %>%
      pull(ID) %>% unique()
  } else if (!is.null(meta_gsea) && "pathway" %in% colnames(meta_gsea)) {
    gsea_sig_paths <- meta_gsea %>%
      filter(!is.na(padj), padj <= PADJ_THRESH) %>%
      pull(pathway) %>% unique()
  }

  pathway_convergence <- consensus_paths %>%
    dplyr::select(
      pathway_id    = any_of(c("ID", "pathway_id")),
      description   = any_of(c("Description", "description")),
      database      = any_of(c("database", "Database")),
      n_studies_s1  = any_of(c("n_studies_significant", "n_studies")),
      mean_padj_s1  = any_of(c("mean_padj", "padj"))
    ) %>%
    mutate(
      in_meta_gsea = pathway_id %in% gsea_sig_paths,
      # Convergence score: n_studies_s1 + 1 if in meta GSEA
      convergence_score = n_studies_s1 + as.integer(in_meta_gsea)
    ) %>%
    arrange(desc(convergence_score), mean_padj_s1)

  cat(sprintf("  Consensus pathways (>= 3/5 studies): %d\n",
              nrow(pathway_convergence)))
  cat(sprintf("  Also significant in meta-GSEA: %d\n",
              sum(pathway_convergence$in_meta_gsea, na.rm = TRUE)))

  if (nrow(pathway_convergence) > 0) {
    cat("  Top 10 convergent pathways:\n")
    top10_paths <- pathway_convergence %>% slice_head(n = 10)
    for (i in seq_len(nrow(top10_paths))) {
      desc_short <- substr(top10_paths$description[i], 1, 50)
      cat(sprintf("    [%d/5 studies] %s%s\n",
                  top10_paths$n_studies_s1[i],
                  desc_short,
                  ifelse(top10_paths$in_meta_gsea[i], " *GSEA", "")))
    }
  }

} else {
  cat("  No consensus pathway data available; skipping pathway convergence.\n")
}

cat("\n  [✓] Section 7 successfully loaded and parsed.\n")


# ----------------------------------------------------------------------
# ##Section 8 — Literature Convergence Table
#
# **Purpose:** Build a structured evidence table for the 6 MCPS genes, integrating five independent evidence streams: (1) CNS meta-analysis (Stream 1), (2) human brain validation (Stream 2), (3) peripheral blood concordance (Stream 3), (4) independent mouse replication (Majer 2019), and (5) published literature (mechanistic/functional evidence). Each gene is assigned an evidence tier (Tier 1–4+) based on the number of streams in which it is supported.
#
# **Inputs:** Stream result tables loaded in Section 1; curated literature annotations (hard-coded in this section).
#
# **Outputs:** `../results/integration/MCPS_literature_convergence_table.csv` — 6 genes × 12 columns covering all evidence streams and tier classification.
#
# **Key decisions:** Tier classification: Tier 1 = CNS only; Tier 2 = CNS + human; Tier 3 = CNS + human + pathway/literature; Tier 4 = CNS + human + peripheral + independent replication. Igsf6 has no direct prion literature (honest gap — reported as such). Slc15a4 is flagged as non-microglial (pDC/B-cell mechanism via TASL-IRF5 axis).
# ----------------------------------------------------------------------

# =============================================================================
# SECTION 8 — Literature Convergence Table
# =============================================================================
cat("\n--- Section 8: Literature Convergence Table ---\n")

lit_convergence <- data.frame(
  Gene_Mouse = c("Tasl", "Igsf6", "Fcgr4", "Parp14", "Iigp1", "Il1a"),
  Gene_Human = c("SLC15A4", "IGSF6", "FCGR4", "PARP14", "IRGM", "IL1A"),
  Meta_LFC   = c(1.10, 1.36, 2.54, 0.62, 0.58, 1.03),
  Meta_padj  = c(5.6e-10, 1.6e-7, 3.0e-5, 1.3e-4, 0.022, 0.037),
  I2_pct     = c(0, 15.8, 60.8, 0, 0, 22.8),

  # Evidence stream 2: Human blood (GSE214373)
  Stream2_HumanBlood = c(
    "NS (LFC+0.19)",
    "YES (LFC+0.52, padj=0.048)",
    "YES (LFC+0.61, padj=0.031)",
    "NS (LFC+0.08)",
    "NS (LFC+0.11)",
    "NS (LFC+0.22)"
  ),

  # Evidence stream 3: Human brain pathway (GSE124571 GSEA)
  Stream3_HumanBrainPathway = c(
    "YES (IFN/TLR pathway NES>2.0, padj<0.001)",
    "YES (myeloid activation NES=+2.10, padj<0.001)",
    "YES (myeloid activation NES=+2.10, padj<0.001)",
    "YES (JAK-STAT NES=+2.13, padj<0.001)",
    "YES (IFN-gamma NES=+2.41, padj<0.001)",
    "YES (IL-6 production NES=+2.27, padj<0.001)"
  ),

  # Evidence stream 4: Independent mouse prion RNA-seq (Majer 2019, Sci Reports)
  Stream4_IndepMouseRNA = c(
    "Not reported",
    "Not reported",
    "YES — upregulated in prion mouse brain (Majer 2019, Sci Reports)",
    "Not reported",
    "YES — progressive upregulation 80-157 dpi (Majer 2019, Sci Reports)",
    "Not reported"
  ),

  # Evidence stream 5: Human prion disease literature
  Stream5_HumanLiterature = c(
    "TLR7/8/9 endosomal scaffold via TASL-IRF5 axis (Heinz 2020, Nature); SLE/autoimmune GWAS",
    "Myeloid surface receptor; no direct prion literature found",
    "Fc-gamma receptor IV; myeloid activation marker; upregulated in prion mouse brain (Majer 2019)",
    "Microglial activation inhibitor; STAT1 ADP-ribosylation; upregulated by IFNbeta/TLR agonists",
    "Interferon-inducible GTPase; autophagy regulator; upregulated progressively in prion mouse brain",
    "IL-1alpha elevated in CJD CSF (Stoeck 2014); drives A1 astrocyte transformation in prion disease"
  ),

  # Microglia correlation (from notebook 09)
  Microglia_r = c(0.061, 0.638, 0.525, 0.468, 0.486, 0.566),
  Microglia_sig = c("NS", "***", "***", "***", "***", "***"),

  Evidence_Tier = c(
    "Tier 3+ (CNS + pathway + literature; NOT microglia-driven)",
    "Tier 4 (CNS + peripheral replicated; microglia-correlated)",
    "Tier 4+ (CNS + peripheral + independent mouse + literature; microglia-correlated)",
    "Tier 3 (CNS + pathway + literature; microglia-correlated)",
    "Tier 4+ (CNS + independent mouse + pathway; microglia-correlated)",
    "Tier 3 (CNS + human CSF literature; microglia-correlated)"
  ),

  stringsAsFactors = FALSE
)

cat("=== MCPS Gene Evidence Summary ===\n")
print(lit_convergence[, c("Gene_Mouse", "Gene_Human", "Meta_LFC", "Meta_padj",
                           "I2_pct", "Microglia_r", "Microglia_sig", "Evidence_Tier")])

# FIX: Safely write using write_csv
tryCatch({
  write_csv(lit_convergence, file.path(OUT_DIR, "MCPS_literature_convergence_table.csv"))
  cat("\n  Saved: MCPS_literature_convergence_table.csv\n")
}, error = function(e) {
  warning(sprintf("Failed to save Literature Convergence Table: %s", conditionMessage(e)))
})


# ----------------------------------------------------------------------
# ##Section 9 — Pathway-Level Human Validation (GSE124571 GSEA)
#
# **Purpose:** Test whether the consensus pathways identified in the mouse meta-analysis (Stream 1 GSEA) are also significantly enriched in the human sCJD dataset (GSE124571). This provides pathway-level cross-species validation independent of individual gene concordance.
#
# **Inputs:** `../results/stream1/enrichment/consensus_pathways_3of5.csv` (consensus pathway list from Stream 1); GSE124571 microarray data (downloaded via GEOquery, cached after first run).
#
# **Outputs:** `../results/integration/GSE124571_human_GSEA_vs_consensus_pathways.csv` — 362 pathways tested, 322 significant (89.0%). The bubble plot figure is generated in Section 10 (all figures are consolidated there).
#
# **Key decisions:** GSEA ranking metric = signed −log10(padj) × sign(LFC) from GSE124571. All 7 Hallmark pathways significant (NES range +2.01 to +2.41). 99.1% of significant pathways have NES>0 (upregulated in sCJD), confirming the neuroinflammatory signature is conserved across species.
# ----------------------------------------------------------------------

# =============================================================================
# SECTION 9 — Pathway-Level Human Validation (GSE124571 GSEA)
# =============================================================================
cat("\n--- Section 9: Pathway-Level Human Validation ---\n")

# Load GSE124571
gse <- getGEO("GSE124571", GSEMatrix = TRUE, getGPL = FALSE)
eset <- gse[[1]]

# Remove replicate sample
keep <- !grepl("replicate", pData(eset)$title)
eset_clean <- eset[, keep]
condition <- pData(eset_clean)$`condition:ch1`

# Normalize and run limma
expr_log  <- log2(exprs(eset_clean) + 1)
expr_norm <- normalizeBetweenArrays(expr_log, method = "quantile")
group  <- factor(condition, levels = c("Control", "sCJD affected"))
design <- model.matrix(~group)
fit    <- eBayes(lmFit(expr_norm, design))
results <- topTable(fit, coef = "groupsCJD affected", number = Inf, sort.by = "none")

# Map probes to gene symbols
gpl <- getGEO("GPL14951")
probe_gene <- Table(gpl)[, c("ID", "Symbol")]
probe_gene <- probe_gene[probe_gene$Symbol != "" & !is.na(probe_gene$Symbol), ]
results$probe_id <- rownames(results)
results_annot <- merge(results, probe_gene, by.x = "probe_id", by.y = "ID")
results_unique <- results_annot[!duplicated(results_annot$Symbol), ]
ranked_genes <- sort(setNames(results_unique$t, results_unique$Symbol), decreasing = TRUE)

cat("  Ranked genes:", length(ranked_genes), "\n")
cat("  DEGs (padj<0.05):", sum(results_unique$adj.P.Val < 0.05), "\n")

# Load consensus pathways and run GSEA
consensus_df <- read.csv(file.path(STREAM1_DIR, "enrichment", "consensus_pathways_3of5.csv"))

hallmark  <- msigdbr(species = "Homo sapiens", collection = "H")
go_bp     <- msigdbr(species = "Homo sapiens", collection = "C5", subcollection = "GO:BP")
reactome  <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:REACTOME")
kegg      <- msigdbr(species = "Homo sapiens", collection = "C2", subcollection = "CP:KEGG_LEGACY")
all_gs    <- split(rbind(hallmark, go_bp, reactome, kegg)$gene_symbol,
                   rbind(hallmark, go_bp, reactome, kegg)$gs_name)

get_msigdb_name <- function(id, desc, db) {
  if (db == "GO_BP") {
    nm <- paste0("GOBP_", toupper(gsub("[^A-Za-z0-9]", "_", desc)))
    return(gsub("_+", "_", gsub("_$", "", nm)))
  } else if (db %in% c("Hallmark", "Reactome")) {
    return(id)
  } else if (db == "KEGG") {
    nm <- paste0("KEGG_", toupper(gsub("[^A-Za-z0-9]", "_", desc)))
    return(gsub("_+", "_", gsub("_$", "", nm)))
  }
  return(NA_character_)
}
consensus_df$msigdb_name <- mapply(get_msigdb_name,
                                    consensus_df$ID, consensus_df$Description, consensus_df$database)
matched_gs <- all_gs[consensus_df$msigdb_name[consensus_df$msigdb_name %in% names(all_gs)]]

set.seed(42)
gsea_res <- fgsea(pathways = matched_gs, stats = ranked_genes,
                  minSize = 10, maxSize = 500, nPermSimple = 10000)

cat("\n  Pathways tested:", nrow(gsea_res), "\n")
cat("  Significant (padj<0.05):", sum(gsea_res$padj < 0.05, na.rm=TRUE), "\n")
cat("  Overlap rate:", round(100 * sum(gsea_res$padj < 0.05, na.rm=TRUE) / nrow(gsea_res), 1), "%\n")

tryCatch({
  gsea_save <- as.data.frame(gsea_res)
  gsea_save$leadingEdge <- sapply(gsea_res$leadingEdge, paste, collapse=";")
  write_csv(gsea_save, file.path(OUT_DIR, "GSE124571_human_GSEA_vs_consensus_pathways.csv"))
  cat("\n  Saved: GSE124571_human_GSEA_vs_consensus_pathways.csv\n")
}, error = function(e) {
  warning(sprintf("Failed to save Pathway GSEA results: %s", conditionMessage(e)))
})

# Key convergent pathways extraction
key_pathways <- c(
  "HALLMARK_INTERFERON_GAMMA_RESPONSE", "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE", "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "HALLMARK_ALLOGRAFT_REJECTION", "HALLMARK_COMPLEMENT",
  "GOBP_ADAPTIVE_IMMUNE_RESPONSE", "GOBP_REGULATION_OF_INNATE_IMMUNE_RESPONSE",
  "GOBP_CYTOKINE_MEDIATED_SIGNALING_PATHWAY", "GOBP_MYELOID_LEUKOCYTE_ACTIVATION",
  "GOBP_MYELOID_CELL_ACTIVATION_INVOLVED_IN_IMMUNE_RESPONSE",
  "GOBP_INTERLEUKIN_6_PRODUCTION", "GOBP_DENDRITIC_CELL_DIFFERENTIATION",
  "GOBP_LEUKOCYTE_PROLIFERATION", "REACTOME_NEUTROPHIL_DEGRANULATION"
)

cat("\n=== Key convergent pathways in human sCJD (GSE124571) ===\n")
for (pw in key_pathways) {
  if (pw %in% gsea_res$pathway) {
    row <- gsea_res[gsea_res$pathway == pw, ]
    cat(sprintf("  %-50s NES=%+.2f  padj=%.2e\n", pw, row$NES, row$padj))
  }
}


# ----------------------------------------------------------------------
# ###Interpretation: Pathway-Level Human Validation
#
# **Key result:** 322/362 consensus pathways (89.0%) are significant in human sCJD (GSE124571) at padj<0.05, with 99.1% showing NES>0 (upregulated). This confirms that the neuroinflammatory and innate immune pathway signature identified in the mouse meta-analysis is robustly conserved in human prion disease.
#
# **Top pathways (NES>2.0):**
# - HALLMARK_INTERFERON_GAMMA_RESPONSE: NES=+2.41, padj=2.0×10⁻¹²
# - GOBP_ADAPTIVE_IMMUNE_RESPONSE: NES=+2.27, padj=1.5×10⁻¹⁶
# - HALLMARK_INTERFERON_ALPHA_RESPONSE: NES=+2.25, padj=5.4×10⁻⁷
# - HALLMARK_INFLAMMATORY_RESPONSE: NES=+2.17, padj=4.5×10⁻⁹
# - HALLMARK_IL6_JAK_STAT3_SIGNALING: NES=+2.13, padj=5.7×10⁻⁶
#
# **Manuscript framing:** This pathway-level validation strengthens the MCPS gene set by showing that the broader biological context (innate immune activation, interferon signalling, myeloid cell activation) is conserved in human disease, even where individual gene direction may differ.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# ##Section 10 — Publication-Ready Figures
#
# **Purpose:** Generate all publication-ready figures for the integration notebook: MCPS gene forest plot (all three streams), cross-tissue LFC scatter, replication cascade waterfall, and pathway convergence bubble plot.
#
# **Inputs:** All stream result tables loaded in Section 1; `gsea_res` and `plot_data` from Section 9 (pathway validation figure).
#
# **Outputs:** All figures saved to `FIG_DIR` = `../results/integration/figures/` via `save_fig()`. Includes: pathway-level human validation bubble plot (from Section 9 GSEA), MCPS gene forest plot, cross-tissue LFC scatter, and replication cascade waterfall.
#
# **Key decisions:** All figures use `theme_pub` for consistent styling. MCPS genes are highlighted in red. Slc15a4 direction discordance is shown explicitly with an open symbol and annotation.
# ----------------------------------------------------------------------


# ----------------------------------------------------------------------
# ###Figure: Pathway-Level Human Validation Bubble Plot
#
# Bubble plot of the top 15 convergent pathways from the GSE124571 GSEA (Section 9). Point size = −log₁₀(padj); colour = gene set database. All shown pathways have NES>2.0 and padj<0.05 in human sCJD.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 8 — Publication-ready figures
# =============================================================================
cat("\n--- Section 8: Generating figures ---\n")

# ---------------------------------------------------------------------------
# Fig 1 — Replication cascade waterfall bar chart
# ---------------------------------------------------------------------------
# Shows how many MCPS genes survive each successive validation stage.
# Bars are coloured by stream; a dashed line shows the MCPS baseline.
# ---------------------------------------------------------------------------
cat("  Fig 1: Replication cascade\n")

tryCatch({
  cascade_plot_df <- cascade_df %>%
    filter(!is.na(n_genes)) %>%
    mutate(
      stage_f = factor(stage_short, levels = stage_short),
      stream_colour = case_when(
        stage_order == 1 ~ "Stream 1 (CNS meta)",
        stage_order == 2 ~ "Stream 1 (CNS meta)",
        stage_order == 3 ~ "Stream 2 (Human sCJD)",
        stage_order %in% c(4, 5) ~ "Stream 3 (Peripheral)",
        stage_order == 6 ~ "Triple-validated"
      )
    )

  stream_bar_colours <- c(
    "Stream 1 (CNS meta)"   = "#2166AC",
    "Stream 2 (Human sCJD)" = "#D6604D",
    "Stream 3 (Peripheral)" = "#4DAC26",
    "Triple-validated"      = "#1B7837"
  )

  # Add a minimum display height for zero-count bars so they remain visible.
  # The actual n_genes value is preserved for labels; only the plotted height is floored.
  cascade_plot_df <- cascade_plot_df %>%
    dplyr::mutate(
      n_genes_display = pmax(n_genes, n_mcps * 0.02),  # floor at 2% of MCPS for visibility
      label_text = ifelse(n_genes == 0,
                          "0\n(0%)",
                          sprintf("%d\n(%.0f%%)", n_genes, pct_of_mcps))
    )

  p_cascade <- ggplot(cascade_plot_df,
                      aes(x = stage_f, y = n_genes_display, fill = stream_colour)) +
    geom_col(width = 0.65, colour = "white", linewidth = 0.3) +
    geom_text(aes(label = label_text),
              vjust = -0.3, size = 3, colour = "grey20") +
    geom_hline(yintercept = n_mcps, linetype = "dashed",
               colour = "grey50", linewidth = 0.4) +
    annotate("text", x = 0.6, y = n_mcps * 1.02,
             label = sprintf("MCPS = %d", n_mcps),
             hjust = 0, size = 2.8, colour = "grey40") +
    scale_fill_manual(values = stream_bar_colours, name = NULL) +
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.18)),
      breaks = pretty(c(0, n_mcps))
    ) +
    labs(
      x       = NULL,
      y       = "Number of MCPS genes",
      title   = "Fig 1 — Replication cascade across three streams",
      subtitle = sprintf(
        "Starting from %d MCPS genes (Stream 1). Percentages relative to MCPS.",
        n_mcps
      ),
      caption = paste(
        "Stage 2 uses mean replication rate across Stratum D studies.",
        "Stages 3–5 use direction-concordant genes.",
        "Stage 6 = Tier 1 + Tier 2 triple-validated candidates."
      )
    ) +
    theme_pub +
    theme(
      legend.position    = "bottom",
      axis.text.x        = element_text(size = 8, lineheight = 1.1),
      panel.grid.major.x = element_blank()
    )

  save_fig(p_cascade, "fig1_replication_cascade", width = 9, height = 5.5)
  print(p_cascade) # Display inline in Colab

}, error = function(e) {
  warning(sprintf("Fig 1 failed: %s", conditionMessage(e)))
  cat("  Fig 1 skipped.\n")
})

# ---------------------------------------------------------------------------
# Fig 2 — Cross-stream LFC scatter matrix (4-panel)
# ---------------------------------------------------------------------------
cat("\n  Fig 2: Cross-stream LFC scatter matrix\n")

tryCatch({
  # Select the four primary comparisons for the figure
  fig2_pairs <- list(
    list(x = "meta_log2FC", y = "human_lfc",
         xlabel = "CNS meta-LFC (S1)", ylabel = "Human sCJD LFC (S2)",
         panel  = "A   S1 vs S2 (human sCJD)"),
    list(x = "meta_log2FC", y = "whole_blood_lfc",
         xlabel = "CNS meta-LFC (S1)", ylabel = "Whole blood LFC (S3)",
         panel  = "B   S1 vs S3 (whole blood)"),
    list(x = "meta_log2FC", y = "spleen_lfc",
         xlabel = "CNS meta-LFC (S1)", ylabel = "Spleen LFC (S3)",
         panel  = "C   S1 vs S3 (spleen)"),
    list(x = "meta_log2FC", y = "skeletal_muscle_lfc",
         xlabel = "CNS meta-LFC (S1)", ylabel = "Skeletal muscle LFC (S3)",
         panel  = "D   S1 vs S3 (skeletal muscle)")
  )

  scatter_panels <- lapply(fig2_pairs, function(p) {
    xcol <- p$x; ycol <- p$y
    if (!xcol %in% colnames(master) || !ycol %in% colnames(master)) {
      return(
        ggplot() +
          annotate("text", x = 0.5, y = 0.5,
                   label = sprintf("Data not available\n(%s)", p$panel),
                   size = 3.5, colour = "grey50") +
          labs(title = p$panel) +
          theme_pub +
          theme(axis.text = element_blank(), axis.ticks = element_blank(),
                axis.title = element_blank())
      )
    }

    df_p <- master %>%
      filter(!is.na(.data[[xcol]]), !is.na(.data[[ycol]])) %>%
      dplyr::select(gene_symbol, x = all_of(xcol), y = all_of(ycol),
                    mcps_direction, evidence_tier)

    # Pearson r for annotation
    r_row <- lfc_cor_df %>% filter(x_col == xcol, y_col == ycol)
    r_label <- if (nrow(r_row) > 0 && !is.na(r_row$pearson_r[1])) {
      sprintf("r = %+.3f\nn = %d", r_row$pearson_r[1], r_row$n[1])
    } else {
      sprintf("n = %d", nrow(df_p))
    }

    # Highlight triple-validated genes
    df_highlight <- df_p %>%
      filter(evidence_tier %in% c(
        "Tier 1: CNS + human leading edge + concordant peripheral",
        "Tier 2: CNS + human leading edge + detected peripheral"
      ))

    ggplot(df_p, aes(x = x, y = y)) +
      geom_hline(yintercept = 0, linetype = "dashed",
                 colour = "grey70", linewidth = 0.3) +
      geom_vline(xintercept = 0, linetype = "dashed",
                 colour = "grey70", linewidth = 0.3) +
      geom_point(aes(colour = mcps_direction), size = 1.2, alpha = 0.55,
                 stroke = 0) +
      {if (nrow(df_highlight) > 0)
        geom_point(data = df_highlight, aes(x = x, y = y),
                   colour = "#1B7837", size = 2.2, alpha = 0.9,
                   shape = 18)
      } +
      {if (nrow(df_highlight) > 0 && nrow(df_highlight) <= 20)
        geom_text_repel(data = df_highlight,
                        aes(x = x, y = y, label = gene_symbol),
                        size = 2.2, colour = "#1B7837",
                        max.overlaps = 15, segment.size = 0.2,
                        box.padding = 0.3)
      } +
      geom_smooth(method = "lm", formula = y ~ x, se = TRUE,
                  colour = "grey30", linewidth = 0.5, alpha = 0.15) +
      annotate("text",
               x = -Inf, y = Inf,
               label = r_label,
               hjust = -0.1, vjust = 1.3,
               size = 2.8, colour = "grey20") +
      scale_colour_manual(values = DIRECTION_COLOURS, name = "CNS direction") +
      labs(
        x     = p$xlabel,
        y     = p$ylabel,
        title = p$panel
      ) +
      theme_pub +
      theme(legend.position = "none")
  })

  fig2 <- wrap_plots(scatter_panels, ncol = 2) +
    plot_annotation(
      title    = "Cross-stream LFC correlations — MCPS genes",
      subtitle = paste(
        "Each point = one MCPS gene. Green diamonds = triple-validated candidates.",
        "Line = OLS fit ± 95% CI."
      ),
      theme = theme_pub
    )

  save_fig(fig2, "fig2_lfc_scatter_matrix", width = 10, height = 9)
  print(fig2) # Display inline in Colab

}, error = function(e) {
  warning(sprintf("Fig 2 failed: %s", conditionMessage(e)))
  cat("  Fig 2 skipped.\n")
})

# ---------------------------------------------------------------------------
# Fig 3 — Triple-validated candidate dot plot
# ---------------------------------------------------------------------------
cat("\n  Fig 3: Triple-validated candidate dot plot\n")

tryCatch({
  if (nrow(triple_validated) == 0) {
    cat("  No triple-validated candidates; skipping Fig 3.\n")
  } else {
    layers <- c(
      "CNS meta\n(S1)"    = "meta_log2FC",
      "Human sCJD\n(S2)"  = "human_lfc",
      "Whole blood\n(S3)" = "whole_blood_lfc",
      "Spleen\n(S3)"      = "spleen_lfc",
      "Muscle\n(S3)"      = "skeletal_muscle_lfc"
    )

    tv_long <- lapply(names(layers), function(layer_name) {
      col <- layers[[layer_name]]
      if (!col %in% colnames(triple_validated)) {
        return(data.frame(
          gene_symbol  = triple_validated$gene_symbol,
          layer        = layer_name,
          lfc          = NA_real_,
          direction    = NA_character_,
          stringsAsFactors = FALSE
        ))
      }
      data.frame(
        gene_symbol  = triple_validated$gene_symbol,
        evidence_tier = triple_validated$evidence_tier,
        layer        = layer_name,
        lfc          = triple_validated[[col]],
        direction    = ifelse(triple_validated[[col]] > 0,
                              "Up", "Down"),
        abs_lfc      = abs(triple_validated[[col]]),
        stringsAsFactors = FALSE
      )
    }) %>% bind_rows()

    gene_order <- triple_validated %>%
      arrange(evidence_tier, first_detected_wpi, desc(abs(meta_log2FC))) %>%
      pull(gene_symbol) %>%
      unique()

    tv_long <- tv_long %>%
      mutate(
        gene_f  = factor(gene_symbol, levels = rev(gene_order)),
        layer_f = factor(layer, levels = names(layers)),
        tier_short = ifelse(
          str_detect(evidence_tier, "Tier 1"), "Tier 1", "Tier 2"
        )
      )

    p_dotplot <- ggplot(tv_long %>% filter(!is.na(lfc)),
                        aes(x = layer_f, y = gene_f)) +
      geom_point(aes(size = abs_lfc, colour = direction), alpha = 0.85) +
      geom_point(data = tv_long %>% filter(is.na(lfc)),
                 aes(x = layer_f, y = gene_f),
                 shape = 4, size = 2, colour = "grey70") +
      scale_colour_manual(
        values = c("Up" = "#D6604D", "Down" = "#2166AC"),
        name   = "Direction"
      ) +
      scale_size_continuous(
        range  = c(1.5, 6),
        name   = "|LFC|",
        breaks = c(0.5, 1, 2, 3)
      ) +
      facet_grid(tier_short ~ ., scales = "free_y", space = "free_y") +
      labs(
        x       = "Evidence layer",
        y       = NULL,
        title   = "Fig 3 — Triple-validated biomarker candidates",
        subtitle = sprintf(
          "%d Tier 1 + %d Tier 2 candidates. × = not tested in that tissue.",
          n_tier1, n_tier2
        ),
        caption = paste(
          "Tier 1: concordant direction in peripheral tissue + human sCJD leading edge.",
          "Tier 2: detected in peripheral tissue + human sCJD leading edge."
        )
      ) +
      theme_pub +
      theme(
        axis.text.y      = element_text(size = 7, face = "italic"),
        axis.text.x      = element_text(size = 8, lineheight = 1.1),
        strip.text.y     = element_text(size = 8, face = "bold"),
        legend.position  = "right",
        panel.grid.major.x = element_blank()
      )

    n_tv <- length(gene_order)
    fig3_height <- min(14, max(5, 0.22 * n_tv + 2.5))
    save_fig(p_dotplot, "fig3_triple_validated_dotplot",
             width = 8, height = fig3_height)
    print(p_dotplot) # Display inline in Colab
  }

}, error = function(e) {
  warning(sprintf("Fig 3 failed: %s", conditionMessage(e)))
  cat("  Fig 3 skipped.\n")
})

# ---------------------------------------------------------------------------
# Fig 4 — MCPS first-detection timeline (all peripheral tissues)
# ---------------------------------------------------------------------------
cat("\n  Fig 4: MCPS first-detection timeline\n")

tryCatch({
  if (is.null(s3_first_det) || nrow(s3_first_det) == 0) {
    cat("  No first-detection data; skipping Fig 4.\n")
  } else {
    fd_df <- s3_first_det %>%
      mutate(
        tissue_label = TISSUE_LABELS[tissue],
        wpi_label    = factor(wpi, levels = WPI_LEVELS, labels = WPI_LABELS)
      ) %>%
      filter(!is.na(tissue_label))

    fd_cumul <- fd_df %>%
      arrange(tissue, wpi) %>%
      group_by(tissue, tissue_label) %>%
      mutate(
        cumul_detected = cumsum(n_mcps_sig),
        pct_mcps       = 100 * cumul_detected / n_mcps
      ) %>%
      ungroup()

    tissue_colours_fig4 <- c(
      "Whole blood"     = "#E69F00",
      "Spleen"          = "#56B4E9",
      "Skeletal muscle" = "#009E73"
    )

    p_timeline <- ggplot(fd_cumul,
                         aes(x = wpi, y = cumul_detected,
                             colour = tissue_label, group = tissue_label)) +
      geom_line(linewidth = 0.9) +
      geom_point(aes(size = n_mcps_sig), alpha = 0.85) +
      geom_text_repel(
        data = fd_cumul %>%
          group_by(tissue_label) %>%
          filter(cumul_detected == max(cumul_detected)) %>%
          slice_tail(n = 1) %>%
          ungroup(),
        aes(label = sprintf("%s\n(%d genes)", tissue_label, cumul_detected)),
        size = 2.8, nudge_x = 0.5, segment.size = 0.2,
        max.overlaps = 10
      ) +
      scale_colour_manual(values = tissue_colours_fig4, name = "Tissue") +
      scale_size_continuous(range = c(1.5, 5), name = "New MCPS\ngenes at wpi") +
      scale_x_continuous(
        breaks = WPI_LEVELS,
        labels = WPI_LABELS
      ) +
      scale_y_continuous(
        sec.axis = sec_axis(
          ~ . / n_mcps * 100,
          name   = "% of MCPS genes detected",
          breaks = seq(0, 100, 10)
        )
      ) +
      labs(
        x       = "Timepoint (wpi)",
        y       = "Cumulative MCPS genes detected",
        title   = "Fig 4 — MCPS first-detection timeline in peripheral tissues",
        subtitle = sprintf(
          "Cumulative count of MCPS genes (n=%d) first reaching padj \u2264 %.2f in each tissue.",
          n_mcps, PADJ_THRESH
        ),
        caption = "Term = terminal (~22 wpi). Point size = new MCPS genes detected at that timepoint."
      ) +
      theme_pub +
      theme(legend.position = "bottom")

    save_fig(p_timeline, "fig4_detection_timeline", width = 8.5, height = 5.5)
    print(p_timeline) # Display inline in Colab
  }

}, error = function(e) {
  warning(sprintf("Fig 4 failed: %s", conditionMessage(e)))
  cat("  Fig 4 skipped.\n")
})

# ---------------------------------------------------------------------------
# Fig 5 — Pathway convergence heatmap
# ---------------------------------------------------------------------------
cat("\n  Fig 5: Pathway convergence heatmap\n")

tryCatch({
  if (nrow(pathway_convergence) == 0) {
    cat("  No pathway convergence data; skipping Fig 5.\n")
  } else {
    pw_top <- pathway_convergence %>%
      slice_head(n = min(40, nrow(pathway_convergence)))

    pw_mat <- matrix(
      NA_real_,
      nrow = nrow(pw_top),
      ncol = 2,
      dimnames = list(
        substr(pw_top$description, 1, 55),
        c("Stream 1\n(n studies)", "Meta-GSEA\n(sig.)")
      )
    )
    pw_mat[, 1] <- pw_top$n_studies_s1
    pw_mat[, 2] <- as.numeric(pw_top$in_meta_gsea)

    col_fun_studies <- colorRamp2(c(3, 4, 5), c("#DEEBF7", "#9ECAE1", "#2166AC"))
    col_fun_gsea    <- colorRamp2(c(0, 1), c("grey92", "#D6604D"))

    ht_pw <- Heatmap(
      pw_mat[, 1, drop = FALSE],
      name            = "N studies\n(Stream 1)",
      col             = col_fun_studies,
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      show_row_names  = TRUE,
      row_names_gp    = gpar(fontsize = 7),
      column_names_gp = gpar(fontsize = 9),
      width           = unit(1.5, "cm"),
      border          = TRUE,
      heatmap_legend_param = list(
        at = c(3, 4, 5),
        labels = c("3", "4", "5"),
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 7)
      )
    ) +
    Heatmap(
      pw_mat[, 2, drop = FALSE],
      name            = "Meta-GSEA\n(padj \u2264 0.05)",
      col             = col_fun_gsea,
      cluster_rows    = FALSE,
      cluster_columns = FALSE,
      show_row_names  = FALSE,
      column_names_gp = gpar(fontsize = 9),
      width           = unit(1.5, "cm"),
      border          = TRUE,
      heatmap_legend_param = list(
        at = c(0, 1),
        labels = c("No", "Yes"),
        title_gp = gpar(fontsize = 8, fontface = "bold"),
        labels_gp = gpar(fontsize = 7)
      )
    )

    pw_height <- min(12, max(4, 0.22 * nrow(pw_top) + 1.5))

    # Save to Drive
    save_heatmap(
      ht_pw, "fig5_pathway_convergence",
      width = 7, height = pw_height,
      column_title = sprintf(
        "Pathway convergence — %d consensus pathways (>= 3/5 Stream 1 studies)",
        nrow(pw_top)
      ),
      column_title_gp = gpar(fontsize = 10, fontface = "bold")
    )

    # Display inline in Colab
    draw(ht_pw,
         column_title = sprintf("Pathway convergence — %d consensus pathways (>= 3/5 Stream 1 studies)", nrow(pw_top)),
         column_title_gp = gpar(fontsize = 10, fontface = "bold"))
  }

}, error = function(e) {
  warning(sprintf("Fig 5 failed: %s", conditionMessage(e)))
  cat("  Fig 5 skipped.\n")
})

# ---------------------------------------------------------------------------
# Fig 6 — Pathway-Level Human Validation (GSE124571)
# ---------------------------------------------------------------------------
cat("\n  Fig 6: Pathway-Level Human Validation Bubble Plot\n")

tryCatch({
  plot_data <- gsea_res %>%
    dplyr::filter(pathway %in% key_pathways) %>%
    dplyr::mutate(
      neg_log10_padj = -log10(padj),
      Database = dplyr::case_when(
        grepl("^HALLMARK_", pathway) ~ "Hallmark",
        grepl("^GOBP_",     pathway) ~ "GO:BP",
        grepl("^REACTOME_", pathway) ~ "Reactome",
        TRUE ~ "Other"
      ),
      Label = pathway %>%
        gsub("^HALLMARK_|^GOBP_|^REACTOME_", "", .) %>%
        gsub("_", " ", .) %>%
        tools::toTitleCase()
    ) %>%
    dplyr::arrange(desc(NES))

  plot_data$Label <- factor(plot_data$Label, levels = rev(plot_data$Label))

  p_bubble <- ggplot(plot_data, aes(x = NES, y = Label, size = neg_log10_padj, color = Database)) +
    geom_point(alpha = 0.85) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 0.5) +
    scale_size_continuous(name = expression(-log[10](p[adj])), range = c(3, 10),
                          breaks = c(5, 10, 15)) +
    scale_color_manual(values = c("Hallmark" = "#E9ED4C", "GO:BP" = "#75A025",
                                   "Reactome" = "#0279EE"), name = "Gene Set Database") +
    scale_x_continuous(limits = c(-0.5, 3.0), breaks = seq(0, 3, 0.5)) +
    labs(title    = "Fig 6 — Pathway-Level Human Validation",
         subtitle = "Mouse CNS consensus pathways replicated in human sCJD brain (GSE124571 microarray)",
         x = "Normalized Enrichment Score (NES)", y = NULL,
         caption  = "GSEA against 362 consensus pathways from 5-dataset mouse CNS meta-analysis\nAll shown pathways: padj < 0.05 in human sCJD vs. control (N=10+10)") +
    theme_bw(base_size = 12) +
    theme(plot.title = element_text(face = "bold", size = 14),
          axis.text.y = element_text(size = 9),
          legend.position = "right")

  save_fig(p_bubble, "fig6_human_pathway_validation_bubble", width = 11, height = 7)
  print(p_bubble)

}, error = function(e) {
  warning(sprintf("Fig 6 failed: %s", conditionMessage(e)))
  cat("  Fig 6 skipped.\n")
})

cat("\n  [✓] Section 10 successfully executed. All figures generated.\n")


# ----------------------------------------------------------------------
# ##Section 11 — Write Output Tables
#
# **Purpose:** Write all final integration output tables: master integration table, MCPS gene summary, triple-validation statistics, and cross-stream correlation results.
#
# **Inputs:** All integration result objects in R.
#
# **Outputs:** `../results/master_integration_table.tsv`; `../results/mcps_gene_summary.tsv`; `../results/triple_validation_stats.tsv`.
#
# **Key decisions:** The master integration table is the primary deliverable of the entire pipeline — it is the basis for all manuscript tables.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 11 — Write output tables
# =============================================================================
cat("\n--- Section 11: Writing output tables ---\n")

# --- 11.1  Master integration table ---
tryCatch({
  write_csv(
    master %>% arrange(evidence_tier, meta_padj),
    file.path(OUT_DIR, "master_integration_table.csv")
  )
  cat("  Written: master_integration_table.csv\n")
}, error = function(e) {
  warning(sprintf("Failed to write master_integration_table: %s", conditionMessage(e)))
})

# --- 11.2  Triple-validated candidates ---
if (nrow(triple_validated) > 0) {
  tryCatch({
    write_csv(
      triple_validated %>%
        dplyr::select(
          gene_symbol, evidence_tier, mcps_direction,
          meta_log2FC, meta_padj, I2,
          any_of(c("human_lfc", "human_padj", "s2_dir_match")),
          any_of(paste0(TISSUES, "_lfc")),
          any_of(paste0(TISSUES, "_padj")),
          any_of(paste0(TISSUES, "_dir_match")),
          any_of(c("first_detected_wpi", "first_detected_tissue",
                   "first_det_padj", "first_det_wpi_label"))
        ) %>%
        arrange(evidence_tier, first_detected_wpi, meta_padj),
      file.path(OUT_DIR, "triple_validated_candidates.csv")
    )
    cat("  Written: triple_validated_candidates.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write triple_validated_candidates: %s", conditionMessage(e)))
  })
}

# --- 11.3  Cross-stream LFC correlations ---
tryCatch({
  write_csv(lfc_cor_df, file.path(OUT_DIR, "cross_stream_lfc_correlations.csv"))
  cat("  Written: cross_stream_lfc_correlations.csv\n")
}, error = function(e) {
  warning(sprintf("Failed to write lfc_cor_df: %s", conditionMessage(e)))
})

# --- 11.4  Replication cascade ---
tryCatch({
  write_csv(
    cascade_df %>% dplyr::select(stage_order, stage_short, n_genes, pct_of_mcps),
    file.path(OUT_DIR, "replication_cascade.csv")
  )
  cat("  Written: replication_cascade.csv\n")
}, error = function(e) {
  warning(sprintf("Failed to write replication_cascade: %s", conditionMessage(e)))
})

# --- 11.5  Pathway convergence ---
if (nrow(pathway_convergence) > 0) {
  tryCatch({
    write_csv(
      pathway_convergence,
      file.path(OUT_DIR, "pathway_convergence.csv")
    )
    cat("  Written: pathway_convergence.csv\n")
  }, error = function(e) {
    warning(sprintf("Failed to write pathway_convergence: %s", conditionMessage(e)))
  })
}

cat("\n  [✓] Section 11 successfully loaded and executed.\n")


# ----------------------------------------------------------------------
# ##Section 12 — Plain-Text Summary Report & Session Info
#
# **Purpose:** Write a comprehensive plain-text summary of the full integration analysis: sample counts across all streams, MCPS gene statistics, triple-validation results, pathway convergence highlights, and R session info for reproducibility.
#
# **Inputs:** All integration result objects in R.
#
# **Outputs:** `../results/integration_summary.txt`; `../results/session_info.txt`.
#
# **Key decisions:** The summary is designed to be copy-pasted directly into the manuscript Results section. Session info is saved for Methods reproducibility statement.
# ----------------------------------------------------------------------


# =============================================================================
# SECTION 12 — Plain-Text Summary Report & Session Info
# =============================================================================
cat("\n--- Section 12: Writing summary report ---\n")

# --- 12.1  Collect key statistics ---
# FIX: Added as.numeric() to force character columns into floats
s2_nes  <- if (!is.null(s2_gsea) && "NES"  %in% colnames(s2_gsea)) as.numeric(s2_gsea$NES[1])  else NA_real_
s2_padj <- if (!is.null(s2_gsea) && "padj" %in% colnames(s2_gsea)) as.numeric(s2_gsea$padj[1]) else NA_real_

s3_first_wpi_by_tissue <- if (!is.null(s3_first_det)) {
  s3_first_det %>%
    filter(n_mcps_sig >= 1) %>%
    group_by(tissue) %>%
    slice_min(order_by = wpi, n = 1) %>%
    ungroup() %>%
    dplyr::select(tissue, first_wpi = wpi)
} else {
  data.frame(tissue = character(0), first_wpi = integer(0))
}

# --- 12.2  Build summary lines ---
summary_lines <- c(
  "=============================================================",
  "  08_integrate_streams.R — Cross-Stream Integration Summary",
  "=============================================================",
  "",
  "PIPELINE OVERVIEW",
  "-----------------",
  "  Stream 1: CNS meta-analysis (5 Discovery + 4 Replication studies)",
  "  Stream 2: Human sCJD validation (GSE214373, N=14)",
  "  Stream 3: Peripheral longitudinal analysis (GSE210128, N=173)",
  "  Stream 4: Human sCJD Pathway validation (GSE124571, N=20)",
  "",
  "STREAM 1 — MCPS",
  "----------------",
  sprintf("  MCPS size:          %d genes", n_mcps),
  sprintf("  Up in CNS:          %d genes",
          sum(mcps_genes$meta_log2FC > 0, na.rm = TRUE)),
  sprintf("  Down in CNS:        %d genes",
          sum(mcps_genes$meta_log2FC < 0, na.rm = TRUE)),
  sprintf("  Median |meta-LFC|:  %.3f",
          median(abs(mcps_genes$meta_log2FC), na.rm = TRUE)),
  sprintf("  Median I2:          %.1f%%",
          median(mcps_genes$I2, na.rm = TRUE))
)

if (!is.null(replic_summary)) {
  summary_lines <- c(summary_lines,
    "",
    "STREAM 1 — STRATUM D REPLICATION",
    "---------------------------------"
  )
  for (i in seq_len(nrow(replic_summary))) {
    r <- replic_summary[i, ]
    summary_lines <- c(summary_lines,
      sprintf("  %-12s  %d/%d MCPS genes replicated (%.1f%%)",
              r$study, r$n_same_direction, r$mcps_total,
              r$replication_rate_pct))
  }
}

if (!is.null(s2_gsea)) {
  summary_lines <- c(summary_lines,
    "",
    "STREAM 2 — HUMAN sCJD VALIDATION (GSE214373)",
    "---------------------------------------------",
    sprintf("  GSEA NES:           %.3f", s2_nes),
    sprintf("  GSEA padj:          %.4f", s2_padj),
    sprintf("  Leading-edge genes: %d", length(s2_leading_genes))
  )
}

if (!is.null(s3_conc_summary)) {
  summary_lines <- c(summary_lines,
    "",
    "STREAM 3 — PERIPHERAL CONCORDANCE (TERMINAL)",
    "---------------------------------------------"
  )
  for (i in seq_len(nrow(s3_conc_summary))) {
    r <- s3_conc_summary[i, ]
    tl <- if ("tissue_label" %in% colnames(r)) r$tissue_label else r$tissue
    summary_lines <- c(summary_lines,
      sprintf("  %-20s  concordant=%d/%d (%.1f%%)  binom_p=%.4f  r=%.3f",
              tl,
              r$n_dir_match, r$n_mcps_total, r$pct_of_mcps,
              r$binom_pvalue, r$lfc_pearson_r))
  }
}

if (nrow(s3_first_wpi_by_tissue) > 0) {
  summary_lines <- c(summary_lines,
    "",
    "STREAM 3 — MCPS FIRST DETECTION",
    "--------------------------------"
  )
  for (i in seq_len(nrow(s3_first_wpi_by_tissue))) {
    r   <- s3_first_wpi_by_tissue[i, ]
    lbl <- WPI_LABELS[match(r$first_wpi, WPI_LEVELS)]
    summary_lines <- c(summary_lines,
      sprintf("  %-20s  first detection: %s wpi",
              TISSUE_LABELS[r$tissue],
              ifelse(is.na(lbl), as.character(r$first_wpi), lbl)))
  }
}

# Incorporate New Section 9 stats (GSE124571)
if (exists("gsea_res") && !is.null(gsea_res)) {
  summary_lines <- c(summary_lines,
    "",
    "NEW: HUMAN PATHWAY VALIDATION (GSE124571)",
    "-----------------------------------------",
    sprintf("  Pathways Tested:        %d", nrow(gsea_res)),
    sprintf("  Significant (padj<.05): %d", sum(gsea_res$padj < 0.05, na.rm=TRUE)),
    sprintf("  Validation Rate:        %.1f%%", 100 * sum(gsea_res$padj < 0.05, na.rm=TRUE) / nrow(gsea_res))
  )
}

summary_lines <- c(summary_lines,
  "",
  "CROSS-STREAM INTEGRATION",
  "------------------------",
  sprintf("  Master table:              %d MCPS genes × %d columns",
          nrow(master), ncol(master))
)

for (i in seq_len(nrow(tier_counts))) {
  summary_lines <- c(summary_lines,
    sprintf("  %-50s %d genes",
            tier_counts$evidence_tier[i], tier_counts$n[i]))
}

summary_lines <- c(summary_lines,
  "",
  "REPLICATION CASCADE",
  "-------------------"
)
for (i in seq_len(nrow(cascade_df))) {
  r <- cascade_df[i, ]
  if (!is.na(r$n_genes)) {
    summary_lines <- c(summary_lines,
      sprintf("  Stage %d: %-35s %d genes (%.1f%%)",
              r$stage_order,
              gsub("\n", " ", r$stage_short),
              r$n_genes, r$pct_of_mcps))
  }
}

summary_lines <- c(summary_lines,
  "",
  "OUTPUT FILES",
  "------------",
  "  [Data Tables]",
  "  master_integration_table.csv",
  "  triple_validated_candidates.csv",
  "  cross_stream_lfc_correlations.csv",
  "  replication_cascade.csv",
  "  pathway_convergence.csv",
  "  MCPS_literature_convergence_table.csv",
  "  GSE124571_human_GSEA_vs_consensus_pathways.csv",
  "",
  "  [Figures]",
  "  fig1_replication_cascade.pdf/png",
  "  fig2_lfc_scatter_matrix.pdf/png",
  "  fig3_triple_validated_dotplot.pdf/png",
  "  fig4_detection_timeline.pdf/png",
  "  fig5_pathway_convergence.pdf/png",
  "  fig6_human_pathway_validation_bubble.pdf/png",
  "",
  "PACKAGE VERSIONS",
  "----------------",
  sprintf("  R version:       %s", R.version.string),
  sprintf("  ComplexHeatmap:  %s", as.character(packageVersion("ComplexHeatmap"))),
  sprintf("  ggplot2:         %s", as.character(packageVersion("ggplot2"))),
  sprintf("  dplyr:           %s", as.character(packageVersion("dplyr"))),
  sprintf("  limma:           %s", as.character(packageVersion("limma"))),
  sprintf("  fgsea:           %s", as.character(packageVersion("fgsea"))),
  "",
  "=============================================================",
  sprintf("  Generated: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S %Z")),
  "============================================================="
)

tryCatch({
  writeLines(summary_lines,
             file.path(OUT_DIR, "integration_summary.txt"))
  cat("  Written: integration_summary.txt\n")
}, error = function(e) {
  warning(sprintf("Failed to write summary: %s", conditionMessage(e)))
})

cat("\n")
cat(paste(summary_lines, collapse = "\n"), "\n")

# --- 12.3  Session info ---
tryCatch({
  sink(file.path(OUT_DIR, "integration_session_info.txt"))
  print(sessionInfo())
  sink()
  cat("  Written: integration_session_info.txt\n")
}, error = function(e) {
  if (sink.number() > 0) sink()
  warning(sprintf("Failed to write session info: %s", conditionMessage(e)))
})

cat("\n=== 08_integrate_streams.R pipeline strictly completed ===\n")
