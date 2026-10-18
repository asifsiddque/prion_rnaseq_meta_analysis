## =============================================================
##  renv_setup.R
##  Prion Disease RNA-Seq Meta-Analysis Pipeline
##  R Reproducibility Environment Bootstrap
## =============================================================
##
##  Run this script ONCE before running any analysis script (04–08).
##  It installs renv (if needed) and restores the exact package
##  versions recorded in renv.lock.
##
##  Usage
##  -----
##  MobaXterm / any terminal:
##    Rscript scripts/renv_setup.R
##
##  RStudio:
##    source("scripts/renv_setup.R")
##
##  Kaggle Notebook (paste into a code cell):
##    source("/kaggle/working/scripts/renv_setup.R")
##
##  After setup, start every R session with:
##    renv::activate()   # or open the project in RStudio (auto-activates)
##
##  Requirements
##  ------------
##  - R >= 4.3
##  - Bioconductor 3.18 (installed automatically by renv::restore())
##  - Internet access for first-time package installation
##  - ~3–5 GB disk space for the project library
##
##  Notes
##  -----
##  - renv installs packages into a project-local library (renv/library/).
##    Your system R library is NOT modified.
##  - If renv::restore() fails for a specific package, run:
##      renv::install("<package_name>")
##    then re-run renv::restore().
##  - To update the lockfile after adding new packages:
##      renv::snapshot()
## =============================================================

cat("=============================================================\n")
cat("  renv_setup.R — Prion RNA-Seq Pipeline\n")
cat("=============================================================\n\n")

## ── Step 1: Install renv if not present ──────────────────────
if (!requireNamespace("renv", quietly = TRUE)) {
  cat("  [1/3] Installing renv from CRAN...\n")
  install.packages("renv", repos = "https://cloud.r-project.org")
} else {
  cat(sprintf("  [1/3] renv %s already installed.\n",
              as.character(packageVersion("renv"))))
}

## ── Step 2: Initialise renv project (if not already done) ────
##  renv::init() is idempotent — safe to call on an already-
##  initialised project. We use bare = TRUE to avoid copying
##  the current library (we want a clean restore from lockfile).
cat("  [2/3] Initialising renv project...\n")
if (!file.exists("renv.lock")) {
  stop(paste(
    "renv.lock not found in the current working directory.",
    "Make sure you are running this script from the project root:",
    "  setwd('/path/to/prion_rnaseq_meta_analysis')",
    "  source('scripts/renv_setup.R')",
    sep = "\n"
  ))
}

## Activate renv for this session (creates renv/ folder if absent)
renv::activate()

## ── Step 3: Restore exact package versions from renv.lock ────
cat("  [3/3] Restoring packages from renv.lock...\n")
cat("        This may take 20–60 minutes on first run.\n\n")

tryCatch({
  renv::restore(prompt = FALSE)
  cat("\n  renv::restore() completed successfully.\n")
}, error = function(e) {
  cat(sprintf("\n  ERROR during renv::restore(): %s\n", conditionMessage(e)))
  cat("  Try running renv::restore() manually in an interactive R session.\n")
  cat("  If a specific package fails, install it with:\n")
  cat("    renv::install('<package_name>')\n")
  stop(e)
})

## ── Done ─────────────────────────────────────────────────────
cat("\n=============================================================\n")
cat("  Setup complete. Package versions are now locked to renv.lock.\n")
cat("\n  To activate this environment in future R sessions:\n")
cat("    renv::activate()   # or open the .Rproj in RStudio\n")
cat("\n  To run the pipeline:\n")
cat("    Rscript scripts/04_build_coldata.R\n")
cat("    Rscript scripts/05_stream1_cns_meta_analysis.R\n")
cat("    Rscript scripts/06_stream2_human_validation.R\n")
cat("    Rscript scripts/07_stream3_peripheral_analysis.R\n")
cat("    Rscript scripts/08_integrate_streams.R\n")
cat("=============================================================\n")
