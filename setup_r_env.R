# setup_r_env.R
# ─────────────────────────────────────────────────────────────────────────────
# Run this ONCE after creating the conda environment to install all R packages
# and capture an renv.lock for reproducibility.
#
# Usage:
#   conda activate scrnaseq
#   Rscript setup_r_env.R
# ─────────────────────────────────────────────────────────────────────────────

# ── 1. Bootstrap renv ────────────────────────────────────────────────────────
# renv isolates R packages per project (like Python venvs)
if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

# Initialise renv in the project root (creates renv/ and renv.lock)
# bare = TRUE means don't scan project files yet — we'll snapshot manually
renv::init(bare = TRUE)

# ── 2. Install CRAN packages ─────────────────────────────────────────────────
cran_pkgs <- c(
  # Seurat ecosystem
  "Seurat",          # v5 — single-cell analysis
  "SeuratObject",
  # Data manipulation
  "dplyr",
  "data.table",
  "Matrix",
  # Visualisation
  "ggplot2",
  "patchwork",
  # File utilities
  "R.utils",         # gunzip()
  # Jupyter kernel (must be installed into renv library, not just conda)
  "IRkernel",
  # Package management
  "BiocManager",
  "remotes"
)

install.packages(cran_pkgs, repos = "https://cloud.r-project.org")

# ── 3. Install Bioconductor packages ────────────────────────────────────────
# Note: glmGamPoi (used by SCTransform) is excluded — C++14 compile issue on
# macOS ARM. This pipeline uses NormalizeData+ScaleData, not SCTransform.
bioc_pkgs <- c(
  "BiocGenerics"
)
BiocManager::install(bioc_pkgs, ask = FALSE)

# ── 4. Register IRkernel so papermill finds the 'ir' kernel ─────────────────
# Requires jupyter to be on PATH — set it explicitly for robustness
conda_bin <- Sys.getenv("CONDA_PREFIX")
if (nchar(conda_bin) > 0) {
  old_path <- Sys.getenv("PATH")
  Sys.setenv(PATH = paste0(conda_bin, "/bin:", old_path))
}
IRkernel::installspec(user = TRUE)
cat("ir kernel registered.\n")

# ── 5. Snapshot — writes renv.lock with exact versions of everything installed
renv::snapshot()
cat("\nrenv.lock written. Commit this file to reproduce the environment.\n")

# ── 6. Print session info for the log ────────────────────────────────────────
cat("\n── Session Info ──────────────────────────────────────────────────────\n")
sessionInfo()
