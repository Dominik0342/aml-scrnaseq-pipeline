# =============================================================
# install_packages.R
# Run once before using the pipeline:
#   Rscript install_packages.R
# =============================================================
#
# NOTE: on a fresh Linux VM, you may need to install some system
# libraries first. On Ubuntu/Debian, run from the terminal:
#
#   sudo apt-get update
#   sudo apt-get install -y libcurl4-openssl-dev libssl-dev \
#                            libxml2-dev libfontconfig1-dev \
#                            libharfbuzz-dev libfribidi-dev \
#                            libfreetype6-dev libpng-dev \
#                            libtiff5-dev libjpeg-dev
# =============================================================

# CRAN packages
cran_pkgs <- c("Seurat", "dplyr", "ggplot2", "patchwork",
               "viridis", "Matrix")
for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager", repos = "https://cloud.r-project.org")
}
bioc_pkgs <- c("SingleCellExperiment", "scDblFinder")
for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    BiocManager::install(pkg, update = FALSE, ask = FALSE)
  }
}

# BoneMarrowMap and Symphony (from GitHub)
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools", repos = "https://cloud.r-project.org")
}
if (!requireNamespace("symphony", quietly = TRUE)) {
  devtools::install_github("immunogenomics/symphony", upgrade = "never")
}
if (!requireNamespace("BoneMarrowMap", quietly = TRUE)) {
  devtools::install_github("andygxzeng/BoneMarrowMap", upgrade = "never")
}

cat("\n>>> Installation complete.\n")
cat(">>> You can now run: Rscript pipeline.R\n\n")
