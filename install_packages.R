# Install all required packages
lib <- "C:/Users/a3036/R/win-library/4.6"
if (!dir.exists(lib)) dir.create(lib, recursive = TRUE)
.libPaths(c(lib, .libPaths()))

cat("R version:", R.version.string, "\n")
cat("Library path:", lib, "\n\n")

# CRAN packages one by one
cran_pkgs <- c("here", "Rtsne", "ggrepel", "pheatmap", "randomForest",
               "glmnet", "survival", "survminer", "WGCNA", "caret",
               "patchwork", "tidyverse")

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE, lib.loc = lib)) {
    cat(sprintf("Installing %s...\n", pkg))
    tryCatch(
      install.packages(pkg, lib = lib, repos = "https://cloud.r-project.org", quiet = TRUE),
      error = function(e) cat(sprintf("  FAILED: %s\n", e$message))
    )
  } else {
    cat(sprintf("%s already installed.\n", pkg))
  }
}

cat("\n--- CRAN packages done ---\n\n")

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE, lib.loc = lib)) {
  cat("Installing BiocManager...\n")
  install.packages("BiocManager", lib = lib, repos = "https://cloud.r-project.org")
}

cat("Installing Bioconductor packages...\n")
BiocManager::install(c("DESeq2", "gprofiler2"), lib = lib, ask = FALSE, update = FALSE)

cat("\n--- All installations attempted ---\n")
installed <- .packages(all.available = TRUE, lib.loc = lib)
cat(sprintf("Total installed: %d packages\n", length(installed)))
