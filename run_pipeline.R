# Quick setup + run script
.libPaths(c("C:/Users/a3036/R/win-library/4.6", .libPaths()))
setwd("G:/CloudeCoding/Workplace/r-work")

cat("Setting up and running pipeline...\n")
cat("Library paths:\n")
cat(paste(.libPaths(), collapse = "\n"), "\n\n")

# Source each script in order
scripts <- c(
  "code/01_data_preprocessing.R",
  "code/02_diff_expression.R",
  "code/03_classification.R",
  "code/04_clustering.R",
  "code/05_wgcna.R",
  "code/06_survival_analysis.R",
  "code/07_mutation_analysis.R",
  "code/08_visualization.R"
)

for (s in scripts) {
  cat(sprintf("\n=== Running %s ===\n", s))
  t0 <- Sys.time()
  tryCatch({
    source(s, local = FALSE)
    elapsed <- round(as.numeric(difftime(Sys.time(), t0, units = "mins")), 1)
    cat(sprintf("  Done in %.1f min\n", elapsed))
  }, error = function(e) {
    cat(sprintf("  ERROR: %s\n", e$message))
  })
}

cat("\n=== Pipeline Complete ===\n")
cat(sprintf("End time: %s\n", Sys.time()))
