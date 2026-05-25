lib <- "C:/Users/a3036/R/win-library/4.6"
.libPaths(c(lib, .libPaths()))
setwd("G:/CloudeCoding/Workplace/r-work")

cat("Running WGCNA + Visualization...\n")

tryCatch({
  cat("\n=== WGCNA ===\n")
  source("code/05_wgcna.R", local = FALSE)
}, error = function(e) cat(sprintf("WGCNA ERROR: %s\n", e$message)))

tryCatch({
  cat("\n=== Visualization ===\n")
  source("code/08_visualization.R", local = FALSE)
}, error = function(e) cat(sprintf("Visualization ERROR: %s\n", e$message)))

cat("\nDone!\n")
