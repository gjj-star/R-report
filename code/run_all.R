# ===========================================================================
# Master Script: Run All Analysis Steps
# 项目：TCGA-BRCA多组学数据挖掘
# 说明：按顺序执行所有分析脚本
# ===========================================================================

# User library path (winget-installed R requires this)
.user_lib <- "C:/Users/a3036/R/win-library/4.6"
if (dir.exists(.user_lib)) .libPaths(c(.user_lib, .libPaths()))

cat("\n")
cat("================================================================\n")
cat("  TCGA-BRCA Multi-Omics Data Mining Pipeline\n")
cat("  All-in-One Execution Script\n")
cat("================================================================\n")
cat(sprintf("  Start time: %s\n", Sys.time()))
cat(sprintf("  R version: %s\n", R.version.string))
cat(sprintf("  Library paths: %s\n", paste(.libPaths(), collapse = "; ")))
cat("================================================================\n\n")

# Set working directory to project root
BASE_DIR <- "G:/CloudeCoding/Workplace/r-work"
if (!dir.exists(BASE_DIR)) {
  BASE_DIR <- tryCatch(here::here(), error = function(e) getwd())
}

setwd(BASE_DIR)
cat(sprintf("  Working directory: %s\n\n", BASE_DIR))

# Script execution log
log_file <- file.path("results", "run_log.txt")
dir.create("results", showWarnings = FALSE, recursive = TRUE)

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

results <- data.frame(
  Script = character(), Status = character(), Time_min = numeric(),
  stringsAsFactors = FALSE
)

for (script in scripts) {
  script_name <- basename(script)
  cat(sprintf("\n{'='*60}\n"))
  cat(sprintf("Running: %s\n", script_name))
  cat(sprintf("{'='*60}\n"))

  t0 <- Sys.time()

  success <- tryCatch({
    source(script, local = FALSE)
    TRUE
  }, error = function(e) {
    cat(sprintf("\n  ERROR in %s: %s\n", script_name, e$message))
    FALSE
  })

  elapsed <- as.numeric(difftime(Sys.time(), t0, units = "mins"))
  status <- ifelse(success, "OK", "FAILED")

  results <- rbind(results, data.frame(
    Script = script_name, Status = status, Time_min = round(elapsed, 2),
    stringsAsFactors = FALSE
  ))

  cat(sprintf("\n  %s: %s (%.1f min)\n", script_name, status, elapsed))
}

# Final summary
cat("\n\n")
cat("================================================================\n")
cat("  Pipeline Execution Summary\n")
cat("================================================================\n")
print(results)
cat(sprintf("\n  Total time: %.1f minutes\n", sum(results$Time_min)))
cat(sprintf("  End time: %s\n", Sys.time()))
cat(sprintf("  Results: %s\n", file.path(BASE_DIR, "results")))
cat("================================================================\n")

write.csv(results, log_file, row.names = FALSE)
