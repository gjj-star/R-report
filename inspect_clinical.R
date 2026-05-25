lib <- "C:/Users/a3036/R/win-library/4.6"
.libPaths(c(lib, .libPaths()))

e <- new.env()
load("G:/CloudeCoding/Workplace/r-work/TCGA六种癌症转录本数据及临床数据/TCGA六种癌症转录本数据及临床数据/brca_clinical.Rdata", envir = e)
clin <- e$brca_clinical

cat("=== All Clinical Columns ===\n")
for (i in seq_along(colnames(clin))) {
  cat(sprintf("  %3d: %s\n", i, colnames(clin)[i]))
}

cat("\n=== Searching for key columns ===\n")
for (p in c("estrogen", "progesterone", "her2", "stage", "er_status", "pr_status", "lymph", "neoplasm")) {
  idx <- grep(p, colnames(clin), ignore.case = TRUE)
  if (length(idx) > 0) {
    cat(sprintf("  '%s': %s\n", p, paste(colnames(clin)[idx], collapse = ", ")))
  }
}

cat("\n=== Sample values for potential ER/PR/HER2 columns ===\n")
for (col in colnames(clin)) {
  if (grepl("estrogen|progesterone|her2|er_status|pr_status|neoplasm", col, ignore.case = TRUE)) {
    vals <- unique(as.character(clin[[col]]))
    vals <- vals[!is.na(vals)]
    cat(sprintf("  %s: %s\n", col, paste(head(vals, 5), collapse = " | ")))
  }
}
