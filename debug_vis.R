lib <- "C:/Users/a3036/R/win-library/4.6"
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages(library(tidyverse))

si <- read.csv("G:/CloudeCoding/Workplace/r-work/data/processed/sample_info.csv", stringsAsFactors = FALSE)
cat("sample_info cols:", paste(colnames(si), collapse = ", "), "\n")
cat("sample_type values:", paste(unique(si$sample_type), collapse = ", "), "\n")

clin <- readRDS("G:/CloudeCoding/Workplace/r-work/data/processed/brca_clinical_clean.rds")
cat("clinical class:", class(clin), "\n")
cat("has stage_simple:", "stage_simple" %in% colnames(clin), "\n")
cat("has molecular_subtype:", "molecular_subtype" %in% colnames(clin), "\n")
cat("Subtype values:\n")
if ("molecular_subtype" %in% colnames(clin)) print(table(clin$molecular_subtype, useNA = "always"))

# Test ggplot
tryCatch({
  type_counts <- si %>% count(sample_type) %>% arrange(desc(n))
  cat("\ntype_counts:\n")
  print(type_counts)
  cat("\nCreating plot...\n")
  p <- ggplot(type_counts, aes(x = reorder(sample_type, -n), y = n, fill = sample_type)) +
    geom_bar(stat = "identity")
  cat("ggplot created successfully\n")
}, error = function(e) cat(sprintf("ggplot error: %s\n", e$message)))
