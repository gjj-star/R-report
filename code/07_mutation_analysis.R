# ===========================================================================
# Script 07: Somatic Mutation Analysis
# 项目：TCGA-BRCA多组学数据挖掘
# 方法：maftools突变景观分析
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("\n========== Script 07: Mutation Analysis ==========\n\n")

BASE_DIR  <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
FIG_DIR   <- file.path(BASE_DIR, "results", "figures")
TBL_DIR   <- file.path(BASE_DIR, "results", "tables")

# ---- 1. Check for maftools and mutation data ----
cat("Step 1: Checking for mutation data and maftools...\n")

if (!requireNamespace("maftools", quietly = TRUE)) {
  cat("  maftools not installed. Skipping mutation analysis.\n")
  cat("  Install with: BiocManager::install('maftools')\n")
  cat("  Note: TCGA mutation data (MAF files) need to be downloaded separately.\n")
  cat("  The analysis will proceed without mutation data.\n")
} else {
  library(maftools)

  # Check if mutation data exists locally
  maf_file <- file.path(BASE_DIR, "data", "processed", "brca_mutation.maf")
  maf_rds  <- file.path(BASE_DIR, "data", "processed", "brca_maf.rds")

  maf_data <- NULL

  if (file.exists(maf_rds)) {
    cat("  Loading cached MAF data...\n")
    maf_data <- readRDS(maf_rds)
  } else if (file.exists(maf_file)) {
    cat("  Reading MAF file...\n")
    maf_data <- read.maf(maf = maf_file)
  } else {
    # Try to download from TCGA using TCGAbiolinks
    cat("  No local mutation data found.\n")
    cat("  Attempting to download TCGA-BRCA mutation data via TCGAbiolinks...\n")

    if (requireNamespace("TCGAbiolinks", quietly = TRUE)) {
      library(TCGAbiolinks)

      tryCatch({
        query_maf <- GDCquery(
          project = "TCGA-BRCA",
          data.category = "Simple Nucleotide Variation",
          data.type = "Masked Somatic Mutation",
          access = "open",
          workflow.type = "Aliquot Ensemble Somatic Variant Merging and Masking"
        )

        maf_dir <- file.path(BASE_DIR, "data", "raw", "maf")
        dir.create(maf_dir, showWarnings = FALSE, recursive = TRUE)
        GDCdownload(query_maf, directory = maf_dir)

        maf_files <- GDCprepare_clinic(query_maf, clinical.info = "patient")
        maf_data <- read.maf(maf = query_maf)
        saveRDS(maf_data, maf_rds)
        cat("  Mutation data downloaded and cached.\n")
      }, error = function(e) {
        cat(sprintf("  Download failed: %s\n", e$message))
      })
    }
  }

  # ---- 2. Mutation analysis ----
  if (!is.null(maf_data)) {
    cat("\nStep 2: Mutation landscape analysis...\n")

    # Mutation summary
    cat(sprintf("  Samples: %d\n", nrow(maf_data@summary$summary)))
    cat(sprintf("  Genes: %d\n", maf_data@summary$total))
    cat(sprintf("  Variants: %d\n", maf_data@variants.per.sample$Variants))

    # Oncoplot (Top 20 genes)
    cat("\nStep 3: Generating oncoplot...\n")
    pdf(file.path(FIG_DIR, "oncoplot_top20.pdf"), width = 12, height = 8)
    tryCatch({
      oncoplot(maf = maf_data, top = 20, fontSize = 0.7)
    }, error = function(e) {
      cat(sprintf("  Oncoplot failed: %s\n", e$message))
    })
    dev.off()

    # Variant type distribution
    cat("\nStep 4: Variant type distribution...\n")
    vcs <- maf_data@summary$variant.type.summary
    write.csv(vcs, file.path(TBL_DIR, "mutated_genes_summary.csv"), row.names = FALSE)

    pdf(file.path(FIG_DIR, "mutation_types.pdf"), width = 8, height = 5)
    tryCatch({
      plotmafSummary(maf = maf_data, rmOutlier = TRUE, addStat = 'median')
    }, error = function(e) {
      cat(sprintf("  Mutation summary plot failed: %s\n", e$message))
    })
    dev.off()

    # Mutually exclusive / co-occurring genes
    cat("\nStep 5: Pairwise mutation analysis...\n")
    tryCatch({
      co_occur <- somaticInteractions(maf = maf_data, top = 25, pvalue = c(0.05, 0.01))
      if (!is.null(co_occur)) {
        write.csv(as.data.frame(co_occur), file.path(TBL_DIR, "mutation_interactions.csv"),
                  row.names = FALSE)
      }
    }, error = function(e) {
      cat(sprintf("  Interaction analysis failed: %s\n", e$message))
    })

    # Lollipop plot for top gene
    cat("\nStep 6: Lollipop plot for top gene...\n")
    top_gene <- maf_data@summary$gene.summary$Hugo_Symbol[1]
    if (!is.na(top_gene)) {
      pdf(file.path(FIG_DIR, "lollipop_top_gene.pdf"), width = 10, height = 5)
      tryCatch({
        lollipopPlot(maf = maf_data, gene = top_gene, AACol = "HGVSp_Short",
                     showMutationRate = TRUE)
      }, error = function(e) {
        cat(sprintf("  Lollipop plot failed: %s\n", e$message))
      })
      dev.off()
    }

    cat("\n========== Mutation Analysis Complete ==========\n")
  } else {
    cat("\n  No mutation data available. Mutation analysis skipped.\n")
    cat("  To include mutation analysis:\n")
    cat("  1. Download TCGA-BRCA MAF files from GDC portal\n")
    cat("  2. Save as 'brca_mutation.maf' in data/processed/\n")
    cat("  3. Re-run this script\n")
  }
}

cat("\n========== Script 07 Done ==========\n")
