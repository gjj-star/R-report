# ===========================================================================
# Script 02: Differential Expression Analysis (DESeq2) & Enrichment
# 项目：TCGA-BRCA多组学数据挖掘
# 方法：DESeq2差异表达分析 + g:Profiler功能富集分析
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(DESeq2)
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
})

cat("\n========== Script 02: Differential Expression Analysis ==========\n\n")

BASE_DIR   <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
INPUT_DIR  <- file.path(BASE_DIR, "data", "processed")
FIG_DIR    <- file.path(BASE_DIR, "results", "figures")
TBL_DIR    <- file.path(BASE_DIR, "results", "tables")

dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. Load data ----
cat("Step 1: Loading preprocessed data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
counts_tumor <- brca$counts
clinical     <- brca$clinical

# Load normal counts
normal_file <- file.path(INPUT_DIR, "brca_normal_counts.rds")
has_normal <- file.exists(normal_file)
if (has_normal) {
  counts_normal <- readRDS(normal_file)
  cat(sprintf("  Tumor: %d genes x %d samples\n", nrow(counts_tumor), ncol(counts_tumor)))
  cat(sprintf("  Normal: %d genes x %d samples\n", nrow(counts_normal), ncol(counts_normal)))
} else {
  stop("Normal tissue counts not found. Cannot perform Tumor vs Normal analysis.")
}

# ---- 2. Build Tumor vs Normal DESeq2 dataset ----
cat("\nStep 2: Building Tumor vs Normal comparison...\n")

common_genes <- intersect(rownames(counts_tumor), rownames(counts_normal))
counts_combined <- cbind(
  counts_tumor[common_genes, , drop = FALSE],
  counts_normal[common_genes, , drop = FALSE]
)

col_data <- data.frame(
  sample_id = colnames(counts_combined),
  condition = c(rep("Tumor", ncol(counts_tumor)), rep("Normal", ncol(counts_normal))),
  stringsAsFactors = FALSE
)
rownames(col_data) <- col_data$sample_id
col_data$condition <- factor(col_data$condition, levels = c("Normal", "Tumor"))

cat(sprintf("  Combined: %d genes x %d samples (%d Tumor + %d Normal)\n",
            nrow(counts_combined), ncol(counts_combined),
            ncol(counts_tumor), ncol(counts_normal)))

# ---- 3. Run DESeq2 ----
cat("\nStep 3: Running DESeq2 (this may take several minutes)...\n")

dds <- DESeqDataSetFromMatrix(
  countData = round(counts_combined),
  colData   = col_data,
  design    = ~ condition
)

# Independent filtering
keep <- rowSums(counts(dds) >= 10) >= min(10, ncol(dds) / 10)
dds <- dds[keep, ]
cat(sprintf("  After filtering: %d genes\n", sum(keep)))

dds <- DESeq(dds)
res <- results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = 0.05)
res <- res[order(res$padj), ]

n_up   <- sum(res$log2FoldChange > 1 & res$padj < 0.05, na.rm = TRUE)
n_down <- sum(res$log2FoldChange < -1 & res$padj < 0.05, na.rm = TRUE)
cat(sprintf("  DEGs (|log2FC|>1, padj<0.05): %d up, %d down, %d total\n",
            n_up, n_down, n_up + n_down))

# ---- 4. Save DEG results ----
cat("\nStep 4: Saving DEG results...\n")

deg_table <- data.frame(
  gene     = rownames(res),
  baseMean = res$baseMean,
  log2FC   = res$log2FoldChange,
  lfcSE    = res$lfcSE,
  stat     = res$stat,
  pvalue   = res$pvalue,
  padj     = res$padj,
  stringsAsFactors = FALSE
)
deg_table$regulation <- ifelse(deg_table$padj < 0.05 & deg_table$log2FC > 1, "Up",
                       ifelse(deg_table$padj < 0.05 & deg_table$log2FC < -1, "Down", "NS"))

write.csv(deg_table, file.path(TBL_DIR, "brca_degs_deseq2.csv"), row.names = FALSE)

top_degs <- deg_table %>% filter(!is.na(padj), padj < 0.05) %>% arrange(padj) %>% head(100)
write.csv(top_degs, file.path(TBL_DIR, "brca_degs_top100.csv"), row.names = FALSE)

# ---- 5. Volcano Plot ----
cat("\nStep 5: Generating volcano plot...\n")

deg_plot <- deg_table %>%
  mutate(
    log10p = -log10(pvalue),
    sig = case_when(
      padj < 0.05 & log2FC > 1  ~ "Up",
      padj < 0.05 & log2FC < -1 ~ "Down",
      TRUE ~ "NS"
    ),
    log2FC_cap = pmax(pmin(log2FC, 8), -8),
    log10p_cap = pmin(log10p, 50)
  )

# Label top genes
top_labels <- deg_plot %>%
  filter(sig != "NS") %>%
  arrange(padj) %>%
  head(15)

cols <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

p <- ggplot(deg_plot, aes(x = log2FC_cap, y = log10p_cap, color = sig)) +
  geom_point(size = 0.5, alpha = 0.6) +
  scale_color_manual(values = cols, name = "Regulation") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_labels, aes(label = gene),
                  size = 3, max.overlaps = 20, color = "black") +
  labs(title = "BRCA: Tumor vs Normal Differential Expression",
       subtitle = sprintf("%d DEGs (|log2FC|>1, padj<0.05)", n_up + n_down),
       x = "log2 Fold Change", y = "-log10(p-value)") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  annotate("text", x = 5, y = 2, label = paste("Up:", n_up), color = "#E41A1C", size = 5) +
  annotate("text", x = -5, y = 2, label = paste("Down:", n_down), color = "#377EB8", size = 5)

ggsave(file.path(FIG_DIR, "volcano_brca.pdf"), p, width = 8, height = 7)
cat("  Volcano plot saved.\n")

# ---- 6. MA Plot ----
cat("\nStep 6: Generating MA plot...\n")
pdf(file.path(FIG_DIR, "ma_plot_brca.pdf"), width = 8, height = 7)
plotMA(res, ylim = c(-5, 5), alpha = 0.05, main = "BRCA: MA Plot (Tumor vs Normal)")
dev.off()
cat("  MA plot saved.\n")

# ---- 7. Top DEG Heatmap ----
cat("\nStep 7: Generating top DEG heatmap...\n")

top50 <- rownames(res)[order(res$padj)][1:min(50, nrow(res))]
top50 <- top50[!is.na(top50)]

if (length(top50) >= 20) {
  vsd <- vst(dds, blind = FALSE)
  vsd_mat <- assay(vsd)[top50, ]

  # Z-score normalization for heatmap
  vsd_scaled <- t(scale(t(vsd_mat)))

  annot_col <- data.frame(Condition = col_data$condition)
  rownames(annot_col) <- colnames(vsd_scaled)
  ann_colors <- list(Condition = c(Tumor = "#E41A1C", Normal = "#377EB8"))

  pdf(file.path(FIG_DIR, "deg_heatmap_top50.pdf"), width = 12, height = 10)
  pheatmap(vsd_scaled,
           scale = "none",
           annotation_col = annot_col,
           annotation_colors = ann_colors,
           show_rownames = TRUE,
           show_colnames = FALSE,
           cluster_rows = TRUE,
           cluster_cols = TRUE,
           fontsize_row = 6,
           main = "BRCA Top 50 DEGs: Tumor vs Normal")
  dev.off()
  cat("  Heatmap saved.\n")
}

# ---- 8. Functional Enrichment Analysis ----
cat("\nStep 8: Functional enrichment analysis (g:Profiler)...\n")

if (requireNamespace("gprofiler2", quietly = TRUE)) {
  library(gprofiler2)

  # All DEGs
  all_degs <- deg_table %>%
    filter(regulation != "NS") %>%
    pull(gene)
  up_degs <- deg_table %>% filter(regulation == "Up") %>% pull(gene)
  down_degs <- deg_table %>% filter(regulation == "Down") %>% pull(gene)

  cat(sprintf("  Running enrichment for: All=%d, Up=%d, Down=%d\n",
              length(all_degs), length(up_degs), length(down_degs)))

  # Enrichment for each group
  run_enrichment <- function(gene_list, label) {
    if (length(gene_list) < 5) return(NULL)
    tryCatch({
      gost_result <- gost(
        query = gene_list,
        organism = "hsapiens",
        significant = TRUE,
        evcodes = TRUE,
        correction_method = "fdr"
      )
      if (!is.null(gost_result$result)) {
        df <- gost_result$result %>% arrange(p_value)
        saveRDS(gost_result, file.path(TBL_DIR, paste0("gprofiler_", label, ".rds")))
        write.csv(df, file.path(TBL_DIR, paste0("enrichment_", label, ".csv")), row.names = FALSE)
        cat(sprintf("    %s: %d enriched terms\n", label, nrow(df)))
        return(df)
      }
    }, error = function(e) {
      cat(sprintf("    %s enrichment failed: %s\n", label, e$message))
    })
    return(NULL)
  }

  enr_all  <- run_enrichment(all_degs, "all_degs")
  enr_up   <- run_enrichment(up_degs, "upregulated")
  enr_down <- run_enrichment(down_degs, "downregulated")

  # Enrichment visualization
  plot_enrichment <- function(enr_df, title, fname) {
    if (is.null(enr_df) || nrow(enr_df) == 0) return()

    top_terms <- enr_df %>%
      head(20) %>%
      mutate(
        term_name_short = ifelse(nchar(term_name) > 50,
                                 paste0(substr(term_name, 1, 47), "..."),
                                 term_name),
        neg_log10p = -log10(p_value)
      )

    p <- ggplot(top_terms, aes(x = reorder(term_name_short, neg_log10p), y = neg_log10p,
                                fill = source)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      labs(title = title, x = "", y = "-log10(p-value)") +
      theme_bw(base_size = 12) +
      theme(legend.position = "bottom")

    ggsave(file.path(FIG_DIR, fname), p, width = 10, height = 8)
  }

  plot_enrichment(enr_up, "Upregulated Genes: Top Enriched Terms", "enrichment_upregulated.pdf")
  plot_enrichment(enr_down, "Downregulated Genes: Top Enriched Terms", "enrichment_downregulated.pdf")
  cat("  Enrichment plots saved.\n")
} else {
  cat("  gprofiler2 not installed. Skipping enrichment analysis.\n")
  cat("  Install with: BiocManager::install('gprofiler2')\n")
}

# ---- 9. Summary ----
cat("\n========== Differential Expression Analysis Complete ==========\n")
cat(sprintf("  Total DEGs: %d (Up: %d, Down: %d)\n", n_up + n_down, n_up, n_down))
cat(sprintf("  Genes tested: %d\n", nrow(res)))
cat("=============================================================\n")
