# ===========================================================================
# Script 08: Comprehensive Visualization (Publication-Quality Figures)
# 项目：TCGA-BRCA多组学数据挖掘
# 说明：汇总各分析模块的关键可视化，生成论文级图表
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
  library(patchwork)
})

cat("\n========== Script 08: Comprehensive Visualization ==========\n\n")

BASE_DIR  <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
INPUT_DIR <- file.path(BASE_DIR, "data", "processed")
FIG_DIR   <- file.path(BASE_DIR, "results", "figures")
TBL_DIR   <- file.path(BASE_DIR, "results", "tables")
PUB_DIR   <- file.path(BASE_DIR, "results", "figures_pub")

dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)

# Color palette (NPG style)
npg_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
              "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# ---- 1. Load data ----
cat("Step 1: Loading analysis results...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat  <- brca$tpm
clinical <- brca$clinical

# ---- 2. Fig1: Data overview ----
cat("\nStep 2: Figure 1 - Data Overview (sample distribution)...\n")

sample_info_file <- file.path(INPUT_DIR, "sample_info.csv")
if (file.exists(sample_info_file)) {
  sample_info <- read.csv(sample_info_file, stringsAsFactors = FALSE)

  # Sample type distribution
  type_counts <- sample_info %>% count(sample_type) %>% arrange(desc(n))

  p1a <- ggplot(type_counts, aes(x = reorder(sample_type, -n), y = n, fill = sample_type)) +
    geom_bar(stat = "identity", width = 0.6) +
    geom_text(aes(label = n), vjust = -0.5, size = 5) +
    scale_fill_manual(values = npg_cols) +
    labs(title = "TCGA-BRCA Sample Distribution",
         x = "", y = "Number of Samples") +
    theme_bw(base_size = 14) +
    theme(legend.position = "none")

  # Stage distribution
  if ("stage_simple" %in% colnames(clinical)) {
    stage_counts <- clinical %>% filter(!is.na(stage_simple)) %>% count(stage_simple)

    p1b <- ggplot(stage_counts, aes(x = stage_simple, y = n, fill = stage_simple)) +
      geom_bar(stat = "identity", width = 0.6) +
      geom_text(aes(label = n), vjust = -0.5, size = 5) +
      scale_fill_manual(values = npg_cols[1:nrow(stage_counts)]) +
      labs(title = "Pathologic Stage Distribution", x = "", y = "Count") +
      theme_bw(base_size = 14) +
      theme(legend.position = "none")

    ggsave(file.path(PUB_DIR, "fig1_sample_overview.pdf"),
           p1a + p1b, width = 12, height = 5)
  } else {
    ggsave(file.path(PUB_DIR, "fig1_sample_overview.pdf"), p1a, width = 6, height = 5)
  }
}

# ---- 3. Fig2: Volcano plot (enhanced) ----
cat("\nStep 3: Figure 2 - Enhanced Volcano Plot...\n")

deg_file <- file.path(TBL_DIR, "brca_degs_deseq2.csv")
if (file.exists(deg_file)) {
  deg <- read.csv(deg_file, stringsAsFactors = FALSE)

  n_up   <- sum(deg$regulation == "Up")
  n_down <- sum(deg$regulation == "Down")

  top_up <- deg %>% filter(regulation == "Up") %>% arrange(padj) %>% head(10)
  top_dn <- deg %>% filter(regulation == "Down") %>% arrange(padj) %>% head(10)
  top_labels <- bind_rows(top_up, top_dn)

  deg_plot <- deg %>%
    mutate(log2FC_cap = pmax(pmin(log2FC, 8), -8),
           log10p = pmin(-log10(pvalue), 50))

  cols <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey75")

  p2 <- ggplot(deg_plot, aes(x = log2FC_cap, y = log10p, color = regulation)) +
    geom_point(size = 0.5, alpha = 0.5) +
    scale_color_manual(values = cols, name = "") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = top_labels, aes(label = gene),
                    size = 3, max.overlaps = 20, color = "black",
                    segment.color = "grey50") +
    labs(title = "Differential Expression: Tumor vs Normal",
         subtitle = sprintf("Up: %d | Down: %d | Total DEGs: %d", n_up, n_down, n_up + n_down),
         x = expression(log[2]~"Fold Change"),
         y = expression(-log[10]~"p-value")) +
    theme_bw(base_size = 14) +
    theme(legend.position = "top")

  ggsave(file.path(PUB_DIR, "fig2_volcano_enhanced.pdf"), p2, width = 9, height = 7)
}

# ---- 4. Fig3: KM Curves by Stage ----
cat("\nStep 4: Figure 3 - KM Curves...\n")
# (Already generated in Script 06; copy if needed)
km_file <- file.path(FIG_DIR, "km_curves_stage.pdf")
if (file.exists(km_file)) file.copy(km_file, file.path(PUB_DIR, "fig3_km_stage.pdf"), overwrite = TRUE)

# ---- 5. Fig4: PCA Plot ----
cat("\nStep 5: Figure 4 - PCA Plot...\n")
pca_file <- file.path(FIG_DIR, "pca_subtype.pdf")
if (file.exists(pca_file)) file.copy(pca_file, file.path(PUB_DIR, "fig4_pca_subtype.pdf"), overwrite = TRUE)

# ---- 6. Fig5: DEG Heatmap ----
cat("\nStep 6: Figure 5 - DEG Heatmap...\n")
heatmap_file <- file.path(FIG_DIR, "deg_heatmap_top50.pdf")
if (file.exists(heatmap_file)) file.copy(heatmap_file, file.path(PUB_DIR, "fig5_deg_heatmap.pdf"), overwrite = TRUE)

# ---- 7. Fig6: Classification Model Comparison ----
cat("\nStep 7: Figure 6 - Classification Comparison...\n")
class_file <- file.path(FIG_DIR, "classification_comparison.pdf")
if (file.exists(class_file)) file.copy(class_file, file.path(PUB_DIR, "fig6_classification.pdf"), overwrite = TRUE)

# ---- 8. Fig7: WGCNA Module-Trait ----
cat("\nStep 8: Figure 7 - WGCNA Module-Trait...\n")
wgcna_file <- file.path(FIG_DIR, "wgcna_module_trait.pdf")
if (file.exists(wgcna_file)) file.copy(wgcna_file, file.path(PUB_DIR, "fig7_wgcna_module_trait.pdf"), overwrite = TRUE)

# ---- 9. Fig8: WGCNA Soft Power ----
cat("\nStep 9: Figure 8 - WGCNA Soft Power...\n")
sft_file <- file.path(FIG_DIR, "wgcna_soft_power.pdf")
if (file.exists(sft_file)) file.copy(sft_file, file.path(PUB_DIR, "fig8_wgcna_soft_power.pdf"), overwrite = TRUE)

# ---- 10. Fig9: Cox Forest Plot ----
cat("\nStep 10: Figure 9 - Cox Forest Plot...\n")
cox_file <- file.path(FIG_DIR, "cox_forest_plot.pdf")
if (file.exists(cox_file)) file.copy(cox_file, file.path(PUB_DIR, "fig9_cox_forest.pdf"), overwrite = TRUE)

# ---- 11. Fig10: Comprehensive Summary Figure ----
cat("\nStep 11: Figure 10 - Analysis Summary...\n")

summary_info <- list()

# DEG stats
if (exists("deg")) {
  summary_info$degs <- sprintf("DEGs: %d (Up: %d, Down: %d)",
                               sum(deg$regulation != "NS"),
                               sum(deg$regulation == "Up"),
                               sum(deg$regulation == "Down"))
}

# Classification
class_file <- file.path(TBL_DIR, "classification_metrics.csv")
if (file.exists(class_file)) {
  metrics <- read.csv(class_file, stringsAsFactors = FALSE)
  best <- metrics[which.max(metrics$Accuracy), ]
  summary_info$classification <- sprintf("Best classifier: %s (%.1f%% accuracy)",
                                          best$Model, best$Accuracy * 100)
}

# WGCNA
wgcna_mod_file <- file.path(TBL_DIR, "wgcna_modules.csv")
if (file.exists(wgcna_mod_file)) {
  mods <- read.csv(wgcna_mod_file, stringsAsFactors = FALSE)
  summary_info$wgcna <- sprintf("WGCNA modules: %d", length(unique(mods$ModuleColor)))
}

# Survival
cox_res_file <- file.path(TBL_DIR, "cox_regression.csv")
if (file.exists(cox_res_file)) {
  cox <- read.csv(cox_res_file, stringsAsFactors = FALSE)
  summary_info$survival <- sprintf("Cox model: %d variables", nrow(cox))
}

# Write summary
summary_text <- paste(unlist(summary_info), collapse = "\n")
writeLines(summary_text, file.path(TBL_DIR, "analysis_summary.txt"))

# ---- 12. List all generated figures ----
cat("\n========== Generated Publication Figures ==========\n")
pub_figs <- list.files(PUB_DIR, pattern = "\\.(pdf|png)$", full.names = FALSE)
for (f in pub_figs) {
  cat(sprintf("  %s\n", f))
}
cat(sprintf("\n  Total: %d publication figures\n", length(pub_figs)))
cat(sprintf("  Output directory: %s\n", PUB_DIR))

cat("\n========== Visualization Complete ==========\n")
