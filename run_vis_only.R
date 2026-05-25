lib <- "C:/Users/a3036/R/win-library/4.6"
.libPaths(c(lib, .libPaths()))
setwd("G:/CloudeCoding/Workplace/r-work")

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
  library(ggrepel)
})

TBL_DIR <- "results/tables"
FIG_DIR <- "results/figures"
PUB_DIR <- "results/figures_pub"
dir.create(PUB_DIR, showWarnings = FALSE, recursive = TRUE)

npg_cols <- c("#E64B35", "#4DBBD5", "#00A087", "#3C5488", "#F39B7F",
              "#8491B4", "#91D1C2", "#DC0000", "#7E6148", "#B09C85")

# Fig 1: Sample overview
cat("Fig 1: Sample overview...\n")
si <- read.csv(file.path("data/processed", "sample_info.csv"), stringsAsFactors = FALSE)
clin <- readRDS(file.path("data/processed", "brca_clinical_clean.rds"))

tp <- table(si$sample_type)
tc <- data.frame(sample_type = names(tp), n = as.integer(tp), stringsAsFactors = FALSE)
tc <- tc[order(-tc$n), ]

p1a <- ggplot(tc, aes(x = reorder(sample_type, -n), y = n, fill = sample_type)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  scale_fill_manual(values = npg_cols) +
  labs(title = "TCGA-BRCA Sample Distribution", x = "", y = "Number of Samples") +
  theme_bw(base_size = 14) + theme(legend.position = "none")

stg <- table(clin$stage_simple[!is.na(clin$stage_simple)])
sc <- data.frame(stage = names(stg), n = as.integer(stg), stringsAsFactors = FALSE)
p1b <- ggplot(sc, aes(x = stage, y = n, fill = stage)) +
  geom_bar(stat = "identity", width = 0.6) +
  geom_text(aes(label = n), vjust = -0.5, size = 5) +
  scale_fill_manual(values = npg_cols[1:nrow(sc)]) +
  labs(title = "Pathologic Stage Distribution", x = "", y = "Count") +
  theme_bw(base_size = 14) + theme(legend.position = "none")

library(patchwork)
ggsave(file.path(PUB_DIR, "fig1_sample_overview.pdf"), p1a + p1b, width = 12, height = 5)
cat("  Done.\n")

# Fig 2: Volcano
cat("Fig 2: Volcano plot...\n")
deg <- read.csv(file.path(TBL_DIR, "brca_degs_deseq2.csv"), stringsAsFactors = FALSE)
n_up <- sum(deg$regulation == "Up")
n_down <- sum(deg$regulation == "Down")

top_labels <- deg %>% filter(regulation != "NS") %>% arrange(padj) %>% head(15)

deg_plot <- deg %>%
  mutate(log2FC_cap = pmax(pmin(log2FC, 8), -8),
         log10p = pmin(-log10(pvalue), 50))

cols <- c("Up" = "#E64B35", "Down" = "#4DBBD5", "NS" = "grey75")

p2 <- ggplot(deg_plot, aes(x = log2FC_cap, y = log10p, color = regulation)) +
  geom_point(size = 0.5, alpha = 0.5) +
  scale_color_manual(values = cols, name = "") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "grey40") +
  geom_text_repel(data = top_labels, aes(x = pmax(pmin(log2FC, 8), -8),
                  y = pmin(-log10(pvalue), 50), label = gene),
                  size = 3, max.overlaps = 20, color = "black") +
  labs(title = "Differential Expression: Tumor vs Normal",
       subtitle = sprintf("Up: %d | Down: %d | Total DEGs: %d", n_up, n_down, n_up + n_down),
       x = expression(log[2]~"Fold Change"),
       y = expression(-log[10]~"p-value")) +
  theme_bw(base_size = 14) + theme(legend.position = "top")

ggsave(file.path(PUB_DIR, "fig2_volcano_enhanced.pdf"), p2, width = 9, height = 7)
cat("  Done.\n")

# Copy remaining figures from figures/ to figures_pub/
cat("Copying remaining figures...\n")
copies <- list(
  c("km_curves_stage.pdf", "fig3_km_stage.pdf"),
  c("pca_subtype.pdf", "fig4_pca_subtype.pdf"),
  c("deg_heatmap_top50.pdf", "fig5_deg_heatmap.pdf"),
  c("classification_comparison.pdf", "fig6_classification.pdf"),
  c("wgcna_module_trait.pdf", "fig7_wgcna_module_trait.pdf"),
  c("wgcna_soft_power.pdf", "fig8_wgcna_soft_power.pdf"),
  c("cox_forest_plot.pdf", "fig9_cox_forest.pdf"),
  c("km_curves_signature.pdf", "fig10_km_signature.pdf"),
  c("hclust_heatmap_top500.pdf", "fig11_hclust_heatmap.pdf"),
  c("silhouette_scores.pdf", "fig12_silhouette.pdf"),
  c("tsne_brca.pdf", "fig13_tsne.pdf"),
  c("classification_roc.pdf", "fig14_roc.pdf"),
  c("classification_cm.pdf", "fig15_confusion_matrix.pdf"),
  c("wgcna_modules_dendrogram.pdf", "fig16_wgcna_dendrogram.pdf"),
  c("wgcna_network.pdf", "fig17_wgcna_network.pdf")
)

for (cp in copies) {
  src <- file.path(FIG_DIR, cp[1])
  dst <- file.path(PUB_DIR, cp[2])
  if (file.exists(src)) file.copy(src, dst, overwrite = TRUE)
}

pub_figs <- list.files(PUB_DIR, pattern = "\\.pdf$")
cat(sprintf("\nTotal publication figures: %d\n", length(pub_figs)))
cat("Visualization complete!\n")
