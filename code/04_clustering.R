# ===========================================================================
# Script 04: Clustering Analysis (PCA, t-SNE, K-means, Consensus Clustering)
# 项目：TCGA-BRCA多组学数据挖掘
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(pheatmap)
})

cat("\n========== Script 04: Clustering Analysis ==========\n\n")

BASE_DIR  <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
INPUT_DIR <- file.path(BASE_DIR, "data", "processed")
FIG_DIR   <- file.path(BASE_DIR, "results", "figures")
TBL_DIR   <- file.path(BASE_DIR, "results", "tables")

# ---- 1. Load data ----
cat("Step 1: Loading data...\n")

brca     <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat  <- brca$tpm
clinical <- brca$clinical

log_tpm  <- log2(tpm_mat + 1)
gene_var <- apply(log_tpm, 1, var, na.rm = TRUE)
top2000  <- names(sort(gene_var, decreasing = TRUE))[1:min(2000, nrow(log_tpm))]
log_tpm_top <- log_tpm[top2000, ]

cat(sprintf("  Using top %d variable genes for clustering\n", length(top2000)))

# Map clinical annotations
expr_patients <- substr(colnames(log_tpm_top), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)
annot_sub  <- clinical$molecular_subtype[match_idx]
annot_stage <- clinical$stage_simple[match_idx]

subtype_colors <- c("Luminal A" = "#1B9E77", "Luminal B" = "#D95F02",
                    "HER2-enriched" = "#7570B3", "Triple Negative" = "#E7298A")
stage_colors <- c("Stage I" = "#1B9E77", "Stage II" = "#D95F02",
                   "Stage III" = "#7570B3", "Stage IV" = "#E7298A")

# ---- 2. PCA ----
cat("\nStep 2: PCA analysis...\n")

pca <- prcomp(t(log_tpm_top), center = TRUE, scale. = TRUE)
pca_var <- summary(pca)$importance[2, ] * 100

pca_df <- data.frame(
  PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
  sample_id = rownames(pca$x),
  subtype = annot_sub,
  stage   = annot_stage,
  stringsAsFactors = FALSE
)

cat(sprintf("  PC1: %.1f%% | PC2: %.1f%% | PC3: %.1f%% variance\n",
            pca_var[1], pca_var[2], pca_var[3]))

# PCA by subtype
p_sub <- ggplot(pca_df %>% filter(!is.na(subtype)),
                aes(x = PC1, y = PC2, color = subtype)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = subtype_colors, name = "Molecular Subtype") +
  labs(title = "BRCA PCA: Colored by Molecular Subtype",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_bw(base_size = 14) +
  stat_ellipse(level = 0.95, linewidth = 1)
ggsave(file.path(FIG_DIR, "pca_subtype.pdf"), p_sub, width = 9, height = 7)

# PCA by stage
p_stage <- ggplot(pca_df %>% filter(!is.na(stage)),
                  aes(x = PC1, y = PC2, color = stage)) +
  geom_point(size = 2, alpha = 0.7) +
  scale_color_manual(values = stage_colors, name = "Stage") +
  labs(title = "BRCA PCA: Colored by Pathologic Stage",
       x = sprintf("PC1 (%.1f%%)", pca_var[1]),
       y = sprintf("PC2 (%.1f%%)", pca_var[2])) +
  theme_bw(base_size = 14)
ggsave(file.path(FIG_DIR, "pca_stage.pdf"), p_stage, width = 9, height = 7)
cat("  PCA plots saved.\n")

# ---- 3. t-SNE ----
cat("\nStep 3: t-SNE dimensionality reduction...\n")

if (requireNamespace("Rtsne", quietly = TRUE)) {
  library(Rtsne)
  set.seed(42)
  tsne_res <- Rtsne(t(log_tpm_top), perplexity = min(30, ncol(log_tpm_top) / 4),
                     max_iter = 1000, check_duplicates = FALSE)

  tsne_df <- data.frame(
    tSNE1 = tsne_res$Y[, 1], tSNE2 = tsne_res$Y[, 2],
    sample_id = colnames(log_tpm_top),
    subtype = annot_sub,
    stringsAsFactors = FALSE
  )

  p_tsne <- ggplot(tsne_df %>% filter(!is.na(subtype)),
                   aes(x = tSNE1, y = tSNE2, color = subtype)) +
    geom_point(size = 2, alpha = 0.7) +
    scale_color_manual(values = subtype_colors, name = "Molecular Subtype") +
    labs(title = "BRCA t-SNE: Top 2000 Variable Genes",
         x = "t-SNE 1", y = "t-SNE 2") +
    theme_bw(base_size = 14)
  ggsave(file.path(FIG_DIR, "tsne_brca.pdf"), p_tsne, width = 9, height = 7)
  cat("  t-SNE plot saved.\n")
} else {
  cat("  Rtsne not installed. Skipping t-SNE.\n")
}

# ---- 4. Hierarchical Clustering Heatmap ----
cat("\nStep 4: Hierarchical clustering heatmap...\n")

top500 <- names(sort(gene_var, decreasing = TRUE))[1:min(500, length(gene_var))]
log_tpm_500 <- log_tpm[top500, ]

annot_col <- data.frame(
  Subtype = ifelse(is.na(annot_sub), "Unknown", annot_sub),
  Stage   = ifelse(is.na(annot_stage), "Unknown", annot_stage),
  stringsAsFactors = FALSE
)
rownames(annot_col) <- colnames(log_tpm_500)

ann_colors_hm <- list(
  Subtype = c(subtype_colors, "Unknown" = "grey80"),
  Stage   = c(stage_colors, "Unknown" = "grey80")
)

pdf(file.path(FIG_DIR, "hclust_heatmap_top500.pdf"), width = 14, height = 10)
pheatmap(log_tpm_500,
         scale = "row",
         annotation_col = annot_col,
         annotation_colors = ann_colors_hm,
         show_rownames = FALSE,
         show_colnames = FALSE,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation",
         clustering_method = "ward.D2",
         main = "BRCA Hierarchical Clustering: Top 500 Variable Genes")
dev.off()
cat("  Hierarchical clustering heatmap saved.\n")

# ---- 5. K-means Silhouette Analysis ----
cat("\nStep 5: K-means silhouette analysis...\n")

sil_scores <- numeric(9)
names(sil_scores) <- as.character(2:10)

set.seed(42)
for (k in 2:10) {
  km <- kmeans(t(log_tpm_top), centers = k, nstart = 25)
  ss <- cluster::silhouette(km$cluster, dist(t(log_tpm_top)))
  sil_scores[as.character(k)] <- mean(ss[, 3])
}

sil_df <- data.frame(K = 2:10, Silhouette = sil_scores, stringsAsFactors = FALSE)
write.csv(sil_df, file.path(TBL_DIR, "clustering_silhouette.csv"), row.names = FALSE)

p_sil <- ggplot(sil_df, aes(x = K, y = Silhouette)) +
  geom_line(color = "#377EB8", linewidth = 1) +
  geom_point(size = 3, color = "#E41A1C") +
  labs(title = "BRCA K-means Silhouette Analysis",
       x = "Number of Clusters (K)", y = "Average Silhouette Width") +
  theme_bw(base_size = 14) +
  scale_x_continuous(breaks = 2:10)
ggsave(file.path(FIG_DIR, "silhouette_scores.pdf"), p_sil, width = 7, height = 5)

best_k <- which.max(sil_scores) + 1
cat(sprintf("  Best K by silhouette: %d (score=%.4f)\n", best_k, max(sil_scores)))

# ---- 6. Save K-means clusters ----
set.seed(42)
km_best <- kmeans(t(log_tpm_top), centers = best_k, nstart = 25)

cluster_df <- data.frame(
  sample_id     = colnames(log_tpm_top),
  kmeans_cluster = km_best$cluster,
  subtype        = annot_sub,
  stringsAsFactors = FALSE
)
write.csv(cluster_df, file.path(TBL_DIR, "clustering_kmeans.csv"), row.names = FALSE)

# Cluster vs Subtype contingency
if (!all(is.na(annot_sub))) {
  ct <- table(cluster_df$kmeans_cluster, cluster_df$subtype)
  write.csv(ct, file.path(TBL_DIR, "clustering_vs_subtype.csv"))
  cat("  Cluster vs Subtype table saved.\n")
}

# ---- 7. Consensus Clustering (optional) ----
cat("\nStep 7: Consensus clustering (ConsensusClusterPlus)...\n")

if (requireNamespace("ConsensusClusterPlus", quietly = TRUE)) {
  library(ConsensusClusterPlus)

  cc_dir <- file.path(TBL_DIR, "consensus_clustering")
  dir.create(cc_dir, showWarnings = FALSE, recursive = TRUE)

  cc_results <- ConsensusClusterPlus(
    d = as.matrix(log_tpm_top),
    maxK = 8,
    reps = 1000,
    pItem = 0.8,
    pFeature = 1,
    clusterAlg = "hc",
    distance = "pearson",
    innerLinkage = "ward.D2",
    finalLinkage = "ward.D2",
    seed = 42,
    title = cc_dir,
    plot = "pdf"
  )

  # Save consensus cluster assignments
  consensus_df <- data.frame(
    sample_id = colnames(log_tpm_top),
    K2 = cc_results[[2]]$consensusClass,
    K3 = cc_results[[3]]$consensusClass,
    K4 = cc_results[[4]]$consensusClass,
    stringsAsFactors = FALSE
  )
  write.csv(consensus_df, file.path(cc_dir, "consensus_clusters.csv"), row.names = FALSE)
  cat(sprintf("  Consensus clustering completed (K=2..8).\n"))
} else {
  cat("  ConsensusClusterPlus not installed. Skipping.\n")
}

cat("\n========== Clustering Analysis Complete ==========\n")
