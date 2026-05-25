# ===========================================================================
# Script 05: WGCNA Co-expression Network Analysis
# 项目：TCGA-BRCA多组学数据挖掘
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

# Disable all parallel backends to avoid .doSnowGlobals error
suppressWarnings({
  if (requireNamespace("BiocParallel", quietly = TRUE)) {
    BiocParallel::register(BiocParallel::SerialParam())
  }
  if (requireNamespace("doParallel", quietly = TRUE)) {
    foreach::registerDoSEQ()
  }
})

suppressPackageStartupMessages({
  library(WGCNA)
  library(tidyverse)
  library(pheatmap)
})

# Configure WGCNA threading (minimum 2 required by blockwiseModules)
allowWGCNAThreads(2)

cat("\n========== Script 05: WGCNA Co-expression Network ==========\n\n")

BASE_DIR  <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
INPUT_DIR <- file.path(BASE_DIR, "data", "processed")
FIG_DIR   <- file.path(BASE_DIR, "results", "figures")
TBL_DIR   <- file.path(BASE_DIR, "results", "tables")

# ---- 1. Load and prepare data ----
cat("Step 1: Loading data...\n")

brca     <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat  <- brca$tpm
clinical <- brca$clinical

log_tpm  <- log2(tpm_mat + 1)
gene_var <- apply(log_tpm, 1, var, na.rm = TRUE)
n_genes  <- min(5000, nrow(log_tpm))
top_genes <- names(sort(gene_var, decreasing = TRUE))[1:n_genes]

datExpr <- t(log_tpm[top_genes, ])
cat(sprintf("  WGCNA input: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

# ---- 2. Sample quality check ----
cat("\nStep 2: Sample quality check...\n")

gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat(sprintf("  After filtering: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))
} else {
  cat("  All samples and genes pass QC.\n")
}

# Sample dendrogram
sampleTree <- hclust(dist(datExpr), method = "average")
pdf(file.path(FIG_DIR, "wgcna_sample_dendrogram.pdf"), width = 12, height = 6)
par(cex = 0.5)
plot(sampleTree, main = "BRCA Sample Clustering for WGCNA", sub = "", xlab = "")
dev.off()

# ---- 3. Soft-threshold power selection ----
cat("\nStep 3: Selecting soft-threshold power...\n")

powers <- c(1:20)
sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5, networkType = "signed")

pdf(file.path(FIG_DIR, "wgcna_soft_power.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
cex1 <- 0.9

plot(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     xlab = "Soft Threshold (power)", ylab = "Scale Free Topology Model Fit, signed R^2",
     type = "n", main = "Scale Independence")
text(sft$fitIndices[, 1], -sign(sft$fitIndices[, 3]) * sft$fitIndices[, 2],
     labels = powers, cex = cex1, col = "red")
abline(h = 0.8, col = "blue", lty = 2)

plot(sft$fitIndices[, 1], sft$fitIndices[, 5],
     xlab = "Soft Threshold (power)", ylab = "Mean Connectivity",
     type = "n", main = "Mean Connectivity")
text(sft$fitIndices[, 1], sft$fitIndices[, 5], labels = powers, cex = cex1, col = "red")
par(mfrow = c(1, 1))
dev.off()

soft_power <- sft$powerEstimate
if (is.na(soft_power) || soft_power < 3) soft_power <- 6
cat(sprintf("  Soft power: %d\n", soft_power))

# ---- 4. Build WGCNA network ----
cat("\nStep 4: Building WGCNA network (this may take several minutes)...\n")

net <- blockwiseModules(
  datExpr,
  power = soft_power,
  TOMType = "signed",
  minModuleSize = 30,
  reassignThreshold = 0,
  mergeCutHeight = 0.25,
  numericLabels = TRUE,
  pamRespectsDendro = FALSE,
  saveTOMs = FALSE,
  verbose = 2,
  maxBlockSize = n_genes + 100,
  nThreads = 2
)

module_colors <- labels2colors(net$colors)
cat(sprintf("  Modules identified: %d\n", length(unique(module_colors))))

# Module assignments
module_df <- data.frame(
  Gene = colnames(datExpr),
  Module = net$colors,
  ModuleColor = module_colors,
  stringsAsFactors = FALSE
)
write.csv(module_df, file.path(TBL_DIR, "wgcna_modules.csv"), row.names = FALSE)

mod_sizes <- table(module_colors)
cat("  Module sizes:\n")
for (mc in names(sort(mod_sizes, decreasing = TRUE))) {
  cat(sprintf("    %s: %d genes\n", mc, mod_sizes[mc]))
}

# Module dendrogram
pdf(file.path(FIG_DIR, "wgcna_modules_dendrogram.pdf"), width = 12, height = 6)
plotDendroAndColors(net$dendrograms[[1]],
                    module_colors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE,
                    hang = 0.03,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "WGCNA Module Dendrogram")
dev.off()

# ---- 5. Module-trait correlation ----
cat("\nStep 5: Module-trait correlation analysis...\n")

expr_patients <- substr(rownames(datExpr), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

# Build trait matrix
traits <- data.frame(
  Age = as.numeric(clinical$age_at_diagnosis[match_idx]),
  stringsAsFactors = FALSE
)

if ("positive_lymph_nodes" %in% colnames(clinical))
  traits$LymphNodes <- as.numeric(clinical$positive_lymph_nodes[match_idx])
if ("os_time" %in% colnames(clinical))
  traits$OS_Time <- as.numeric(clinical$os_time[match_idx])
if ("vital_status" %in% colnames(clinical))
  traits$OS_Status <- as.numeric(clinical$vital_status[match_idx])

# Encode subtypes
if ("molecular_subtype" %in% colnames(clinical)) {
  for (st in unique(na.omit(clinical$molecular_subtype[match_idx]))) {
    col_name <- paste0("Subtype_", gsub("[ -]", "_", st))
    traits[[col_name]] <- as.numeric(clinical$molecular_subtype[match_idx] == st)
  }
}

# Encode stages
if ("stage_simple" %in% colnames(clinical)) {
  for (stg in unique(na.omit(clinical$stage_simple[match_idx]))) {
    col_name <- paste0("Stage_", gsub(" ", "_", stg))
    traits[[col_name]] <- as.numeric(clinical$stage_simple[match_idx] == stg)
  }
}

# Impute NAs
for (i in seq_along(traits)) {
  nas <- is.na(traits[[i]])
  if (any(nas)) traits[[i]][nas] <- median(traits[[i]], na.rm = TRUE)
}

# Module eigengenes
MEs <- net$MEs
colnames(MEs) <- gsub("^ME", "M", colnames(MEs))

# Correlation
module_trait_cor  <- cor(MEs, traits, use = "p")
module_trait_pval <- corPvalueStudent(module_trait_cor, nrow(datExpr))

# Heatmap
text_matrix <- matrix(sprintf("%.3f", module_trait_cor),
                      nrow = nrow(module_trait_cor),
                      dimnames = dimnames(module_trait_cor))

pdf(file.path(FIG_DIR, "wgcna_module_trait.pdf"), width = 14, height = 10)
pheatmap(module_trait_cor,
         main = "WGCNA Module-Trait Correlations",
         display_numbers = text_matrix, fontsize_number = 6,
         fontsize = 8,
         color = colorRampPalette(c("#377EB8", "white", "#E41A1C"))(100),
         breaks = seq(-1, 1, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE)
dev.off()
cat("  Module-trait heatmap saved.\n")

# ---- 6. Hub gene identification ----
cat("\nStep 6: Identifying hub genes in key modules...\n")

# Find modules most correlated with survival or stage
if ("OS_Status" %in% colnames(module_trait_cor)) {
  cor_target <- abs(module_trait_cor[, "OS_Status"])
} else {
  cor_target <- apply(abs(module_trait_cor), 1, max)
}

sig_modules <- names(sort(cor_target, decreasing = TRUE))[1:min(3, length(cor_target))]
cat(sprintf("  Key modules: %s\n", paste(sig_modules, collapse = ", ")))

hub_genes_list <- list()
color_to_num <- setNames(unique(net$colors), labels2colors(unique(net$colors)))

for (mod_name in sig_modules) {
  mod_num <- as.numeric(gsub("M", "", mod_name))
  mod_genes <- colnames(datExpr)[net$colors == mod_num]
  me_col <- mod_name

  if (length(mod_genes) >= 20 && me_col %in% colnames(MEs)) {
    kME <- cor(datExpr[, mod_genes, drop = FALSE], MEs[, me_col, drop = FALSE], use = "p")
    colnames(kME) <- "MM"

    hub_df <- data.frame(Gene = rownames(kME), MM = kME[, 1], stringsAsFactors = FALSE) %>%
      arrange(desc(abs(MM))) %>% head(20)
    hub_genes_list[[mod_name]] <- hub_df

    cat(sprintf("    %s: %d genes, top hub=%s (MM=%.3f)\n",
                mod_name, length(mod_genes), hub_df$Gene[1], abs(hub_df$MM[1])))
  }
}

hub_all <- bind_rows(hub_genes_list, .id = "Module")
write.csv(hub_all, file.path(TBL_DIR, "wgcna_hub_genes.csv"), row.names = FALSE)

# ---- 7. Module size distribution plot ----
cat("\nStep 7: Module size distribution...\n")

mod_size_df <- data.frame(
  Module = names(mod_sizes),
  Size = as.numeric(mod_sizes),
  stringsAsFactors = FALSE
) %>% filter(Module != "grey") %>% arrange(desc(Size))

p_mod <- ggplot(mod_size_df, aes(x = reorder(Module, -Size), y = Size, fill = Module)) +
  geom_bar(stat = "identity") +
  scale_fill_identity() +
  labs(title = "WGCNA Module Size Distribution",
       x = "Module", y = "Number of Genes") +
  theme_bw(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(FIG_DIR, "wgcna_module_sizes.pdf"), p_mod, width = 8, height = 5)

# ---- 8. Gene network heatmap for top module ----
cat("\nStep 8: Gene network visualization for top module...\n")

top_mod_name <- sig_modules[1]
top_mod_num  <- as.numeric(gsub("M", "", top_mod_name))
top_mod_genes <- colnames(datExpr)[net$colors == top_mod_num]

if (length(top_mod_genes) >= 30) {
  me_col <- top_mod_name
  kME_all <- cor(datExpr[, top_mod_genes], MEs[, me_col, drop = FALSE], use = "p")
  top_hub_idx <- order(abs(kME_all[, 1]), decreasing = TRUE)[1:min(50, length(top_mod_genes))]
  top_hub <- top_mod_genes[top_hub_idx]

  adj <- abs(cor(datExpr[, top_hub], use = "p")) ^ soft_power
  diag(adj) <- 0

  pdf(file.path(FIG_DIR, "wgcna_network.pdf"), width = 10, height = 10)
  pheatmap(adj,
           main = sprintf("WGCNA %s Module: Gene Network", top_mod_name),
           show_rownames = (ncol(adj) <= 40),
           show_colnames = (ncol(adj) <= 40),
           color = colorRampPalette(c("white", "#FFA500", "#E41A1C"))(50))
  dev.off()
  cat("  Network heatmap saved.\n")
}

cat("\n========== WGCNA Analysis Complete ==========\n")
cat(sprintf("  Modules: %d | Soft power: %d\n", length(unique(module_colors)), soft_power))
