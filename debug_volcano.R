lib <- "C:/Users/a3036/R/win-library/4.6"
.libPaths(c(lib, .libPaths()))
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggrepel)
})

TBL_DIR <- "G:/CloudeCoding/Workplace/r-work/results/tables"
FIG_DIR <- "G:/CloudeCoding/Workplace/r-work/results/figures"

deg_file <- file.path(TBL_DIR, "brca_degs_deseq2.csv")
deg <- read.csv(deg_file, stringsAsFactors = FALSE)

cat("DEG columns:", paste(colnames(deg), collapse = ", "), "\n")
cat("DEG rows:", nrow(deg), "\n")
cat("regulation values:\n")
print(table(deg$regulation))

# Try creating volcano
deg_plot <- deg %>%
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

cat("\nsig values:\n")
print(table(deg_plot$sig))

n_up   <- sum(deg_plot$sig == "Up")
n_down <- sum(deg_plot$sig == "Down")

top_labels <- deg_plot %>%
  filter(sig != "NS") %>%
  arrange(padj) %>%
  head(15)

cat("\nTop labels:\n")
print(head(top_labels[, c("gene", "log2FC", "padj", "sig")]))

cols <- c("Up" = "#E41A1C", "Down" = "#377EB8", "NS" = "grey70")

cat("\nCreating volcano plot...\n")
tryCatch({
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
    theme(legend.position = "top")

  ggsave(file.path(FIG_DIR, "volcano_brca_v2.pdf"), p, width = 8, height = 7)
  cat("Volcano plot saved!\n")
}, error = function(e) {
  cat(sprintf("ERROR: %s\n", e$message))
  cat(sprintf("Class of deg_plot$log2FC_cap: %s\n", class(deg_plot$log2FC_cap)))
  cat(sprintf("Class of deg_plot$log10p_cap: %s\n", class(deg_plot$log10p_cap)))
  cat(sprintf("Class of deg_plot$sig: %s\n", class(deg_plot$sig)))
  cat(sprintf("Any NA in log2FC_cap: %s\n", any(is.na(deg_plot$log2FC_cap))))
  cat(sprintf("Any Inf in log10p_cap: %s\n", any(is.infinite(deg_plot$log10p_cap))))
})
