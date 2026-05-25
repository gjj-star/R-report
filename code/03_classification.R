# ===========================================================================
# Script 03: Molecular Subtype Classification
# 项目：TCGA-BRCA多组学数据挖掘
# 方法：Random Forest + LASSO Multinomial + XGBoost
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(tidyverse)
  library(caret)
  library(randomForest)
  library(glmnet)
})

cat("\n========== Script 03: Molecular Subtype Classification ==========\n\n")

BASE_DIR  <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
INPUT_DIR <- file.path(BASE_DIR, "data", "processed")
FIG_DIR   <- file.path(BASE_DIR, "results", "figures")
TBL_DIR   <- file.path(BASE_DIR, "results", "tables")

# ---- 1. Load and prepare data ----
cat("Step 1: Loading and preparing classification data...\n")

brca <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
tpm_mat  <- brca$tpm
clinical <- brca$clinical

cat(sprintf("  TPM: %d genes x %d samples\n", nrow(tpm_mat), ncol(tpm_mat)))

# Map expression samples to clinical subtypes
expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

# Use molecular_subtype if available, otherwise try stage_simple
if ("molecular_subtype" %in% colnames(clinical)) {
  labels_raw <- clinical$molecular_subtype[match_idx]
  label_name <- "Molecular Subtype"
} else if ("stage_simple" %in% colnames(clinical)) {
  labels_raw <- clinical$stage_simple[match_idx]
  label_name <- "Pathologic Stage"
  cat("  Using pathologic stage as classification target.\n")
} else {
  stop("No suitable classification target found in clinical data.")
}

names(labels_raw) <- colnames(tpm_mat)

# Remove NA labels
valid_idx <- !is.na(labels_raw) & labels_raw != ""
cat(sprintf("  Samples with known %s: %d / %d\n", label_name, sum(valid_idx), length(valid_idx)))

tpm_sub   <- tpm_mat[, valid_idx]
labels_sub <- labels_raw[valid_idx]

cat("  Distribution:\n")
for (lv in names(sort(table(labels_sub), decreasing = TRUE))) {
  cat(sprintf("    %s: %d\n", lv, sum(labels_sub == lv)))
}

# Filter small classes
min_samples <- 15
label_counts <- table(labels_sub)
keep_labels  <- names(label_counts[label_counts >= min_samples])
keep_idx     <- labels_sub %in% keep_labels

tpm_sub  <- tpm_sub[, keep_idx]
labels   <- factor(labels_sub[keep_idx])

cat(sprintf("  After filtering (>=%d samples): %d samples, %d classes\n",
            min_samples, ncol(tpm_sub), length(levels(labels))))

# ---- 2. Feature selection ----
cat("\nStep 2: Feature selection (top 500 variable genes)...\n")

log_tpm <- log2(tpm_sub + 1)
gene_var <- apply(log_tpm, 1, var, na.rm = TRUE)
top500 <- names(sort(gene_var, decreasing = TRUE))[1:min(500, length(gene_var))]

X <- t(log_tpm[top500, ])
y <- labels

cat(sprintf("  Feature matrix: %d samples x %d genes\n", nrow(X), ncol(X)))

# ---- 3. Train/Test split ----
cat("\nStep 3: Train/Test split (70/30 stratified)...\n")

set.seed(42)
train_idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[train_idx, ]
X_test  <- X[-train_idx, ]
y_train <- y[train_idx]
y_test  <- y[-train_idx]

cat(sprintf("  Train: %d | Test: %d\n", length(y_train), length(y_test)))

# ---- 4. Random Forest ----
cat("\nStep 4: Training Random Forest (ntree=500)...\n")

set.seed(42)
rf_model <- randomForest(x = X_train, y = y_train, ntree = 500, importance = TRUE)

rf_pred <- predict(rf_model, X_test)
rf_cm   <- confusionMatrix(rf_pred, y_test)
cat(sprintf("  RF Accuracy: %.4f | Kappa: %.4f\n",
            rf_cm$overall["Accuracy"], rf_cm$overall["Kappa"]))

rf_imp <- importance(rf_model)
rf_imp_df <- data.frame(
  Gene = rownames(rf_imp),
  MeanDecreaseGini = rf_imp[, "MeanDecreaseGini"],
  stringsAsFactors = FALSE
) %>% arrange(desc(MeanDecreaseGini))
write.csv(rf_imp_df, file.path(TBL_DIR, "classification_rf_features.csv"), row.names = FALSE)

# ---- 5. LASSO Multinomial ----
cat("\nStep 5: Training LASSO Multinomial classifier...\n")

set.seed(42)
cv_lasso <- cv.glmnet(x = X_train, y = y_train, family = "multinomial", alpha = 1, nfolds = 5)

lasso_pred <- predict(cv_lasso, X_test, s = "lambda.min", type = "class")
lasso_pred <- factor(as.vector(lasso_pred), levels = levels(y_train))
lasso_cm   <- confusionMatrix(lasso_pred, y_test)
cat(sprintf("  LASSO Accuracy: %.4f | Kappa: %.4f\n",
            lasso_cm$overall["Accuracy"], lasso_cm$overall["Kappa"]))

# Extract selected features
lasso_coef <- coef(cv_lasso, s = "lambda.min")
lasso_genes <- unique(unlist(lapply(lasso_coef, function(x) {
  nz <- which(x[-1] != 0)
  if (length(nz) > 0) colnames(X_train)[nz]
})))
cat(sprintf("  LASSO selected %d genes from %d candidates\n", length(lasso_genes), ncol(X_train)))
write.csv(data.frame(Gene = lasso_genes), file.path(TBL_DIR, "classification_lasso_features.csv"), row.names = FALSE)

# ---- 6. XGBoost ----
cat("\nStep 6: Training XGBoost classifier...\n")

has_xgboost <- requireNamespace("xgboost", quietly = TRUE)
xgboost_acc <- NA

if (has_xgboost) {
  library(xgboost)

  y_train_num <- as.numeric(y_train) - 1
  y_test_num  <- as.numeric(y_test) - 1

  dtrain <- xgb.DMatrix(data = X_train, label = y_train_num)
  dtest  <- xgb.DMatrix(data = X_test,  label = y_test_num)

  params <- list(
    objective      = "multi:softmax",
    num_class      = length(levels(y_train)),
    max_depth      = 6,
    eta            = 0.1,
    subsample      = 0.8,
    colsample_bytree = 0.8
  )

  set.seed(42)
  xgb_model <- xgb.train(params = params, data = dtrain, nrounds = 500,
                          watchlist = list(train = dtrain, test = dtest), verbose = 0)

  xgb_pred_raw   <- predict(xgb_model, dtest)
  xgb_pred_class <- factor(levels(y_train)[xgb_pred_raw + 1], levels = levels(y_train))
  xgb_cm         <- confusionMatrix(xgb_pred_class, y_test)
  xgboost_acc    <- xgb_cm$overall["Accuracy"]
  cat(sprintf("  XGBoost Accuracy: %.4f | Kappa: %.4f\n",
              xgb_cm$overall["Accuracy"], xgb_cm$overall["Kappa"]))

  xgb_imp <- xgb.importance(model = xgb_model, feature_names = colnames(X_train))
  write.csv(xgb_imp, file.path(TBL_DIR, "classification_xgb_features.csv"), row.names = FALSE)
} else {
  cat("  xgboost not installed, skipping.\n")
}

# ---- 7. Model comparison ----
cat("\nStep 7: Model comparison...\n")

model_metrics <- data.frame(
  Model    = c("Random Forest", "LASSO"),
  Accuracy = c(rf_cm$overall["Accuracy"], lasso_cm$overall["Accuracy"]),
  Kappa    = c(rf_cm$overall["Kappa"], lasso_cm$overall["Kappa"]),
  F1_Macro = c(mean(rf_cm$byClass[, "F1"], na.rm = TRUE),
               mean(lasso_cm$byClass[, "F1"], na.rm = TRUE)),
  stringsAsFactors = FALSE
)

if (has_xgboost) {
  model_metrics <- rbind(model_metrics, data.frame(
    Model = "XGBoost", Accuracy = xgb_cm$overall["Accuracy"],
    Kappa = xgb_cm$overall["Kappa"],
    F1_Macro = mean(xgb_cm$byClass[, "F1"], na.rm = TRUE),
    stringsAsFactors = FALSE
  ))
}

write.csv(model_metrics, file.path(TBL_DIR, "classification_metrics.csv"), row.names = FALSE)
cat("  Model metrics:\n")
print(model_metrics)

# Model comparison bar chart
metrics_long <- model_metrics %>%
  pivot_longer(-Model, names_to = "Metric", values_to = "Value")

p <- ggplot(metrics_long, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.3f", Value)), position = position_dodge(width = 0.8),
            vjust = -0.5, size = 3.5) +
  scale_fill_manual(values = c("Accuracy" = "#E41A1C", "Kappa" = "#377EB8", "F1_Macro" = "#4DAF4A")) +
  labs(title = paste("Classification Model Comparison:", label_name),
       y = "Score", x = "") +
  theme_bw(base_size = 14) +
  theme(legend.position = "top") +
  ylim(0, 1)

ggsave(file.path(FIG_DIR, "classification_comparison.pdf"), p, width = 9, height = 6)
cat("  Model comparison plot saved.\n")

# ---- 8. ROC Curves (multi-class) ----
cat("\nStep 8: Generating multi-class ROC curves...\n")

if (requireNamespace("pROC", quietly = TRUE)) {
  rf_prob <- predict(rf_model, X_test, type = "prob")
  roc_colors <- c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00")

  pdf(file.path(FIG_DIR, "classification_roc.pdf"), width = 10, height = 8)
  par(mfrow = c(2, ceiling(length(levels(y_train)) / 2)))

  for (i in seq_along(levels(y_train))) {
    class_name <- levels(y_train)[i]
    binary_truth <- (y_test == class_name)
    roc_obj <- pROC::roc(binary_truth, rf_prob[, class_name], quiet = TRUE)
    pROC::plot.roc(roc_obj, main = paste("ROC:", class_name),
                   col = roc_colors[i], lwd = 2, legacy.axes = TRUE)
    legend("bottomright", sprintf("AUC = %.3f", pROC::auc(roc_obj)), bty = "n", cex = 1.2)
  }

  par(mfrow = c(1, 1))
  dev.off()
  cat("  ROC curves saved.\n")
}

# ---- 9. Confusion Matrix Heatmap ----
cat("\nStep 9: Generating confusion matrix heatmap...\n")

pdf(file.path(FIG_DIR, "classification_cm.pdf"), width = 14, height = 6)
par(mfrow = c(1, 2))

# RF confusion matrix
rf_cm_table <- as.data.frame.matrix(rf_cm$table)
rf_cm_prop  <- sweep(rf_cm_table, 1, rowSums(rf_cm_table), "/")
pheatmap(as.matrix(rf_cm_prop),
         main = "Random Forest: Confusion Matrix",
         display_numbers = TRUE, number_format = "%.2f",
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#377EB8"))(50))

# LASSO confusion matrix
lasso_cm_table <- as.data.frame.matrix(lasso_cm$table)
lasso_cm_prop  <- sweep(lasso_cm_table, 1, rowSums(lasso_cm_table), "/")
pheatmap(as.matrix(lasso_cm_prop),
         main = "LASSO: Confusion Matrix",
         display_numbers = TRUE, number_format = "%.2f",
         cluster_rows = FALSE, cluster_cols = FALSE,
         color = colorRampPalette(c("white", "#E41A1C"))(50))

par(mfrow = c(1, 1))
dev.off()
cat("  Confusion matrices saved.\n")

cat("\n========== Classification Analysis Complete ==========\n")
cat(sprintf("  Best model: %s (Accuracy=%.3f)\n",
            model_metrics$Model[which.max(model_metrics$Accuracy)],
            max(model_metrics$Accuracy)))
