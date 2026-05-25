# ===========================================================================
# Script 06: Survival Analysis (KM, Cox Regression, Prognostic Signature)
# 项目：TCGA-BRCA多组学数据挖掘
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(tidyverse)
  library(glmnet)
})

cat("\n========== Script 06: Survival Analysis ==========\n\n")

BASE_DIR  <- tryCatch(here::here(), error = function(e) normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")))
INPUT_DIR <- file.path(BASE_DIR, "data", "processed")
FIG_DIR   <- file.path(BASE_DIR, "results", "figures")
TBL_DIR   <- file.path(BASE_DIR, "results", "tables")

# ---- 1. Load data ----
cat("Step 1: Loading data...\n")

brca     <- readRDS(file.path(INPUT_DIR, "brca_tumor_processed.rds"))
clinical <- brca$clinical
tpm_mat  <- brca$tpm

clinical <- clinical %>%
  mutate(
    os_time_years = os_time / 365.25,
    os_status = vital_status
  )

cat(sprintf("  Patients: %d | Events (death): %d\n",
            nrow(clinical), sum(clinical$os_status == 1, na.rm = TRUE)))

# ---- 2. KM by molecular subtype ----
cat("\nStep 2: KM curves by molecular subtype...\n")

if ("molecular_subtype" %in% colnames(clinical)) {
  df_sub <- clinical %>%
    filter(!is.na(molecular_subtype), molecular_subtype != "",
           !is.na(os_time_years), os_time_years > 0)

  fit_sub <- survfit(Surv(os_time_years, os_status) ~ molecular_subtype, data = df_sub)

  subtype_colors <- c("Luminal A" = "#1B9E77", "Luminal B" = "#D95F02",
                      "HER2-enriched" = "#7570B3", "Triple Negative" = "#E7298A")
  # Only include colors that exist in data
  subtype_colors <- subtype_colors[names(subtype_colors) %in% df_sub$molecular_subtype]

  pdf(file.path(FIG_DIR, "km_curves_subtype.pdf"), width = 9, height = 7)
  p <- ggsurvplot(fit_sub, data = df_sub, pval = TRUE,
                  palette = subtype_colors,
                  legend.title = "Molecular Subtype",
                  xlab = "Time (Years)", ylab = "Overall Survival Probability",
                  risk.table = TRUE, risk.table.height = 0.25,
                  ggtheme = theme_bw(base_size = 14))
  print(p)
  dev.off()

  lr_sub <- survdiff(Surv(os_time_years, os_status) ~ molecular_subtype, data = df_sub)
  cat(sprintf("  Subtype log-rank p = %.4e\n", lr_sub$pvalue))
}

# ---- 3. KM by stage ----
cat("\nStep 3: KM curves by pathologic stage...\n")

if ("stage_simple" %in% colnames(clinical)) {
  df_stg <- clinical %>%
    filter(!is.na(stage_simple), stage_simple != "",
           !is.na(os_time_years), os_time_years > 0)

  fit_stg <- survfit(Surv(os_time_years, os_status) ~ stage_simple, data = df_stg)

  stage_colors <- c("Stage I" = "#1B9E77", "Stage II" = "#D95F02",
                    "Stage III" = "#7570B3", "Stage IV" = "#E7298A")
  stage_colors <- stage_colors[names(stage_colors) %in% df_stg$stage_simple]

  pdf(file.path(FIG_DIR, "km_curves_stage.pdf"), width = 9, height = 7)
  p <- ggsurvplot(fit_stg, data = df_stg, pval = TRUE,
                  palette = stage_colors,
                  legend.title = "Pathologic Stage",
                  xlab = "Time (Years)", ylab = "Overall Survival Probability",
                  risk.table = TRUE, risk.table.height = 0.25,
                  ggtheme = theme_bw(base_size = 14))
  print(p)
  dev.off()
  cat("  Stage KM curves saved.\n")
}

# ---- 4. Multivariate Cox regression ----
cat("\nStep 4: Multivariate Cox regression...\n")

cox_data <- data.frame(
  os_time   = clinical$os_time_years,
  os_status = clinical$os_status,
  age       = as.numeric(clinical$age_at_diagnosis),
  stringsAsFactors = FALSE
)

if ("positive_lymph_nodes" %in% colnames(clinical))
  cox_data$lymph_nodes <- as.numeric(clinical$positive_lymph_nodes)

if ("stage_simple" %in% colnames(clinical)) {
  cox_data$stage_II  <- as.numeric(clinical$stage_simple == "Stage II")
  cox_data$stage_III <- as.numeric(clinical$stage_simple == "Stage III")
  cox_data$stage_IV  <- as.numeric(clinical$stage_simple == "Stage IV")
}

if ("molecular_subtype" %in% colnames(clinical)) {
  cox_data$subtype_HER2 <- as.numeric(clinical$molecular_subtype == "HER2-enriched")
  cox_data$subtype_TN   <- as.numeric(clinical$molecular_subtype == "Triple Negative")
}

cox_data <- cox_data %>% filter(os_time > 0, !is.na(os_status))
cox_data[is.na(cox_data)] <- 0

# Build formula dynamically
covars <- setdiff(names(cox_data), c("os_time", "os_status"))
formula_str <- paste("Surv(os_time, os_status) ~", paste(covars, collapse = " + "))
cox_model <- coxph(as.formula(formula_str), data = cox_data)

cox_summary <- summary(cox_model)
cat(sprintf("  C-index: %.3f | LR p-value: %.2e\n",
            cox_summary$concordance[1], cox_summary$logtest["pvalue"]))

# Save Cox results
cox_table <- data.frame(
  Variable = rownames(cox_summary$coefficients),
  HR = cox_summary$coefficients[, "exp(coef)"],
  CI_lower = cox_summary$conf.int[, "lower .95"],
  CI_upper = cox_summary$conf.int[, "upper .95"],
  pvalue = cox_summary$coefficients[, "Pr(>|z|)"],
  stringsAsFactors = FALSE
)
write.csv(cox_table, file.path(TBL_DIR, "cox_regression.csv"), row.names = FALSE)

# Forest plot
pdf(file.path(FIG_DIR, "cox_forest_plot.pdf"), width = 10, height = 6)
tryCatch({
  ggforest(cox_model, data = cox_data, main = "BRCA: Multivariate Cox Regression")
}, error = function(e) {
  # Manual forest plot fallback
  conf_int <- cox_summary$conf.int
  plot_data <- data.frame(
    Variable = rownames(conf_int),
    HR = conf_int[, "exp(coef)"],
    Lower = conf_int[, "lower .95"],
    Upper = conf_int[, "upper .95"],
    stringsAsFactors = FALSE
  )
  ggplot(plot_data, aes(x = HR, y = rev(Variable))) +
    geom_point(size = 3, color = "#E41A1C") +
    geom_errorbarh(aes(xmin = Lower, xmax = Upper), height = 0.2) +
    geom_vline(xintercept = 1, linetype = "dashed") +
    labs(title = "BRCA: Multivariate Cox Regression", x = "Hazard Ratio (95% CI)", y = "") +
    theme_bw(base_size = 14) +
    scale_x_log10()
})
dev.off()
cat("  Cox forest plot saved.\n")

# ---- 5. LASSO-Cox Prognostic Signature ----
cat("\nStep 5: Building prognostic gene signature (LASSO-Cox)...\n")

expr_patients <- substr(colnames(tpm_mat), 9, 12)
clinical$pid_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)
match_idx <- match(expr_patients, clinical$pid_short)

surv_ok <- !is.na(clinical$os_time[match_idx]) & clinical$os_time[match_idx] > 0 &
           !is.na(clinical$os_status[match_idx])

cat(sprintf("  Patients with complete survival data: %d\n", sum(surv_ok)))

# Use top DEGs as candidate genes
deg_file <- file.path(TBL_DIR, "brca_degs_deseq2.csv")
if (file.exists(deg_file) && sum(surv_ok) > 100) {
  deg_table <- read.csv(deg_file, stringsAsFactors = FALSE)
  top_degs <- deg_table %>%
    filter(!is.na(padj), padj < 0.01, abs(log2FC) > 2) %>%
    pull(gene)

  degs_in_expr <- intersect(top_degs, rownames(tpm_mat))
  cat(sprintf("  Candidate prognostic genes: %d\n", length(degs_in_expr)))

  if (length(degs_in_expr) > 50) {
    X <- t(log2(tpm_mat[degs_in_expr, surv_ok] + 1))
    y <- Surv(
      clinical$os_time[match_idx][surv_ok] / 365.25,
      clinical$os_status[match_idx][surv_ok]
    )

    set.seed(42)
    cv_fit <- cv.glmnet(x = X, y = y, family = "cox", alpha = 1, nfolds = 5)

    coefs <- coef(cv_fit, s = "lambda.min")
    selected_idx <- which(as.matrix(coefs) != 0)
    selected_genes <- rownames(coefs)[selected_idx]
    selected_coefs <- coefs[selected_idx]

    cat(sprintf("  LASSO-Cox selected %d prognostic genes\n", length(selected_genes)))

    if (length(selected_genes) >= 3) {
      risk_score <- as.vector(X[, selected_genes, drop = FALSE] %*% selected_coefs)
      risk_group <- ifelse(risk_score > median(risk_score), "High Risk", "Low Risk")

      risk_df <- data.frame(
        time  = clinical$os_time[match_idx][surv_ok] / 365.25,
        status = clinical$os_status[match_idx][surv_ok],
        risk_group = risk_group,
        risk_score = risk_score,
        stringsAsFactors = FALSE
      )

      fit_risk <- survfit(Surv(time, status) ~ risk_group, data = risk_df)

      pdf(file.path(FIG_DIR, "km_curves_signature.pdf"), width = 9, height = 7)
      p <- ggsurvplot(fit_risk, data = risk_df, pval = TRUE,
                      palette = c("High Risk" = "#E41A1C", "Low Risk" = "#377EB8"),
                      legend.title = "Risk Group",
                      xlab = "Time (Years)", ylab = "Overall Survival Probability",
                      risk.table = TRUE, risk.table.height = 0.25,
                      ggtheme = theme_bw(base_size = 14))
      print(p)
      dev.off()

      prog_genes <- data.frame(Gene = selected_genes, Coefficient = selected_coefs,
                               stringsAsFactors = FALSE) %>% arrange(desc(abs(Coefficient)))
      write.csv(prog_genes, file.path(TBL_DIR, "prognostic_genes.csv"), row.names = FALSE)
      cat("  Prognostic signature KM curve saved.\n")
    }
  }
} else {
  cat("  Insufficient data for prognostic signature.\n")
}

cat("\n========== Survival Analysis Complete ==========\n")
