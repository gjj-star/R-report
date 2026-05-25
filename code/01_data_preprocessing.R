# ===========================================================================
# Script 01: BRCA数据预处理
# 项目：TCGA-BRCA多组学数据挖掘
# 说明：加载TCGA-BRCA表达和临床数据，执行清洗、标准化、缺失值处理
# 输入：brca_exp.Rdata, brca_clinical.Rdata
# 输出：处理后的RDS文件，供后续分析使用
# ===========================================================================

Sys.setlocale("LC_ALL", "English")
options(stringsAsFactors = FALSE, timeout = 600)

suppressPackageStartupMessages({
  library(tidyverse)
})

cat("\n========== Script 01: BRCA Data Preprocessing ==========\n\n")

# 0. 路径定义 ===============================================================
# Auto-detect project root: try here, then script dir, then working dir
BASE_DIR <- tryCatch(here::here(), error = function(e) {
  tryCatch(normalizePath(file.path(dirname(sys.frame(1)$ofile), "..")), error = function(e2) getwd())
})
DATA_DIR    <- file.path(BASE_DIR, "TCGA六种癌症转录本数据及临床数据",
                         "TCGA六种癌症转录本数据及临床数据")
OUTPUT_DIR  <- file.path(BASE_DIR, "data", "processed")
FIG_DIR     <- file.path(BASE_DIR, "results", "figures")
TBL_DIR     <- file.path(BASE_DIR, "results", "tables")

dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(FIG_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(TBL_DIR, showWarnings = FALSE, recursive = TRUE)

# ---- 1. 加载原始数据 ----
cat("Step 1: Loading TCGA-BRCA data...\n")

e <- new.env()
load(file.path(DATA_DIR, "brca_exp.Rdata"), envir = e)
counts_raw <- e$brca
rm(e); gc()

e <- new.env()
load(file.path(DATA_DIR, "brca_clinical.Rdata"), envir = e)
clinical_raw <- e$brca_clinical
rm(e); gc()

cat(sprintf("  Expression matrix: %d genes x %d samples\n", nrow(counts_raw), ncol(counts_raw)))
cat(sprintf("  Clinical data: %d patients x %d variables\n", nrow(clinical_raw), ncol(clinical_raw)))
cat(sprintf("  Data type: %s (range: %.0f - %.0f)\n", typeof(counts_raw), min(counts_raw), max(counts_raw)))

# ---- 2. TCGA Barcode解析 ----
cat("\nStep 2: Parsing TCGA barcodes (Tumor vs Normal)...\n")

barcodes    <- colnames(counts_raw)
patient_ids <- substr(barcodes, 9, 12)
sample_code <- substr(barcodes, 14, 15)

sample_info <- data.frame(
  barcode     = barcodes,
  patient_id  = patient_ids,
  sample_code = sample_code,
  sample_type = case_when(
    sample_code == "01" ~ "Tumor",
    sample_code == "06" ~ "Metastatic",
    sample_code == "11" ~ "Normal",
    TRUE                ~ "Other"
  ),
  stringsAsFactors = FALSE
)

cat("  Sample distribution:\n")
for (tp in names(table(sample_info$sample_type))) {
  cat(sprintf("    %s: %d\n", tp, table(sample_info$sample_type)[tp]))
}
write.csv(sample_info, file.path(OUTPUT_DIR, "sample_info.csv"), row.names = FALSE)

# ---- 3. Split Tumor/Normal ----
cat("\nStep 3: Splitting tumor and normal expression matrices...\n")

tumor_idx  <- which(sample_info$sample_type == "Tumor")
normal_idx <- which(sample_info$sample_type == "Normal")

counts_tumor  <- counts_raw[, tumor_idx, drop = FALSE]
counts_normal <- counts_raw[, normal_idx, drop = FALSE]

cat(sprintf("  Tumor: %d | Normal: %d\n", ncol(counts_tumor), ncol(counts_normal)))

# ---- 4. Low-expression gene filtering ----
cat("\nStep 4: Filtering low-expression genes...\n")

min_samples <- max(2, ceiling(ncol(counts_tumor) * 0.1))
keep_genes  <- rowSums(counts_tumor >= 10) >= min_samples

counts_tumor  <- counts_tumor[keep_genes, , drop = FALSE]
counts_normal <- counts_normal[keep_genes, , drop = FALSE]

cat(sprintf("  Before: %d genes | After: %d genes (%.1f%% retained)\n",
            length(keep_genes), sum(keep_genes), 100 * mean(keep_genes)))

# ---- 5. Gene identifier mapping ----
cat("\nStep 5: Gene identifier mapping (Ensembl -> Symbol)...\n")

ensembl_ids   <- rownames(counts_tumor)
ensembl_clean <- gsub("\\.\\d+$", "", ensembl_ids)
cat(sprintf("  Total Ensembl IDs: %d\n", length(ensembl_clean)))

gene_map_file <- file.path(OUTPUT_DIR, "ensembl_gene_map.rds")

if (file.exists(gene_map_file)) {
  cat("  Loading cached gene annotation...\n")
  gene_map <- readRDS(gene_map_file)
} else {
  # Try org.Hs.eg.db first (preferred), fall back to Ensembl FTP
  gene_map <- NULL

  if (requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
    cat("  Using org.Hs.eg.db for gene mapping...\n")
    library(org.Hs.eg.db)

    map_result <- tryCatch({
      AnnotationDbi::select(org.Hs.eg.db,
                            keys = ensembl_clean,
                            columns = c("SYMBOL", "GENETYPE"),
                            keytype = "ENSEMBL")
    }, error = function(e) NULL)

    if (!is.null(map_result) && nrow(map_result) > 0) {
      gene_map <- data.frame(
        gene_stable_id = map_result$ENSEMBL,
        gene_name      = map_result$SYMBOL,
        stringsAsFactors = FALSE
      )
      gene_map <- gene_map[!is.na(gene_map$gene_name), ]
      gene_map <- gene_map[!duplicated(gene_map$gene_stable_id), ]
      saveRDS(gene_map, gene_map_file)
      cat(sprintf("  org.Hs.eg.db mapped %d genes\n", nrow(gene_map)))
    }
  }

  # Fallback: Ensembl FTP
  if (is.null(gene_map)) {
    ftp_url <- "https://ftp.ensembl.org/pub/current_tsv/homo_sapiens/Homo_sapiens.GRCh38.113.gene.txt.gz"
    cat(sprintf("  Downloading Ensembl gene annotation from: %s\n", ftp_url))

    tmp_gz <- tempfile(fileext = ".gz")
    success <- tryCatch({
      download.file(ftp_url, tmp_gz, mode = "wb", timeout = 600)
      TRUE
    }, error = function(e) {
      cat(sprintf("  Download failed: %s\n", e$message))
      FALSE
    })

    if (success) {
      gene_data <- read.table(gzfile(tmp_gz), header = TRUE, sep = "\t",
                              quote = "", comment.char = "", stringsAsFactors = FALSE,
                              fill = TRUE, na.strings = "")
      unlink(tmp_gz)

      avail_cols <- intersect(c("gene_stable_id", "gene_name", "gene_biotype",
                                "description", "start_position", "end_position"),
                              colnames(gene_data))
      gene_map <- gene_data[, avail_cols]

      if (all(c("start_position", "end_position") %in% avail_cols)) {
        gene_map$gene_length <- gene_map$end_position - gene_map$start_position + 1
      }

      saveRDS(gene_map, gene_map_file)
      cat(sprintf("  Downloaded %d gene annotations\n", nrow(gene_map)))
    } else {
      cat("  Ensembl download also failed; using Ensembl IDs as symbols\n")
      gene_map <- data.frame(
        gene_stable_id = ensembl_clean,
        gene_name      = ensembl_clean,
        gene_biotype   = "unknown",
        description    = "unknown",
        gene_length    = 1000,
        stringsAsFactors = FALSE
      )
    }
  }
}

# ---- 6. Assign gene symbols and lengths ----
cat("\nStep 6: Assigning gene symbols and lengths...\n")

map_idx <- match(ensembl_clean, gene_map$gene_stable_id)

gene_symbol <- ifelse(!is.na(map_idx) & !is.na(gene_map$gene_name[map_idx]) &
                        gene_map$gene_name[map_idx] != "",
                      gene_map$gene_name[map_idx], ensembl_ids)
gene_symbol[is.na(gene_symbol) | gene_symbol == ""] <- ensembl_ids[is.na(gene_symbol) | gene_symbol == ""]

# Handle duplicate symbols
dup_idx <- which(duplicated(gene_symbol) | duplicated(gene_symbol, fromLast = TRUE))
gene_symbol[dup_idx] <- paste0(gene_symbol[dup_idx], "_", ensembl_clean[dup_idx])

gene_length_vec <- if ("gene_length" %in% colnames(gene_map) && !all(is.na(map_idx))) {
  len <- gene_map$gene_length[map_idx]
  len[is.na(len) | len <= 0] <- median(len[!is.na(len) & len > 0], na.rm = TRUE)
  if (all(is.na(len))) len <- rep(1000, length(len))
  len
} else {
  rep(1000, nrow(counts_tumor))
}

cat(sprintf("  Unique gene symbols: %d\n", length(unique(gene_symbol))))
cat(sprintf("  Median gene length: %.0f bp\n", median(gene_length_vec, na.rm = TRUE)))

rownames(counts_tumor)  <- gene_symbol
rownames(counts_normal) <- gene_symbol

# Gene annotation table
gene_anno <- data.frame(
  ensembl_id    = ensembl_ids,
  ensembl_clean = ensembl_clean,
  gene_symbol   = gene_symbol,
  gene_length   = gene_length_vec,
  stringsAsFactors = FALSE
)

# ---- 7. Counts -> TPM normalization ----
cat("\nStep 7: Counts -> TPM normalization...\n")

counts_to_tpm <- function(counts, gene_lengths) {
  rpk <- counts / (gene_lengths / 1000)
  sweep(rpk, 2, colSums(rpk) / 1e6, "/")
}

tpm_tumor  <- counts_to_tpm(counts_tumor, gene_length_vec)
tpm_normal <- counts_to_tpm(counts_normal, gene_length_vec)

cat(sprintf("  TPM range (log2): Tumor [%.2f, %.2f], Normal [%.2f, %.2f]\n",
            min(log2(tpm_tumor + 1)), max(log2(tpm_tumor + 1)),
            min(log2(tpm_normal + 1)), max(log2(tpm_normal + 1))))

# ---- 8. Clinical data cleaning ----
cat("\nStep 8: Cleaning clinical data...\n")

cat(sprintf("  Raw clinical columns: %s\n", paste(colnames(clinical_raw)[1:min(10, ncol(clinical_raw))], collapse = ", ")))
cat(sprintf("  Total columns: %d\n", ncol(clinical_raw)))

# Identify key clinical variables - flexible matching
col_names <- colnames(clinical_raw)

find_col <- function(patterns) {
  for (p in patterns) {
    idx <- grep(p, col_names, ignore.case = TRUE)
    if (length(idx) > 0) return(col_names[idx[1]])
  }
  return(NA_character_)
}

barcode_col   <- find_col("bcr_patient_barcode|patient_barcode|barcode")
vital_col     <- find_col("vital_status")
death_col     <- find_col("days_to_death")
followup_col  <- find_col("days_to_last_follow|last_follow")
age_col       <- find_col("age_at.*diagnosis|age_at_initial")
gender_col    <- find_col("^gender$")
race_col      <- find_col("race")
stage_col     <- find_col("stage.*pathologic|pathologic_stage|ajcc_pathologic_stage")
er_col        <- find_col("estrogen_receptor|er_status")
pr_col        <- find_col("progesterone_receptor|pr_status")
her2_col      <- find_col("her2_neu|her2_status")
hist_col      <- find_col("histological_type")
lymph_col     <- find_col("lymphnodes_positive|lymph_node")
radiation_col <- find_col("radiation_therapy")

cat(sprintf("  Identified clinical columns:\n"))
var_map <- list(
  barcode = barcode_col, vital = vital_col, death = death_col,
  followup = followup_col, age = age_col, gender = gender_col,
  race = race_col, stage = stage_col, er = er_col,
  pr = pr_col, her2 = her2_col, hist = hist_col,
  lymph = lymph_col, radiation = radiation_col
)
for (nm in names(var_map)) {
  cat(sprintf("    %-12s -> %s\n", nm, ifelse(is.na(var_map[[nm]]), "NOT FOUND", var_map[[nm]])))
}

# Build clean clinical dataframe
clinical <- data.frame(
  patient_id = if (!is.na(barcode_col)) clinical_raw[[barcode_col]] else rownames(clinical_raw),
  stringsAsFactors = FALSE
)

# Extract clinical variables
if (!is.na(vital_col))    clinical$vital_status <- clinical_raw[[vital_col]]
if (!is.na(death_col))    clinical$days_to_death <- as.numeric(clinical_raw[[death_col]])
if (!is.na(followup_col)) clinical$days_to_last_followup <- as.numeric(clinical_raw[[followup_col]])
if (!is.na(age_col))      clinical$age_at_diagnosis <- as.numeric(clinical_raw[[age_col]])
if (!is.na(gender_col))   clinical$gender <- clinical_raw[[gender_col]]
if (!is.na(race_col))     clinical$race <- clinical_raw[[race_col]]
if (!is.na(stage_col))    clinical$pathologic_stage <- clinical_raw[[stage_col]]
if (!is.na(er_col))       clinical$er_status <- clinical_raw[[er_col]]
if (!is.na(pr_col))       clinical$pr_status <- clinical_raw[[pr_col]]
if (!is.na(her2_col))     clinical$her2_status <- clinical_raw[[her2_col]]
if (!is.na(hist_col))     clinical$histological_type <- clinical_raw[[hist_col]]
if (!is.na(lymph_col))    clinical$positive_lymph_nodes <- as.numeric(clinical_raw[[lymph_col]])
if (!is.na(radiation_col)) clinical$radiation_therapy <- clinical_raw[[radiation_col]]

# Compute OS time and event
clinical <- clinical %>%
  mutate(
    vital_status = ifelse(grepl("Dead|DEAD|dead", vital_status), 1, 0),
    os_time = ifelse(vital_status == 1,
                     days_to_death,
                     days_to_last_followup),
    os_time = ifelse(is.na(os_time) | os_time <= 0, 1, os_time)
  )

# Simplified stage
if ("pathologic_stage" %in% colnames(clinical)) {
  clinical <- clinical %>%
    mutate(
      stage_simple = case_when(
        grepl("^(Stage I[A-C]?|I[A-C]?)$", pathologic_stage, ignore.case = TRUE) ~ "Stage I",
        grepl("Stage I$|^I$", pathologic_stage, ignore.case = TRUE) ~ "Stage I",
        grepl("Stage II[A-C]?$|^II[A-C]?$", pathologic_stage, ignore.case = TRUE) ~ "Stage II",
        grepl("Stage III[A-C]?$|^III[A-C]?$", pathologic_stage, ignore.case = TRUE) ~ "Stage III",
        grepl("Stage IV[A-C]?$|^IV[A-C]?$", pathologic_stage, ignore.case = TRUE) ~ "Stage IV",
        TRUE ~ NA_character_
      )
    )
}

# Molecular subtype classification (ER/PR/HER2)
if (all(c("er_status", "pr_status", "her2_status") %in% colnames(clinical))) {
  clinical <- clinical %>%
    mutate(
      er_status_clean = case_when(
        grepl("Positive|pos", er_status, ignore.case = TRUE) ~ "Positive",
        grepl("Negative|neg", er_status, ignore.case = TRUE) ~ "Negative",
        TRUE ~ NA_character_
      ),
      pr_status_clean = case_when(
        grepl("Positive|pos", pr_status, ignore.case = TRUE) ~ "Positive",
        grepl("Negative|neg", pr_status, ignore.case = TRUE) ~ "Negative",
        TRUE ~ NA_character_
      ),
      her2_status_clean = case_when(
        grepl("Positive|pos", her2_status, ignore.case = TRUE) ~ "Positive",
        grepl("Negative|neg", her2_status, ignore.case = TRUE) ~ "Negative",
        TRUE ~ NA_character_
      ),
      molecular_subtype = case_when(
        (er_status_clean == "Positive" | pr_status_clean == "Positive") & her2_status_clean == "Negative" ~ "Luminal A",
        (er_status_clean == "Positive" | pr_status_clean == "Positive") & her2_status_clean == "Positive" ~ "Luminal B",
        er_status_clean == "Negative" & pr_status_clean == "Negative" & her2_status_clean == "Positive" ~ "HER2-enriched",
        er_status_clean == "Negative" & pr_status_clean == "Negative" & her2_status_clean == "Negative" ~ "Triple Negative",
        TRUE ~ NA_character_
      )
    )
}

cat(sprintf("  Clean clinical variables: %d\n", ncol(clinical)))

# ---- 9. Missing value imputation ----
cat("\nStep 9: Missing value imputation...\n")

for (v in names(clinical)) {
  miss_n <- sum(is.na(clinical[[v]]))
  if (miss_n > 0) cat(sprintf("    %-25s: %d (%.1f%%)\n", v, miss_n, 100 * miss_n / nrow(clinical)))
}

# Impute: median for numeric, mode for categorical
impute_targets <- c("age_at_diagnosis", "positive_lymph_nodes", "os_time",
                     "stage_simple", "molecular_subtype", "race")
impute_targets <- intersect(impute_targets, colnames(clinical))

for (v in impute_targets) {
  nas <- is.na(clinical[[v]])
  if (any(nas)) {
    if (is.numeric(clinical[[v]])) {
      clinical[[v]][nas] <- median(clinical[[v]], na.rm = TRUE)
    } else {
      tbl <- table(clinical[[v]])
      if (length(tbl) > 0) clinical[[v]][nas] <- names(tbl)[which.max(tbl)]
    }
  }
}
cat("  Imputation complete.\n")

# ---- 10. Clinical-expression alignment ----
cat("\nStep 10: Aligning clinical and expression data...\n")

expr_patients_id <- substr(colnames(counts_tumor), 9, 12)
clinical$patient_id_short <- gsub("^TCGA-\\w{2}-", "", clinical$patient_id)

match_idx <- match(expr_patients_id, clinical$patient_id_short)
matched   <- !is.na(match_idx)

counts_tumor_final <- counts_tumor[, matched, drop = FALSE]
tpm_tumor_final    <- tpm_tumor[, matched, drop = FALSE]
clinical_final     <- clinical[match_idx[matched], , drop = FALSE]

cat(sprintf("  Matched: %d / %d (%.1f%%)\n",
            ncol(counts_tumor_final), ncol(counts_tumor), 100 * mean(matched)))

# ---- 11. Save processed data ----
cat("\nStep 11: Saving processed data...\n")

saveRDS(counts_tumor_final, file.path(OUTPUT_DIR, "brca_tumor_counts.rds"))
saveRDS(tpm_tumor_final, file.path(OUTPUT_DIR, "brca_tumor_tpm.rds"))
saveRDS(clinical_final, file.path(OUTPUT_DIR, "brca_clinical_clean.rds"))
saveRDS(gene_anno, file.path(OUTPUT_DIR, "gene_annotation.rds"))

if (ncol(counts_normal) > 0) {
  saveRDS(counts_normal, file.path(OUTPUT_DIR, "brca_normal_counts.rds"))
  saveRDS(tpm_normal, file.path(OUTPUT_DIR, "brca_normal_tpm.rds"))
}

# Full bundle
saveRDS(list(
  counts            = counts_tumor_final,
  tpm               = tpm_tumor_final,
  clinical          = clinical_final,
  gene_annotation   = gene_anno,
  sample_info       = sample_info,
  preprocessing_date = Sys.time()
), file.path(OUTPUT_DIR, "brca_tumor_processed.rds"))

# ---- 12. Summary report ----
cat("\n========== Preprocessing Summary ==========\n")
cat(sprintf("  Output: %s\n", OUTPUT_DIR))
cat(sprintf("  Tumor samples: %d (%d genes x %d samples)\n",
            ncol(counts_tumor_final), nrow(counts_tumor_final), ncol(counts_tumor_final)))
cat(sprintf("  Normal samples: %d\n", ncol(counts_normal)))
cat(sprintf("  Clinical variables: %d\n", ncol(clinical_final)))
if ("stage_simple" %in% colnames(clinical_final)) {
  cat(sprintf("  Staging: %s\n",
              paste(names(table(clinical_final$stage_simple)),
                    table(clinical_final$stage_simple), sep = "=", collapse = ", ")))
}
if ("molecular_subtype" %in% colnames(clinical_final)) {
  cat(sprintf("  Subtypes: %s\n",
              paste(names(table(clinical_final$molecular_subtype)),
                    table(clinical_final$molecular_subtype), sep = "=", collapse = ", ")))
}
cat(sprintf("  Survival: Death=%d, Alive=%d\n",
            sum(clinical_final$vital_status == 1, na.rm = TRUE),
            sum(clinical_final$vital_status == 0, na.rm = TRUE)))
cat("=========================================\n")

# Summary table
summary_table <- data.frame(
  Metric = c("Original genes", "Filtered genes", "Tumor samples", "Normal samples",
             "Clinical variables", "Death events", "Alive",
             "Stage I", "Stage II", "Stage III", "Stage IV"),
  Value = c(
    length(keep_genes), nrow(counts_tumor_final), ncol(counts_tumor_final),
    ncol(counts_normal), ncol(clinical_final),
    sum(clinical_final$vital_status == 1, na.rm = TRUE),
    sum(clinical_final$vital_status == 0, na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage I", na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage II", na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage III", na.rm = TRUE),
    sum(clinical_final$stage_simple == "Stage IV", na.rm = TRUE)
  )
)

if ("molecular_subtype" %in% colnames(clinical_final)) {
  for (st in c("Luminal A", "Luminal B", "HER2-enriched", "Triple Negative")) {
    summary_table <- rbind(summary_table, data.frame(
      Metric = st,
      Value = sum(clinical_final$molecular_subtype == st, na.rm = TRUE)
    ))
  }
}

write.csv(summary_table, file.path(OUTPUT_DIR, "preprocessing_summary.csv"), row.names = FALSE)
cat("Preprocessing complete!\n")
