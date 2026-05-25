# TCGA-BRCA Multi-Omics Data Mining

《生物大数据分析》期末报告项目

## 项目结构

```
r-work/
├── TCGA六种癌症转录本数据及临床数据/   # 原始数据 (.Rdata)
│   ├── brca_exp.Rdata
│   ├── brca_clinical.Rdata
│   └── ... (其他癌种数据)
├── code/                                 # R分析脚本
│   ├── 01_data_preprocessing.R          # 数据预处理
│   ├── 02_diff_expression.R             # DESeq2差异表达 + g:Profiler富集
│   ├── 03_classification.R              # RF + LASSO + XGBoost分类
│   ├── 04_clustering.R                  # PCA + t-SNE + K-means + 共识聚类
│   ├── 05_wgcna.R                       # WGCNA共表达网络
│   ├── 06_survival_analysis.R           # KM + Cox + LASSO-Cox
│   ├── 07_mutation_analysis.R           # maftools突变分析
│   ├── 08_visualization.R               # 论文级图表汇总
│   └── run_all.R                        # 一键运行全部脚本
├── data/processed/                       # 处理后数据 (RDS)
├── results/
│   ├── figures/                          # 分析图表 (PDF)
│   ├── figures_pub/                      # 论文用图
│   └── tables/                           # 结果表格 (CSV)
└── report/
    ├── BRCA_Final_Report.md             # 报告 (Markdown)
    ├── BRCA_Paper.tex                    # 报告 (LaTeX)
    └── references.bib                    # 参考文献
```

## 运行方法

### 前置要求

- R >= 4.0
- 依赖包：tidyverse, DESeq2, caret, randomForest, glmnet, WGCNA, survival, survminer, gprofiler2, pheatmap, ggrepel, Rtsne, xgboost (可选), ConsensusClusterPlus (可选), maftools (可选)

### 安装依赖

```r
# CRAN packages
install.packages(c("tidyverse", "caret", "randomForest", "glmnet",
                    "WGCNA", "survival", "survminer", "pheatmap",
                    "ggrepel", "Rtsne", "patchwork"))

# Bioconductor packages
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "gprofiler2", "ConsensusClusterPlus", "maftools"))

# Optional
install.packages("xgboost")
```

### 运行分析

```r
# 一键运行全部
source("code/run_all.R")

# 或逐个运行
source("code/01_data_preprocessing.R")
source("code/02_diff_expression.R")
# ...
```

## 分析方法

| 模块 | 方法 | 评分权重 |
|------|------|----------|
| 数据预处理 | 条形码解析、过滤、TPM标准化、缺失值处理 | 20分 |
| 差异表达 | DESeq2 + g:Profiler功能富集 | 属于数据挖掘(40分) |
| 机器学习分类 | Random Forest + LASSO + XGBoost | 属于数据挖掘(40分) |
| 聚类分析 | PCA + t-SNE + K-means + 共识聚类 | 属于数据挖掘(40分) |
| WGCNA | 共表达网络 + 枢纽基因识别 | 属于数据挖掘(40分) |
| 生存分析 | KM + Cox + LASSO-Cox | 属于数据挖掘(40分) |
| 可视化 | 火山图、热图、PCA、KM、森林图等 | 30分 |
