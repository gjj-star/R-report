# TCGA-BRCA 乳腺癌分子特征的数据挖掘与预后分析

《生物大数据分析》期末报告项目

---

## 项目概述

本项目基于 TCGA-BRCA（乳腺癌）多组学数据，综合运用差异表达分析、机器学习分类、聚类分析、WGCNA 共表达网络和生存分析等多种生物信息学方法，系统挖掘乳腺癌的分子特征并构建预后评估模型。

**最终提交版本：**
- 论文（Word）：`report/乳腺癌分子特征的数据挖掘与预后分析（最终版）.docx`
- 论文（PDF）：`report/乳腺癌分子特征的数据挖掘与预后分析（最终提交版）.pdf`

---

## 分析工作流

```
TCGA-BRCA 原始数据 (.Rdata)
        │
        ▼
┌──────────────────────────────────┐
│ ① 数据预处理                      │
│   · 临床数据清洗、分期标准化        │
│   · 表达矩阵过滤、TPM 标准化        │
│   · 样本条形码解析                  │
└──────────┬───────────────────────┘
           │
           ▼
┌──────────────────────────────────┐
│ ② 差异表达分析 (DESeq2)           │
│   · Tumor vs Normal 差异基因       │
│   · g:Profiler GO/KEGG 富集分析    │
└──────────┬───────────────────────┘
           │
     ┌─────┴─────┬─────────────┐
     ▼           ▼             ▼
┌─────────┐ ┌─────────┐ ┌──────────┐
│③ 分类    │ │④ 聚类    │ │⑤ WGCNA   │
│ RF+LASSO │ │ PCA+tSNE │ │ 共表达网络 │
│ +XGBoost │ │ K-means  │ │ 枢纽基因   │
└────┬─────┘ └────┬─────┘ └────┬─────┘
     │            │            │
     └────────────┼────────────┘
                  ▼
┌──────────────────────────────────┐
│ ⑥ 生存分析                        │
│   · KM 生存曲线（Stage/Subtype/Signature）│
│   · Cox 比例风险模型               │
│   · LASSO-Cox 预后模型             │
│   · Cox 森林图                    │
└──────────┬───────────────────────┘
           │
           ▼
┌──────────────────────────────────┐
│ ⑦ 突变分析 (maftools)             │
│   · 突变全景图                    │
│   · 显著突变基因筛选               │
└──────────┬───────────────────────┘
           │
           ▼
┌──────────────────────────────────┐
│ ⑧ 论文级可视化与报告撰写           │
│   · 18 张论文用图                  │
│   · 中文论文格式排版               │
└──────────────────────────────────┘
```

---

## 分析方法

| 模块 | 方法 | 工具/包 |
|------|------|---------|
| 数据预处理 | 条形码解析、低表达过滤、TPM 标准化、缺失值填补 | tidyverse, edgeR |
| 差异表达 | DESeq2 + g:Profiler GO/KEGG 富集分析 | DESeq2, gprofiler2 |
| 分类 | Random Forest + LASSO + XGBoost 三模型对比 | randomForest, glmnet, xgboost |
| 聚类 | PCA + t-SNE + K-means + 共识聚类 | prcomp, Rtsne, ConsensusClusterPlus |
| WGCNA | 加权基因共表达网络 + 模块-性状关联 + 枢纽基因 | WGCNA |
| 生存分析 | KM 曲线 + Cox 回归 + LASSO-Cox | survival, survminer, glmnet |
| 突变分析 | 突变全景图 + 瀑布图 | maftools |
| 可视化 | 火山图、热图、PCA、t-SNE、KM 曲线、森林图、ROC | ggplot2, pheatmap, ggrepel |

---

## 项目结构

```
r-work/
├── code/                         # R 分析脚本（核心代码）
│   ├── 01_data_preprocessing.R   # 数据预处理
│   ├── 02_diff_expression.R      # DESeq2 差异表达 + 富集分析
│   ├── 03_classification.R       # RF + LASSO + XGBoost 分类
│   ├── 04_clustering.R           # PCA + t-SNE + K-means + 共识聚类
│   ├── 05_wgcna.R                # WGCNA 共表达网络
│   ├── 06_survival_analysis.R    # KM + Cox + LASSO-Cox
│   ├── 07_mutation_analysis.R    # maftools 突变分析
│   ├── 08_visualization.R        # 论文级图表生成
│   └── run_all.R                 # 一键运行全部脚本
├── data/processed/               # 预处理后的数据 (RDS/CSV)
├── results/
│   ├── figures/                  # 分析图表 (PDF)
│   ├── figures_png/              # 分析图表 (PNG)
│   ├── figures_pub/              # 论文用图 (最终版)
│   ├── tables/                   # 结果表格 (CSV)
│   └── enrichment/               # 富集分析结果
├── report/
│   ├── 乳腺癌分子特征的数据挖掘与预后分析（最终版）.docx
│   ├── 乳腺癌分子特征的数据挖掘与预后分析（最终提交版）.pdf
│   ├── references.bib            # 参考文献
│   └── figures/                  # 论文嵌入图片
├── run_pipeline.R                # 便捷运行脚本
└── 《生物大数据分析》期末报告要求.txt
```

---

## 运行方法

### 环境要求

- R >= 4.0
- RStudio（推荐）

### 安装依赖

```r
# CRAN
install.packages(c("tidyverse", "caret", "randomForest", "glmnet",
                    "WGCNA", "survival", "survminer", "pheatmap",
                    "ggrepel", "Rtsne", "patchwork", "xgboost"))

# Bioconductor
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "gprofiler2", "ConsensusClusterPlus",
                        "maftools", "edgeR"))
```

### 执行分析

```r
# 一键运行全部流程
source("code/run_all.R")

# 或分步执行
source("code/01_data_preprocessing.R")
source("code/02_diff_expression.R")
# ... 依此类推
```

---

## 遇到的主要困难与解决方案

### 1. 临床数据列名映射问题

**困难：** TCGA 临床数据中的列名在不同癌种间不一致，`stage_simple` 和 `molecular_subtype` 等关键字段的提取逻辑依赖特定的列名模式，初次运行时因列名不匹配导致下游分析中断。

**解决：** 编写列名自动匹配逻辑，通过模糊搜索关键字段名（如 `ajcc_pathologic_stage`、`brca_subtype` 等），增加防御性检查和 fallback 机制。提交记录：`2650a64`。

### 2. WGCNA 并行计算与线程管理

**困难：** WGCNA 的 `blockwiseModules()` 默认开启多线程，在 Windows 环境下与 `mclapply` 不兼容，导致函数卡死或内存溢出。

**解决：** 显式设置 `nThreads = 1` 和 `allowWGCNAThreads()` 单线程模式，禁用 `mclapply` 回退到 `lapply`。在 `05_wgcna.R` 中添加线程控制开关。

### 3. 生存曲线（KM 曲线）内容错误

**困难：** 初次生成的 KM 曲线中，风险表（risk table）与曲线不对齐，部分分组的 p 值标注位置错误，Cox 森林图的 HR 置信区间计算有误。

**解决：** 重写 `ggsurvplot()` 调用参数，修正 risk table 的 `risk.table.col` 映射，确保每个分组的颜色、曲线和风险表行一致对应。使用 `survminer::ggforest()` 替代手动绘制 Cox 森林图。提交记录：`3c1b89e`。

### 4. 火山图可视化极端值处理

**困难：** 部分基因的 p 值极小（p < 1e-50），log2FC 超过 ±8，直接用原始值绘图导致散点过度压缩、标注重叠严重。

**解决：** 对 log2FC 设置 ±8 的 cap，对 -log10(p) 设置 50 的上限，使主分布区域清晰可见。极端值点用不同形状标注以示区分。见 `debug_volcano.R` 调试过程。

### 5. 中文论文排版迭代

**困难：** 报告要求使用中文撰写、宋体全文、代码附录。从 Markdown → LaTeX → Word 经历多轮转换：
- LaTeX 编译中文需要 `xelatex` + `ctex` 宏包，但本地 LaTeX 发行版缺少部分中文字体
- `pandoc` 转换 Markdown → Word 时，代码块和表格格式丢失
- Word 中图片嵌入后分辨率下降，编号无法自动更新

**解决：** 最终采用 Word 直接排版，全程手动调整图片锚点和交叉引用。18 张图片统一 300 DPI 嵌入。代码段以等宽字体 Consolas 嵌入方法章节。经过 `v1 → v2 → v3 → v4` 四轮格式微调后定稿。提交记录：`38a5920`、`8ad28d0`、`213717c`。

### 6. 数据文件路径跨环境兼容

**困难：** 原始数据 `.Rdata` 文件位于项目根目录的 `TCGA六种癌症转录本数据及临床数据/` 夹中，路径含中文，在不同系统和工作目录下读取不稳定。

**解决：** 在 `run_all.R` 开头统一设置 `setwd()` 并定义所有路径常量，各子脚本使用基于项目根目录的相对路径读取数据。

---

## 关键结果摘要

- **差异表达基因：** 鉴定出 ~2000 个显著 DEGs（|log2FC| > 1, padj < 0.05），富集于细胞周期、PI3K-Akt、MAPK 等癌症相关通路
- **分类性能：** Random Forest 达到最高 AUC > 0.98，XGBoost 次之，LASSO 提供可解释性
- **分子亚型聚类：** t-SNE 清晰分离 LumA/LumB/HER2/Basal/Normal 五类亚型
- **WGCNA 关键模块：** 识别出与肿瘤分期显著相关的 blue/turquoise 模块，筛选出 TOP2A、BIRC5 等枢纽基因
- **预后模型：** LASSO-Cox 筛选出 12 基因 signature，KM 曲线显示高风险组生存率显著降低（p < 0.0001）

---

## 评分标准对应

| 评分项 | 满分 | 对应内容 |
|--------|------|----------|
| 数据背景介绍 | 10 分 | 报告第 1 章 |
| 数据预处理 | 20 分 | `01_data_preprocessing.R` + 报告第 2.1-2.3 节 |
| 数据挖掘算法 | 40 分 | `02-07_*.R`（差异表达、分类、聚类、WGCNA、生存分析、突变分析） |
| 数据可视化 | 30 分 | `08_visualization.R` + 报告中 18 张图表 |
