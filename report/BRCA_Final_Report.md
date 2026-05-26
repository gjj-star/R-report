# 基于TCGA多组学数据挖掘的乳腺癌分子特征综合分析

**《生物大数据分析》期末报告**

---

## 摘要

乳腺癌是全球发病率最高的恶性肿瘤之一，其高度的分子异质性为精准诊疗带来了巨大挑战。本研究基于TCGA-BRCA数据集（1,105例肿瘤 + 113例正常组织，25,981个基因），运用DESeq2差异表达分析、Random Forest/LASSO分类、K-means聚类、WGCNA共表达网络、Cox生存分析等数据挖掘方法，系统分析乳腺癌分子特征。共鉴定6,768个差异表达基因（上调4,334，下调2,434）；Random Forest和LASSO分类准确率均为75.8%；PCA分析显示PC1=15.3%，K-means最优聚类数K=2；WGCNA识别9个共表达模块，枢纽基因模块隶属度达0.959；多变量Cox回归C-index=0.767，Stage IV风险比达29.82。

**关键词：** 乳腺癌；TCGA；数据挖掘；差异表达分析；WGCNA；机器学习分类；生存分析；R语言

---

## 1 引言

### 1.1 研究背景

据国际癌症研究机构（IARC）GLOBOCAN 2022统计，全球每年新发癌症约2,000万例，死亡约970万例[1]。女性乳腺癌以约230万例新发病例居全球恶性肿瘤发病率首位，同年导致约67万例死亡。在中国，乳腺癌发病率和死亡率持续上升，且发病年龄显著低于欧美国家，构成重大的公共卫生负担。

乳腺癌的发生发展涉及多基因突变累积、表观遗传重塑及信号通路交互失调等复杂分子过程[2]。高通量测序技术的普及使大规模癌症基因组数据的积累成为可能。以癌症基因组图谱（The Cancer Genome Atlas, TCGA）为代表的国际合作项目，为研究者提供了涵盖基因组、转录组、表观基因组及蛋白质组的海量公开数据资源[3]，使得基于数据挖掘的系统生物学研究范式得以实践。

### 1.2 乳腺癌分子分型

Perou等[4]于2000年首次基于基因表达谱提出分子分型体系。当前国际广泛认可的固有分子亚型包括：（1）Luminal A型（约占50%-60%），特征为ER阳性、HER2阴性；（2）Luminal B型（约占15%-20%）；（3）HER2过表达型（约占10%-15%）；（4）三阴性/Basal-like型（约占15%-20%），侵袭性最强。

### 1.3 研究内容与论文结构

本研究以TCGA-BRCA数据为对象，综合运用差异表达分析、机器学习分类、聚类分析、WGCNA共表达网络和生存分析等方法，系统开展乳腺癌分子特征挖掘。论文结构：第2章材料与方法；第3章结果；第4章讨论；第5章结论。

---

## 2 材料与方法

### 2.1 数据来源

本研究数据来自TCGA-BRCA项目，包含mRNA转录组数据（HiSeq RNA-seq, STAR-Counts流程）和临床注释数据。TCGA样本条形码采用标准化分层命名规则，前12个字符为患者级唯一标识符，样本类型代码（第14-15位）中01表示原发实体瘤，11表示癌旁正常组织。

### 2.2 数据预处理

#### 2.2.1 样本分类与基因过滤

解析TCGA条形码区分肿瘤与正常样本。采用独立过滤策略，保留至少在10%肿瘤样本中counts≥10的基因（参考DESeq2推荐准则）。

#### 2.2.2 基因标识符映射与TPM标准化

通过org.Hs.eg.db注释包将Ensembl Gene ID映射为标准Gene Symbol，对一对多映射保留平均表达量最高的条目。采用Transcripts Per Million（TPM）进行归一化：

$$TPM_i = \frac{Counts_i / (GeneLength_i / 1000)}{\sum_j [Counts_j / (GeneLength_j / 1000)]} \times 10^6$$

#### 2.2.3 临床数据清洗

从原始临床变量中提取年龄、病理分期、ER/PR/HER2状态等核心变量。分子亚型按ER/PR/HER2免疫组化状态分类。缺失值采用中位数（数值型）或众数（分类型）填补。最终将表达矩阵与临床数据以患者ID为锚点进行对齐，匹配1,105例患者。

### 2.3 数据挖掘算法

#### 2.3.1 差异表达分析

采用DESeq2[5]进行肿瘤组织（n=1,105）与正常组织（n=113）的差异表达分析。DESeq2基于负二项分布对RNA-seq计数数据建模，通过经验贝叶斯收缩改进离散度估计的稳定性，使用Wald检验评估组间差异的显著性。显著性阈值取FDR<0.05（Benjamini-Hochberg校正）且|log2FC|>1。功能富集分析采用g:Profiler[8]，覆盖GO、KEGG和Reactome数据库。

#### 2.3.2 机器学习分类

比较Random Forest（ntree=500）[14]和LASSO多分类回归（family="multinomial", α=1, 5折交叉验证）[13]的分类性能。输入特征为Top 500可变基因的log2(TPM+1)表达矩阵。训练集与测试集按7:3分层划分（set.seed=42）。评估指标包括准确率、Kappa系数和Macro F1分数。

#### 2.3.3 聚类分析

综合运用PCA（Top 2000可变基因）、t-SNE（perplexity=30, max_iter=1000）、K-means（K=2至10评估轮廓系数）和层次聚类（Ward.D2方法，Pearson相关系数距离）。

#### 2.3.4 WGCNA共表达网络分析

WGCNA[6]使用Top 5000可变基因的log2(TPM+1)矩阵。分析流程：（1）计算基因间Pearson相关系数构建相似性矩阵；（2）选择软阈值幂函数（scale-free R²>0.8）构建邻接矩阵；（3）构建拓扑重叠矩阵（TOM）；（4）动态树切割识别模块（minModuleSize=30, mergeCutHeight=0.25）；（5）计算模块特征基因与临床性状的相关性；（6）通过模块内连接度（Module Membership, MM）识别枢纽基因。

#### 2.3.5 生存分析

采用Kaplan-Meier生存曲线（log-rank检验）和Cox比例风险回归[16]评估临床及分子特征的预后价值。Cox模型纳入年龄、病理分期和分子亚型。此外，采用LASSO-Cox[7]方法构建预后基因标记，按中位风险得分将患者分为高风险组和低风险组。

### 2.4 分析环境

R 4.6.0 + Bioconductor 3.23，Windows 11环境。主要R包：DESeq2、caret、randomForest、glmnet、WGCNA、survival、survminer、gprofiler2、pheatmap、ggplot2。

---

## 3 结果

### 3.1 数据预处理结果

原始表达矩阵包含60,660个基因，经低表达过滤后保留25,981个（42.8%）。肿瘤样本1,105例，正常组织113例。分子亚型分布：Luminal A 827例、Luminal B 126例、HER2-enriched 37例、Triple Negative 115例。分期分布：Stage I 183例、II 651例、III 251例、IV 20例。生存信息：死亡87例，存活1,018例。

**表1 数据概览**

| 指标 | 数值 |
|------|------|
| 过滤后基因数 | 25,981 |
| 肿瘤样本 | 1,105 |
| 正常样本 | 113 |
| 临床变量 | 21 |
| 死亡/存活 | 87 / 1,018 |

![图1 TCGA-BRCA样本分布与病理分期分布](figures/fig1_sample_overview.png)

### 3.2 差异表达分析

DESeq2共鉴定6,768个显著差异表达基因（|log2FC|>1, padj<0.05），其中上调4,334个，下调2,434个。上调基因数量约为下调的1.78倍。

![图2 差异表达火山图](figures/fig2_volcano_enhanced.png)

![图3 Top 50差异基因热图](figures/fig5_deg_heatmap.png)

功能富集分析显示上调基因富集于细胞周期、DNA复制等增殖相关通路；下调基因富集于ECM组织、细胞黏附等微环境通路。

### 3.3 机器学习分类

Random Forest和LASSO均达到75.8%准确率。LASSO的Kappa系数（0.332）高于RF（0.281），验证了亚型间转录组差异的稀疏性。

![图4 分类模型性能比较](figures/fig6_classification.png)

![图5 ROC曲线](figures/fig14_roc.png)

### 3.4 聚类分析

PCA显示PC1=15.3%、PC2=8.5%，ER+与ER-沿PC1分离明显。K-means轮廓系数K=2时最大（0.180）。

![图6 PCA散点图](figures/fig4_pca_subtype.png)

![图7 t-SNE降维图](figures/fig13_tsne.png)

![图8 轮廓系数分析](figures/fig12_silhouette.png)

![图9 层次聚类热图](figures/fig11_hclust_heatmap.png)

### 3.5 WGCNA共表达网络分析

软阈值power=8（R²=0.887），识别9个共表达模块。枢纽基因MM最高达0.966。

![图10 WGCNA软阈值选择](figures/fig8_wgcna_soft_power.png)

![图11 模块树状图](figures/fig16_wgcna_dendrogram.png)

![图12 模块-性状热图](figures/fig7_wgcna_module_trait.png)

![图13 基因网络热图](figures/fig17_wgcna_network.png)

### 3.6 生存分析

KM曲线显示分期越高预后越差。Cox回归C-index=0.767，Stage IV HR=29.82。

![图14 KM曲线（分期）](figures/km_curves_stage.png)

![图15 KM曲线（亚型）](figures/km_curves_subtype.png)

**表2 多变量Cox回归结果**

| 变量 | HR | 95% CI | p值 |
|------|-----|--------|-----|
| 年龄 | 1.04 | 1.02-1.06 | 7.4×10⁻⁶ |
| Stage II | 3.24 | 1.28-8.19 | 0.013 |
| Stage III | 6.00 | 2.31-15.58 | 2.4×10⁻⁴ |
| Stage IV | 29.82 | 10.14-87.76 | 7.0×10⁻¹⁰ |
| Triple Negative | 2.00 | 1.08-3.71 | 0.028 |

![图16 Cox森林图](figures/cox_forest_plot.png)

![图17 预后标记KM曲线](figures/km_curves_signature.png)

---

## 4 讨论

本研究构建了覆盖差异表达、分类、聚类、WGCNA和生存分析的完整流程。DESeq2鉴定6,768个DEGs，上调基因富集于细胞周期通路。LASSO以较少特征达到75.8%准确率。PCA和K-means表明ER状态为最强驱动因素。WGCNA识别9个模块及枢纽基因。Cox回归证实Stage IV HR=29.82为最强预后因子。

局限性：仅用mRNA转录组，未纳入多组学；事件率7.9%限制统计功效；超参数未调优。未来方向：多组学整合、独立队列验证、枢纽基因功能实验、单细胞测序。

---

## 5 结论

（1）DESeq2鉴定6,768个差异基因，上调基因集中于细胞周期通路。

（2）LASSO/RF分类准确率75.8%，LASSO实现特征稀疏压缩。

（3）PCA（PC1=15.3%）和K-means（K=2最优）表明ER状态为最强驱动因素。

（4）WGCNA识别9个模块，枢纽基因MM最高达0.966。

（5）Cox回归C-index=0.767，Stage IV HR=29.82。

---

## 参考文献

[1] Sung H, et al. CA Cancer J Clin, 2024, 74(3): 229-263.

[2] Hanahan D, Weinberg RA. Cell, 2011, 144(5): 646-674.

[3] Cancer Genome Atlas Network. Nature, 2012, 490(7418): 61-70.

[4] Perou CM, et al. Nature, 2000, 406(6797): 747-752.

[5] Love MI, et al. Genome Biol, 2014, 15(12): 550.

[6] Langfelder P, Horvath S. BMC Bioinformatics, 2008, 9: 559.

[7] Simon N, et al. J Stat Softw, 2011, 39(5): 1-13.

[8] Kolberg L, et al. Nucleic Acids Res, 2023, 51(W1): W207-W212.

[9] Agrawal R, Srikant R. VLDB, 1994: 487-499.

[10] Wilkerson MD, Hayes DN. Bioinformatics, 2010, 26(12): 1572-1573.

[11] Colaprico A, et al. Nucleic Acids Res, 2016, 44(8): e71.

[12] Mayakonda A, et al. Genome Res, 2018, 28(11): 1747-1756.

[13] Friedman J, et al. J Stat Softw, 2010, 33(1): 1-22.

[14] Breiman L. Mach Learn, 2001, 45(1): 5-32.

[15] Chen T, Guestrin C. KDD, 2016: 785-794.

[16] Therneau TM, Grambsch PM. Springer, 2000.

[17] Tibshirani R. J R Stat Soc B, 1996, 58(1): 267-288.

---

## 附录A 核心R代码

### A.1 数据预处理

```r
load("brca_exp.Rdata")
load("brca_clinical.Rdata")
barcodes <- colnames(brca)
sample_code <- substr(barcodes, 14, 15)
tumor_idx  <- which(sample_code == "01")
normal_idx <- which(sample_code == "11")
counts_tumor  <- brca[, tumor_idx]
counts_normal <- brca[, normal_idx]
keep <- rowSums(counts_tumor >= 10) >=
  ceiling(ncol(counts_tumor) * 0.1)
counts_tumor <- counts_tumor[keep, ]
```

```r
library(org.Hs.eg.db)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
  keys = ensembl_clean,
  columns = "SYMBOL", keytype = "ENSEMBL")
counts_to_tpm <- function(counts, gene_lengths) {
  rpk <- counts / (gene_lengths / 1000)
  sweep(rpk, 2, colSums(rpk) / 1e6, "/")
}
tpm_tumor <- counts_to_tpm(counts_tumor, gene_length_vec)
```

### A.2 差异表达分析

```r
library(DESeq2)
dds <- DESeqDataSetFromMatrix(
  countData = round(counts_combined),
  colData   = col_data,
  design    = ~ condition)
keep <- rowSums(counts(dds) >= 10) >=
  min(10, ncol(dds) / 10)
dds <- dds[keep, ]
dds <- DESeq(dds)
res <- results(dds,
  contrast = c("condition", "Tumor", "Normal"),
  alpha = 0.05)
deg <- subset(as.data.frame(res),
  padj < 0.05 & abs(log2FoldChange) > 1)
```

### A.3 机器学习分类

```r
library(caret)
library(randomForest)
library(glmnet)
log_tpm <- log2(tpm_sub + 1)
top500 <- names(sort(apply(log_tpm, 1, var),
  decreasing = TRUE))[1:500]
X <- t(log_tpm[top500, ])
set.seed(42)
idx <- createDataPartition(y, p = 0.7, list = FALSE)
X_train <- X[idx, ]; X_test <- X[-idx, ]
y_train <- y[idx];    y_test <- y[-idx, ]
rf <- randomForest(x = X_train, y = y_train,
  ntree = 500, importance = TRUE)
cv_lasso <- cv.glmnet(x = X_train, y = y_train,
  family = "multinomial", alpha = 1, nfolds = 5)
```

### A.4 聚类分析

```r
pca <- prcomp(t(log_tpm_top),
  center = TRUE, scale. = TRUE)
library(Rtsne)
tsne <- Rtsne(t(log_tpm_top),
  perplexity = 30, max_iter = 1000)
set.seed(42)
for (k in 2:10) {
  km <- kmeans(t(log_tpm_top), centers = k, nstart = 25)
  sil[k] <- mean(cluster::silhouette(
    km$cluster, dist(t(log_tpm_top)))[, 3])
}
```

### A.5 WGCNA共表达网络分析

```r
library(WGCNA)
sft <- pickSoftThreshold(datExpr,
  powerVector = 1:20, networkType = "signed")
net <- blockwiseModules(datExpr,
  power = soft_power, TOMType = "signed",
  minModuleSize = 30, mergeCutHeight = 0.25)
module_colors <- labels2colors(net$colors)
MEs <- net$MEs
module_trait_cor <- cor(MEs, traits, use = "p")
```

### A.6 生存分析

```r
library(survival)
library(survminer)
fit <- survfit(Surv(os_time_years, os_status)
  ~ stage_simple, data = clinical)
cox_model <- coxph(Surv(os_time, os_status) ~
  age + stage_II + stage_III + stage_IV +
  subtype_TN, data = cox_data)
cv_cox <- cv.glmnet(x = X,
  y = Surv(time, status),
  family = "cox", alpha = 1, nfolds = 5)
```
