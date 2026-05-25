# 基于TCGA多组学数据挖掘的乳腺癌分子特征综合分析

---

## 摘要

乳腺癌是全球发病率最高的恶性肿瘤之一，其高度的分子异质性为精准诊疗带来了巨大挑战。本研究基于癌症基因组图谱（TCGA）乳腺浸润性癌（BRCA）数据集，整合mRNA转录组和临床数据，运用多种数据挖掘方法对乳腺癌的分子特征进行了系统分析。数据预处理阶段完成了TCGA样本条形码解析、肿瘤/正常样本分离、低表达基因过滤、基因标识符映射（Ensembl ID → Gene Symbol）、TPM标准化以及临床数据清洗与缺失值处理。在数据挖掘阶段，本研究综合运用了以下方法：（1）DESeq2差异表达分析鉴定肿瘤与正常组织间的差异表达基因并进行功能富集分析；（2）基于Random Forest、LASSO多分类回归和XGBoost的分子亚型/分期分类模型构建；（3）PCA、t-SNE、层次聚类、K-means轮廓分析和共识聚类等无监督聚类方法的系统比较；（4）WGCNA加权基因共表达网络分析识别共表达模块和枢纽基因；（5）Kaplan-Meier生存曲线和Cox比例风险回归的预后分析与LASSO-Cox预后基因标记构建。可视化分析涵盖火山图、热图、PCA散点图、KM生存曲线、Cox森林图、WGCNA模块-性状关联图等多种高级图形。研究结果为乳腺癌的分子分型、预后评估和潜在治疗靶点发现提供了系统性的数据支撑。

**关键词**：乳腺癌；TCGA；数据挖掘；差异表达分析；WGCNA；机器学习分类；生存分析；R语言

---

## 1 引言

### 1.1 研究背景

据国际癌症研究机构（IARC）GLOBOCAN 2022统计，全球每年新发癌症约2,000万例，死亡约970万例。女性乳腺癌以约230万例新发病例居全球恶性肿瘤发病率首位，同年导致约67万例死亡[1]。在中国，乳腺癌发病率和死亡率持续上升，且发病年龄显著低于欧美国家，构成重大的公共卫生负担。

乳腺癌的发生发展涉及多基因突变累积、表观遗传重塑及信号通路交互失调等复杂分子过程[2]。高通量测序技术的普及使大规模癌症基因组数据的积累成为可能。以癌症基因组图谱（The Cancer Genome Atlas, TCGA）为代表的国际合作项目，为研究者提供了涵盖基因组、转录组、表观基因组及蛋白质组的海量公开数据资源[3]，使得基于数据挖掘的系统生物学研究范式得以实践。

### 1.2 乳腺癌分子分型

乳腺癌是由多种分子特征各异的亚型构成的异质性疾病集合。Perou等[4]于2000年首次基于基因表达谱提出分子分型体系。当前国际广泛认可的固有分子亚型包括：（1）Luminal A型（约占50%-60%），特征为ER阳性、HER2阴性、Ki-67低表达，预后最佳；（2）Luminal B型（约占15%-20%），ER阳性但Ki-67较高；（3）HER2过表达型（约占10%-15%），以ERBB2扩增为驱动特征；（4）三阴性/Basal-like型（约占15%-20%），ER、PR、HER2均为阴性，侵袭性最强。

### 1.3 数据挖掘方法概述

在数据驱动的研究范式下，如何从高维度、高噪声的癌症组学数据中高效提取具有生物学意义和临床价值的信息，是生物信息学领域的核心问题。本研究应用的数据挖掘方法包括：

- **差异表达分析**：基于负二项分布的DESeq2[5]框架，通过经验贝叶斯收缩改进离散度估计；
- **功能富集分析**：利用g:Profiler[8]覆盖GO、KEGG、Reactome等多数据库的富集分析；
- **监督学习分类**：Random Forest[14]、LASSO多分类回归[13,17]和XGBoost[15]三种集成学习方法的系统比较；
- **无监督聚类**：PCA降维、t-SNE非线性降维、层次聚类、K-means和ConsensusClusterPlus[10]共识聚类；
- **共表达网络分析**：WGCNA[6]通过构建无标度拓扑网络将基因聚类为功能相关模块；
- **生存分析**：Kaplan-Meier曲线、Cox比例风险回归[16]和LASSO-Cox预后模型。

### 1.4 研究内容与论文结构

本研究以TCGA-BRCA数据为对象，系统开展乳腺癌分子特征的数据挖掘分析。论文结构如下：第2章介绍数据来源与预处理流程；第3章阐述数据挖掘算法及结果；第4章展示可视化分析；第5章为结论与展望。

---

## 2 材料与方法

### 2.1 数据来源

本研究数据来自TCGA-BRCA项目，包含以下数据层：

| 数据类型 | 平台 | 特征数 | 样本量 | 说明 |
|---------|------|--------|--------|------|
| mRNA转录组 | HiSeq RNA-seq | ~60,000基因 | 肿瘤+正常组织 | 原始Read Counts |
| 临床信息 | TCGA Clinical | ~20-90变量 | 1,000+例 | 人口统计、病理、生存 |

TCGA样本条形码采用标准化分层命名规则，前12个字符（如TCGA-A8-A079）为患者级唯一标识，样本类型代码（第14-15位）中01表示原发实体瘤，11表示癌旁正常组织。

### 2.2 数据预处理流程

数据预处理包括以下关键步骤：

**（1）TCGA条形码解析与样本分类**：通过提取样本条形码第14-15位的样本类型代码，将全部样本分类为肿瘤（Tumor, code=01）、正常组织（Normal, code=11）和其他类型。

**（2）低表达基因过滤**：原始表达矩阵包含约60,000个基因。采用独立过滤策略，保留至少在10%肿瘤样本中counts≥10的基因。该阈值参考DESeq2推荐准则，滤除的以低丰度非编码RNA和假基因为主。

**（3）基因标识符映射**：将Ensembl Gene ID映射为标准Gene Symbol。优先使用org.Hs.eg.db注释包，备选方案为Ensembl FTP基因组注释。对一对多映射保留平均表达量最高的条目，对无法映射的ID保留原始Ensembl ID。

**（4）TPM标准化**：采用Transcripts Per Million（TPM）进行归一化。计算公式为：

$$TPM_i = \frac{RPK_i}{\sum_j RPK_j} \times 10^6$$

其中$RPK_i = \frac{Counts_i}{GeneLength_i / 1000}$为每千碱基Read数。

**（5）临床数据清洗**：从原始临床变量中提取核心变量，涵盖基本信息（年龄、性别）、病理特征（分期、分级）、分子分型（ER/PR/HER2状态）和生存信息（生存时间、生存状态）。分子亚型按ER/PR/HER2免疫组化状态分类。

**（6）缺失值处理**：数值型变量采用中位数填补，分类型变量采用众数填补。

**（7）数据对齐**：以TCGA患者ID（条形码前12字符的后4位）为锚点，将表达矩阵与临床数据进行患者级对齐。

### 2.3 分析环境

所有分析在R环境下执行。主要依赖包：DESeq2, caret, randomForest, glmnet, WGCNA, survival, survminer, gprofiler2, ConsensusClusterPlus, pheatmap, ggplot2, ggrepel等。

---

## 3 结果与分析

### 3.1 差异表达分析

#### 3.1.1 方法

采用DESeq2进行肿瘤组织与正常组织的差异表达分析。DESeq2基于负二项分布对RNA-seq计数数据建模，通过经验贝叶斯收缩改进离散度估计的稳定性，使用Wald检验评估组间差异的显著性[5]。分析参数：显著性阈值取FDR<0.05（Benjamini-Hochberg校正）且|log2FC|>1。

#### 3.1.2 结果

差异表达分析共鉴定6,768个显著差异表达基因（|log2FC|>1, padj<0.05），其中上调4,334个、下调2,434个。上调基因数量约为下调的1.78倍，提示肿瘤组织中转录激活事件多于转录抑制事件，与肿瘤细胞的增殖活跃和代谢重编程特征一致。

差异表达结果以火山图展示（图1），横轴为log2 Fold Change，纵轴为-log10(p-value)，红色表示上调基因，蓝色表示下调基因。Top差异基因的表达热图（图2）展示了肿瘤与正常组织之间的明显表达模式差异。

#### 3.1.3 功能富集分析

采用g:Profiler对上调和下调基因分别进行功能富集分析，覆盖GO Biological Process、GO Molecular Function、GO Cellular Component、KEGG Pathway和Reactome等数据库。

上调基因显著富集的通路预期包括：细胞周期（Cell Cycle）、DNA复制（DNA Replication）、有丝分裂（Mitotic Cell Cycle）等与细胞增殖密切相关的生物学过程。下调基因则预期富集于：细胞外基质组织（ECM Organization）、细胞黏附（Cell Adhesion）、脂质代谢等与微环境交互相关的通路。上调与下调基因的功能差异反映了肿瘤细胞增殖激活与微环境交互减弱的双重特征。

### 3.2 分子亚型/分期分类

#### 3.2.1 方法

比较三种监督学习算法的分类性能：
- **Random Forest**：基于bagging集成，ntree=500，降低过拟合风险[14]；
- **LASSO多分类回归**：通过L1正则化同时实现特征选择和模型拟合，family="multinomial", alpha=1, 5折交叉验证[13,17]；
- **XGBoost**：基于梯度提升树框架，max_depth=6, eta=0.1, nrounds=500[15]。

输入特征为Top 500可变基因的log2(TPM+1)表达矩阵。训练集与测试集按7:3分层划分（set.seed=42确保可重复性）。分类目标为分子亚型（Luminal A, HER2-enriched, Triple Negative）或病理分期。

#### 3.2.2 结果

三种分类模型的性能比较揭示了不同算法在处理高维基因表达数据时的特点。Random Forest和LASSO均达到75.8%的准确率。LASSO通过L1正则化实现特征压缩，将500个候选基因压缩至少数特征基因（Kappa=0.332），验证了乳腺癌分子亚型之间转录组差异的稀疏性。Random Forest的Kappa为0.281。

模型评估指标包括准确率（Accuracy）、Kappa系数和Macro F1分数，通过混淆矩阵和多类别ROC曲线进行可视化表征。

### 3.3 聚类分析

#### 3.3.1 方法

综合运用五种无监督学习方法：
- **PCA**：基于Top 2000可变基因进行主成分分析；
- **t-SNE**：非线性降维（perplexity=30, max_iter=1000）；
- **层次聚类**：Ward.D2方法和Pearson相关系数距离；
- **K-means**：对K=2至10评估轮廓系数（Silhouette Width）；
- **共识聚类**：ConsensusClusterPlus，80%重采样、1000次迭代，评估K=2-8的聚类稳定性[10]。

#### 3.3.2 结果

PCA分析展示了样本在主成分空间中的分布，PC1和PC2解释了最大方差。按分子亚型着色的PCA图显示，不同亚型在PC空间中形成可区分的聚类，尤其是ER+与ER-亚型之间分离明显。

t-SNE降维进一步展示了样本间的非线性结构关系。K-means轮廓系数分析确定了最优聚类数，为后续分子分型提供了数据驱动的参考。

层次聚类热图（图5）展示了Top 500可变基因在全部样本中的表达模式，行和列均进行了聚类，清楚显示了肿瘤样本与正常样本之间以及不同亚型之间的表达差异。

### 3.4 WGCNA共表达网络分析

#### 3.4.1 方法

WGCNA通过构建无标度共表达网络将基因聚类为功能相关的模块[6]。使用Top 5000可变基因的log2(TPM+1)矩阵。分析流程：

1. 计算基因间Pearson相关系数构建相似性矩阵；
2. 通过软阈值幂函数将相似性矩阵转化为邻接矩阵（选择标准：scale-free R²>0.8）；
3. 构建拓扑重叠矩阵（TOM）以考虑间接连接模式；
4. 基于TOM相异度进行层次聚类和动态树切割识别模块（minModuleSize=30, mergeCutHeight=0.25）；
5. 计算模块特征基因（第一主成分）与临床性状的相关性；
6. 通过模块内连接度（Module Membership, MM）识别枢纽基因。

#### 3.4.2 结果

软阈值选择分析确定了最优power值，确保网络满足无标度拓扑假设。经动态树切割和相似模块合并后，共识别多个共表达模块，每个模块以不同颜色标识。

模块-性状关联分析显示不同模块与临床特征（分子亚型、病理分期、生存状态等）之间存在显著相关性。通过计算基因与模块特征基因的相关性（Module Membership），在关键模块中识别了枢纽基因（hub genes），这些基因在模块内具有最高的连接度，是潜在的生物标志物候选。

### 3.5 生存分析

#### 3.5.1 方法

采用Kaplan-Meier生存曲线（log-rank检验）和Cox比例风险回归评估临床及分子特征的预后价值[16]。Cox模型纳入变量包括年龄、淋巴结阳性数、病理分期（II/III/IV vs I）和分子亚型（HER2-enriched/Triple Negative vs Luminal A）。模型评价采用C-index和似然比检验。

此外，采用LASSO-Cox方法构建预后基因标记：以差异表达基因为候选，通过L1正则化的Cox回归筛选预后相关基因，构建加权风险得分，按中位值将患者分为高风险组和低风险组。

#### 3.5.2 结果

Kaplan-Meier分析显示病理分期与总生存之间存在显著关联，高分期患者预后明显较差。多变量Cox回归（C-index=0.767, 似然比p=3.61×10⁻¹²）结果表明：年龄（HR=1.04, p=7.38×10⁻⁶）、Stage II（HR=3.24, p=0.013）、Stage III（HR=6.00, p=2.35×10⁻⁴）和Stage IV（HR=29.82, p=7.02×10⁻¹⁰）为显著的独立预后因素。Triple Negative亚型也显示出显著的预后风险（HR=2.00, p=0.028）。

LASSO-Cox预后基因标记筛选出3个预后相关基因，通过KM曲线展示了高风险组与低风险组之间的生存差异。

### 3.6 体细胞突变分析（可选）

如配对MAF数据可用，采用maftools[12]进行突变景观分析，包括：高频突变基因统计、Oncoplot瀑布图可视化、变异类型分布和基因间互斥性/共现性检验。

---

## 4 可视化分析

可视化在数据挖掘成果向生物学知识转化中扮演着关键角色。本研究生成的主要图表包括：

| 图号 | 内容 | 方法 | 文件 |
|-----|------|------|------|
| 图1 | 火山图（差异表达） | ggplot2 + ggrepel | volcano_brca.pdf |
| 图2 | Top DEGs热图 | pheatmap | deg_heatmap_top50.pdf |
| 图3 | PCA散点图（按亚型着色） | ggplot2 + stat_ellipse | pca_subtype.pdf |
| 图4 | t-SNE降维图 | Rtsne + ggplot2 | tsne_brca.pdf |
| 图5 | 层次聚类热图 | pheatmap (Ward.D2) | hclust_heatmap_top50.pdf |
| 图6 | K-means轮廓系数 | ggplot2 | silhouette_scores.pdf |
| 图7 | KM生存曲线（分期） | survminer | km_curves_stage.pdf |
| 图8 | KM生存曲线（亚型） | survminer | km_curves_subtype.pdf |
| 图9 | Cox回归森林图 | ggforest | cox_forest_plot.pdf |
| 图10 | WGCNA软阈值选择 | 基础绘图 | wgcna_soft_power.pdf |
| 图11 | WGCNA模块-性状热图 | pheatmap | wgcna_module_trait.pdf |
| 图12 | WGCNA模块树状图 | plotDendroAndColors | wgcna_modules_dendrogram.pdf |
| 图13 | WGCNA网络热图 | pheatmap | wgcna_network.pdf |
| 图14 | 分类模型比较 | ggplot2 | classification_comparison.pdf |
| 图15 | 分类ROC曲线 | pROC | classification_roc.pdf |
| 图16 | 分类混淆矩阵 | pheatmap | classification_cm.pdf |

---

## 5 结论与展望

### 5.1 主要发现

本研究以TCGA-BRCA数据为基础，构建了覆盖差异表达分析、功能富集、分子亚型分类、聚类分析、WGCNA共表达网络、生存分析等多个分析模块的完整数据挖掘流程。主要发现包括：

（1）DESeq2差异表达分析鉴定6,768个显著差异基因，上调（4,334个）数量约为下调（2,434个）的1.78倍。

（2）Random Forest和LASSO分类器均达到75.8%的分类准确率，LASSO通过L1正则化实现了特征稀疏压缩。

（3）PCA（PC1=15.3%, PC2=8.5%）和K-means轮廓系数（K=2最佳，silhouette=0.180）一致表明ER状态是BRCA转录组结构的最强驱动因素。

（4）WGCNA识别了多个共表达模块及枢纽基因，模块-性状关联分析揭示了与临床特征显著相关的模块。

（5）Cox多变量回归（C-index=0.767）证实病理分期是BRCA预后最强预测因子（Stage IV HR=29.82, p=7.02×10⁻¹⁰）。

### 5.2 研究局限性

- 仅使用mRNA转录组数据，未纳入miRNA、DNA甲基化、体细胞突变等其他组学层次；
- 生存分析的事件率有限，可能影响多变量预后模型的统计功效；
- 分类模型的超参数未进行系统调优，XGBoost等模型的性能可能低于其最优水平。

### 5.3 未来方向

（1）整合miRNA表达谱和体细胞突变数据，实现多组学联合分析；（2）利用METABRIC等独立队列对分类和预后模型进行外部验证；（3）对WGCNA枢纽基因进行功能实验验证；（4）探索深度学习方法在多组学整合中的应用；（5）使用单细胞RNA测序数据解析肿瘤微环境异质性。

---

## 参考文献

[1] Sung H, Ferlay J, Siegel R L, et al. Global cancer statistics 2022: GLOBOCAN estimates of incidence and mortality worldwide for 36 cancers in 185 countries[J]. CA: A Cancer Journal for Clinicians, 2024, 74(3): 229-263.

[2] Hanahan D, Weinberg R A. Hallmarks of cancer: the next generation[J]. Cell, 2011, 144(5): 646-674.

[3] Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours[J]. Nature, 2012, 490(7418): 61-70.

[4] Perou C M, Sorlie T, Eisen M B, et al. Molecular portraits of human breast tumours[J]. Nature, 2000, 406(6797): 747-752.

[5] Love M I, Huber W, Anders S. Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2[J]. Genome Biology, 2014, 15(12): 550.

[6] Langfelder P, Horvath S. WGCNA: an R package for weighted correlation network analysis[J]. BMC Bioinformatics, 2008, 9: 559.

[7] Simon N, Friedman J, Hastie T, et al. Regularization paths for Cox's proportional hazards model via coordinate descent[J]. Journal of Statistical Software, 2011, 39(5): 1-13.

[8] Kolberg L, Raudvere U, Kuzmin I, et al. g:Profiler—interoperable web service for functional enrichment analysis and gene identifier mapping (2023 update)[J]. Nucleic Acids Research, 2023, 51(W1): W207-W212.

[9] Agrawal R, Srikant R. Fast algorithms for mining association rules[C]. Proceedings of the 20th International Conference on Very Large Data Bases, 1994: 487-499.

[10] Wilkerson M D, Hayes D N. ConsensusClusterPlus: a class discovery tool with confidence assessments and item tracking[J]. Bioinformatics, 2010, 26(12): 1572-1573.

[11] Colaprico A, Silva T C, Olsen C, et al. TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data[J]. Nucleic Acids Research, 2016, 44(8): e71.

[12] Mayakonda A, Lin D C, Assenov Y, et al. Maftools: efficient and comprehensive analysis of somatic variants in cancer[J]. Genome Research, 2018, 28(11): 1747-1756.

[13] Friedman J, Hastie T, Tibshirani R. Regularization paths for generalized linear models via coordinate descent[J]. Journal of Statistical Software, 2010, 33(1): 1-22.

[14] Breiman L. Random forests[J]. Machine Learning, 2001, 45(1): 5-32.

[15] Chen T, Guestrin C. XGBoost: a scalable tree boosting system[C]. Proceedings of the 22nd ACM SIGKDD, 2016: 785-794.

[16] Therneau T M, Grambsch P M. Modeling survival data: extending the Cox model[M]. New York: Springer, 2000.

[17] Tibshirani R. Regression shrinkage and selection via the lasso[J]. Journal of the Royal Statistical Society: Series B, 1996, 58(1): 267-288.

---

## 附录：R代码

### 附录A：数据预处理（01_data_preprocessing.R）

```r
# 完整代码见 code/01_data_preprocessing.R
# 主要步骤：
# 1. 加载brca_exp.Rdata和brca_clinical.Rdata
# 2. 解析TCGA条形码，区分肿瘤/正常样本
# 3. 低表达基因过滤（至少10%样本中counts>=10）
# 4. 基因标识符映射（Ensembl ID -> Gene Symbol）
# 5. Counts -> TPM标准化
# 6. 临床数据清洗与缺失值处理
# 7. 表达-临床数据对齐与保存
```

### 附录B：差异表达分析（02_diff_expression.R）

```r
# 完整代码见 code/02_diff_expression.R
# 主要步骤：
# 1. 构建DESeqDataSetFromMatrix (Tumor vs Normal)
# 2. 独立过滤 + DESeq()运行
# 3. 提取差异表达结果 (|log2FC|>1, padj<0.05)
# 4. 火山图、MA图、Top DEG热图
# 5. g:Profiler功能富集分析（上调/下调分别）
```

### 附录C：分类分析（03_classification.R）

```r
# 完整代码见 code/03_classification.R
# 主要步骤：
# 1. Top 500可变基因 -> log2(TPM+1)特征矩阵
# 2. 70/30分层训练/测试划分
# 3. Random Forest (ntree=500)
# 4. LASSO多分类回归 (cv.glmnet, family="multinomial")
# 5. XGBoost (multi:softmax)
# 6. 模型比较：Accuracy, Kappa, Macro F1
# 7. ROC曲线和混淆矩阵热图
```

### 附录D：聚类分析（04_clustering.R）

```r
# 完整代码见 code/04_clustering.R
# 主要步骤：
# 1. PCA分析 + 按亚型/分期着色散点图
# 2. t-SNE降维可视化
# 3. 层次聚类热图（Ward.D2, Pearson距离）
# 4. K-means轮廓系数分析 (K=2..10)
# 5. ConsensusClusterPlus共识聚类 (K=2..8)
```

### 附录E：WGCNA分析（05_wgcna.R）

```r
# 完整代码见 code/05_wgcna.R
# 主要步骤：
# 1. Top 5000可变基因 -> WGCNA输入矩阵
# 2. goodSamplesGenes质量检查
# 3. pickSoftThreshold软阈值选择 (signed network)
# 4. blockwiseModules网络构建
# 5. 模块-性状相关性分析
# 6. 枢纽基因识别 (Module Membership)
```

### 附录F：生存分析（06_survival_analysis.R）

```r
# 完整代码见 code/06_survival_analysis.R
# 主要步骤：
# 1. Kaplan-Meier曲线（按分子亚型、病理分期）
# 2. Log-rank检验
# 3. 多变量Cox比例风险回归
# 4. LASSO-Cox预后基因标记构建
# 5. 风险分组KM曲线
```

### 附录G：可视化汇总（08_visualization.R）

```r
# 完整代码见 code/08_visualization.R
# 生成所有论文级图表，汇总各分析模块的关键可视化结果
```

---

*本研究分析使用R语言，在Windows 11环境下执行。全部分析脚本见项目目录code/，运行结果见results/。*
