# GSE80655 精神疾病转录组分析

基于 GEO 数据集 GSE80655（Ramaker et al., 2017, Genome Biology）的死后脑组织 RNA-seq 数据，对双相情感障碍（BD）、精神分裂症（SCZ）、重度抑郁症（MDD）和正常对照做转录组分析。

## 数据

- 281 样本 / 96 个体，三个脑区（AnCg、nAcc、DLPFC）
- raw counts，过滤后保留 20274 个基因
- 脑区作为协变量

## 分析流程

1. QC 与 PCA
2. 差异表达（limma + voom；FDR < 0.05，|logFC| > 0.5）
3. 跨疾病比较（韦恩图、logFC 相关性）
4. Hallmark 富集（超几何检验）
5. PPI 网络（STRING，top 200 BD 差异基因）
6. 诊断模型（LASSO / 弹性网 / 随机森林；按个体划分训练测试集）

## 主要结果

| | BD | SCZ | MDD |
|---|---|---|---|
| 显著差异基因 | 582 | 480 | 29 |

- BD 与 SCZ 的 logFC 相关性 r = 0.86，共享 252 个差异基因
- BD 上调：IL6-JAK-STAT3、炎症反应、干扰素、补体
- BD 下调：突触信号、氧化磷酸化
- PPI 网络 189 条边，hub 基因为 HSPA8（degree 20）
- 最优模型：弹性网，个体级 AUC = 0.938

## 运行

```bash
Rscript scripts/01_differential.R
Rscript scripts/02_modeling.R
```

## 依赖

GEOquery、limma、ggplot2、pheatmap、org.Hs.eg.db、msigdbr、igraph、glmnet、pROC、randomForest、patchwork
