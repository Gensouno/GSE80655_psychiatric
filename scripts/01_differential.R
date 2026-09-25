# GSE80655: data processing, PCA, DEG, cross-disorder, enrichment
# three brain regions combined, region as covariate

library(GEOquery)
library(limma)
library(ggplot2)
library(pheatmap)
library(org.Hs.eg.db)
library(msigdbr)
library(patchwork)

setwd("/Users/jschen/Downloads/GSE80655_psychiatric")
dir.create("output", showWarnings = FALSE)

# 1. read raw counts
cat("reading counts...\n")
raw <- read.delim("data/GSE80655_GeneExpressionData_Updated_3-26-2018.txt",
                  check.names = FALSE)
counts <- as.matrix(raw[, -1])
rownames(counts) <- raw$gene_id
cat("raw matrix:", nrow(counts), "genes x", ncol(counts), "samples\n")

# 2. parse phenotype from title
gse <- getGEO("GSE80655", destdir = "data/", GSEMatrix = TRUE)
pd <- pData(gse[[1]])
pd$slid <- sub(".*_(SL[0-9]+)$", "\\1", pd$title)
parts <- strsplit(pd$title, "_")
pd$region <- sapply(parts, `[`, 2)
pd$dcode  <- sapply(parts, `[`, 3)
pd$diagnosis <- c(C = "Control", M = "MDD", B = "BD", S = "SCZ")[pd$dcode]

pd <- pd[pd$slid %in% colnames(counts), ]
counts <- counts[, pd$slid]
cat("all samples:\n"); print(table(pd$region, pd$diagnosis))

# 3. Ensembl ID to gene symbol
sym <- mapIds(org.Hs.eg.db, keys = rownames(counts), column = "SYMBOL",
              keytype = "ENSEMBL", multiVals = "first")
has <- !is.na(sym)
counts <- counts[has, ]; sym <- sym[has]
tot <- rowSums(counts)
ord <- order(sym, -tot)
counts <- counts[ord, ]; sym <- sym[ord]
keep <- !duplicated(sym)
counts <- counts[keep, ]
rownames(counts) <- sym[keep]
cat("after id conversion:", nrow(counts), "genes\n")

# 4. filter low-expression genes
keep_g <- rowSums(counts >= 10) >= 20
counts_f <- counts[keep_g, ]
cat("after low-expression filter:", nrow(counts_f), "genes\n")

region <- factor(pd$region[match(colnames(counts_f), pd$slid)])
group  <- factor(pd$diagnosis[match(colnames(counts_f), pd$slid)],
                 levels = c("Control", "BD", "MDD", "SCZ"))

write.csv(counts_f, "data/expr_counts.csv")
write.csv(data.frame(sample = colnames(counts_f), region = region,
                    diagnosis = group), "data/pheno.csv", row.names = FALSE)

libsize <- colSums(counts_f)
logcpm <- log2(t(t(counts_f) / libsize * 1e6) + 1)

# 5. PCA
do_pca <- function(mat) {
  pca <- prcomp(t(mat), scale. = TRUE)
  pvar <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 1)
  list(df = data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2],
                       region = region, group = group), pvar = pvar)
}
p1 <- do_pca(logcpm)
pa <- ggplot(p1$df, aes(PC1, PC2, color = region)) +
  geom_point(size = 2) +
  labs(x = paste0("PC1 (", p1$pvar[1], "%)"),
       y = paste0("PC2 (", p1$pvar[2], "%)"), title = "colored by brain region",
       color = NULL) + theme_bw() + scale_color_brewer(palette = "Set2")

# remove region effect, then PCA by diagnosis
logcpm_nr <- removeBatchEffect(logcpm, batch = region)
p2 <- do_pca(logcpm_nr)
pb <- ggplot(p2$df, aes(PC1, PC2, color = group)) +
  geom_point(size = 2) +
  stat_ellipse(level = 0.95, show.legend = FALSE) +
  labs(x = paste0("PC1 (", p2$pvar[1], "%)"),
       y = paste0("PC2 (", p2$pvar[2], "%)"), title = "region effect removed",
       color = NULL) + theme_bw() +
  scale_color_manual(values = c(Control = "gray50", BD = "#d95f02",
                                MDD = "#1b9e77", SCZ = "#7570b3"))
ggsave("output/pca_plot.pdf", pa / pb + plot_annotation(title = "PCA"),
       width = 6, height = 9)
cat("PCA saved\n")

# 6. DEG: each disorder vs Control, region as covariate
run_deg <- function(disease) {
  sel <- group %in% c("Control", disease)
  g <- factor(as.character(group[sel]), levels = c("Control", disease))
  rg <- region[sel]
  des <- model.matrix(~ rg + g)
  vv <- voom(counts_f[, sel], des, plot = FALSE)
  fit <- lmFit(vv, des); fit <- eBayes(fit)
  tt <- topTable(fit, coef = ncol(des), number = Inf)
  tt$gene <- rownames(tt)
  tt
}
deg_BD  <- run_deg("BD")
deg_SCZ <- run_deg("SCZ")
deg_MDD <- run_deg("MDD")
write.csv(deg_BD,  "output/deg_result_BD.csv",  row.names = FALSE)
write.csv(deg_SCZ, "output/deg_result_SCZ.csv", row.names = FALSE)
write.csv(deg_MDD, "output/deg_result_MDD.csv", row.names = FALSE)

sig_list <- list()
for (nm in c("BD", "SCZ", "MDD")) {
  d <- get(paste0("deg_", nm))
  sig_list[[nm]] <- d$gene[d$adj.P.Val < 0.05 & abs(d$logFC) > 0.5]
  cat(nm, "significant genes:", length(sig_list[[nm]]), "\n")
}

plot_volcano <- function(d, ttl) {
  d$status <- ifelse(d$adj.P.Val < 0.05 & d$logFC > 0.5, "up",
              ifelse(d$adj.P.Val < 0.05 & d$logFC < -0.5, "down", "ns"))
  d$label <- ifelse(d$gene %in% d$gene[order(d$adj.P.Val)][1:10], d$gene, "")
  ggplot(d, aes(logFC, -log10(adj.P.Val))) +
    geom_point(aes(color = status), size = 1, alpha = 0.7) +
    geom_text(aes(label = label), size = 2.5, vjust = 1.5, check_overlap = TRUE) +
    scale_color_manual(values = c(up = "#d95f02", down = "#1b9e77", ns = "gray80")) +
    labs(title = ttl, x = "log2 fold change", y = "-log10(FDR)") +
    theme_bw() + theme(legend.position = "none")
}
ggsave("output/volcano_BD.pdf",  plot_volcano(deg_BD,  "BD vs Control"),  width = 6, height = 5)
ggsave("output/volcano_SCZ.pdf", plot_volcano(deg_SCZ, "SCZ vs Control"), width = 6, height = 5)
ggsave("output/volcano_MDD.pdf", plot_volcano(deg_MDD, "MDD vs Control"), width = 6, height = 5)

# BD top50 heatmap (BD + Control samples)
sel_bc <- group %in% c("Control", "BD")
top50 <- deg_BD$gene[order(deg_BD$adj.P.Val)][1:50]
ann <- data.frame(region = region[sel_bc], diagnosis = group[sel_bc])
rownames(ann) <- colnames(counts_f)[sel_bc]
pdf("output/deg_heatmap_BD.pdf", width = 9, height = 9)
pheatmap(logcpm[top50, sel_bc],
         annotation_col = ann,
         show_colnames = FALSE, fontsize_row = 7,
         color = colorRampPalette(c("#3b4cc0", "white", "#d95f02"))(100),
         main = "Top 50 BD DEGs (BD vs Control)")
dev.off()
cat("heatmap saved\n")

# 7. cross-disorder comparison
pdf("output/venn_diseases.pdf", width = 6, height = 6)
vc <- vennCounts(data.frame(BD  = rownames(counts_f) %in% sig_list$BD,
                           SCZ = rownames(counts_f) %in% sig_list$SCZ,
                           MDD = rownames(counts_f) %in% sig_list$MDD))
vennDiagram(vc, circle.col = c("#d95f02", "#7570b3", "#1b9e77"),
            main = "Significant DEGs across disorders")
dev.off()

bd_specific <- setdiff(sig_list$BD, union(sig_list$SCZ, sig_list$MDD))
write.csv(data.frame(gene = bd_specific), "output/bd_specific_genes.csv",
          row.names = FALSE)
cat("BD-specific genes:", length(bd_specific), "\n")

lfc <- data.frame(BD = deg_BD$logFC[match(rownames(counts_f), deg_BD$gene)],
                  SCZ = deg_SCZ$logFC[match(rownames(counts_f), deg_SCZ$gene)],
                  MDD = deg_MDD$logFC[match(rownames(counts_f), deg_MDD$gene)])
cmat <- cor(lfc, use = "pairwise.complete.obs")
cordf <- expand.grid(x = colnames(cmat), y = colnames(cmat))
cordf$r <- as.vector(cmat)
p <- ggplot(cordf, aes(x, y, fill = r)) +
  geom_tile(color = "white") + geom_text(aes(label = round(r, 2))) +
  scale_fill_gradient2(low = "#3b4cc0", mid = "white", high = "#d95f02",
                       limits = c(0, 1)) +
  labs(title = "logFC correlation across disorders", x = NULL, y = NULL) +
  theme_bw()
ggsave("output/logfc_correlation_heatmap.pdf", p, width = 5, height = 4.5)

# 8. Hallmark enrichment
hall <- msigdbr(species = "Homo sapiens", category = "H")
universe <- rownames(counts_f); N <- length(universe)
enrich <- function(query) {
  query <- intersect(query, universe)
  sets <- split(hall$gene_symbol, hall$gs_name)
  res <- do.call(rbind, lapply(names(sets), function(s) {
    gs <- intersect(sets[[s]], universe)
    k <- length(intersect(query, gs))
    pval <- phyper(k - 1, length(gs), N - length(gs), length(query),
                   lower.tail = FALSE)
    data.frame(pathway = s, overlap = k, set_size = length(gs),
               query_size = length(query), pval = pval)
  }))
  res$FDR <- p.adjust(res$pval, method = "BH")
  res[order(res$pval), ]
}
up_bd <- deg_BD$gene[deg_BD$adj.P.Val < 0.05 & deg_BD$logFC > 0.5]
dn_bd <- deg_BD$gene[deg_BD$adj.P.Val < 0.05 & deg_BD$logFC < -0.5]
en_up <- enrich(up_bd); en_dn <- enrich(dn_bd)
write.csv(rbind(cbind(direction = "up", en_up), cbind(direction = "down", en_dn)),
          "output/enrichment_result.csv", row.names = FALSE)

plot_bubble <- function(en, ttl) {
  d <- head(en, 15)
  d$pathway <- factor(sub("HALLMARK_", "", d$pathway),
                      levels = rev(sub("HALLMARK_", "", d$pathway)))
  ggplot(d, aes(overlap / set_size, pathway)) +
    geom_point(aes(size = overlap, color = -log10(pval))) +
    scale_color_gradient(low = "#f0c9a8", high = "#b32929") +
    labs(title = ttl, x = "enrichment ratio", y = NULL,
         size = "genes", color = "-log10(P)") +
    theme_bw() + theme(axis.text.y = element_text(size = 8))
}
ggsave("output/enrichment_bubble_up.pdf",
       plot_bubble(en_up, "BD up-regulated Hallmark"), width = 7, height = 6)
ggsave("output/enrichment_bubble_down.pdf",
       plot_bubble(en_dn, "BD down-regulated Hallmark"), width = 7, height = 6)

cat("\n=== 01 differential analysis done ===\n")
