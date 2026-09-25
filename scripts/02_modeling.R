# GSE80655: PPI network, LASSO, classifier comparison
# BD vs Control, individual-level train/test split

library(limma)
library(ggplot2)
library(igraph)
library(glmnet)
library(pROC)
library(randomForest)
library(GEOquery)

setwd("/Users/jschen/Downloads/GSE80655_psychiatric")

counts_f <- as.matrix(read.csv("data/expr_counts.csv", row.names = 1,
                               check.names = FALSE))
pheno <- read.csv("data/pheno.csv")
deg_BD <- read.csv("output/deg_result_BD.csv")

# individual ID from series matrix title (X prefix)
gse <- getGEO("GSE80655", destdir = "data/", GSEMatrix = TRUE)
pd <- pData(gse[[1]])
pd$slid <- sub(".*_(SL[0-9]+)$", "\\1", pd$title)
pd$indiv <- sub("^(X[0-9]+)_.*$", "\\1", pd$title)
indiv <- pd$indiv[match(colnames(counts_f), pd$slid)]

keep_s <- pheno$diagnosis %in% c("Control", "BD")
counts_bc <- counts_f[, keep_s]
grp <- factor(pheno$diagnosis[keep_s], levels = c("Control", "BD"))
indiv_bc <- indiv[keep_s]
cat("BD vs Control:", ncol(counts_bc), "samples,",
    length(unique(indiv_bc)), "individuals\n")

# voom with region covariate, then remove batch effect
rg <- factor(pheno$region[keep_s])
des <- model.matrix(~ rg + grp)
v <- voom(counts_bc, des, plot = FALSE)
expr_clean <- t(removeBatchEffect(v$E, batch = rg))  # samples x genes
y <- as.numeric(grp) - 1   # Control=0, BD=1

# 1. PPI network (top 200 BD DEGs)
cat("building PPI...\n")
deg_sorted <- deg_BD[order(deg_BD$adj.P.Val), ]
top200 <- head(deg_sorted$gene, 200)
identifiers <- URLencode(paste(top200, collapse = "%0d"))
url <- paste0("https://string-db.org/api/tsv/network?identifiers=", identifiers,
              "&species=9606&required_score=400&network_type=functional")
ppi <- try(read.delim(url), silent = TRUE)
if (inherits(ppi, "try-error") || nrow(ppi) == 0) {
  cat("STRING fetch failed, skipping PPI\n")
} else {
  write.table(ppi, "output/ppi_string.tsv", sep = "\t", row.names = FALSE)
  edges <- unique(ppi[, c("preferredName_A", "preferredName_B")])
  g <- graph_from_data_frame(edges, directed = FALSE)
  g <- simplify(g)
  dg <- sort(degree(g), decreasing = TRUE)
  write.csv(data.frame(Gene = names(dg), Degree = as.numeric(dg)),
            "output/hub_genes.csv", row.names = FALSE)
  cat("PPI edges:", ecount(g), " top hub:", names(dg)[1], dg[1], "\n")

  upgenes <- deg_sorted$gene[deg_sorted$logFC > 0]
  V(g)$size <- degree(g) * 1.5 + 4
  V(g)$frame.color <- "gray30"
  set.seed(1)
  pdf("output/ppi_network.pdf", width = 10, height = 9)
  l <- layout_with_fr(g)
  plot(g, layout = l, vertex.label = ifelse(degree(g) >= 4, V(g)$name, NA),
       vertex.label.cex = 0.8, vertex.label.color = "black",
       vertex.color = ifelse(names(V(g)) %in% upgenes, "#f4a261", "#76c893"),
       edge.color = "gray80", main = "PPI network of top BD DEGs")
  dev.off()
}

# 2. individual-level train/test split (70/30, stratified)
indiv_label <- tapply(y, indiv_bc, function(z) z[1])
all_ind <- names(indiv_label)
set.seed(42)
train_ind <- c(sample(all_ind[indiv_label == 0], floor(sum(indiv_label == 0) * 0.7)),
               sample(all_ind[indiv_label == 1], floor(sum(indiv_label == 1) * 0.7)))
test_ind <- setdiff(all_ind, train_ind)
train_idx <- which(indiv_bc %in% train_ind)
test_idx  <- which(indiv_bc %in% test_ind)
cat("train individuals:", length(train_ind), " test individuals:", length(test_ind), "\n")

# candidate features: FDR < 0.05, |logFC| > 0.25
feat <- deg_sorted$gene[deg_sorted$adj.P.Val < 0.05 & abs(deg_sorted$logFC) > 0.25]
feat <- intersect(feat, colnames(expr_clean))
cat("candidate features:", length(feat), "\n")

Xtr <- expr_clean[train_idx, feat, drop = FALSE]
Xte <- expr_clean[test_idx,  feat, drop = FALSE]
ytr <- y[train_idx]; yte <- y[test_idx]

# 3. LASSO
set.seed(42)
cv_lasso <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 1, nfolds = 10)
lam_l <- cv_lasso$lambda.min
co <- as.matrix(coef(cv_lasso, s = lam_l))
nz <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
cat("LASSO non-zero genes:", length(nz), "\n")

pdf("output/lasso_coef_path.pdf", width = 6, height = 5)
plot(glmnet(Xtr, ytr, family = "binomial", alpha = 1), xvar = "lambda",
     main = "LASSO coefficient path")
abline(v = log(lam_l), lty = 2, col = "gray50"); dev.off()
pdf("output/lasso_cv.pdf", width = 6, height = 5)
plot(cv_lasso, main = "LASSO cross-validation"); dev.off()
write.csv(data.frame(gene = nz, coef = co[nz, 1]),
          "output/lasso_nonzero_coef.csv", row.names = FALSE)

prob_l <- as.numeric(predict(cv_lasso, newx = Xte, s = lam_l, type = "response"))

# 4. Elastic Net
set.seed(42)
cv_en <- cv.glmnet(Xtr, ytr, family = "binomial", alpha = 0.5, nfolds = 10)
prob_e <- as.numeric(predict(cv_en, newx = Xte, s = cv_en$lambda.min,
                             type = "response"))

# 5. Random Forest
set.seed(42)
rf <- randomForest(x = Xtr, y = factor(ytr), ntree = 500)
prob_r <- predict(rf, newdata = Xte, type = "prob")[, 2]

# 6. sample-level ROC comparison
roc_l <- roc(yte, prob_l, quiet = TRUE)
roc_e <- roc(yte, prob_e, quiet = TRUE)
roc_r <- roc(yte, prob_r, quiet = TRUE)
auc_l <- as.numeric(auc(roc_l)); auc_e <- as.numeric(auc(roc_e)); auc_r <- as.numeric(auc(roc_r))
cat("sample-level test AUC - LASSO:", round(auc_l,3),
    " Elastic Net:", round(auc_e,3), " RF:", round(auc_r,3), "\n")

pdf("output/classifier_comparison_roc.pdf", width = 6, height = 6)
plot(roc_l, col = "#d95f02", lwd = 2, main = "Test set ROC (individual split)")
plot(roc_e, col = "#1b9e77", lwd = 2, add = TRUE)
plot(roc_r, col = "#7570b3", lwd = 2, add = TRUE)
legend("bottomright",
       legend = c(paste0("LASSO (", round(auc_l,3), ")"),
                  paste0("Elastic Net (", round(auc_e,3), ")"),
                  paste0("Random Forest (", round(auc_r,3), ")")),
       col = c("#d95f02", "#1b9e77", "#7570b3"), lwd = 2)
dev.off()

# 7. individual-level aggregated AUC
test_indiv_id <- indiv_bc[test_idx]
agg <- function(prob) tapply(prob, test_indiv_id, mean)
indiv_true <- tapply(yte, test_indiv_id, function(z) z[1])
auc_il <- as.numeric(auc(roc(indiv_true, agg(prob_l), quiet = TRUE)))
auc_ie <- as.numeric(auc(roc(indiv_true, agg(prob_e), quiet = TRUE)))
auc_ir <- as.numeric(auc(roc(indiv_true, agg(prob_r), quiet = TRUE)))
cat("individual-level test AUC - LASSO:", round(auc_il,3),
    " Elastic Net:", round(auc_ie,3), " RF:", round(auc_ir,3), "\n")

auc_tab <- data.frame(
  Classifier = c("LASSO", "Elastic Net", "Random Forest"),
  Sample_AUC = round(c(auc_l, auc_e, auc_r), 3),
  Individual_AUC = round(c(auc_il, auc_ie, auc_ir), 3))
write.csv(auc_tab, "output/classifier_auc_comparison.csv", row.names = FALSE)
print(auc_tab)

pred_df <- data.frame(sample = rownames(Xte), individual = test_indiv_id,
                      true = ifelse(yte == 1, "BD", "Control"),
                      prob_BD = round(prob_l, 3),
                      predicted = ifelse(prob_l > 0.5, "BD", "Control"))
write.csv(pred_df, "output/prediction_test.csv", row.names = FALSE)

cat("\n=== 02 modeling done ===\n")
