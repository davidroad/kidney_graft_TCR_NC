#!/usr/bin/env Rscript
options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({ library(Seurat); library(limma) })

root <- "/home/26_immune_NC"
obj_file <- "/home/0.collaboration/20_sc_VDJ/merge_4_with_vdj/merged_4_with_vdj_T_cells_filter_10_28_2025.RData"
out <- file.path(root, "figure_revision_2026_07_08/outputs/current_mainline_recompute")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

e <- new.env(parent = emptyenv())
load(obj_file, envir = e)
obj_name <- intersect(ls(e), c("merged_all_T_batchcorrected_filter"))[1]
if (is.na(obj_name)) stop("filtered T-cell object not found")
obj <- get(obj_name, e)

# This is the original Fig. 2B definition: proportions of the 11 retained
# annotated T-cell clusters within each graft/PBMC sample.
cluster_levels <- c("Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+",
  "CXCR6+ effector CD8+", "CCR2+ CD4+", "Naive/Memory-like CD8+",
  "Naive CD4+", "Treg", "CX3CR1+ effector CD8+",
  "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+")
md <- obj@meta.data
if (!all(c("batch", "label") %in% names(md))) stop("batch/label metadata missing")
cl <- as.character(Idents(obj))
md$cluster <- factor(cl, levels = cluster_levels)
md$patient <- sub("^(graft|PBMC)_", "", as.character(md$batch))
md$condition <- ifelse(grepl("^graft", md$batch, ignore.case = TRUE), "Graft", "PBMC")
md <- md[!is.na(md$cluster) & md$condition %in% c("Graft", "PBMC"), , drop = FALSE]

counts <- as.data.frame(table(md$patient, md$condition, md$cluster), stringsAsFactors = FALSE)
names(counts) <- c("patient", "condition", "tcluster", "n")
denom <- aggregate(n ~ patient + condition, counts, sum)
names(denom)[3] <- "denom"
tab <- merge(counts, denom, by = c("patient", "condition"), all.x = TRUE)
tab$proportion_pct <- 100 * tab$n / tab$denom
tab <- tab[order(tab$patient, tab$condition, match(tab$tcluster, cluster_levels)), ]
write.csv(tab, file.path(out, "fig2b_latest_11cluster_sample_level_proportions.csv"), row.names = FALSE)
write.csv(tab, file.path(out, "fig2b_current_sample_level_cluster_proportions.csv"), row.names = FALSE)

# Donor-blocked empirical-Bayes model on logit proportions.
tab$patient <- factor(tab$patient)
tab$condition <- factor(tab$condition, levels = c("PBMC", "Graft"))
tab$y <- qlogis(pmin(pmax(tab$proportion_pct / 100, 0.005), 0.995))
res <- lapply(cluster_levels, function(k) {
  z <- tab[tab$tcluster == k, , drop = FALSE]
  design <- model.matrix(~ patient + condition, data = z)
  fit <- eBayes(lmFit(matrix(z$y, nrow = 1), design), trend = TRUE, robust = TRUE)
  tt <- topTable(fit, coef = "conditionGraft", number = 1, sort.by = "none")
  data.frame(analysis = "Fig2B_latest_11cluster_frequency", feature = k,
    n_pairs = length(unique(z$patient)),
    logit_effect_graft_vs_pbmc = tt$logFC[1], moderated_t = tt$t[1],
    p_value = tt$P.Value[1], FDR = tt$adj.P.Val[1],
    model = "limma lmFit + eBayes; design ~ patient + condition; logit proportion",
    stringsAsFactors = FALSE)
})
stats <- do.call(rbind, res)
write.csv(stats, file.path(out, "fig2b_limma_ebayes_logit_proportions.csv"), row.names = FALSE)
writeLines(c("Source: merged_4_with_vdj_T_cells_filter_10_28_2025.RData",
  "Object: merged_all_T_batchcorrected_filter",
  "Definition: 11 retained annotated T-cell clusters; cluster 8, 11, and 12 were removed by the original script.",
  "Model: limma lmFit/eBayes on donor-level logit-transformed proportions; design ~ patient + condition."),
  file.path(out, "fig2b_latest_provenance.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "fig2b_latest_sessionInfo.txt"))
cat("Wrote latest 11-cluster Fig. 2B tables.\n")
