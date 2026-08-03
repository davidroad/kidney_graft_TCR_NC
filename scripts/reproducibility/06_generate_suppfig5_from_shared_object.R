#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4L) {
  stop("Usage: Rscript 06_generate_suppfig5_from_shared_object.R <shared-clone.RData> <CD8-T-cell.RData> <panel-e-cell-selection.csv.gz> <output-dir>")
}
input_file <- normalizePath(args[[1]], mustWork = TRUE)
cd8_input_file <- normalizePath(args[[2]], mustWork = TRUE)
selection_file <- normalizePath(args[[3]], mustWork = TRUE)
out_dir <- args[[4]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(input_file)
if (!exists("all_clone_test_filter") || !inherits(all_clone_test_filter, "Seurat")) {
  stop("Input must contain the Seurat object 'all_clone_test_filter'.")
}
obj <- all_clone_test_filter
DefaultAssay(obj) <- "RNA"
obj$label <- factor(as.character(obj$label), levels = c("Graft", "PBMC"))

panel_genes <- list(
  a = c("KLF2", "IL7R", "SELL"),
  b = c("TOX", "PDCD1", "HAVCR2"),
  c = c("IFNG", "PRF1", "ID2"),
  d = c("MKI67", "CDKN2A")
)
needed <- unique(c(unlist(panel_genes), "TCF7"))
missing <- setdiff(needed, rownames(obj[["RNA"]]))
if (length(missing)) stop("Missing Supplementary Fig. 5 genes: ", paste(missing, collapse = ", "))

make_feature_group <- function(genes, tag) {
  plots <- lapply(genes, function(gene) {
    FeaturePlot(
      obj, features = gene, reduction = "umap", split.by = "label",
      cols = c("lightgrey", "red"), keep.scale = "feature", order = TRUE,
      pt.size = 0.16, ncol = 2, raster = FALSE
    ) & theme_classic(base_size = 8) &
      theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size = 9))
  })
  p <- wrap_plots(plots, nrow = 1)
  ggsave(file.path(out_dir, paste0("suppfig5", tag, "_shared_clone_featureplots.pdf")), p,
         width = if (length(genes) == 3) 15 else 10, height = 4.2, units = "in", device = cairo_pdf)
}
invisible(mapply(make_feature_group, panel_genes, names(panel_genes), SIMPLIFY = FALSE))

cd8_env <- new.env(parent = emptyenv())
load(cd8_input_file, envir = cd8_env)
if (!exists("merged_all_T_batchcorrected_filter_CD8", envir = cd8_env, inherits = FALSE)) {
  stop("CD8 input must contain 'merged_all_T_batchcorrected_filter_CD8'.")
}
cd8_cells <- intersect(colnames(obj), colnames(cd8_env$merged_all_T_batchcorrected_filter_CD8))
if (!length(cd8_cells)) stop("No shared-clonotype cells were found in the CD8 object.")
rm(cd8_env)
gc()
selected_cells <- read.csv(gzfile(selection_file), stringsAsFactors = FALSE)$cell_id
if (!length(selected_cells)) stop("The panel-e cell selection is empty.")
if (any(!selected_cells %in% cd8_cells)) stop("The panel-e selection contains cells outside the CD8 shared-clonotype intersection.")
cd8_cells <- selected_cells

stats_gene_set <- c(
  "LDHA", "MKI67", "TOP2A", "CDK1", "CCNB1", "CDC20", "PLK1", "AURKB", "MCM5", "PCNA",
  "CDKN1A", "CDKN1B", "CDKN2A", "RB1", "TP53", "BCL2", "MCL1", "BAX", "BAK1", "BCL2L11",
  "BAD", "BID", "CASP3", "TNFRSF1A"
)
stats_gene_set <- intersect(stats_gene_set, rownames(obj[["RNA"]]))
expr <- FetchData(obj, vars = unique(c("TCF7", stats_gene_set)))[cd8_cells, , drop = FALSE]
expr$group <- ifelse(expr$TCF7 > 0, "TCF7hi", "TCF7-")
p_raw <- vapply(stats_gene_set, function(gene) {
  suppressWarnings(wilcox.test(expr[[gene]] ~ expr$group, exact = FALSE)$p.value)
}, numeric(1))
p_fdr <- p.adjust(p_raw, method = "BH")
stats <- data.frame(gene = stats_gene_set, wilcoxon_p = p_raw, wilcoxon_fdr = p_fdr)
write.csv(stats, file.path(out_dir, "suppfig5e_cell_level_wilcoxon_statistics.csv"), row.names = FALSE)

vdat <- rbind(
  data.frame(group = expr$group, gene = "MKI67", expression = expr$MKI67),
  data.frame(group = expr$group, gene = "CDKN2A", expression = expr$CDKN2A)
)
vdat$group <- factor(vdat$group, levels = c("TCF7hi", "TCF7-"))
vdat$gene <- factor(vdat$gene, levels = c("MKI67", "CDKN2A"))
ann <- stats[match(levels(vdat$gene), stats$gene), ]
ann$gene <- factor(ann$gene, levels = levels(vdat$gene))
ann$label <- ifelse(ann$wilcoxon_fdr < 1e-4, "<0.0001", sprintf("%.4f", ann$wilcoxon_fdr))
ann$y <- vapply(levels(vdat$gene), function(g) max(vdat$expression[vdat$gene == g], na.rm = TRUE) * 1.08 + 0.04, numeric(1))

p5e <- ggplot(vdat, aes(group, expression, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  geom_segment(data = ann, aes(x = 1, xend = 2, y = y, yend = y), inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = ann, aes(x = 1.5, y = y, label = label), inherit.aes = FALSE, vjust = -0.5, size = 3) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("TCF7hi" = "#4C78A8", "TCF7-" = "#B33A3A")) +
  scale_x_discrete(labels = c("TCF7hi" = "TCF7hi", "TCF7-" = "TCF7−")) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.18))) +
  labs(x = NULL, y = "Log-normalized RNA expression") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", strip.background = element_blank(), strip.text = element_text(face = "italic"))
ggsave(file.path(out_dir, "suppfig5e_mki67_cdkn2a_by_tcf7_status.pdf"), p5e,
       width = 5.8, height = 4.2, units = "in", device = cairo_pdf)

umap <- Embeddings(obj, "umap")
metric_expr <- FetchData(obj, vars = needed)
patient_map <- c("0701" = "P1", "0708" = "P2", "0712" = "P3", "0825" = "P4")
metric <- data.frame(
  cell_id = sprintf("shared_%06d", seq_len(nrow(umap))), UMAP_1 = umap[, 1], UMAP_2 = umap[, 2],
  tissue = as.character(obj$label), patient = unname(patient_map[as.character(obj$patient)]),
  included_in_cd8_violin = rownames(umap) %in% cd8_cells,
  metric_expr[rownames(umap), needed, drop = FALSE], check.names = FALSE
)
write.csv(metric, gzfile(file.path(out_dir, "suppfig5_shared_clone_umap_expression.csv.gz")), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "suppfig5_sessionInfo.txt"))
