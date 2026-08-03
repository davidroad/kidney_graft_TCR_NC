#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 04_generate_suppfig1_2_from_objects.R <all-cell.RData> <T-cell.RData> <output-dir>")
}
all_cell_file <- normalizePath(args[[1]], mustWork = TRUE)
t_cell_file <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

batch_levels <- c(
  "graft_0701", "graft_0708", "graft_0712", "graft_0825",
  "PBMC_0701", "PBMC_0708", "PBMC_0712", "PBMC_0825"
)
batch_labels <- c(
  "Graft 07012022", "Graft 07082022", "Graft 07122021", "Graft 08252021",
  "PBMC 07012022", "PBMC 07082022", "PBMC 07122021", "PBMC 08252021"
)

all_markers <- c(
  "S100A9", "S100A8", "LYZ", "IL7R", "LTB", "AQP3", "GZMK", "CCL5", "CD8B",
  "TCL1A", "CD79A", "MS4A1", "GNLY", "NKG7", "FGFBP2", "CXCL13", "LAG3", "CCL4L2",
  "VCAN", "CST3", "CCR7", "LEF1", "TCF7", "TRBV13", "TRBV12-4", "TRBV5-4", "FOXP3",
  "TNFRSF4", "IL32", "TRDV1", "LINC02446", "ANK3", "TSHZ2", "BCL2", "HLA-DRA", "XCL2",
  "XCL1", "SMIM25", "SERPINA1", "MS4A7", "IGLV3-21", "JCHAIN", "IGHA1", "APOE", "C1QA",
  "FN1", "NAMPT", "STMN1", "TUBA1B", "TYMS", "AFF3", "RALGPS2", "BACH2", "IGLV1-51",
  "IGHV1-24", "IGLV3-25", "IGHM", "BANK1", "HLA-DPA1", "FCER1A", "PLD4", "ITM2C",
  "IRF8", "ALAS2", "SLC25A37", "HBQ1", "CLC", "GATA2", "HDC", "IGKV1-39", "IGKC"
)
t_markers <- c(
  "FOS", "IL7R", "TCF7", "CCL4L2", "CCL4", "GZMB", "TRBV4-2", "TRBV13", "TRBV27",
  "KLRB1", "LTB", "GZMK", "AOAH", "GZMA", "CCR7", "FHIT", "LEF1", "FOXP3", "RTKN2",
  "IKZF2", "GNLY", "NKG7", "GZMH", "ZNF683", "LINC02446", "XCL1", "CXCL13", "TOX2",
  "TNFRSF4", "ISG15", "MX1", "IFI6"
)

load(all_cell_file)
if (!exists("merged.all.batchcorrected") || !inherits(merged.all.batchcorrected, "Seurat")) {
  stop("The all-cell RData file must contain the Seurat object 'merged.all.batchcorrected'.")
}
all_obj <- merged.all.batchcorrected
rm(merged.all.batchcorrected)
if (exists("merged_combined_regressed")) rm(merged_combined_regressed)
gc()

all_obj$plot_cluster <- factor(as.character(all_obj$integrated_snn_res.0.5), levels = as.character(0:27))
all_obj$plot_batch <- factor(as.character(all_obj$batch), levels = batch_levels, labels = batch_labels)
Idents(all_obj) <- "plot_cluster"

p_s1a <- DimPlot(
  all_obj, reduction = "umap", group.by = "plot_cluster", split.by = "plot_batch",
  ncol = 4, pt.size = 0.12, shuffle = TRUE, seed = 1, raster = FALSE
)
p_s1a <- (p_s1a + patchwork::plot_annotation(title = "")) &
  labs(x = "UMAP_1", y = "UMAP_2", colour = "Cluster")
ggsave(file.path(out_dir, "suppfig1a_all_cd45_umap_by_sample.pdf"), p_s1a,
       width = 15.5, height = 7.5, units = "in", device = cairo_pdf)

DefaultAssay(all_obj) <- "integrated"
missing_all <- setdiff(all_markers, rownames(all_obj[["integrated"]]))
if (length(missing_all)) stop("Missing Supplementary Fig. 1 markers: ", paste(missing_all, collapse = ", "))
set.seed(1)
all_small <- subset(all_obj, downsample = 200)
p_s1b <- DoHeatmap(
  all_small, features = all_markers, group.by = "plot_cluster", slot = "scale.data",
  size = 2.4, label = FALSE, raster = TRUE, angle = 0
) + labs(x = "Cluster", y = NULL)
ggsave(file.path(out_dir, "suppfig1b_all_cd45_marker_heatmap.pdf"), p_s1b,
       width = 15.5, height = 11, units = "in", device = cairo_pdf)
rm(all_obj, all_small, p_s1a, p_s1b)
gc()

load(t_cell_file)
if (!exists("merged_all_T_batchcorrected_filter") || !inherits(merged_all_T_batchcorrected_filter, "Seurat")) {
  stop("The T-cell RData file must contain 'merged_all_T_batchcorrected_filter'.")
}
t_obj <- merged_all_T_batchcorrected_filter
rm(merged_all_T_batchcorrected_filter)
if (exists("combined")) rm(combined)
gc()

annotation_order <- c(
  "Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+", "CXCR6+ effector CD8+",
  "CCR2+ CD4+", "Naive/Memory-like CD8+", "Naive CD4+", "Treg",
  "CX3CR1+ effector CD8+", "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+"
)
t_obj$plot_cluster <- factor(
  match(as.character(t_obj$cell_type_annotated_granularity), annotation_order) - 1L,
  levels = 0:10
)
t_obj$plot_batch <- factor(as.character(t_obj$batch), levels = batch_levels, labels = batch_labels)
if (anyNA(t_obj$plot_cluster)) stop("An unrecognized T-cell annotation was found.")
Idents(t_obj) <- "plot_cluster"

p_s2a <- DimPlot(
  t_obj, reduction = "umap", group.by = "plot_cluster", split.by = "plot_batch",
  ncol = 4, pt.size = 0.22, shuffle = TRUE, seed = 1, raster = FALSE
)
p_s2a <- (p_s2a + patchwork::plot_annotation(title = "")) &
  labs(x = "UMAP_1", y = "UMAP_2", colour = "Cluster")
ggsave(file.path(out_dir, "suppfig2a_tcell_umap_by_sample.pdf"), p_s2a,
       width = 15.5, height = 7.5, units = "in", device = cairo_pdf)

DefaultAssay(t_obj) <- "integrated"
missing_t <- setdiff(t_markers, rownames(t_obj[["integrated"]]))
if (length(missing_t)) stop("Missing Supplementary Fig. 2 markers: ", paste(missing_t, collapse = ", "))
set.seed(1)
t_small <- subset(t_obj, downsample = 200)
p_s2b <- DoHeatmap(
  t_small, features = t_markers, group.by = "plot_cluster", slot = "scale.data",
  size = 2.8, label = FALSE, raster = TRUE, angle = 0
) + labs(x = "Cluster", y = NULL)
ggsave(file.path(out_dir, "suppfig2b_tcell_marker_heatmap.pdf"), p_s2b,
       width = 13.5, height = 8.5, units = "in", device = cairo_pdf)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "suppfig1_2_sessionInfo.txt"))
