#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop("Usage: Rscript 04_plot_suppfig1_5_from_metrics.R <plotting-input-dir> <table-dir> <output-dir>")
}
input_dir <- normalizePath(args[[1]], mustWork = TRUE)
table_dir <- normalizePath(args[[2]], mustWork = TRUE)
out_dir <- args[[3]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

save_panel <- function(plot, subdir, filename, width, height) {
  target <- file.path(out_dir, subdir)
  dir.create(target, recursive = TRUE, showWarnings = FALSE)
  ggsave(file.path(target, filename), plot, width = width, height = height,
         units = "in", device = cairo_pdf)
}

umap_panel <- function(data, cluster_levels, filename, subdir, width = 15.5) {
  data$cluster_id <- factor(as.character(data$cluster_id), levels = cluster_levels)
  data$sample <- factor(as.character(data$sample), levels = unique(as.character(data$sample)))
  palette <- setNames(hcl.colors(length(cluster_levels), "Dynamic"), cluster_levels)
  p <- ggplot(data, aes(UMAP_1, UMAP_2, colour = cluster_id)) +
    geom_point(size = 0.11, alpha = 0.9) +
    facet_wrap(~sample, ncol = 4) +
    scale_colour_manual(values = palette, drop = FALSE) +
    labs(x = "UMAP_1", y = "UMAP_2", colour = "Cluster") +
    theme_classic(base_size = 8) +
    theme(strip.background = element_blank(), strip.text = element_text(face = "bold"))
  save_panel(p, subdir, filename, width, 7.5)
}

heatmap_panel <- function(metric_file, filename, subdir, width, height) {
  data <- read.csv(gzfile(metric_file), check.names = FALSE)
  genes <- setdiff(names(data), c("cell_id", "cluster_id"))
  cluster_values <- as.character(data$cluster_id)
  cluster_levels <- if (all(grepl("^[0-9]+$", cluster_values))) {
    as.character(sort(unique(as.integer(cluster_values))))
  } else {
    unique(cluster_values)
  }
  data$cluster_id <- factor(cluster_values, levels = cluster_levels)
  data <- data[order(data$cluster_id, data$cell_id), , drop = FALSE]
  data$x <- seq_len(nrow(data))
  matrix_values <- as.matrix(data[, genes, drop = FALSE])
  long <- data.frame(
    x = rep(data$x, times = length(genes)),
    gene = factor(rep(genes, each = nrow(data)), levels = rev(genes)),
    expression = as.vector(matrix_values)
  )
  breaks <- cumsum(as.integer(table(data$cluster_id)))
  breaks <- breaks[breaks < nrow(data)] + 0.5
  cluster_cols <- setNames(hcl.colors(nlevels(data$cluster_id), "Dynamic"), levels(data$cluster_id))
  strip <- ggplot(data, aes(x, 1, fill = cluster_id)) +
    geom_raster(key_glyph = "polygon") +
    geom_vline(xintercept = breaks, colour = "white", linewidth = 0.65) +
    annotate("rect", xmin = 0.5, xmax = nrow(data) + 0.5,
             ymin = 0.5, ymax = 1.5, fill = NA, colour = "black", linewidth = 0.25) +
    scale_fill_manual(values = cluster_cols, drop = FALSE) +
    theme_void(base_size = 7) +
    theme(
      legend.key = element_rect(fill = "white", colour = "grey30", linewidth = 0.35),
      legend.key.spacing.x = grid::unit(2.5, "mm"),
      legend.key.spacing.y = grid::unit(1.5, "mm"),
      legend.text = element_text(margin = margin(r = 2))
    ) +
    guides(fill = guide_legend(
      title = "Cluster", nrow = 2, byrow = TRUE,
      keywidth = grid::unit(5, "mm"), keyheight = grid::unit(3.5, "mm"),
      override.aes = list(colour = "grey30", linewidth = 0.35)
    ))
  heat <- ggplot(long, aes(x, gene, fill = expression)) +
    geom_raster() +
    geom_vline(xintercept = breaks, colour = "white", linewidth = 0.18) +
    scale_fill_gradient2(low = "#D800D8", mid = "black", high = "yellow", midpoint = 0,
                         limits = c(-2.5, 2.5), oob = scales::squish, name = "Expression") +
    labs(x = "Cells ordered by cluster", y = NULL) +
    theme_classic(base_size = 7) +
    theme(axis.ticks.y = element_blank())
  p <- strip / heat + plot_layout(heights = c(0.32, 10), guides = "collect")
  save_panel(p, subdir, filename, width, height)
}

fig1_umap <- read.csv(gzfile(file.path(input_dir, "fig1_all_cd45_umap_annotations.csv.gz")))
fig2_umap <- read.csv(gzfile(file.path(input_dir, "fig2_tcell_umap_annotations.csv.gz")))
annotation_order <- c(
  "Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+", "CXCR6+ effector CD8+",
  "CCR2+ CD4+", "Naive/Memory-like CD8+", "Naive CD4+", "Treg",
  "CX3CR1+ effector CD8+", "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+"
)
fig2_umap$cluster_id <- match(fig2_umap$cluster, annotation_order) - 1L
if (anyNA(fig2_umap$cluster_id)) stop("Unrecognized T-cell annotation in Fig. 2 metric.")

umap_panel(fig1_umap, as.character(0:27), "suppfig1a_all_cd45_umap_by_sample.pdf",
           "supplementary_figure_1")
heatmap_panel(file.path(input_dir, "suppfig1b_all_cd45_marker_heatmap_metrics.csv.gz"),
              "suppfig1b_all_cd45_marker_heatmap.pdf", "supplementary_figure_1", 15.5, 11)
umap_panel(fig2_umap, as.character(0:10), "suppfig2a_tcell_umap_by_sample.pdf",
           "supplementary_figure_2")
heatmap_panel(file.path(input_dir, "suppfig2b_tcell_marker_heatmap_metrics.csv.gz"),
              "suppfig2b_tcell_marker_heatmap.pdf", "supplementary_figure_2", 13.5, 8.5)

occ <- read.csv(file.path(table_dir, "supplementary_figure_3", "suppfig3a_clonal_occupancy_counts.csv"))
occ$cluster <- factor(occ$cluster, levels = unique(occ$cluster))
occ$clone_size <- factor(occ$clone_size, levels = c(
  "Single (<=1)", "Small (1-5)", "Medium (5-20)", "Large (20-100)", "Hyperexpanded (100-500)"
))
p3a <- ggplot(occ, aes(cluster, cell_number, fill = clone_size)) +
  geom_col(width = 0.92, colour = "black", linewidth = 0.12) +
  geom_text(aes(label = ifelse(cell_number >= 10, cell_number, "")),
            position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = c("#111111", "#542057", "#D41473", "#F08A12", "#FFF37A"),
                    drop = FALSE) +
  labs(x = NULL, y = "Cell number", fill = "Clone size") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1))
save_panel(p3a, "supplementary_figure_3", "suppfig3a_clonal_occupancy_by_tcell_cluster.pdf", 12, 7.5)

overlap <- read.csv(file.path(table_dir, "supplementary_figure_3", "suppfig3b_clonotype_overlap_values.csv"))
sample_order <- unique(c(overlap$x, overlap$y))
overlap$x <- factor(overlap$x, levels = sample_order)
overlap$y <- factor(overlap$y, levels = rev(sample_order))
heat_theme <- theme_minimal(base_size = 7) +
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1))
p_j <- ggplot(overlap, aes(x, y, fill = jaccard)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = ifelse(jaccard >= 0.001, sprintf("%.3f", jaccard), "0")), size = 2) +
  scale_fill_viridis_c(option = "B", name = "Jaccard") + coord_equal() + heat_theme
p_r <- ggplot(overlap, aes(x, y, fill = raw)) +
  geom_tile(colour = "white", linewidth = 0.25) +
  geom_text(aes(label = raw), size = 2) +
  scale_fill_viridis_c(option = "B", name = "Raw") + coord_equal() + heat_theme
save_panel(p_j + p_r, "supplementary_figure_3", "suppfig3b_tcr_clonotype_overlap.pdf", 13.5, 6.5)

scatter <- read.csv(file.path(table_dir, "supplementary_figure_4", "suppfig4_clonotype_expansion_plot_data.csv"))
scatter$patient <- sprintf("%04d", as.integer(scatter$patient))
make_scatter <- function(patient_id) {
  data <- scatter[scatter$patient == patient_id, , drop = FALSE]
  ggplot(data, aes(prop_pbmc, prop_graft, size = total_n, colour = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_point(alpha = 0.82) +
    scale_size_area(max_size = 5) +
    labs(x = "Proportion in PBMC", y = "Proportion in Graft", colour = "Class",
         size = "Total n", title = paste("Patient", patient_id)) +
    theme_classic(base_size = 8) +
    theme(legend.position = "right")
}
p4 <- wrap_plots(lapply(sort(unique(scatter$patient)), make_scatter), ncol = 2)
save_panel(p4, "supplementary_figure_4", "suppfig4_clonotype_expansion_graft_vs_pbmc.pdf", 13.5, 11)

shared <- read.csv(gzfile(file.path(input_dir, "suppfig5_shared_clone_umap_expression.csv.gz")),
                   check.names = FALSE)
panel_genes <- list(
  a = c("KLF2", "IL7R", "SELL"),
  b = c("TOX", "PDCD1", "HAVCR2"),
  c = c("IFNG", "PRF1", "ID2"),
  d = c("MKI67", "CDKN2A")
)
make_feature_group <- function(genes, tag) {
  plots <- lapply(genes, function(gene) {
    data <- shared[, c("UMAP_1", "UMAP_2", "tissue", gene)]
    names(data)[4] <- "expression"
    ggplot(data[order(data$expression), ], aes(UMAP_1, UMAP_2, colour = expression)) +
      geom_point(size = 0.14) +
      facet_wrap(~tissue, nrow = 1) +
      scale_colour_gradient(low = "lightgrey", high = "red") +
      labs(x = "UMAP_1", y = "UMAP_2", colour = gene, title = gene) +
      theme_classic(base_size = 8) +
      theme(legend.position = "none", plot.title = element_text(face = "italic", hjust = 0.5),
            strip.background = element_blank())
  })
  save_panel(wrap_plots(plots, nrow = 1), "supplementary_figure_5",
             paste0("suppfig5", tag, "_shared_clone_featureplots.pdf"),
             if (length(genes) == 3) 15 else 10, 4.2)
}
invisible(mapply(make_feature_group, panel_genes, names(panel_genes), SIMPLIFY = FALSE))

stats5 <- read.csv(file.path(table_dir, "supplementary_figure_5",
                             "suppfig5e_cell_level_wilcoxon_statistics.csv"))
selected <- shared[shared$included_in_cd8_violin, , drop = FALSE]
vdat <- rbind(
  data.frame(group = ifelse(selected$TCF7 > 0, "TCF7hi", "TCF7-"),
             gene = "MKI67", expression = selected$MKI67),
  data.frame(group = ifelse(selected$TCF7 > 0, "TCF7hi", "TCF7-"),
             gene = "CDKN2A", expression = selected$CDKN2A)
)
vdat$group <- factor(vdat$group, levels = c("TCF7hi", "TCF7-"))
vdat$gene <- factor(vdat$gene, levels = c("MKI67", "CDKN2A"))
ann <- stats5[match(levels(vdat$gene), stats5$gene), ]
ann$gene <- factor(ann$gene, levels = levels(vdat$gene))
ann$label <- ifelse(ann$wilcoxon_fdr < 1e-4, "<0.0001", sprintf("%.4f", ann$wilcoxon_fdr))
ann$y <- vapply(levels(vdat$gene), function(g) {
  max(vdat$expression[vdat$gene == g], na.rm = TRUE) * 1.08 + 0.04
}, numeric(1))
p5e <- ggplot(vdat, aes(group, expression, fill = group)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.25) +
  geom_segment(data = ann, aes(x = 1, xend = 2, y = y, yend = y),
               inherit.aes = FALSE, linewidth = 0.4) +
  geom_text(data = ann, aes(x = 1.5, y = y, label = label),
            inherit.aes = FALSE, vjust = -0.5, size = 3) +
  facet_wrap(~gene, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = c("TCF7hi" = "#4C78A8", "TCF7-" = "#B33A3A")) +
  scale_x_discrete(labels = c("TCF7hi" = "TCF7hi", "TCF7-" = "TCF7−")) +
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.18))) +
  labs(x = NULL, y = "Log-normalized RNA expression") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none", strip.background = element_blank(),
        strip.text = element_text(face = "italic"))
save_panel(p5e, "supplementary_figure_5", "suppfig5e_mki67_cdkn2a_by_tcf7_status.pdf", 5.8, 4.2)

writeLines(capture.output(sessionInfo()), file.path(out_dir, "suppfig1_5_metrics_sessionInfo.txt"))
