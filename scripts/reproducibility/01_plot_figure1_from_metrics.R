#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "inputs/plotting_inputs"
output_dir <- if (length(args) >= 2) args[[2]] else "figures/figure_1"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_metric <- function(name) read.csv(file.path(input_dir, name), check.names = FALSE)
save_pdf <- function(plot, name, width, height) {
  ggsave(file.path(output_dir, name), plot, width = width, height = height, units = "in", device = cairo_pdf)
}

umap <- read_metric("fig1_all_cd45_umap_annotations.csv.gz")
props <- read_metric("fig1_all_cd45_sample_cluster_proportions.csv")
dot <- read_metric("fig1d_marker_dotplot_metrics.csv")

umap_theme <- theme_void(base_size = 9) + theme(
  legend.title = element_blank(),
  legend.text = element_text(size = 6),
  legend.key.height = unit(0.28, "cm"),
  plot.margin = margin(4, 4, 4, 4)
)

cluster_levels <- unique(umap$cluster[order(as.numeric(umap$cluster_id))])
umap$cluster <- factor(umap$cluster, levels = cluster_levels)
cluster_colors <- setNames(scales::hue_pal()(length(cluster_levels)), cluster_levels)

p_b <- ggplot(umap, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.08, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = cluster_colors, na.value = "grey75") +
  coord_equal() + umap_theme
save_pdf(p_b, "fig1b_all_cd45_cluster_umap.pdf", 9.0, 6.2)

sample_levels <- unique(umap$sample)
p_c <- ggplot(umap, aes(UMAP_1, UMAP_2, color = factor(sample, levels = sample_levels))) +
  geom_point(size = 0.08, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = setNames(scales::hue_pal()(length(sample_levels)), sample_levels)) +
  coord_equal() + umap_theme
save_pdf(p_c, "fig1c_all_cd45_sample_umap.pdf", 7.6, 6.2)

dot$cluster <- factor(dot$cluster, levels = rev(cluster_levels))
dot$gene <- factor(dot$gene, levels = unique(dot$gene))
p_d <- ggplot(dot, aes(gene, cluster)) +
  geom_point(aes(size = percent_expressing, color = scaled_average_expression)) +
  scale_color_gradient(low = "grey90", high = "#B2182B", name = "Scaled average\nexpression") +
  scale_size(range = c(0.2, 4.0), name = "% expressing") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1), legend.position = "right")
save_pdf(p_d, "fig1d_all_cd45_marker_dotplot.pdf", 11.0, 8.5)

tissue_colors <- c(PBMC = "#2166AC", Graft = "#B2182B")
p_e <- ggplot(umap, aes(UMAP_1, UMAP_2, color = tissue)) +
  geom_point(size = 0.08, alpha = 0.8, stroke = 0) +
  scale_color_manual(values = tissue_colors, na.value = "grey75") +
  coord_equal() + umap_theme
save_pdf(p_e, "fig1e_all_cd45_tissue_umap.pdf", 6.8, 6.2)

pie <- aggregate(cbind(n_cells, sample_total_cells) ~ tissue + cluster, props, sum)
pie$proportion <- pie$n_cells / ave(pie$n_cells, pie$tissue, FUN = sum)
pie$cluster <- factor(pie$cluster, levels = cluster_levels)
p_f <- ggplot(pie, aes(x = 1, y = proportion, fill = cluster)) +
  geom_col(width = 1, color = NA) +
  coord_polar(theta = "y") +
  facet_wrap(~ tissue, nrow = 1) +
  scale_fill_manual(values = cluster_colors, na.value = "grey75") +
  theme_void(base_size = 9) +
  theme(strip.text = element_text(face = "bold"), legend.title = element_blank(), legend.text = element_text(size = 6))
save_pdf(p_f, "fig1f_all_cd45_cluster_pies.pdf", 9.5, 4.8)

message("Figure 1 metric-only panels written to ", normalizePath(output_dir))
