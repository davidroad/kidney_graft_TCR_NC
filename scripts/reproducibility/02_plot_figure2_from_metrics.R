#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1) args[[1]] else "inputs/plotting_inputs"
output_dir <- if (length(args) >= 2) args[[2]] else "figures/figure_2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

read_metric <- function(name) read.csv(file.path(input_dir, name), check.names = FALSE)
save_pdf <- function(plot, name, width, height) {
  ggsave(file.path(output_dir, name), plot, width = width, height = height, units = "in", device = cairo_pdf)
}

tcell <- read_metric("fig2_tcell_umap_annotations.csv.gz")
shared <- read_metric("fig2_shared_clone_umap_expression.csv.gz")
de <- read_metric("fig2f_shared_clone_graft_vs_pbmc_de.csv")

umap_theme <- theme_void(base_size = 9) + theme(
  legend.title = element_blank(), legend.text = element_text(size = 7),
  plot.margin = margin(4, 4, 4, 4)
)

cluster_levels <- unique(tcell$cluster[order(as.numeric(tcell$cluster_id))])
cluster_colors <- setNames(scales::hue_pal()(length(cluster_levels)), cluster_levels)
tcell$cluster <- factor(tcell$cluster, levels = cluster_levels)
shared$cluster <- factor(shared$cluster, levels = cluster_levels)

p_a <- ggplot(tcell, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.12, alpha = 0.85, stroke = 0) +
  scale_color_manual(values = cluster_colors, na.value = "grey75") +
  coord_equal() + umap_theme
save_pdf(p_a, "fig2a_tcell_cluster_umap.pdf", 8.5, 6.0)

p_d <- ggplot(shared, aes(UMAP_1, UMAP_2, color = cluster)) +
  geom_point(size = 0.28, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = cluster_colors, na.value = "grey75") +
  coord_equal() + umap_theme
save_pdf(p_d, "fig2d_shared_clone_cluster_umap.pdf", 8.5, 6.0)
save_pdf(p_d, "fig2h_shared_clone_cluster_umap.pdf", 8.5, 6.0)

p_e <- ggplot(shared, aes(UMAP_1, UMAP_2, color = tissue)) +
  geom_point(size = 0.28, alpha = 0.9, stroke = 0) +
  scale_color_manual(values = c(PBMC = "#B2182B", Graft = "#2166AC"), na.value = "grey75") +
  coord_equal() + umap_theme
save_pdf(p_e, "fig2e_shared_clone_tissue_umap.pdf", 6.6, 5.8)

de$p_for_plot <- pmax(de$p_val_adj, 1e-300)
de$direction <- ifelse(de$avg_log2FC > 0 & de$p_val_adj < 0.05, "Graft", ifelse(de$avg_log2FC < 0 & de$p_val_adj < 0.05, "PBMC", "Not significant"))
label_genes <- c("KLF2", "TCF7", "CISH", "JUN", "CXCR6", "IFNG", "GZMB")
p_f <- ggplot(de, aes(avg_log2FC, -log10(p_for_plot), color = direction)) +
  geom_point(size = 0.5, alpha = 0.65) +
  geom_text(data = de[de$gene %in% label_genes, ], aes(label = gene), color = "black", size = 2.5, check_overlap = TRUE) +
  scale_color_manual(values = c(PBMC = "#2166AC", `Not significant` = "grey75", Graft = "#B2182B")) +
  labs(x = "Log2 fold change", y = "-Log10 adjusted P value") +
  theme_classic(base_size = 9) + theme(legend.title = element_blank())
save_pdf(p_f, "fig2f_shared_clone_graft_vs_pbmc_volcano.pdf", 6.6, 5.4)

feature_genes <- c("TCF7", "CXCR6", "GZMB")
feature_plots <- lapply(feature_genes, function(gene) {
  ggplot(shared, aes(UMAP_1, UMAP_2, color = .data[[gene]])) +
    geom_point(size = 0.25, alpha = 0.9, stroke = 0) +
    facet_wrap(~ tissue, nrow = 1) +
    scale_color_gradient(low = "grey88", high = "#B2182B", name = gene) +
    coord_equal() + umap_theme + ggtitle(gene) +
    theme(strip.text = element_text(face = "bold"), plot.title = element_text(face = "italic", hjust = 0.5))
})
p_g <- wrap_plots(feature_plots, ncol = 1, guides = "keep")
save_pdf(p_g, "fig2g_shared_clone_featureplots_split_tissue.pdf", 7.2, 12.0)

message("Figure 2 metric-only panels written to ", normalizePath(output_dir))
message("Fig. 2C is generated from filtered V(D)J contig annotations with scripts/reproducibility/02b_generate_fig2c_screpertoire_alluvial.R.")
