#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
input_dir <- if (length(args) >= 1L) args[[1]] else "inputs/plotting_inputs"
table_dir <- if (length(args) >= 2L) args[[2]] else "tables"
output_dir <- if (length(args) >= 3L) args[[3]] else "figures/supplementary_figure_13"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

score <- read.csv(
  gzfile(file.path(input_dir, "suppfig13a_cd8_glycolysis_score_umap_metrics.csv.gz")),
  check.names = FALSE
)
paired <- read.delim(file.path(table_dir, "paired_glycolysis_exact_datapoints.tsv"), check.names = FALSE)
gene <- read.delim(
  file.path(table_dir, "suppfig13c_paired_pseudobulk_glycolysis_gene_statistics.tsv"),
  check.names = FALSE
)

theme_umap <- theme_void(base_size = 9) +
  theme(legend.position = "right", plot.margin = margin(4, 4, 4, 4))

score <- score[order(score$glycolysis_score), ]
p_a <- ggplot(score, aes(UMAP_1, UMAP_2, colour = glycolysis_score)) +
  geom_point(size = 0.25, alpha = 0.9, stroke = 0) +
  scale_colour_gradientn(
    colours = c("#D8D8D8", "#F8D4CC", "#E75A4B", "firebrick3"),
    name = "Glycolysis score"
  ) +
  coord_equal() +
  theme_umap

paired$condition <- factor(paired$condition, levels = c("PBMC", "Graft"))
p_b <- ggplot(paired, aes(condition, mean_score, group = patient)) +
  geom_line(colour = "grey65", size = 0.45) +
  geom_point(aes(fill = condition), shape = 21, size = 2.2, colour = "black", stroke = 0.3) +
  scale_fill_manual(values = c(PBMC = "#C9D5E3", Graft = "#E7B4AE")) +
  labs(x = NULL, y = "Mean glycolysis score") +
  theme_classic(base_size = 9) +
  theme(legend.position = "none")

gene$gene <- factor(gene$gene, levels = rev(gene$gene[order(gene$display_order)]))
gene$significant <- gene$FDR < 0.05
p_c <- ggplot(gene, aes(logFC, gene, colour = significant)) +
  geom_vline(xintercept = 0, size = 0.3, colour = "grey70") +
  geom_segment(aes(x = 0, xend = logFC, yend = gene), size = 0.55) +
  geom_point(aes(size = -log10(pmax(FDR, .Machine$double.xmin))), shape = 16) +
  scale_colour_manual(values = c(`FALSE` = "grey55", `TRUE` = "#B2182B")) +
  labs(x = "Log2 fold change, graft versus PBMC", y = NULL, size = "-log10 FDR") +
  theme_classic(base_size = 9) +
  theme(legend.title = element_text(size = 8), legend.text = element_text(size = 7))

ggsave(
  file.path(output_dir, "suppfig13a_cd8_glycolysis_score_umap.pdf"),
  p_a, width = 6.2, height = 5.2, units = "in", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "suppfig13b_paired_glycolysis_score.pdf"),
  p_b, width = 3.2, height = 3.2, units = "in", device = cairo_pdf
)
ggsave(
  file.path(output_dir, "suppfig13c_paired_pseudobulk_glycolysis_genes_lollipop.pdf"),
  p_c, width = 5.5, height = 4.4, units = "in", device = cairo_pdf
)

message("Supplementary Fig. 13A-C outputs written to ", normalizePath(output_dir))
