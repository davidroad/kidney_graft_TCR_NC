#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript 05_plot_suppfig6a_from_metrics.R ",
    "<cxcr6_umap_metrics.csv.gz> <cxcr6_frequency_datapoints.csv> <output_directory>"
  )
}

umap_file <- normalizePath(args[[1]], mustWork = TRUE)
frequency_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

umap <- read.csv(gzfile(umap_file), stringsAsFactors = FALSE)
frequency <- read.csv(frequency_file, stringsAsFactors = FALSE)

required_umap <- c("UMAP_1", "UMAP_2", "condition", "dataset", "CXCR6_relative")
required_frequency <- c("source", "sample_id", "figure_group", "cxcr6_frequency_pct")
if (!all(required_umap %in% colnames(umap))) {
  stop("UMAP metrics are missing: ", paste(setdiff(required_umap, colnames(umap)), collapse = ", "))
}
if (!all(required_frequency %in% colnames(frequency))) {
  stop("Frequency table is missing: ", paste(setdiff(required_frequency, colnames(frequency)), collapse = ", "))
}

conditions <- c("PBMC:TOT", "PBMC:STA", "PBMC:BPAR", "PBMC:Rejection", "Graft:Rejection")
condition_titles <- c(
  "PBMC:TOT" = "PBMC;TOT",
  "PBMC:STA" = "PBMC;STA",
  "PBMC:BPAR" = "PBMC;BPAR",
  "PBMC:Rejection" = "PBMC",
  "Graft:Rejection" = "Graft"
)
umap$condition <- factor(umap$condition, levels = conditions)

dataset_limits <- lapply(split(umap, umap$dataset), function(data) {
  x_range <- range(data$UMAP_1, finite = TRUE)
  y_range <- range(data$UMAP_2, finite = TRUE)
  span <- max(diff(x_range), diff(y_range))
  list(
    x = mean(x_range) + c(-0.53, 0.53) * span,
    y = mean(y_range) + c(-0.53, 0.53) * span
  )
})

make_umap_panel <- function(condition) {
  dataset_name <- if (condition %in% c("PBMC:TOT", "PBMC:STA", "PBMC:BPAR")) {
    "GSE224445"
  } else {
    "GSE319007"
  }
  data <- umap[umap$condition == condition & umap$dataset == dataset_name, , drop = FALSE]
  data <- data[order(data$CXCR6_relative), ]
  limits <- dataset_limits[[dataset_name]]

  ggplot(data, aes(UMAP_1, UMAP_2, colour = CXCR6_relative)) +
    geom_point(size = 0.22, alpha = 0.9, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#D9D9D9", "#F5B39D", "#D7553A", "#7F0000"),
      values = c(0, 0.08, 0.42, 1), limits = c(0, 1), guide = "none"
    ) +
    coord_fixed(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = unname(condition_titles[[condition]])) +
    theme_void(base_size = 8) +
    theme(
      plot.title = element_text(size = 8, hjust = 0.5),
      panel.border = element_rect(colour = "#B7B7B7", fill = NA, linewidth = 0.25),
      plot.margin = margin(1, 1, 1, 1)
    )
}

gse_umaps <- wrap_plots(lapply(conditions[1:3], make_umap_panel), nrow = 1) +
  plot_annotation(title = "Bae et al., 2023\n(GSE224445)") &
  theme(plot.title = element_text(size = 8.2, hjust = 0.5, margin = margin(b = 2)))

primary_umaps <- wrap_plots(lapply(conditions[4:5], make_umap_panel), nrow = 1) +
  plot_annotation(title = "This study\n(GSE319007)") &
  theme(plot.title = element_text(size = 8.2, hjust = 0.5, margin = margin(b = 2)))

frequency$figure_group <- factor(frequency$figure_group, levels = conditions)
gse_frequency <- frequency[frequency$source == "GSE224445", , drop = FALSE]
primary_frequency <- frequency[frequency$source == "This study", , drop = FALSE]

mean_frame <- function(data) {
  aggregate(cxcr6_frequency_pct ~ figure_group, data = data, FUN = mean)
}

gse_means <- mean_frame(gse_frequency)
gse_plot <- ggplot(gse_frequency, aes(figure_group, cxcr6_frequency_pct, colour = figure_group)) +
  geom_col(
    data = gse_means,
    aes(figure_group, cxcr6_frequency_pct, colour = figure_group),
    fill = NA, width = 0.58, linewidth = 0.65, inherit.aes = FALSE
  ) +
  geom_point(position = position_jitter(width = 0.08, height = 0, seed = 17), size = 1.5) +
  scale_colour_manual(values = c("PBMC:TOT" = "#222222", "PBMC:STA" = "#4267A9", "PBMC:BPAR" = "#D94B59")) +
  scale_x_discrete(labels = c("PBMC;TOT", "PBMC;STA", "PBMC;BPAR"), drop = TRUE) +
  scale_y_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80), expand = c(0, 0)) +
  labs(title = "Bae et al., 2023\n(GSE224445)", x = NULL, y = "% CXCR6+ among CD8+") +
  theme_classic(base_size = 7.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8, hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(4, 2, 4, 4)
  )

primary_means <- mean_frame(primary_frequency)
primary_plot <- ggplot(primary_frequency, aes(figure_group, cxcr6_frequency_pct)) +
  geom_col(
    data = primary_means,
    aes(figure_group, cxcr6_frequency_pct, colour = figure_group),
    fill = NA, width = 0.58, linewidth = 0.65, inherit.aes = FALSE
  ) +
  geom_line(aes(group = sample_id), colour = "#8D8D8D", linewidth = 0.35) +
  geom_point(aes(colour = figure_group), size = 1.6) +
  scale_colour_manual(values = c("PBMC:Rejection" = "#4267A9", "Graft:Rejection" = "#D94B59")) +
  scale_x_discrete(labels = c("PBMC", "Graft"), drop = TRUE) +
  scale_y_continuous(limits = c(0, 80), breaks = c(0, 20, 40, 60, 80), expand = c(0, 0)) +
  labs(title = "This study\n(GSE319007)", x = NULL, y = "% CXCR6+ among CD8+") +
  theme_classic(base_size = 7.5) +
  theme(
    legend.position = "none",
    plot.title = element_text(size = 8, hjust = 0.5),
    axis.text.x = element_text(angle = 25, hjust = 1),
    plot.margin = margin(4, 2, 4, 4)
  )

figure <- wrap_elements(full = gse_umaps) +
  wrap_elements(full = primary_umaps) +
  gse_plot +
  primary_plot +
  plot_layout(widths = c(3.1, 2.1, 1.15, 1.15)) +
  plot_annotation(tag_levels = list(c("a", "", "", ""))) &
  theme(plot.tag = element_text(face = "bold", size = 11))

ggsave(
  file.path(output_dir, "suppfig6a_cxcr6_expression.pdf"),
  figure,
  width = 12.2,
  height = 3.1,
  units = "in",
  device = cairo_pdf
)

message("Supplementary Fig. 6A written to ", normalizePath(output_dir))
