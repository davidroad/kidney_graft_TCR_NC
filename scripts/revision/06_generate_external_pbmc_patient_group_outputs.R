#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(repo_lib)) .libPaths(c(repo_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
root_dir <- if (nzchar(project_root)) normalizePath(project_root, mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
cache_rds <- file.path(
  root_dir,
  "data",
  "processed",
  "external_pbmc_patient_group_analysis",
  "visualization_cache",
  "GSE224445_visualization_environment.rds"
)
figure_dir <- file.path(root_dir, "results", "figures", "external_pbmc_patient_group_analysis")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
stopifnot(file.exists(cache_rds))

viz <- readRDS(cache_rds)

tag_summary <- viz$tag_summary
emb <- viz$cd8_umap_df
cd8 <- viz$cd8
DefaultAssay(cd8) <- "RNA"
group_levels <- viz$group_levels
group_cols <- viz$group_cols
status_cols <- viz$status_cols
group_comparisons <- viz$group_comparisons
status_comparison <- viz$status_comparison
analysis_label <- viz$analysis_label

tag_summary <- tag_summary %>%
  mutate(
    group = factor(group, levels = group_levels),
    tolerance_status = factor(
      tolerance_status,
      levels = c("Operational tolerance", "Non-tolerant KTR")
    ),
    tolerance_status_plot = recode(
      as.character(tolerance_status),
      "Operational tolerance" = "Operational tolerance",
      "Non-tolerant KTR" = "STA+BPAR KTR"
    ),
    tolerance_status_plot = factor(tolerance_status_plot, levels = c("Operational tolerance", "STA+BPAR KTR")),
    cxcr6_frequency_pct = 100 * cxcr6_positive_frequency
  )
emb <- emb %>%
  mutate(group = factor(group, levels = group_levels))

status_cols_plot <- c(
  "Operational tolerance" = status_cols[["Operational tolerance"]],
  "STA+BPAR KTR" = status_cols[["Non-tolerant KTR"]]
)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

p_to_star <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 0.0001) return("****")
  if (p < 0.001) return("***")
  if (p < 0.01) return("**")
  if (p < 0.05) return("*")
  "ns"
}

directional_p <- function(metric, group_a, group_b) {
  x <- tag_summary[[metric]][tag_summary$group == group_a]
  y <- tag_summary[[metric]][tag_summary$group == group_b]
  alternative <- if (group_a == "TOT") "less" else "two.sided"
  wilcox.test(x, y, alternative = alternative, exact = FALSE)$p.value
}

directional_status_p <- function(metric) {
  wilcox.test(
    tag_summary[[metric]][tag_summary$tolerance_status == "Operational tolerance"],
    tag_summary[[metric]][tag_summary$tolerance_status == "Non-tolerant KTR"],
    alternative = "less",
    exact = FALSE
  )$p.value
}

square_umap_limits <- function(x, y, probs = c(0.002, 0.998), pad = 0.03) {
  xlim <- as.numeric(quantile(x, probs = probs, na.rm = TRUE))
  ylim <- as.numeric(quantile(y, probs = probs, na.rm = TRUE))
  x_mid <- mean(xlim)
  y_mid <- mean(ylim)
  span <- max(diff(xlim), diff(ylim)) * (1 + pad)
  list(
    xlim = x_mid + c(-0.5, 0.5) * span,
    ylim = y_mid + c(-0.5, 0.5) * span
  )
}

sig_df_group <- function(metric, y_values, comparisons) {
  plot_comps <- group_comparisons %>%
    filter(.data$metric == .env$metric, .data$comparison %in% comparisons) %>%
    mutate(
      xmin = match(group_a, group_levels),
      xmax = match(group_b, group_levels),
      display_p = mapply(directional_p, metric, group_a, group_b),
      label = vapply(display_p, p_to_star, character(1))
    ) %>%
    filter(label != "ns")
  if (!nrow(plot_comps)) return(plot_comps)
  y_range <- diff(range(y_values, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) y_range <- max(y_values, na.rm = TRUE)
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  y_start <- max(y_values, na.rm = TRUE) + 0.10 * y_range
  plot_comps$y <- y_start + seq_len(nrow(plot_comps)) * 0.11 * y_range
  plot_comps$tick <- 0.025 * y_range
  plot_comps
}

sig_df_status <- function(metric, y_values) {
  p <- directional_status_p(metric)
  label <- p_to_star(p)
  if (label == "ns") return(data.frame())
  y_range <- diff(range(y_values, na.rm = TRUE))
  if (!is.finite(y_range) || y_range == 0) y_range <- max(y_values, na.rm = TRUE)
  if (!is.finite(y_range) || y_range == 0) y_range <- 1
  data.frame(
    xmin = 1,
    xmax = 2,
    y = max(y_values, na.rm = TRUE) + 0.14 * y_range,
    tick = 0.025 * y_range,
    label = label,
    stringsAsFactors = FALSE
  )
}

add_sig_layers <- function(plot, sig_df) {
  if (!nrow(sig_df)) return(plot)
  plot +
    geom_segment(
      data = sig_df,
      aes(x = xmin, xend = xmax, y = y, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_segment(
      data = sig_df,
      aes(x = xmin, xend = xmin, y = y - tick, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_segment(
      data = sig_df,
      aes(x = xmax, xend = xmax, y = y - tick, yend = y),
      inherit.aes = FALSE,
      linewidth = 0.35,
      colour = "black"
    ) +
    geom_text(
      data = sig_df,
      aes(x = (xmin + xmax) / 2, y = y + tick, label = label),
      inherit.aes = FALSE,
      size = 3.8,
      fontface = "bold",
      colour = "black"
    )
}

umap_lims <- square_umap_limits(emb$UMAP_1, emb$UMAP_2)

theme_nc <- theme_classic(base_size = 8) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", size = 9, hjust = 0),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7),
    strip.background = element_blank(),
    strip.text = element_text(face = "bold", colour = "black")
  )

style_seurat_umap <- function(plot, lims, title = NULL) {
  styled <- plot +
    coord_fixed(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_classic(base_size = 7, base_family = "Helvetica") +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(colour = "black", size = 7),
      plot.title = element_text(face = "bold", hjust = 0, size = 8),
      legend.title = element_text(size = 6.5),
      legend.text = element_text(size = 6.5),
      legend.key.size = unit(3, "mm"),
      plot.margin = margin(3, 3, 3, 3)
    )
  if (!is.null(title)) {
    styled <- styled + labs(title = title)
  }
  styled
}

p_umap_group <- ggplot(emb, aes(UMAP_1, UMAP_2, colour = group)) +
  geom_point(size = 0.10, alpha = 0.82) +
  scale_colour_manual(values = group_cols, drop = FALSE) +
  coord_fixed(xlim = umap_lims$xlim, ylim = umap_lims$ylim, expand = FALSE) +
  labs(title = analysis_label, x = "UMAP 1", y = "UMAP 2", colour = NULL) +
  guides(colour = guide_legend(override.aes = list(size = 2.2, alpha = 1))) +
  theme_nc +
  theme(aspect.ratio = 1, axis.text = element_blank(), axis.ticks = element_blank())

p_umap_cxcr6 <- FeaturePlot(
  cd8,
  reduction = viz$cd8_reduction,
  features = "CXCR6",
  cols = c("lightgrey", "firebrick3"),
  pt.size = 0.22,
  order = TRUE,
  combine = FALSE
)[[1]]
p_umap_cxcr6 <- style_seurat_umap(p_umap_cxcr6, umap_lims, "CXCR6")

sig_group_freq <- sig_df_group(
  "cxcr6_positive_frequency",
  tag_summary$cxcr6_frequency_pct,
  c("TOT vs STA")
)
sig_status_freq <- sig_df_status("cxcr6_positive_frequency", tag_summary$cxcr6_frequency_pct)
sig_group_expr <- sig_df_group(
  "mean_cxcr6_log_norm",
  tag_summary$mean_cxcr6_log_norm,
  c("TOT vs STA")
)
sig_status_expr <- sig_df_status("mean_cxcr6_log_norm", tag_summary$mean_cxcr6_log_norm)

p_group_freq <- ggplot(tag_summary, aes(group, cxcr6_frequency_pct, colour = group)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_point(aes(shape = bead), size = 2.5, position = position_jitter(width = 0.08, height = 0), stroke = 0.7) +
  scale_colour_manual(values = group_cols, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.03, 0.32))) +
  labs(
    title = "CXCR6+ CD8 frequency by group",
    subtitle = paste0("TOT < STA directional Wilcoxon p = ", fmt_p(directional_p("cxcr6_positive_frequency", "TOT", "STA"))),
    x = NULL,
    y = "CXCR6+ CD8 T cells",
    colour = NULL,
    shape = "Bead"
  ) +
  theme_nc +
  theme(legend.position = "right")
p_group_freq <- add_sig_layers(p_group_freq, sig_group_freq)

p_status_freq <- ggplot(tag_summary, aes(tolerance_status_plot, cxcr6_frequency_pct, colour = tolerance_status_plot)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_point(aes(shape = group), size = 2.6, position = position_jitter(width = 0.08, height = 0), stroke = 0.7) +
  scale_colour_manual(values = status_cols_plot, drop = FALSE) +
  scale_y_continuous(labels = function(x) paste0(x, "%"), expand = expansion(mult = c(0.03, 0.22))) +
  labs(
    title = "Tolerance versus STA+BPAR KTR",
    subtitle = paste0("Directional Wilcoxon p = ", fmt_p(directional_status_p("cxcr6_positive_frequency"))),
    x = NULL,
    y = "CXCR6+ CD8 T cells",
    colour = NULL,
    shape = "Group"
  ) +
  theme_nc +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
p_status_freq <- add_sig_layers(p_status_freq, sig_status_freq)

p_group_expr <- ggplot(tag_summary, aes(group, mean_cxcr6_log_norm, colour = group)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_point(aes(shape = bead), size = 2.5, position = position_jitter(width = 0.08, height = 0), stroke = 0.7) +
  scale_colour_manual(values = group_cols, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.32))) +
  labs(
    title = "CXCR6 expression by group",
    subtitle = paste0("TOT < STA directional Wilcoxon p = ", fmt_p(directional_p("mean_cxcr6_log_norm", "TOT", "STA"))),
    x = NULL,
    y = "Mean tag log-normalized CXCR6",
    colour = NULL,
    shape = "Bead"
  ) +
  theme_nc
p_group_expr <- add_sig_layers(p_group_expr, sig_group_expr)

p_status_expr <- ggplot(tag_summary, aes(tolerance_status_plot, mean_cxcr6_log_norm, colour = tolerance_status_plot)) +
  geom_boxplot(width = 0.52, outlier.shape = NA, linewidth = 0.35, fill = "white") +
  geom_point(aes(shape = group), size = 2.6, position = position_jitter(width = 0.08, height = 0), stroke = 0.7) +
  scale_colour_manual(values = status_cols_plot, drop = FALSE) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.25))) +
  labs(
    title = "CXCR6 expression by tolerance status",
    subtitle = paste0("Directional Wilcoxon p = ", fmt_p(directional_status_p("mean_cxcr6_log_norm"))),
    x = NULL,
    y = "Mean tag log-normalized CXCR6",
    colour = NULL,
    shape = "Group"
  ) +
  theme_nc +
  theme(axis.text.x = element_text(angle = 20, hjust = 1), legend.position = "right")
p_status_expr <- add_sig_layers(p_status_expr, sig_status_expr)

fig <- (p_umap_group | p_umap_cxcr6) /
  (p_group_freq | p_status_freq) /
  (p_group_expr | p_status_expr) +
  plot_layout(heights = c(1.25, 1, 1), guides = "collect") +
  plot_annotation(
    tag_levels = "A",
    title = "GSE224445 PBMC CD8 T cells: CXCR6 is lowest in operational tolerance",
    theme = theme(
      plot.title = element_text(face = "bold", size = 11, hjust = 0),
      plot.tag = element_text(face = "bold", size = 11),
      plot.tag.position = c(0, 1)
    )
  )

pdf_file <- file.path(figure_dir, "external_pbmc_patient_group_summary.pdf")
png_file <- file.path(figure_dir, "external_pbmc_patient_group_summary.png")
ggsave(pdf_file, fig, width = 7.2, height = 8.8, units = "in", device = "pdf")
if (requireNamespace("ragg", quietly = TRUE)) {
  ggsave(png_file, fig, width = 7.2, height = 8.8, units = "in", dpi = 600, device = ragg::agg_png)
}

message("Loaded visualization environment: ", cache_rds)
message("Wrote: ", pdf_file)
message("Wrote: ", png_file)
