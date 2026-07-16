#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) default else args[[pos + 1L]]
}

root_dir <- normalizePath(arg_value("--repo", getwd()), mustWork = TRUE)
out_dir <- arg_value("--out", file.path(root_dir, "results", "article_humess_alignment"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
overview_width <- 13.2
overview_height <- 6.8

scores_file <- file.path(root_dir, "results", "reviewer2_metabolic_profile", "input", "cd8_metabolic_scores.tsv.gz")
logcpm_file <- file.path(root_dir, "results", "reviewer2_metabolic_profile", "input", "pseudobulk_logcpm_by_sample.tsv.gz")
counts_file <- file.path(root_dir, "results", "reviewer2_metabolic_profile", "input", "pseudobulk_counts_by_sample.tsv.gz")
edge_file <- file.path(root_dir, "results", "reviewer2_metabolic_profile", "edgeR", "paired_graft_vs_pbmc_edgeR.tsv")
reporter_file <- file.path(
  root_dir,
  "results", "reviewer2_metabolic_profile", "humess", "humess_run",
  "comparisons", "PBMC__vs__Graft", "ReporterMetabolites",
  "reporter_metabolites_Graft_vs_PBMC_all.tsv"
)
humess_model_dir <- file.path(root_dir, "results", "reviewer2_metabolic_profile", "humess", "humess_run", "models")

required_inputs <- c(scores_file, logcpm_file, counts_file, edge_file, reporter_file)
missing_inputs <- required_inputs[!file.exists(required_inputs)]
if (length(missing_inputs)) {
  stop(
    "Required HUMESS-derived input files are missing:\n- ",
    paste(missing_inputs, collapse = "\n- "),
    "\nSee README.md for the expected input layout."
  )
}

glycolysis_score_genes <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
  "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB", "SLC2A1",
  "G6PC3", "AKR1A1"
)

scores <- read_tsv(scores_file, show_col_types = FALSE)
patient_levels <- sort(unique(as.character(scores$patient)))
if (length(patient_levels) != 4L) stop("Expected four patient identifiers in the metabolic-score input.")
patient_map <- setNames(paste0("P", seq_along(patient_levels)), patient_levels)
public_patient <- function(x) unname(patient_map[as.character(x)])
public_sample_name <- function(x) {
  out <- x
  for (id in names(patient_map)) out <- str_replace_all(out, fixed(id), patient_map[[id]])
  out
}

scores <- scores %>%
  mutate(
    condition = factor(condition, levels = c("PBMC", "Graft")),
    patient = factor(public_patient(patient), levels = paste0("P", 1:4)),
    stem_like_annotation = if_else(stem_like_annotation, "TCF7hi stem-like", "Other CD8")
  )

edge <- read_tsv(edge_file, show_col_types = FALSE)
reporter <- read_tsv(reporter_file, show_col_types = FALSE)
logcpm_source <- read_tsv(logcpm_file, show_col_types = FALSE)
counts_source <- read_tsv(counts_file, show_col_types = FALSE)
names(logcpm_source) <- public_sample_name(names(logcpm_source))
names(counts_source) <- public_sample_name(names(counts_source))

pal_condition <- c(PBMC = "#3B5A7A", Graft = "#B33A3A")
pal_fill <- c(PBMC = "#C9D5E3", Graft = "#E7B4AE")
pal_celltype <- c(
  "CXCL13+CXCR6+ effector CD8+" = "#B33A3A",
  "CXCR6+ effector CD8+" = "#D95A50",
  "ZNF683+TCF7hi CD8+" = "#8E6BBE",
  "CX3CR1+ effector CD8+" = "#3B5A7A",
  "Naive/Memory-like CD8+" = "#4E8A5B",
  "ISG-high CD8+" = "#E2A23A"
)

theme_nm <- function(base_size = 7) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.line = element_line(linewidth = 0.28, colour = "black"),
      axis.ticks = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks.length = unit(1.5, "mm"),
      axis.text = element_text(colour = "black", size = base_size),
      axis.title = element_text(colour = "black", size = base_size + 0.5),
      plot.title = element_text(face = "bold", size = base_size + 1.5, hjust = 0),
      plot.subtitle = element_text(size = base_size, colour = "grey25"),
      legend.title = element_text(size = base_size),
      legend.text = element_text(size = base_size),
      legend.key.size = unit(3.2, "mm"),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      plot.margin = margin(3, 4, 3, 4)
    )
}

save_png <- function(plot, filename, width, height, dpi) {
  if (file.exists(filename)) {
    unlink(filename)
  }
  try(
    ggsave(filename, plot, width = width, height = height, units = "in", dpi = dpi),
    silent = TRUE
  )
  file.exists(filename) && file.size(filename) > 0
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

render_pdf_with_pymupdf <- function(pdf_file, png_file, dpi) {
  code <- paste(
    "import fitz, sys",
    "pdf, png, dpi = sys.argv[1], sys.argv[2], float(sys.argv[3])",
    "doc = fitz.open(pdf)",
    "pix = doc[0].get_pixmap(matrix=fitz.Matrix(dpi/72, dpi/72), alpha=False)",
    "pix.save(png)",
    sep = "; "
  )
  status <- system2("python3", c("-c", shQuote(code), shQuote(pdf_file), shQuote(png_file), dpi))
  status == 0 && file.exists(png_file) && file.size(png_file) > 0
}

paired_scores <- scores %>%
  group_by(patient, condition) %>%
  summarise(mean_score = mean(glycolysis_score, na.rm = TRUE), n_cells = n(), .groups = "drop")

paired_wide <- paired_scores %>%
  select(patient, condition, mean_score) %>%
  pivot_wider(names_from = condition, values_from = mean_score)
paired_t_p <- t.test(
  paired_wide$Graft,
  paired_wide$PBMC,
  paired = TRUE,
  alternative = "two.sided"
)$p.value
paired_wilcox_p <- suppressWarnings(wilcox.test(
  paired_wide$Graft,
  paired_wide$PBMC,
  paired = TRUE,
  alternative = "two.sided",
  exact = FALSE
)$p.value)
paired_p_label <- if (paired_t_p < 0.001) {
  "Overall paired Graft > PBMC: p < 0.001"
} else {
  paste0("Overall paired Graft > PBMC: p = ", sprintf("%.3f", paired_t_p))
}
paired_star <- if (paired_t_p < 0.0001) "****" else if (paired_t_p < 0.001) "***" else if (paired_t_p < 0.01) "**" else if (paired_t_p < 0.05) "*" else "ns"
paired_y_range <- diff(range(paired_scores$mean_score, na.rm = TRUE))
if (!is.finite(paired_y_range) || paired_y_range == 0) {
  paired_y_range <- max(abs(paired_scores$mean_score), na.rm = TRUE)
}
if (!is.finite(paired_y_range) || paired_y_range == 0) paired_y_range <- 1
paired_bracket_y <- max(paired_scores$mean_score, na.rm = TRUE) + 0.10 * paired_y_range
paired_star_y <- max(paired_scores$mean_score, na.rm = TRUE) + 0.15 * paired_y_range
paired_label_y <- max(paired_scores$mean_score, na.rm = TRUE) + 0.25 * paired_y_range
paired_label_short <- if (paired_t_p < 0.001) {
  paste0("p = ", sprintf("%.4f", paired_t_p))
} else {
  paste0("p = ", sprintf("%.4f", paired_t_p))
}
paired_label_short <- paste0(paired_label_short, "; paired two-sided t-test")
write_tsv(paired_scores %>% arrange(patient, condition), file.path(out_dir, "paired_glycolysis_exact_datapoints.tsv"))
write_tsv(paired_wide %>% arrange(patient), file.path(out_dir, "paired_glycolysis_exact_datapoints_wide.tsv"))
write_tsv(logcpm_source, file.path(out_dir, "panel_d_exact_pseudobulk_logcpm_datapoints.tsv"))
pseudobulk_metric <- logcpm_source %>%
  filter(gene %in% glycolysis_score_genes) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "pseudobulk_logCPM") %>%
  mutate(patient = sub("_(Graft|PBMC)$", "", sample_id),
         condition = sub("^.*_(Graft|PBMC)$", "\\1", sample_id)) %>%
  group_by(patient, condition) %>%
  summarise(mean_glycolysis_logCPM = mean(pseudobulk_logCPM, na.rm = TRUE), .groups = "drop")
write_tsv(pseudobulk_metric %>% arrange(patient, condition), file.path(out_dir, "panel_e_exact_pseudobulk_glycolysis_metric.tsv"))

state_scores <- scores %>%
  group_by(condition, celltype_granular) %>%
  summarise(mean_score = mean(glycolysis_score, na.rm = TRUE), n_cells = n(), .groups = "drop") %>%
  group_by(celltype_granular) %>%
  mutate(max_mean = max(mean_score)) %>%
  ungroup() %>%
  arrange(max_mean) %>%
  mutate(celltype_granular = factor(celltype_granular, levels = unique(celltype_granular)))

major_cd8_states <- scores %>%
  count(celltype_granular, name = "total_cells") %>%
  filter(total_cells >= 100) %>%
  pull(celltype_granular)

subtype_box_df <- scores %>%
  filter(celltype_granular %in% major_cd8_states) %>%
  group_by(celltype_granular) %>%
  mutate(state_median = median(glycolysis_score, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(celltype_granular = factor(celltype_granular, levels = unique(celltype_granular[order(state_median)])))

state_patient_scores <- scores %>%
  filter(celltype_granular %in% major_cd8_states) %>%
  group_by(celltype_granular, patient, condition) %>%
  summarise(mean_score = mean(glycolysis_score, na.rm = TRUE), n_cells = n(), .groups = "drop")

state_paired_tests <- state_patient_scores %>%
  select(celltype_granular, patient, condition, mean_score) %>%
  pivot_wider(names_from = condition, values_from = mean_score) %>%
  group_by(celltype_granular) %>%
  summarise(
    n_pairs = sum(!is.na(PBMC) & !is.na(Graft)),
    mean_delta_graft_minus_pbmc = mean(Graft - PBMC, na.rm = TRUE),
    graft_higher_in_all_available_pairs = all(Graft > PBMC, na.rm = TRUE),
    paired_t_p_two_sided = ifelse(
      n_pairs >= 3,
      suppressWarnings(t.test(Graft, PBMC, paired = TRUE, alternative = "two.sided")$p.value),
      NA_real_
    ),
    paired_wilcox_p_two_sided = ifelse(
      n_pairs >= 3,
      suppressWarnings(wilcox.test(Graft, PBMC, paired = TRUE, alternative = "two.sided", exact = FALSE)$p.value),
      NA_real_
    ),
    .groups = "drop"
  )

glycolysis_plot_genes <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
  "PGAM1", "ENO1", "PKM", "LDHA", "LDHB", "SLC2A1"
)
edge_plot_df <- edge %>%
  filter(gene %in% glycolysis_plot_genes) %>%
  mutate(
    neglog10_fdr = -log10(pmax(FDR, .Machine$double.xmin)),
    direction = case_when(
      FDR < 0.05 & logFC > 0 ~ "Up in graft, FDR < 0.05",
      FDR < 0.05 & logFC < 0 ~ "Down in graft, FDR < 0.05",
      TRUE ~ "Nominal or trend"
    )
  ) %>%
  arrange(logFC) %>%
  mutate(gene = factor(gene, levels = gene))

count_lines <- function(path) {
  length(readLines(path, warn = FALSE))
}

model_counts <- tibble(
  condition = rep(c("PBMC", "Graft"), each = 3),
  feature = rep(c("Genes", "Reactions", "Metabolites"), times = 2),
  count = c(
    count_lines(file.path(humess_model_dir, "PBMC", "genelist.tsv")),
    count_lines(file.path(humess_model_dir, "PBMC", "stats", "carveme.reactions.list")),
    count_lines(file.path(humess_model_dir, "PBMC", "stats", "carveme.metabolites.list")),
    count_lines(file.path(humess_model_dir, "Graft", "genelist.tsv")),
    count_lines(file.path(humess_model_dir, "Graft", "stats", "carveme.reactions.list")),
    count_lines(file.path(humess_model_dir, "Graft", "stats", "carveme.metabolites.list"))
  )
) %>%
  mutate(
    condition = factor(condition, levels = c("PBMC", "Graft")),
    feature = factor(feature, levels = c("Genes", "Reactions", "Metabolites")),
    feature_short = recode(feature, Genes = "Genes", Reactions = "Rxns", Metabolites = "Mets")
  )

reporter_plot_df <- reporter %>%
  mutate(
    pvalue = `P-VALUE`,
    neglog10_p = -log10(pmax(pvalue, .Machine$double.xmin)),
    label = case_when(
      ID == "f6p_c" ~ "Fructose 6-phosphate",
      ID == "fru_e" ~ "Fructose",
      TRUE ~ NAME
    ),
    highlight = if_else(ID %in% c("f6p_c", "fru_e"), "glycolysis/fructose node", "top reporter")
  ) %>%
  arrange(pvalue) %>%
  slice_head(n = 10) %>%
  mutate(label = str_trunc(label, width = 28), label = factor(label, levels = rev(label)))

umap_lims <- square_umap_limits(scores$UMAP_1, scores$UMAP_2)
scores_umap_plot <- scores %>%
  arrange(glycolysis_score)

celltype_umap_labels <- scores %>%
  filter(celltype_granular %in% major_cd8_states) %>%
  group_by(celltype_granular) %>%
  summarise(
    UMAP_1 = median(UMAP_1, na.rm = TRUE),
    UMAP_2 = median(UMAP_2, na.rm = TRUE),
    mean_score = mean(glycolysis_score, na.rm = TRUE),
    n_cells = n(),
    .groups = "drop"
  ) %>%
  mutate(
    short_label = recode(
      celltype_granular,
      "CXCL13+CXCR6+ effector CD8+" = "CXCL13+CXCR6+\neffector",
      "CXCR6+ effector CD8+" = "CXCR6+\neffector",
      "ZNF683+TCF7hi CD8+" = "ZNF683+TCF7hi",
      "CX3CR1+ effector CD8+" = "CX3CR1+\neffector",
      "Naive/Memory-like CD8+" = "Naive/Memory-like"
    ),
    score_group = case_when(
      mean_score == max(mean_score, na.rm = TRUE) ~ "highest",
      mean_score == min(mean_score, na.rm = TRUE) ~ "lowest",
      TRUE ~ "intermediate"
    ),
    overview_label = case_when(
      score_group == "highest" ~ "CXCL13+CXCR6+\neffector highest",
      score_group == "lowest" ~ "Naive/Memory-like\nlowest",
      TRUE ~ short_label
    ),
    standalone_label = paste0(short_label, "\nmean=", sprintf("%.2f", mean_score)),
    overview_label_x = case_when(
      score_group == "highest" ~ -4.1,
      score_group == "lowest" ~ 4.15,
      TRUE ~ UMAP_1
    ),
    overview_label_y = case_when(
      score_group == "highest" ~ -5.65,
      score_group == "lowest" ~ -2.2,
      TRUE ~ UMAP_2
    )
  )

celltype_umap_overview_labels <- celltype_umap_labels %>%
  filter(score_group %in% c("highest", "lowest"))

p_umap <- ggplot(scores_umap_plot, aes(UMAP_1, UMAP_2, colour = glycolysis_score)) +
  geom_point(size = 0.48, alpha = 0.85, stroke = 0) +
  scale_colour_gradientn(
    colours = c("#D8D8D8", "#F8D4CC", "#E75A4B", "firebrick3"),
    values = scales::rescale(c(0.55, 0.68, 0.84, 1.05), from = c(0.55, 1.05)),
    limits = c(0.55, 1.05),
    breaks = c(0.55, 0.70, 0.85, 1.05),
    labels = c("<=0.55", "0.70", "0.85", ">=1.05"),
    oob = scales::squish,
    name = "Score"
  ) +
  coord_fixed(xlim = umap_lims$xlim, ylim = umap_lims$ylim, expand = FALSE) +
  labs(title = "CD8 glycolysis score by state", x = "UMAP 1", y = "UMAP 2") +
  theme_nm() +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  ) +
  facet_wrap(~condition, nrow = 1)

p_umap_score_only <- ggplot(scores_umap_plot, aes(UMAP_1, UMAP_2, colour = glycolysis_score)) +
  geom_point(size = 0.48, alpha = 0.85, stroke = 0) +
  scale_colour_gradientn(
    colours = c("#D8D8D8", "#F8D4CC", "#E75A4B", "firebrick3"),
    values = scales::rescale(c(0.55, 0.68, 0.84, 1.05), from = c(0.55, 1.05)),
    limits = c(0.55, 1.05),
    breaks = c(0.55, 0.70, 0.85, 1.05),
    labels = c("<=0.55", "0.70", "0.85", ">=1.05"),
    oob = scales::squish,
    name = "Score"
  ) +
  coord_fixed(xlim = umap_lims$xlim, ylim = umap_lims$ylim, expand = FALSE) +
  labs(title = "CD8 glycolysis score", x = "UMAP 1", y = "UMAP 2") +
  theme_nm() +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

p_umap_celltype <- ggplot(scores, aes(UMAP_1, UMAP_2, colour = celltype_granular)) +
  geom_point(size = 0.42, alpha = 0.82, stroke = 0) +
  geom_label_repel(
    data = celltype_umap_labels,
    aes(UMAP_1, UMAP_2, label = standalone_label, fill = score_group),
    inherit.aes = FALSE,
    colour = "black",
    fontface = "bold",
    size = 2.5,
    label.size = 0.2,
    label.padding = unit(0.14, "lines"),
    min.segment.length = 0,
    segment.size = 0.25,
    segment.colour = "grey35",
    box.padding = 0.35,
    point.padding = 0.25,
    max.overlaps = Inf,
    seed = 7
  ) +
  scale_colour_manual(values = pal_celltype, name = NULL) +
  scale_fill_manual(
    values = c(highest = "#F6B7AE", intermediate = "white", lowest = "#DDE7F1"),
    guide = "none"
  ) +
  coord_fixed(xlim = umap_lims$xlim, ylim = umap_lims$ylim, expand = FALSE) +
  labs(title = "CD8 state annotation on UMAP", x = "UMAP 1", y = "UMAP 2") +
  theme_nm(base_size = 8) +
  guides(colour = guide_legend(override.aes = list(size = 2.6, alpha = 1))) +
  theme(
    aspect.ratio = 1,
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    legend.text = element_text(size = 6.5)
  )

p_paired <- ggplot(paired_scores, aes(condition, mean_score, group = patient)) +
  geom_line(colour = "grey55", linewidth = 0.35) +
  geom_point(aes(fill = condition), shape = 21, size = 2.1, stroke = 0.25, colour = "black") +
  annotate(
    "segment",
    x = 1,
    xend = 2,
    y = paired_bracket_y,
    yend = paired_bracket_y,
    linewidth = 0.3
  ) +
  annotate(
    "text",
    x = 1.5,
    y = paired_star_y,
    label = paired_star,
    fontface = "bold",
    size = 3.4,
    colour = "#B33A3A"
  ) +
  annotate(
    "label",
    x = 1.5,
    y = paired_label_y,
    label = paired_label_short,
    fill = "white",
    colour = "#B33A3A",
    fontface = "bold",
    size = 2.35
  ) +
  scale_fill_manual(values = pal_fill) +
  scale_x_discrete(expand = expansion(add = 0.25)) +
  scale_y_continuous(expand = expansion(mult = c(0.04, 0.30))) +
  labs(title = "Paired patient-level signal", subtitle = paired_p_label, x = NULL, y = "Mean glycolysis score") +
  theme_nm() +
  theme(legend.position = "none")

p_state_base <- ggplot(subtype_box_df, aes(glycolysis_score, celltype_granular, fill = condition)) +
  geom_boxplot(
    width = 0.56,
    outlier.shape = NA,
    linewidth = 0.26,
    position = position_dodge2(width = 0.62, preserve = "single")
  ) +
  scale_fill_manual(values = pal_fill, name = NULL) +
  labs(
    title = "Glycolysis score across CD8 states",
    subtitle = paired_p_label,
    x = "Glycolysis score",
    y = NULL
  ) +
  theme_nm() +
  theme(
    axis.text.y = element_text(size = 6.2),
    legend.position = "top",
    plot.margin = margin(3, 5, 3, 5)
  )

p_paired_no_subtitle <- p_paired +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())

p_state_base_no_subtitle <- p_state_base +
  labs(subtitle = NULL) +
  theme(plot.subtitle = element_blank())

p_edge <- ggplot(edge_plot_df, aes(logFC, gene)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_segment(aes(x = 0, xend = logFC, y = gene, yend = gene), linewidth = 0.35, colour = "grey35") +
  geom_point(aes(fill = direction, size = neglog10_fdr), shape = 21, colour = "black", stroke = 0.22) +
  scale_x_continuous(
    breaks = c(-0.6, 0, 0.6),
    expand = expansion(mult = c(0.02, 0.06))
  ) +
  coord_cartesian(xlim = c(-0.72, 0.62), clip = "off") +
  scale_fill_manual(
    values = c(
      "Up in graft, FDR < 0.05" = "#B33A3A",
      "Down in graft, FDR < 0.05" = "#3B5A7A",
      "Nominal or trend" = "#D9D9D9"
    ),
    name = NULL
  ) +
  scale_size_continuous(range = c(1.4, 3.5), name = "-log10 FDR") +
  labs(title = "Paired pseudobulk glycolysis genes", subtitle = "Metric: edgeR-normalized pseudobulk logCPM; exact patient values in source table", x = "log2FC (Graft / PBMC)", y = NULL) +
  theme_nm() +
  theme(
    axis.text.x = element_text(size = 6.7),
    axis.text.y = element_text(size = 6.5),
    legend.position = "right",
    legend.text = element_text(size = 5.9),
    legend.title = element_text(size = 6.2),
    legend.key.size = unit(2.8, "mm"),
    plot.margin = margin(3, 5, 3, 5)
  )

panel_e_labels <- pseudobulk_metric %>%
  mutate(
    x_base = as.numeric(factor(condition, levels = c("PBMC", "Graft"))),
    x_num = x_base + ifelse(condition == "PBMC", -0.12, 0.12),
    y_label = mean_glycolysis_logCPM + c("P1" = 0.10, "P2" = -0.10, "P3" = 0.10, "P4" = -0.10)[patient],
    hjust_label = ifelse(condition == "PBMC", 1, 0)
  )
panel_e_lines <- pseudobulk_metric %>%
  select(patient, condition, mean_glycolysis_logCPM) %>%
  pivot_wider(names_from = condition, values_from = mean_glycolysis_logCPM) %>%
  filter(!is.na(PBMC), !is.na(Graft))
p_model <- ggplot(panel_e_labels, aes(x_base, mean_glycolysis_logCPM, group = patient)) +
  geom_segment(data = panel_e_lines, aes(x = 1, xend = 2, y = PBMC, yend = Graft), inherit.aes = FALSE, colour = "grey55", linewidth = 0.4) +
  geom_point(aes(fill = condition), shape = 21, size = 2.8, colour = "black", stroke = 0.25) +
  geom_segment(data = panel_e_labels, aes(x = x_base, y = mean_glycolysis_logCPM, xend = x_num, yend = y_label), linewidth = 0.22, colour = "grey45", inherit.aes = FALSE) +
  geom_text(data = panel_e_labels, aes(x = x_num, y = y_label, label = sprintf("%.4f", mean_glycolysis_logCPM), hjust = hjust_label), size = 2.0, colour = "black", inherit.aes = FALSE) +
  scale_fill_manual(values = pal_fill, name = NULL) +
  scale_x_continuous(breaks = c(1, 2), labels = c("PBMC", "Graft"), limits = c(0.72, 2.28)) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.22))) +
  labs(title = NULL, x = NULL, y = "Mean glycolysis metric") +
  theme_nm() +
  theme(legend.position = "top")

p_reporter <- ggplot(reporter_plot_df, aes(neglog10_p, label)) +
  geom_segment(aes(x = 0, xend = neglog10_p, y = label, yend = label), linewidth = 0.3, colour = "grey55") +
  geom_point(aes(fill = highlight, size = `Z-SCORE`), shape = 21, colour = "black", stroke = 0.22) +
  scale_fill_manual(values = c("glycolysis/fructose node" = "#B33A3A", "top reporter" = "#D0D0D0"), name = NULL) +
  scale_size_continuous(range = c(1.8, 3.8), name = "Z-score") +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.08))) +
  labs(title = "Reporter metabolites", x = "-log10 P value", y = NULL) +
  theme_nm() +
  guides(
    size = guide_legend(order = 1, nrow = 1),
    fill = guide_legend(order = 2, nrow = 1)
  ) +
  theme(
    axis.text.y = element_text(size = 6.1),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.text = element_text(size = 5.8),
    legend.title = element_text(size = 6.0),
    legend.key.size = unit(2.8, "mm"),
    legend.margin = margin(0, 0, 0, 0),
    plot.margin = margin(3, 5, 3, 5)
  )

strip_panel_title <- function(plot) {
  plot +
    labs(title = NULL, subtitle = NULL) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )
}

p_umap_score_only_no_title <- strip_panel_title(p_umap_score_only)
p_paired_no_title <- strip_panel_title(p_paired)
p_state_base_no_title <- strip_panel_title(p_state_base)
p_edge_no_title <- strip_panel_title(p_edge)
p_model_no_title <- strip_panel_title(p_model)
p_reporter_no_title <- strip_panel_title(p_reporter)

column_one <- plot_grid(
  p_umap,
  p_edge,
  ncol = 1,
  labels = c("A", "D"),
  label_size = 9,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

column_two <- plot_grid(
  p_paired,
  p_model,
  ncol = 1,
  labels = c("B", "E"),
  label_size = 9,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

column_three <- plot_grid(
  p_state_base,
  p_reporter,
  ncol = 1,
  labels = c("C", "F"),
  label_size = 9,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

overview_grid <- plot_grid(
  column_one,
  column_two,
  column_three,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(1.55, 0.95, 2.15)
)

overview_core <- plot_grid(
  ggdraw() +
    draw_label(
      "Human CD8 metabolic profiling supports graft-associated glycolytic remodeling",
      x = 0.025,
      y = 0.58,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 10,
      fontfamily = "Helvetica"
    ),
  overview_grid,
  ncol = 1,
  rel_heights = c(0.075, 1)
)

overview <- plot_grid(
  NULL,
  overview_core,
  NULL,
  ncol = 3,
  rel_widths = c(0.02, 1, 0.02)
)

ggsave(
  filename = file.path(out_dir, "nm_humess_support_overview_v2.pdf"),
  plot = overview,
  width = overview_width,
  height = overview_height,
  units = "in",
  device = "pdf"
)

overview_pdf <- file.path(out_dir, "nm_humess_support_overview_v2.pdf")
overview_png <- file.path(out_dir, "nm_humess_support_overview_v2.png")
overview_preview <- file.path(out_dir, "nm_humess_support_overview_v2_preview.png")

if (!save_png(overview, overview_png, overview_width, overview_height, 450)) {
  invisible(render_pdf_with_pymupdf(overview_pdf, overview_png, 450))
}

column_one_no_subtitle <- plot_grid(
  p_umap_score_only_no_title,
  p_edge_no_title,
  ncol = 1,
  labels = c("A", "D"),
  label_size = 9,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

column_two_no_subtitle <- plot_grid(
  p_paired_no_title,
  p_model_no_title,
  ncol = 1,
  labels = c("B", "E"),
  label_size = 9,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

column_three_no_subtitle <- plot_grid(
  p_state_base_no_title,
  p_reporter_no_title,
  ncol = 1,
  labels = c("C", "F"),
  label_size = 9,
  label_fontface = "bold",
  label_fontfamily = "Helvetica",
  label_x = 0.01,
  label_y = 0.99,
  hjust = 0,
  vjust = 1,
  align = "v",
  axis = "lr",
  rel_heights = c(1, 1)
)

overview_grid_no_subtitle <- plot_grid(
  column_one_no_subtitle,
  column_two_no_subtitle,
  column_three_no_subtitle,
  nrow = 1,
  align = "h",
  axis = "tb",
  rel_widths = c(1.55, 0.95, 2.15)
)

overview_core_no_subtitle <- plot_grid(
  ggdraw() +
    draw_label(
      "Human CD8 metabolic profiling supports graft-associated glycolytic remodeling",
      x = 0.025,
      y = 0.58,
      hjust = 0,
      vjust = 0.5,
      fontface = "bold",
      size = 10,
      fontfamily = "Helvetica"
    ),
  overview_grid_no_subtitle,
  ncol = 1,
  rel_heights = c(0.075, 1)
)

overview_no_subtitle <- plot_grid(
  NULL,
  overview_core_no_subtitle,
  NULL,
  ncol = 3,
  rel_widths = c(0.02, 1, 0.02)
)

overview_no_subtitle_pdf <- file.path(out_dir, "nm_humess_support_overview_v2_no_subtitle.pdf")
overview_no_subtitle_png <- file.path(out_dir, "nm_humess_support_overview_v2_no_subtitle.png")
overview_no_subtitle_preview <- file.path(out_dir, "nm_humess_support_overview_v2_no_subtitle_preview.png")
overview_no_panel_titles_pdf <- file.path(out_dir, "nm_humess_support_overview_v2_no_panel_titles.pdf")
overview_no_panel_titles_png <- file.path(out_dir, "nm_humess_support_overview_v2_no_panel_titles.png")
overview_no_panel_titles_preview <- file.path(out_dir, "nm_humess_support_overview_v2_no_panel_titles_preview.png")

ggsave(
  filename = overview_no_subtitle_pdf,
  plot = overview_no_subtitle,
  width = overview_width,
  height = overview_height,
  units = "in",
  device = "pdf"
)

if (!save_png(overview_no_subtitle, overview_no_subtitle_png, overview_width, overview_height, 450)) {
  invisible(render_pdf_with_pymupdf(overview_no_subtitle_pdf, overview_no_subtitle_png, 450))
}

ggsave(
  filename = overview_no_panel_titles_pdf,
  plot = overview_no_subtitle,
  width = overview_width,
  height = overview_height,
  units = "in",
  device = "pdf"
)

if (!save_png(overview_no_subtitle, overview_no_panel_titles_png, overview_width, overview_height, 450)) {
  invisible(render_pdf_with_pymupdf(overview_no_panel_titles_pdf, overview_no_panel_titles_png, 450))
}

overview_caption <- paste(
  "Figure. Human CD8 metabolic profiling supports graft-associated glycolytic remodeling.",
  "(A) UMAP of human CD8 T cells colored by a glycolysis-associated transcriptional score,",
  "with major CD8 states annotated and highest/lowest mean-score states indicated;",
  "the score was computed as the mean log-normalized expression of available glycolysis/glucose-utilization genes.",
  "(B) Paired patient-level mean glycolysis scores in PBMC and graft CD8 T cells;",
  paste0("paired two-sided t-test, ", paired_label_short, "."),
  "(C) Distribution of single-cell glycolysis scores across major CD8 states, stratified by PBMC versus graft.",
  "(D) Paired pseudobulk differential expression of glycolysis-associated genes in graft versus PBMC CD8 T cells;",
  "point size indicates -log10 FDR and color indicates FDR-significant direction.",
  "(E) Patient-level HUMESS-supported glycolysis metric shown for four matched patients; values are the mean CD8 pseudobulk logCPM across the predefined glycolysis gene set, not direct metabolic flux.",
  "(F) Top reporter metabolites from the HUMESS graft-versus-PBMC comparison,",
  "with glycolysis/fructose-linked nodes highlighted. These data support a graft-associated metabolic",
  "transcriptional and model-level shift, rather than direct measurement of metabolic flux.",
  sep = " "
)
writeLines(
  str_wrap(overview_caption, width = 110),
  file.path(out_dir, "nm_humess_support_overview_v2_caption.txt")
)

overview_no_subtitle_caption <- paste(
  "Figure. Human CD8 metabolic profiling supports graft-associated glycolytic remodeling.",
  "(A) UMAP of human CD8 T cells colored by a glycolysis-associated transcriptional score;",
  "the score was computed as the mean log-normalized expression of available glycolysis/glucose-utilization genes.",
  "(B) Paired patient-level mean glycolysis scores in PBMC and graft CD8 T cells;",
  paste0("paired two-sided t-test, ", paired_label_short, "."),
  "(C) Distribution of single-cell glycolysis scores across major CD8 states, stratified by PBMC versus graft.",
  "(D) Paired pseudobulk differential expression of glycolysis-associated genes in graft versus PBMC CD8 T cells;",
  "point size indicates -log10 FDR and color indicates FDR-significant direction.",
  "(E) Patient-level HUMESS-supported glycolysis metric shown for four matched patients; values are the mean CD8 pseudobulk logCPM across the predefined glycolysis gene set, not direct metabolic flux.",
  "(F) Top reporter metabolites from the HUMESS graft-versus-PBMC comparison,",
  "with glycolysis/fructose-linked nodes highlighted. This version removes panel subtitles and panel titles,",
  "leaving the panel interiors cleaner for manual annotation.",
  sep = " "
)
writeLines(
  str_wrap(overview_no_subtitle_caption, width = 110),
  file.path(out_dir, "nm_humess_support_overview_v2_no_subtitle_caption.txt")
)
writeLines(
  str_wrap(overview_no_subtitle_caption, width = 110),
  file.path(out_dir, "nm_humess_support_overview_v2_no_panel_titles_caption.txt")
)

ggsave(
  filename = file.path(out_dir, "nm_glycolysis_by_cd8_state_boxplot.pdf"),
  plot = p_state_base,
  width = 5.4,
  height = 3.2,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, "nm_cd8_state_annotated_umap.pdf"),
  plot = p_umap_celltype,
  width = 5.6,
  height = 4.4,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, "nm_cd8_glycolysis_score_umap.pdf"),
  plot = p_umap_score_only,
  width = 5.0,
  height = 4.2,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, "nm_cd8_glycolysis_score_umap_no_title.pdf"),
  plot = p_umap_score_only_no_title,
  width = 5.0,
  height = 4.2,
  units = "in",
  device = "pdf"
)

state_boxplot_png <- file.path(out_dir, "nm_glycolysis_by_cd8_state_boxplot.png")
if (!save_png(p_state_base, state_boxplot_png, 5.4, 3.2, 450)) {
  invisible(render_pdf_with_pymupdf(file.path(out_dir, "nm_glycolysis_by_cd8_state_boxplot.pdf"), state_boxplot_png, 450))
}

celltype_umap_png <- file.path(out_dir, "nm_cd8_state_annotated_umap.png")
if (!save_png(p_umap_celltype, celltype_umap_png, 5.6, 4.4, 450)) {
  invisible(render_pdf_with_pymupdf(file.path(out_dir, "nm_cd8_state_annotated_umap.pdf"), celltype_umap_png, 450))
}

score_umap_png <- file.path(out_dir, "nm_cd8_glycolysis_score_umap.png")
if (!save_png(p_umap_score_only, score_umap_png, 5.0, 4.2, 450)) {
  invisible(render_pdf_with_pymupdf(file.path(out_dir, "nm_cd8_glycolysis_score_umap.pdf"), score_umap_png, 450))
}

score_umap_no_title_png <- file.path(out_dir, "nm_cd8_glycolysis_score_umap_no_title.png")
if (!save_png(p_umap_score_only_no_title, score_umap_no_title_png, 5.0, 4.2, 450)) {
  invisible(render_pdf_with_pymupdf(file.path(out_dir, "nm_cd8_glycolysis_score_umap_no_title.pdf"), score_umap_no_title_png, 450))
}

ggsave(
  filename = file.path(out_dir, "nm_paired_pseudobulk_glycolysis_genes_lollipop.pdf"),
  plot = p_edge,
  width = 5.1,
  height = 3.6,
  units = "in",
  device = "pdf"
)

ggsave(
  filename = file.path(out_dir, "nm_humess_reporter_metabolites_lollipop.pdf"),
  plot = p_reporter,
  width = 6.5,
  height = 3.6,
  units = "in",
  device = "pdf"
)

edge_png <- file.path(out_dir, "nm_paired_pseudobulk_glycolysis_genes_lollipop.png")
reporter_png <- file.path(out_dir, "nm_humess_reporter_metabolites_lollipop.png")
if (!save_png(p_edge, edge_png, 5.1, 3.6, 450)) {
  invisible(render_pdf_with_pymupdf(file.path(out_dir, "nm_paired_pseudobulk_glycolysis_genes_lollipop.pdf"), edge_png, 450))
}
if (!save_png(p_reporter, reporter_png, 6.5, 3.6, 450)) {
  invisible(render_pdf_with_pymupdf(file.path(out_dir, "nm_humess_reporter_metabolites_lollipop.pdf"), reporter_png, 450))
}

write_tsv(paired_scores, file.path(out_dir, "paired_glycolysis_patient_means.tsv"))
write_tsv(
  tibble(
    comparison = "Graft vs PBMC",
    n_pairs = nrow(paired_wide),
    mean_delta_graft_minus_pbmc = mean(paired_wide$Graft - paired_wide$PBMC, na.rm = TRUE),
    graft_higher_pairs = sum(paired_wide$Graft > paired_wide$PBMC, na.rm = TRUE),
    paired_t_p_two_sided = paired_t_p,
    paired_wilcox_p_two_sided = paired_wilcox_p
  ),
  file.path(out_dir, "paired_glycolysis_patient_level_test.tsv")
)
write_tsv(state_scores, file.path(out_dir, "cd8_state_glycolysis_means.tsv"))
write_tsv(
  celltype_umap_labels %>%
    mutate(
      across(c(short_label, overview_label, standalone_label), ~ str_replace_all(.x, "\n", " "))
    ),
  file.path(out_dir, "cd8_state_umap_annotation_labels.tsv")
)
write_tsv(state_paired_tests, file.path(out_dir, "cd8_state_glycolysis_patient_level_tests.tsv"))
write_tsv(
  subtype_box_df %>%
    count(condition, celltype_granular, name = "n_cells") %>%
    left_join(
      subtype_box_df %>%
        group_by(condition, celltype_granular) %>%
        summarise(
          median_score = median(glycolysis_score, na.rm = TRUE),
          mean_score = mean(glycolysis_score, na.rm = TRUE),
          .groups = "drop"
        ),
      by = c("condition", "celltype_granular")
    ),
  file.path(out_dir, "cd8_subtype_glycolysis_boxplot_summary.tsv")
)
write_tsv(model_counts, file.path(out_dir, "humess_model_feature_counts.tsv"))
write_tsv(reporter_plot_df, file.path(out_dir, "humess_top_reporter_metabolites_for_plot.tsv"))

glycolysis_gene_set_source <- tibble(
  gene_order = seq_along(glycolysis_score_genes),
  pathway_name = "Glycolysis-associated transcriptional score",
  source_database = "MSigDB Hallmark collection",
  source_gene_set = "HALLMARK_GLYCOLYSIS",
  source_url = "https://www.gsea-msigdb.org/gsea/msigdb/human/geneset/HALLMARK_GLYCOLYSIS.html",
  gene = glycolysis_score_genes,
  used_in_single_cell_score = gene %in% logcpm_source$gene,
  present_in_pseudobulk_logcpm = gene %in% logcpm_source$gene,
  present_in_pseudobulk_counts = gene %in% counts_source$gene,
  score_definition = "Mean log-normalized single-cell RNA expression across available genes",
  curation_rule = "Canonical glycolysis/glucose-utilization genes selected from HALLMARK_GLYCOLYSIS and present in the CD8 RNA assay",
  pathway_source_note = "This is a transcriptional score, not a direct metabolic-flux measurement",
  citation_note = "Subramanian et al., PNAS 2005; Liberzon et al., Bioinformatics 2011; Liberzon et al., Cell Systems 2015"
)

glycolysis_logcpm_wide <- logcpm_source %>%
  filter(gene %in% glycolysis_score_genes) %>%
  mutate(gene_order = match(gene, glycolysis_score_genes)) %>%
  arrange(gene_order) %>%
  select(gene_order, gene, everything())

glycolysis_counts_wide <- counts_source %>%
  filter(gene %in% glycolysis_score_genes) %>%
  mutate(gene_order = match(gene, glycolysis_score_genes)) %>%
  arrange(gene_order) %>%
  select(gene_order, gene, everything())

glycolysis_logcpm_long <- glycolysis_logcpm_wide %>%
  select(-gene_order) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "pseudobulk_logCPM")

glycolysis_counts_long <- glycolysis_counts_wide %>%
  select(-gene_order) %>%
  pivot_longer(-gene, names_to = "sample_id", values_to = "pseudobulk_UMI_count")

glycolysis_expression_source <- glycolysis_logcpm_long %>%
  full_join(glycolysis_counts_long, by = c("gene", "sample_id")) %>%
  mutate(
    gene_order = match(gene, glycolysis_score_genes),
    patient = str_replace(sample_id, "_.*$", ""),
    condition = str_replace(sample_id, "^.*_", ""),
    condition = factor(condition, levels = c("PBMC", "Graft"))
  ) %>%
  arrange(patient, condition, gene_order) %>%
  select(gene_order, gene, sample_id, patient, condition, pseudobulk_logCPM, pseudobulk_UMI_count)

write_tsv(glycolysis_gene_set_source, file.path(out_dir, "glycolysis_score_gene_set.tsv"))
write_tsv(glycolysis_logcpm_wide, file.path(out_dir, "glycolysis_score_gene_expression_logcpm_wide.tsv"))
write_tsv(glycolysis_expression_source, file.path(out_dir, "glycolysis_score_gene_expression_by_sample.tsv"))

if (!save_png(overview, overview_preview, overview_width, overview_height, 150)) {
  invisible(render_pdf_with_pymupdf(overview_pdf, overview_preview, 150))
}

if (!save_png(overview_no_subtitle, overview_no_subtitle_preview, overview_width, overview_height, 150)) {
  invisible(render_pdf_with_pymupdf(overview_no_subtitle_pdf, overview_no_subtitle_preview, 150))
}

if (!save_png(overview_no_subtitle, overview_no_panel_titles_preview, overview_width, overview_height, 150)) {
  invisible(render_pdf_with_pymupdf(overview_no_panel_titles_pdf, overview_no_panel_titles_preview, 150))
}

message("Wrote: ", overview_pdf)
message("Wrote: ", overview_no_subtitle_pdf)
message("Wrote: ", overview_no_panel_titles_pdf)
