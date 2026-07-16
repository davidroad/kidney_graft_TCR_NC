#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(repo_lib)) .libPaths(c(repo_lib, .libPaths()))

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(stringr)
  library(patchwork)
  library(cowplot)
})

root_dir <- normalizePath(getwd(), mustWork = TRUE)
out_dir <- file.path(root_dir, "results", "human_cd8_metabolic_outputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
overview_width <- 9.2
overview_height <- 6.1

scores_file <- file.path(root_dir, "results", "human_cd8_metabolic_analysis", "input", "cd8_metabolic_scores.tsv.gz")
logcpm_file <- file.path(root_dir, "results", "human_cd8_metabolic_analysis", "input", "pseudobulk_logcpm_by_sample.tsv.gz")
counts_file <- file.path(root_dir, "results", "human_cd8_metabolic_analysis", "input", "pseudobulk_counts_by_sample.tsv.gz")
edge_file <- file.path(root_dir, "results", "human_cd8_metabolic_analysis", "edgeR", "paired_graft_vs_pbmc_edgeR.tsv")
reporter_file <- file.path(
  root_dir,
  "results", "human_cd8_metabolic_analysis", "humess", "humess_run",
  "comparisons", "PBMC__vs__Graft", "ReporterMetabolites",
  "reporter_metabolites_Graft_vs_PBMC_all.tsv"
)
humess_model_dir <- file.path(root_dir, "results", "human_cd8_metabolic_analysis", "humess", "humess_run", "models")

glycolysis_score_genes <- c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
  "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB", "SLC2A1",
  "G6PC3", "AKR1A1"
)

scores <- read_tsv(scores_file, show_col_types = FALSE) %>%
  mutate(
    condition = factor(condition, levels = c("PBMC", "Graft")),
    patient = factor(patient),
    stem_like_annotation = if_else(stem_like_annotation, "TCF7hi stem-like", "Other CD8")
  )

edge <- read_tsv(edge_file, show_col_types = FALSE)
reporter <- read_tsv(reporter_file, show_col_types = FALSE)
logcpm_source <- read_tsv(logcpm_file, show_col_types = FALSE)
counts_source <- read_tsv(counts_file, show_col_types = FALSE)

pal_condition <- c(PBMC = "#3B5A7A", Graft = "#B33A3A")
pal_fill <- c(PBMC = "#C9D5E3", Graft = "#E7B4AE")

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

panel_with_label <- function(plot, tag) {
  ggdraw(plot) +
    draw_label(
      tag,
      x = 0.025,
      y = 0.965,
      hjust = 0,
      vjust = 1,
      fontface = "bold",
      size = 9,
      fontfamily = "Helvetica"
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
  "p < 0.001"
} else {
  paste0("p = ", sprintf("%.3f", paired_t_p))
}
if (paired_t_p < 0.05) paired_label_short <- paste0(paired_label_short, " (*)")

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

p_umap <- ggplot(scores, aes(UMAP_1, UMAP_2, colour = glycolysis_score)) +
  geom_point(size = 0.12, alpha = 0.85, stroke = 0) +
  scale_colour_gradientn(
    colours = c("#222222", "#7EA6C8", "#F0C36D", "#B33A3A"),
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
  theme(axis.text.y = element_text(size = 6.2), legend.position = "top")

p_edge <- ggplot(edge_plot_df, aes(logFC, gene)) +
  geom_vline(xintercept = 0, linewidth = 0.25, colour = "grey60") +
  geom_segment(aes(x = 0, xend = logFC, y = gene, yend = gene), linewidth = 0.35, colour = "grey35") +
  geom_point(aes(fill = direction, size = neglog10_fdr), shape = 21, colour = "black", stroke = 0.22) +
  scale_fill_manual(
    values = c(
      "Up in graft, FDR < 0.05" = "#B33A3A",
      "Down in graft, FDR < 0.05" = "#3B5A7A",
      "Nominal or trend" = "#D9D9D9"
    ),
    name = NULL
  ) +
  scale_size_continuous(range = c(1.4, 3.5), name = "-log10 FDR") +
  labs(title = "Paired pseudobulk glycolysis genes", x = "log2FC (Graft / PBMC)", y = NULL) +
  theme_nm() +
  theme(axis.text.y = element_text(size = 6.5), legend.position = "right")

p_model <- ggplot(model_counts, aes(feature_short, count, group = condition)) +
  geom_line(aes(colour = condition), linewidth = 0.35, position = position_dodge(width = 0.35)) +
  geom_point(aes(fill = condition), shape = 21, size = 2.6, colour = "black", stroke = 0.25, position = position_dodge(width = 0.35)) +
  scale_colour_manual(values = pal_condition, guide = "none") +
  scale_fill_manual(values = pal_fill, name = NULL) +
  scale_y_continuous(labels = scales::comma) +
  labs(title = "HUMESS model features", x = NULL, y = "Feature count") +
  theme_nm() +
  theme(legend.position = "top")

p_reporter <- ggplot(reporter_plot_df, aes(neglog10_p, label)) +
  geom_segment(aes(x = 0, xend = neglog10_p, y = label, yend = label), linewidth = 0.3, colour = "grey55") +
  geom_point(aes(fill = highlight, size = `Z-SCORE`), shape = 21, colour = "black", stroke = 0.22) +
  scale_fill_manual(values = c("glycolysis/fructose node" = "#B33A3A", "top reporter" = "#D0D0D0"), name = NULL) +
  scale_size_continuous(range = c(1.8, 3.8), name = "Z-score") +
  labs(title = "Reporter metabolites", x = "-log10 P value", y = NULL) +
  theme_nm() +
  theme(axis.text.y = element_text(size = 6.1), legend.position = "right")

overview_grid <- plot_grid(
  panel_with_label(p_umap, "A"),
  panel_with_label(p_paired, "B"),
  panel_with_label(p_state_base, "C"),
  panel_with_label(p_edge, "D"),
  panel_with_label(p_model, "E"),
  panel_with_label(p_reporter, "F"),
  nrow = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1.0, 0.85, 1.45),
  rel_heights = c(1, 1.05)
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

ggsave(
  filename = file.path(out_dir, "nm_glycolysis_by_cd8_state_boxplot.pdf"),
  plot = p_state_base,
  width = 4.8,
  height = 3.0,
  units = "in",
  device = "pdf"
)

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

message("Wrote: ", overview_pdf)
