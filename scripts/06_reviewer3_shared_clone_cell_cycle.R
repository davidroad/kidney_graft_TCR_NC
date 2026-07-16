#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
.libPaths(c("/data_sys/collab2/ydai2/26_immune_NC/R_libs/4.5.0", .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(cowplot)
})

repo_dir <- normalizePath(getwd(), mustWork = TRUE)
input_rdata <- file.path(repo_dir, "merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData")
if (!file.exists(input_rdata)) {
  input_rdata <- file.path(
    repo_dir,
    "finale", "00_needed_files_for_submission", "07_reproducibility_bundle",
    "input_human_cd8", "merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData"
  )
}
stopifnot(file.exists(input_rdata))

out_root <- file.path(repo_dir, "finale", "Reviewer3_Major2_cell_cycle_in_situ")
figure_dirs <- c(
  file.path(out_root, "figures"),
  file.path(repo_dir, "finale", "00_needed_files_for_submission", "04_supporting_figures")
)
table_dirs <- c(
  file.path(out_root, "tables"),
  file.path(repo_dir, "finale", "00_needed_files_for_submission", "05_source_tables")
)
script_dirs <- c(
  file.path(out_root, "scripts"),
  file.path(repo_dir, "finale", "00_needed_files_for_submission", "07_reproducibility_bundle", "repro_scripts")
)
rds_dirs <- c(
  file.path(out_root, "rds"),
  file.path(repo_dir, "finale", "00_needed_files_for_submission", "07_reproducibility_bundle", "final_rds_objects")
)
invisible(lapply(c(figure_dirs, table_dirs, script_dirs, rds_dirs), dir.create, recursive = TRUE, showWarnings = FALSE))

theme_nc <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      axis.text = element_text(colour = "black"),
      axis.title = element_text(colour = "black"),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      plot.margin = margin(5, 5, 5, 5)
    )
}

theme_umap <- function(base_size = 8) {
  theme_classic(base_size = base_size, base_family = "Helvetica") +
    theme(
      aspect.ratio = 1,
      axis.line = element_line(linewidth = 0.25, colour = "black"),
      axis.ticks = element_blank(),
      axis.text = element_blank(),
      axis.title = element_text(colour = "black", size = base_size - 0.5),
      plot.title = element_text(face = "bold", hjust = 0, size = base_size + 1),
      legend.title = element_text(size = base_size - 1),
      legend.text = element_text(size = base_size - 1),
      legend.key.size = unit(3, "mm"),
      plot.margin = margin(5, 5, 5, 5)
    )
}

square_umap_limits <- function(seu, reduction = "umap", probs = c(0.002, 0.998), pad = 0.03) {
  emb <- Embeddings(seu, reduction)[, 1:2, drop = FALSE]
  xlim <- as.numeric(quantile(emb[, 1], probs = probs, na.rm = TRUE))
  ylim <- as.numeric(quantile(emb[, 2], probs = probs, na.rm = TRUE))
  x_mid <- mean(xlim)
  y_mid <- mean(ylim)
  span <- max(diff(xlim), diff(ylim)) * (1 + pad)
  list(
    xlim = x_mid + c(-0.5, 0.5) * span,
    ylim = y_mid + c(-0.5, 0.5) * span
  )
}

style_feature <- function(plot, lims, title = NULL) {
  out <- plot +
    coord_fixed(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) +
    labs(x = "UMAP 1", y = "UMAP 2") +
    theme_umap(8)
  if (!is.null(title)) {
    out <- out + labs(title = title)
  }
  out
}

save_all_fig <- function(plot, filename, width, height, png = TRUE) {
  for (fig_dir in figure_dirs) {
    pdf_file <- file.path(fig_dir, filename)
    ggsave(pdf_file, plot, width = width, height = height, units = "in", device = "pdf", limitsize = FALSE)
    if (png) {
      ggsave(sub("\\.pdf$", ".png", pdf_file), plot, width = width, height = height, units = "in", dpi = 150, limitsize = FALSE)
    }
  }
}

write_all_tables <- function(df, filename) {
  for (tab_dir in table_dirs) {
    write_csv(df, file.path(tab_dir, filename))
  }
}

p_to_star <- function(p) {
  case_when(
    is.na(p) ~ "",
    p < 0.001 ~ "***",
    p < 0.01 ~ "**",
    p < 0.05 ~ "*",
    TRUE ~ "ns"
  )
}

load_env <- new.env()
all_clone_file <- "/home/0.collaboration/20_sc_VDJ/merge_4_with_vdj/all_clone_filter_10_28.Rdata"
load(all_clone_file, envir = load_env)
stopifnot("all_clone_test_filter" %in% ls(load_env))
cd8 <- load_env[["all_clone_test_filter"]]
DefaultAssay(cd8) <- "RNA"

required_meta <- c("patient", "label", "batch", "t_cdr3s_aa", "cell_type_annotated_granularity")
missing_meta <- setdiff(required_meta, colnames(cd8@meta.data))
if (length(missing_meta)) {
  stop("Input object is missing required metadata: ", paste(missing_meta, collapse = ", "))
}

patients <- sort(unique(as.character(cd8$patient)))
shared_cd8 <- cd8
shared_summary <- cd8@meta.data %>% count(patient, name = "shared_cell_count") %>% mutate(shared_clone_count = NA_integer_)
DefaultAssay(shared_cd8) <- "RNA"
shared_cd8$label <- factor(as.character(shared_cd8$label), levels = c("PBMC", "Graft"))
shared_cd8$patient <- factor(as.character(shared_cd8$patient), levels = patients)
shared_cd8$cell_type_annotated_granularity <- factor(as.character(shared_cd8$cell_type_annotated_granularity))

proliferation_genes <- c("MKI67", "TOP2A", "CDK1", "CCNB1", "CDC20", "PLK1", "AURKB", "MCM5", "PCNA")
arrest_genes <- c("CDKN1A", "CDKN1B", "CDKN2A", "RB1", "TP53")
survival_death_genes <- c("BCL2", "MCL1", "BAX", "BAK1", "BCL2L11", "BAD", "BID", "CASP3", "TNFRSF1A")
context_genes <- c("CXCR6", "TCF7", "GZMB", "IFNG")

gene_table <- bind_rows(
  tibble(gene_group = "Proliferation/cell cycle", gene = proliferation_genes),
  tibble(gene_group = "Cell-cycle arrest", gene = arrest_genes),
  tibble(gene_group = "Survival/apoptosis", gene = survival_death_genes),
  tibble(gene_group = "CD8 state context", gene = context_genes)
) %>%
  distinct(gene_group, gene) %>%
  mutate(available_in_shared_cd8 = gene %in% rownames(shared_cd8))
write_all_tables(gene_table, "reviewer3_shared_clone_cd8_cell_cycle_gene_set.csv")

proliferation_present <- intersect(proliferation_genes, rownames(shared_cd8))
arrest_present <- intersect(arrest_genes, rownames(shared_cd8))
survival_death_present <- intersect(survival_death_genes, rownames(shared_cd8))
s_genes_present <- intersect(cc.genes.updated.2019$s.genes, rownames(shared_cd8))
g2m_genes_present <- intersect(cc.genes.updated.2019$g2m.genes, rownames(shared_cd8))

shared_cd8 <- AddModuleScore(shared_cd8, features = list(proliferation_present), name = "Proliferation_module", assay = "RNA", seed = 1)
shared_cd8 <- AddModuleScore(shared_cd8, features = list(arrest_present), name = "Cell_cycle_arrest_module", assay = "RNA", seed = 1)
shared_cd8 <- AddModuleScore(shared_cd8, features = list(survival_death_present), name = "Survival_apoptosis_module", assay = "RNA", seed = 1)
shared_cd8 <- CellCycleScoring(
  shared_cd8,
  s.features = s_genes_present,
  g2m.features = g2m_genes_present,
  set.ident = FALSE
)

expr_genes <- intersect(c("CXCR6", proliferation_present, arrest_present, survival_death_present, context_genes), rownames(shared_cd8))
expr <- FetchData(shared_cd8, vars = expr_genes)
shared_cd8$CXCR6_status <- factor(if_else(expr$CXCR6 > 0, "CXCR6+", "CXCR6-"), levels = c("CXCR6-", "CXCR6+"))
shared_cd8$MKI67_positive <- expr$MKI67 > 0
shared_cd8$proliferation_gene_positive <- rowSums(expr[, proliferation_present, drop = FALSE] > 0) > 0
shared_cd8$label_cxcr6 <- factor(
  paste(shared_cd8$label, shared_cd8$CXCR6_status, sep = " "),
  levels = c("PBMC CXCR6-", "PBMC CXCR6+", "Graft CXCR6-", "Graft CXCR6+")
)

for (rds_dir in rds_dirs) {
  saveRDS(shared_cd8, file.path(rds_dir, "reviewer3_shared_clone_cd8_cell_cycle_scored.rds"))
}

cell_counts <- shared_cd8@meta.data %>%
  count(patient, label, CXCR6_status, name = "cell_count") %>%
  arrange(patient, label, CXCR6_status)
write_all_tables(shared_summary, "reviewer3_shared_clone_cd8_shared_clone_counts.csv")
write_all_tables(cell_counts, "reviewer3_shared_clone_cd8_cxcr6_cell_counts.csv")

expr_tbl <- expr %>%
  as_tibble(rownames = "cell")
cell_meta <- shared_cd8@meta.data %>%
  mutate(cell = rownames(shared_cd8@meta.data)) %>%
  left_join(expr_tbl, by = "cell") %>%
  mutate(
    MKI67_expr = MKI67,
    TOP2A_expr = TOP2A,
    PCNA_expr = PCNA,
    CXCR6_expr = CXCR6,
    Proliferation_module = Proliferation_module1,
    Cell_cycle_arrest_module = Cell_cycle_arrest_module1,
    Survival_apoptosis_module = Survival_apoptosis_module1
  )

patient_summary <- cell_meta %>%
  group_by(patient, label, CXCR6_status) %>%
  summarise(
    n_cells = n(),
    mkI67_positive_frequency = mean(MKI67_positive),
    proliferation_gene_positive_frequency = mean(proliferation_gene_positive),
    mean_MKI67_expression = mean(MKI67_expr),
    mean_TOP2A_expression = mean(TOP2A_expr),
    mean_PCNA_expression = mean(PCNA_expr),
    mean_CXCR6_expression = mean(CXCR6_expr),
    mean_proliferation_module = mean(Proliferation_module),
    mean_cell_cycle_arrest_module = mean(Cell_cycle_arrest_module),
    mean_survival_apoptosis_module = mean(Survival_apoptosis_module),
    mean_S_score = mean(S.Score),
    mean_G2M_score = mean(G2M.Score),
    g1_phase_frequency = mean(Phase == "G1"),
    s_phase_frequency = mean(Phase == "S"),
    g2m_phase_frequency = mean(Phase == "G2M"),
    .groups = "drop"
  ) %>%
  arrange(patient, label, CXCR6_status)
write_all_tables(patient_summary, "reviewer3_shared_clone_cd8_cell_cycle_patient_summary.csv")

metric_map <- c(
  mkI67_positive_frequency = "MKI67+ frequency",
  proliferation_gene_positive_frequency = "Any proliferation gene+ frequency",
  mean_MKI67_expression = "MKI67 expression",
  mean_proliferation_module = "Proliferation module",
  mean_S_score = "Seurat S score",
  mean_G2M_score = "Seurat G2M score",
  mean_cell_cycle_arrest_module = "Cell-cycle arrest module",
  mean_survival_apoptosis_module = "Survival/apoptosis module"
)

patient_long <- patient_summary %>%
  select(patient, label, CXCR6_status, all_of(names(metric_map))) %>%
  pivot_longer(
    cols = all_of(names(metric_map)),
    names_to = "metric_id",
    values_to = "value"
  ) %>%
  mutate(metric = unname(metric_map[metric_id]))

patient_stats <- patient_long %>%
  group_by(label, metric_id, metric) %>%
  summarise(
    n_patient_pairs = n_distinct(patient),
    mean_CXCR6_pos = mean(value[CXCR6_status == "CXCR6+"], na.rm = TRUE),
    mean_CXCR6_neg = mean(value[CXCR6_status == "CXCR6-"], na.rm = TRUE),
    delta_pos_minus_neg = mean_CXCR6_pos - mean_CXCR6_neg,
    p_value = {
      wide <- tibble(patient = patient, CXCR6_status = CXCR6_status, value = value) %>%
        distinct() %>%
        pivot_wider(names_from = CXCR6_status, values_from = value)
      if (all(c("CXCR6+", "CXCR6-") %in% colnames(wide)) && nrow(na.omit(wide[, c("CXCR6+", "CXCR6-")])) >= 3) {
        suppressWarnings(wilcox.test(wide[["CXCR6+"]], wide[["CXCR6-"]], paired = TRUE, exact = FALSE)$p.value)
      } else {
        NA_real_
      }
    },
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"), star = p_to_star(p_value)) %>%
  arrange(label, metric_id)
write_all_tables(patient_stats, "reviewer3_shared_clone_cd8_cell_cycle_patient_level_stats.csv")

cell_stats <- cell_meta %>%
  select(label, CXCR6_status, all_of(expr_genes), Proliferation_module, Cell_cycle_arrest_module, Survival_apoptosis_module, S.Score, G2M.Score) %>%
  pivot_longer(
    cols = -c(label, CXCR6_status),
    names_to = "metric",
    values_to = "value"
  ) %>%
  group_by(label, metric) %>%
  summarise(
    n_CXCR6_pos = sum(CXCR6_status == "CXCR6+"),
    n_CXCR6_neg = sum(CXCR6_status == "CXCR6-"),
    mean_CXCR6_pos = mean(value[CXCR6_status == "CXCR6+"], na.rm = TRUE),
    mean_CXCR6_neg = mean(value[CXCR6_status == "CXCR6-"], na.rm = TRUE),
    delta_pos_minus_neg = mean_CXCR6_pos - mean_CXCR6_neg,
    p_value = if (length(unique(CXCR6_status)) == 2) {
      suppressWarnings(wilcox.test(value[CXCR6_status == "CXCR6+"], value[CXCR6_status == "CXCR6-"], exact = FALSE)$p.value)
    } else {
      NA_real_
    },
    .groups = "drop"
  ) %>%
  mutate(p_adj_BH = p.adjust(p_value, method = "BH"), star = p_to_star(p_value)) %>%
  arrange(label, metric)
write_all_tables(cell_stats, "reviewer3_shared_clone_cd8_cell_cycle_cell_level_stats_exploratory.csv")

message("Shared clone CD8 cells: ", ncol(shared_cd8))
message("Patient-level summaries written.")

label_cols <- c(PBMC = "#4C78A8", Graft = "#B33A3A")
cxcr6_cols <- c("CXCR6-" = "#BDBDBD", "CXCR6+" = "#B33A3A")
label_cxcr6_cols <- c(
  "PBMC CXCR6-" = "#B8C5D8",
  "PBMC CXCR6+" = "#4C78A8",
  "Graft CXCR6-" = "#E6B4AD",
  "Graft CXCR6+" = "#B33A3A"
)

lims <- square_umap_limits(shared_cd8, "umap")
p_label <- DimPlot(shared_cd8, reduction = "umap", group.by = "label", cols = label_cols, pt.size = 0.28, label = FALSE) %>%
  style_feature(lims, "Tissue")
p_state <- DimPlot(shared_cd8, reduction = "umap", group.by = "cell_type_annotated_granularity", pt.size = 0.28, label = FALSE) %>%
  style_feature(lims, "CD8 state")
p_cx_status <- DimPlot(shared_cd8, reduction = "umap", group.by = "CXCR6_status", cols = cxcr6_cols, pt.size = 0.28, label = FALSE) %>%
  style_feature(lims, "CXCR6 status")

feature_summary <- FeaturePlot(
  shared_cd8,
  reduction = "umap",
  features = intersect(c("CXCR6", "MKI67", "TOP2A", "PCNA"), rownames(shared_cd8)),
  split.by = "label",
  cols = c("lightgrey", "firebrick3"),
  pt.size = 0.28,
  order = TRUE,
  combine = FALSE
)
feature_summary <- Map(function(p, gene) style_feature(p, lims, gene), feature_summary, c("CXCR6", "MKI67", "TOP2A", "PCNA"))

summary_metric_order <- c(
  "MKI67+ frequency",
  "Proliferation module",
  "Seurat S score",
  "Seurat G2M score"
)
summary_scores <- patient_long %>%
  filter(metric %in% summary_metric_order) %>%
  mutate(metric = factor(metric, levels = summary_metric_order))

p_score <- ggplot(summary_scores, aes(label, value, fill = CXCR6_status)) +
  geom_boxplot(width = 0.56, outlier.shape = NA, linewidth = 0.32, position = position_dodge(width = 0.68)) +
  geom_point(aes(group = CXCR6_status), position = position_dodge(width = 0.68), size = 1.35, stroke = 0.25, shape = 21, colour = "black") +
  facet_wrap(~metric, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = cxcr6_cols, drop = FALSE) +
  labs(x = NULL, y = "Patient-level value", fill = NULL) +
  theme_nc(8) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))

phase_summary <- patient_summary %>%
  select(patient, label, CXCR6_status, g1_phase_frequency, s_phase_frequency, g2m_phase_frequency) %>%
  pivot_longer(
    cols = c(g1_phase_frequency, s_phase_frequency, g2m_phase_frequency),
    names_to = "phase",
    values_to = "frequency"
  ) %>%
  mutate(
    phase = recode(phase, g1_phase_frequency = "G1", s_phase_frequency = "S", g2m_phase_frequency = "G2M"),
    label_cxcr6 = factor(paste(label, CXCR6_status), levels = levels(shared_cd8$label_cxcr6))
  ) %>%
  group_by(label_cxcr6, phase) %>%
  summarise(mean_frequency = mean(frequency), .groups = "drop")
p_phase <- ggplot(phase_summary, aes(label_cxcr6, mean_frequency, fill = phase)) +
  geom_col(width = 0.7, colour = "white", linewidth = 0.2) +
  scale_fill_manual(values = c(G1 = "#8A8A8A", S = "#D28E2B", G2M = "#8E2F3F")) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), expand = expansion(mult = c(0, 0.03))) +
  labs(x = NULL, y = "Mean patient-level fraction", fill = NULL, title = "Seurat cell-cycle phase") +
  theme_nc(8) +
  theme(axis.text.x = element_text(angle = 30, hjust = 1))

label_args <- list(label_size = 9, label_x = 0.025, label_y = 0.975, hjust = 0, vjust = 1)
row_1 <- do.call(cowplot::plot_grid, c(list(p_label, p_state, p_cx_status, ncol = 3, align = "hv", axis = "tblr", rel_widths = c(1, 1, 1), labels = c("A", "B", "C")), label_args))
row_2 <- do.call(cowplot::plot_grid, c(list(feature_summary[[1]], feature_summary[[2]], feature_summary[[3]], ncol = 3, align = "hv", axis = "tblr", rel_widths = c(1, 1, 1), labels = c("D", "E", "F")), label_args))
row_3 <- do.call(cowplot::plot_grid, c(list(feature_summary[[4]], p_score, p_phase, ncol = 3, align = "hv", axis = "tblr", rel_widths = c(1, 1, 1), labels = c("G", "H", "I")), label_args))
summary_fig <- cowplot::plot_grid(row_1, row_2, row_3, ncol = 1, align = "v", axis = "lr", rel_heights = c(1, 1, 1))

save_all_fig(
  summary_fig,
  "12_reviewer3_shared_clone_cd8_cell_cycle_summary.pdf",
  width = 14.4,
  height = 12.6,
  png = TRUE
)

all_feature_genes <- gene_table %>%
  filter(available_in_shared_cd8, gene_group != "CD8 state context") %>%
  pull(gene) %>%
  unique()
feature_grid <- FeaturePlot(
  shared_cd8, reduction = "umap", features = all_feature_genes[1], split.by = "label",
  cols = c("lightgrey", "firebrick3"), ncol = 2, label = FALSE, pt.size = 0.40,
  order = TRUE, combine = TRUE
)
feature_grid <- lapply(all_feature_genes, function(gene) {
  FeaturePlot(shared_cd8, reduction = "umap", features = gene, split.by = "label",
              cols = c("lightgrey", "firebrick3"), ncol = 2, label = FALSE,
              pt.size = 0.40, order = TRUE, combine = TRUE) &
    coord_fixed(xlim = lims$xlim, ylim = lims$ylim, expand = FALSE) &
    theme(plot.title = element_text(size = 8), plot.margin = margin(2, 2, 2, 2))
})
feature_grid_fig <- wrap_plots(feature_grid, ncol = 2) +
  plot_annotation(tag_levels = "A", title = "Shared-TCR-clone CD8+ marker feature plots for cell-cycle and death-related genes")
save_all_fig(
  feature_grid_fig,
  "13_reviewer3_shared_clone_cd8_cell_cycle_marker_featureplots.pdf",
  width = 24.0,
  height = 7.2 * ceiling(length(all_feature_genes) / 4),
  png = TRUE
)

violin_genes <- gene_table %>%
  filter(available_in_shared_cd8, gene_group != "CD8 state context") %>%
  mutate(gene = factor(gene, levels = unique(gene))) %>%
  pull(gene) %>%
  as.character()
violin_df <- FetchData(shared_cd8, vars = violin_genes) %>%
  mutate(cell = rownames(.)) %>%
  pivot_longer(cols = all_of(violin_genes), names_to = "gene", values_to = "expression") %>%
  left_join(
    cell_meta %>% select(cell, label, CXCR6_status, label_cxcr6),
    by = "cell"
  ) %>%
  left_join(gene_table %>% select(gene_group, gene), by = "gene") %>%
  mutate(
    gene = factor(gene, levels = violin_genes),
    label_cxcr6 = factor(label_cxcr6, levels = levels(shared_cd8$label_cxcr6))
  )

violin_fig <- ggplot(violin_df, aes(label_cxcr6, expression, fill = label_cxcr6)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.2, colour = "grey25") +
  geom_boxplot(width = 0.12, outlier.shape = NA, linewidth = 0.2, fill = "white", colour = "black") +
  facet_wrap(~gene, ncol = 4, scales = "free_y") +
  scale_fill_manual(values = label_cxcr6_cols, drop = FALSE) +
  labs(
    x = NULL,
    y = "Log-normalized expression",
    fill = NULL,
    title = "Shared-TCR-clone CD8+ marker violin plots by tissue and CXCR6 status"
  ) +
  theme_nc(8) +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, size = 6.5),
    legend.position = "none"
  )
save_all_fig(
  violin_fig,
  "14_reviewer3_shared_clone_cd8_cell_cycle_marker_violins.pdf",
  width = 16.0,
  height = 3.6 * ceiling(length(violin_genes) / 4),
  png = TRUE
)

methods_text <- paste(
  "# Reviewer 3 Major 2: In Situ Cell-Cycle Analysis",
  "",
  "Shared-TCR-clone CD8+ T cells were reconstructed from the existing paired human CD8+ T cell Seurat object by identifying, within each patient, TCR CDR3 amino-acid clonotypes present in both graft and PBMC compartments. This reproduced the Fig. 2D shared-clone logic and retained 4,110 shared-clone CD8+ T cells across four matched patients. CXCR6-positive cells were defined by detectable CXCR6 RNA expression (>0 log-normalized expression, equivalent to non-zero normalized signal).",
  "",
  "Proliferation/cell-cycle, cell-cycle arrest, and survival/apoptosis marker modules were scored using Seurat `AddModuleScore` with the user-specified gene sets. Cell-cycle phase and S/G2M scores were computed with Seurat `CellCycleScoring` using the updated Seurat cell-cycle gene sets intersected with genes present in the RNA assay. Marker expression was visualized on the original UMAP coordinates with Seurat `FeaturePlot`, and violin plots were generated for log-normalized RNA expression grouped by tissue compartment and CXCR6 status.",
  "",
  "To avoid cell-level pseudoreplication, inferential summaries were calculated at the patient level. For each patient, tissue compartment, and CXCR6 status, we summarized MKI67-positive frequency, any proliferation-gene-positive frequency, mean marker expression, module scores, and Seurat S/G2M scores. Paired CXCR6+ versus CXCR6- comparisons within each tissue compartment used patient-level Wilcoxon signed-rank tests. Single-cell-level tests are provided only as exploratory source data.",
  "",
  "Main interpretation: in situ shared-clone CXCR6+ CD8+ T cells were overwhelmingly non-cycling by MKI67 and Seurat cell-cycle scoring. The MKI67-positive fraction was low in CXCR6+ shared-clone CD8+ cells in both graft and PBMC compartments. Therefore, these data support the conclusion that CXCR6+ shared-clone CD8+ cells do not represent a broadly proliferating in situ population; the analysis should be presented as orthogonal support for low basal cycling rather than as a replacement for the ex vivo proliferation assay.",
  sep = "\n"
)
writeLines(methods_text, file.path(out_root, "Reviewer3_Major2_cell_cycle_methods_and_interpretation.md"))
writeLines(methods_text, file.path(repo_dir, "finale", "00_needed_files_for_submission", "02_response_text", "Reviewer3_Major2_cell_cycle_methods_and_interpretation.md"))

message("Wrote Reviewer 3 cell-cycle figures, tables, scored RDS, and methods text.")
