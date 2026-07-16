#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(repo_lib)) .libPaths(c(repo_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(readr)
})

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
root_dir <- if (nzchar(project_root)) normalizePath(project_root, mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
processed_dir <- file.path(root_dir, "data", "processed", "external_pbmc_patient_group_analysis")
table_dir <- file.path(root_dir, "results", "tables", "external_pbmc_patient_group_analysis")
highres_table_dir <- file.path(root_dir, "results", "tables", "external_pbmc_patient_group_analysis_refined_cd8")
cache_dir <- file.path(processed_dir, "visualization_cache")
dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

cache_rds <- file.path(cache_dir, "GSE224445_visualization_environment.rds")

group_levels <- c("TOT", "STA", "BPAR")
group_cols <- c(TOT = "#4C78A8", STA = "#D28E2B", BPAR = "#8E2F3F")
status_cols <- c("Operational tolerance" = "#4C78A8", "Non-tolerant KTR" = "#6F6F6F")
annotation_cols <- c(
  "Naive/Memory-like CD4+" = "#8F6BB1",
  "Effector-like CD4+" = "#A978B8",
  "FOXP3+ Treg" = "#7B4EA3",
  "Naive/Memory-like CD8+" = "#6BAF45",
  "CX3CR1+ effector CD8+" = "#CC6677",
  "CXCR6+ effector CD8+" = "#B33A3A",
  "CXCL13+CXCR6+ effector CD8+" = "#A52A2A",
  "ZNF683+TCF7hi CD8+" = "#6A9D3D",
  "Cycling cell" = "#E77D72",
  "γδ T cell" = "#4C78A8",
  "Classical monocyte" = "#4E9F9C",
  "Resting monocyte" = "#63B8B0",
  "Naïve B cell" = "#C05AA8",
  "Plasma cell" = "#D37ABA",
  "Low-confidence mixed" = "#8A8A8A"
)

fmt_p <- function(p) {
  ifelse(is.na(p), "NA", ifelse(p < 0.001, "<0.001", sprintf("%.3f", p)))
}

cliff_delta <- function(x, y) {
  cmp <- outer(x, y, "-")
  (sum(cmp > 0, na.rm = TRUE) - sum(cmp < 0, na.rm = TRUE)) / length(cmp)
}

make_comparison <- function(tag_summary, a, b, metric) {
  x <- tag_summary[[metric]][tag_summary$group == a]
  y <- tag_summary[[metric]][tag_summary$group == b]
  data.frame(
    comparison = paste(a, "vs", b),
    metric = metric,
    group_a = a,
    group_b = b,
    n_a = length(x),
    n_b = length(y),
    mean_a = mean(x, na.rm = TRUE),
    mean_b = mean(y, na.rm = TRUE),
    median_a = median(x, na.rm = TRUE),
    median_b = median(y, na.rm = TRUE),
    wilcox_p = wilcox.test(x, y, exact = FALSE)$p.value,
    cliff_delta_a_vs_b = cliff_delta(x, y),
    stringsAsFactors = FALSE
  )
}

make_umap_df <- function(seu, reduction, meta_cols, genes = character()) {
  emb <- as.data.frame(Embeddings(seu, reduction)[, 1:2])
  names(emb) <- c("UMAP_1", "UMAP_2")
  meta <- seu@meta.data[, intersect(meta_cols, colnames(seu@meta.data)), drop = FALSE]
  df <- bind_cols(emb, meta)
  genes <- intersect(genes, rownames(seu))
  if (length(genes)) {
    df <- bind_cols(df, FetchData(seu, vars = genes))
  }
  df
}

integrated_rds <- file.path(processed_dir, "GSE224445_doubletfinder_rpca_highres1_annotated_seurat.rds")
if (!file.exists(integrated_rds)) {
  integrated_rds <- file.path(processed_dir, "GSE224445_doubletfinder_rpca_integrated_seurat.rds")
}
cd8_rds <- file.path(processed_dir, "GSE224445_doubletfinder_rpca_highres1_strict_cd8_seurat.rds")
tag_file <- file.path(highres_table_dir, "GSE224445_highres1_strict_cd8_cxcr6_by_tag.csv")
analysis_label <- "High-resolution strict CD8 T cells"
if (!file.exists(cd8_rds) || !file.exists(tag_file)) {
  cd8_rds <- file.path(processed_dir, "GSE224445_doubletfinder_rpca_cd8_cluster_based_seurat.rds")
  tag_file <- file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cxcr6_by_tag.csv")
  analysis_label <- "Doublet-filtered CD8 T cells"
}
stopifnot(file.exists(integrated_rds), file.exists(cd8_rds), file.exists(tag_file))

obj <- readRDS(integrated_rds)
cd8 <- readRDS(cd8_rds)
tag_summary <- read_csv(tag_file, show_col_types = FALSE) %>%
  mutate(
    group = factor(group, levels = group_levels),
    tolerance_status = if_else(group == "TOT", "Operational tolerance", "Non-tolerant KTR"),
    tolerance_status = factor(tolerance_status, levels = c("Operational tolerance", "Non-tolerant KTR")),
    cxcr6_frequency_pct = 100 * cxcr6_positive_frequency
  )

group_comparisons <- bind_rows(
  make_comparison(tag_summary, "TOT", "STA", "cxcr6_positive_frequency"),
  make_comparison(tag_summary, "TOT", "BPAR", "cxcr6_positive_frequency"),
  make_comparison(tag_summary, "STA", "BPAR", "cxcr6_positive_frequency"),
  make_comparison(tag_summary, "TOT", "STA", "mean_cxcr6_log_norm"),
  make_comparison(tag_summary, "TOT", "BPAR", "mean_cxcr6_log_norm"),
  make_comparison(tag_summary, "STA", "BPAR", "mean_cxcr6_log_norm")
)

status_summary <- tag_summary %>%
  group_by(tolerance_status) %>%
  summarise(
    n_tags = n(),
    total_cd8 = sum(n_cd8),
    total_cxcr6_positive = sum(n_cxcr6_positive),
    pooled_cxcr6_positive_frequency = total_cxcr6_positive / total_cd8,
    mean_tag_cxcr6_positive_frequency = mean(cxcr6_positive_frequency),
    sd_tag_cxcr6_positive_frequency = sd(cxcr6_positive_frequency),
    mean_tag_cxcr6_log_norm = mean(mean_cxcr6_log_norm),
    sd_tag_cxcr6_log_norm = sd(mean_cxcr6_log_norm),
    .groups = "drop"
  )

status_comparison <- data.frame(
  comparison = c("TOT_vs_non_tolerant", "TOT_vs_non_tolerant"),
  metric = c("cxcr6_positive_frequency", "mean_cxcr6_log_norm"),
  n_tot = sum(tag_summary$tolerance_status == "Operational tolerance"),
  n_non_tolerant = sum(tag_summary$tolerance_status == "Non-tolerant KTR"),
  mean_tot = c(
    mean(tag_summary$cxcr6_positive_frequency[tag_summary$tolerance_status == "Operational tolerance"]),
    mean(tag_summary$mean_cxcr6_log_norm[tag_summary$tolerance_status == "Operational tolerance"])
  ),
  mean_non_tolerant = c(
    mean(tag_summary$cxcr6_positive_frequency[tag_summary$tolerance_status == "Non-tolerant KTR"]),
    mean(tag_summary$mean_cxcr6_log_norm[tag_summary$tolerance_status == "Non-tolerant KTR"])
  ),
  wilcox_p = c(
    wilcox.test(
      tag_summary$cxcr6_positive_frequency[tag_summary$tolerance_status == "Operational tolerance"],
      tag_summary$cxcr6_positive_frequency[tag_summary$tolerance_status == "Non-tolerant KTR"],
      exact = FALSE
    )$p.value,
    wilcox.test(
      tag_summary$mean_cxcr6_log_norm[tag_summary$tolerance_status == "Operational tolerance"],
      tag_summary$mean_cxcr6_log_norm[tag_summary$tolerance_status == "Non-tolerant KTR"],
      exact = FALSE
    )$p.value
  ),
  stringsAsFactors = FALSE
)

cluster_scores_file <- file.path(table_dir, "GSE224445_doubletfinder_rpca_cluster_annotation_scores.csv")
cluster_audit_file <- file.path(highres_table_dir, "GSE224445_highres1_cluster_annotation_audit_merge8_markers.csv")
strict_cd8_qc_file <- file.path(highres_table_dir, "GSE224445_highres1_strict_cd8_vs_old_cd8_qc.csv")

marker_genes <- intersect(
  c(
    "CD3D", "TRAC", "CD4", "CD8A", "CD8B", "TCF7", "CCR7", "IL7R",
    "NKG7", "GNLY", "TRDC", "KLRB1", "CXCR6", "GZMK", "GZMB", "IFNG",
    "MS4A1", "CD79A", "JCHAIN", "CD14", "FCGR3A", "C1QA", "MKI67",
    "CX3CR1", "KLRG1", "CXCL13", "ZNF683"
  ),
  rownames(obj)
)
obj_reduction <- if ("umap" %in% Reductions(obj)) "umap" else grep("umap", Reductions(obj), value = TRUE)[1]
cd8_reduction <- if ("cd8_umap_highres1_strict" %in% Reductions(cd8)) {
  "cd8_umap_highres1_strict"
} else if ("cd8_umap_doubletfinder_rpca" %in% Reductions(cd8)) {
  "cd8_umap_doubletfinder_rpca"
} else {
  grep("umap", Reductions(cd8), value = TRUE)[1]
}

heatmap_inputs_file <- file.path(table_dir, "GSE224445_doubletfinder_rpca_annotation_heatmap_inputs_article_names.rds")
heatmap_inputs <- if (file.exists(heatmap_inputs_file)) readRDS(heatmap_inputs_file) else NULL

visualization_env <- list(
  created_at = as.character(Sys.time()),
  r_version = R.version.string,
  seurat_version = as.character(packageVersion("Seurat")),
  source_files = list(
    integrated_rds = integrated_rds,
    cd8_rds = cd8_rds,
    tag_file = tag_file,
    cluster_scores_file = cluster_scores_file,
    cluster_audit_file = cluster_audit_file,
    heatmap_inputs_file = heatmap_inputs_file
  ),
  analysis_label = analysis_label,
  group_levels = group_levels,
  group_cols = group_cols,
  status_cols = status_cols,
  annotation_cols = annotation_cols,
  marker_genes = marker_genes,
  obj_reduction = obj_reduction,
  cd8_reduction = cd8_reduction,
  obj = obj,
  cd8 = cd8,
  obj_umap_df = make_umap_df(
    obj,
    obj_reduction,
    c("group", "seurat_clusters", "highres_cluster", "cluster_annotation_article", "article_annotation_highres", "CXCR6_log_norm"),
    marker_genes
  ),
  cd8_umap_df = make_umap_df(
    cd8,
    cd8_reduction,
    c("group", "tag_id", "bead", "Sample_Tag", "article_annotation_highres", "cluster_annotation_article", "CXCR6_log_norm", "CXCR6_positive"),
    intersect(c("CXCR6", "CD8A", "CD8B", "TRAC", "CD3D"), rownames(cd8))
  ),
  tag_summary = tag_summary,
  group_comparisons = group_comparisons,
  status_summary = status_summary,
  status_comparison = status_comparison,
  cluster_scores = if (file.exists(cluster_scores_file)) read_csv(cluster_scores_file, show_col_types = FALSE) else NULL,
  cluster_audit = if (file.exists(cluster_audit_file)) read_csv(cluster_audit_file, show_col_types = FALSE) else NULL,
  strict_cd8_qc = if (file.exists(strict_cd8_qc_file)) read_csv(strict_cd8_qc_file, show_col_types = FALSE) else NULL,
  heatmap_inputs = heatmap_inputs
)

saveRDS(visualization_env, cache_rds, compress = "gzip")

message("Saved visualization environment RDS: ", cache_rds)
message("Use in R: viz <- readRDS(\"", cache_rds, "\")")
