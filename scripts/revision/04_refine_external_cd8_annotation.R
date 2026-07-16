#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(repo_lib)) .libPaths(c(repo_lib, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

set.seed(224445)

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
root_dir <- if (nzchar(project_root)) normalizePath(project_root, mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
processed_dir <- file.path(root_dir, "data", "processed", "external_pbmc_patient_group_analysis")
figure_dir <- file.path(root_dir, "results", "figures", "external_pbmc_patient_group_analysis_refined_cd8")
table_dir <- file.path(root_dir, "results", "tables", "external_pbmc_patient_group_analysis_refined_cd8")
dir.create(figure_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(table_dir, recursive = TRUE, showWarnings = FALSE)

integrated_rds <- file.path(processed_dir, "GSE224445_doubletfinder_rpca_integrated_seurat.rds")
stopifnot(file.exists(integrated_rds))
obj <- readRDS(integrated_rds)

resolution_use <- 1.0
DefaultAssay(obj) <- "integrated"
obj <- FindClusters(
  obj,
  graph.name = "integrated_snn",
  resolution = resolution_use,
  verbose = FALSE
)
highres_col <- paste0("integrated_snn_res.", resolution_use)
if (!highres_col %in% colnames(obj@meta.data)) {
  highres_col <- tail(grep("integrated_snn_res", colnames(obj@meta.data), value = TRUE), 1)
}
obj$highres_cluster <- factor(obj@meta.data[[highres_col]])

get_assay_layer <- function(seu, assay, layer) {
  tryCatch(
    GetAssayData(seu, assay = assay, layer = layer),
    error = function(e) GetAssayData(seu, assay = assay, slot = layer)
  )
}

safe_mean <- function(mat, genes, cells) {
  genes <- intersect(genes, rownames(mat))
  if (!length(genes) || !length(cells)) return(0)
  mean(Matrix::colMeans(mat[genes, cells, drop = FALSE]))
}

safe_pct <- function(mat, genes, cells) {
  genes <- intersect(genes, rownames(mat))
  if (!length(genes) || !length(cells)) return(0)
  mean(Matrix::colSums(mat[genes, cells, drop = FALSE] > 0) > 0)
}

safe_feature_vec <- function(mat, gene) {
  vec <- if (gene %in% rownames(mat)) as.numeric(mat[gene, ]) else rep(0, ncol(mat))
  names(vec) <- colnames(mat)
  vec
}

rna_data <- get_assay_layer(obj, "RNA", "data")
rna_counts <- get_assay_layer(obj, "RNA", "counts")
adt_data <- if ("ADT" %in% Assays(obj)) get_assay_layer(obj, "ADT", "data") else NULL

merge8_marker_sets <- list(
  t_lineage = c("CD3D", "CD3E", "TRAC"),
  cd8_core = c("CD8A", "CD8B"),
  cd4_core = c("CD4"),
  naive_memory_cd8 = c("TCF7", "CCR7", "IL7R", "SELL", "BACH2"),
  cxcr6_effector_cd8 = c("CXCR6", "GZMK", "KLRB1", "KLRG1"),
  cxcl13_cxcr6_cd8 = c("CXCL13", "CXCR6"),
  cx3cr1_effector_cd8 = c("CX3CR1", "KLRG1", "GZMB", "GNLY", "NKG7", "PRF1"),
  znf683_tcf7_cd8 = c("ZNF683", "TCF7"),
  t_activation_exhaustion = c("PDCD1", "HAVCR2", "TBX21", "IFNG"),
  b_lineage = c("MS4A1", "CD79A", "JCHAIN"),
  myeloid = c("CD14", "FCGR3A", "C1QA", "LYZ"),
  cycling = c("MKI67")
)

clusters <- levels(obj$highres_cluster)
cluster_audit <- bind_rows(lapply(clusters, function(cl) {
  cells <- colnames(obj)[obj$highres_cluster == cl]
  out <- data.frame(
    highres_cluster = cl,
    n_cells = length(cells),
    pct_TRAC_or_CD3D_count = safe_pct(rna_counts, c("TRAC", "CD3D"), cells),
    pct_CD8A_count = safe_pct(rna_counts, "CD8A", cells),
    pct_CD8B_count = safe_pct(rna_counts, "CD8B", cells),
    pct_CD8A_or_CD8B_count = safe_pct(rna_counts, c("CD8A", "CD8B"), cells),
    pct_CD4_count = safe_pct(rna_counts, "CD4", cells),
    pct_CXCR6_count = safe_pct(rna_counts, "CXCR6", cells),
    stringsAsFactors = FALSE
  )
  for (nm in names(merge8_marker_sets)) {
    out[[paste0("score_", nm)]] <- safe_mean(rna_data, merge8_marker_sets[[nm]], cells)
  }
  for (gene in intersect(
    c("CXCR6", "CXCL13", "CX3CR1", "KLRG1", "GZMB", "GZMK", "TCF7", "IL7R", "ZNF683", "MKI67"),
    rownames(rna_data)
  )) {
    out[[paste0("avg_", gene)]] <- safe_mean(rna_data, gene, cells)
  }
  if (!is.null(adt_data)) {
    out$adt_CD3 <- safe_mean(adt_data, "CD3:SK7", cells)
    out$adt_CD4 <- safe_mean(adt_data, "CD4:SK3", cells)
    out$adt_CD8 <- safe_mean(adt_data, "CD8:RPA-T8", cells)
    out$adt_TCRgd <- safe_mean(adt_data, "TCR-gamma-delta:B1", cells)
    out$adt_CD19 <- safe_mean(adt_data, "CD19:SJ25C1", cells)
  }
  out
}))

cluster_audit <- cluster_audit %>%
  mutate(
    candidate_cd8_cluster = score_t_lineage > 2 &
      score_cd8_core > 0.5 &
      score_cd8_core >= score_cd4_core &
      score_b_lineage < 1.5 &
      score_myeloid < 1.5,
    candidate_cd4_cluster = score_t_lineage > 2 &
      score_cd4_core > score_cd8_core &
      score_cd4_core > 0.5 &
      score_b_lineage < 1.5 &
      score_myeloid < 1.5,
    candidate_b_cluster = score_b_lineage >= pmax(score_t_lineage, score_myeloid, score_cd8_core, score_cd4_core),
    candidate_myeloid_cluster = score_myeloid >= pmax(score_t_lineage, score_b_lineage, score_cd8_core, score_cd4_core),
    article_annotation_highres = case_when(
      candidate_cd8_cluster & avg_CXCL13 > 0.10 & avg_CXCR6 > 0.20 ~ "CXCL13+CXCR6+ effector CD8+",
      candidate_cd8_cluster & avg_CXCR6 > 0.20 ~ "CXCR6+ effector CD8+",
      candidate_cd8_cluster & avg_MKI67 > 0.15 ~ "Cycling cell",
      candidate_cd8_cluster & (avg_CX3CR1 > 0.50 | avg_KLRG1 > 0.50 | avg_GZMB > 1.00) ~ "CX3CR1+ effector CD8+",
      candidate_cd8_cluster & avg_ZNF683 > 0.20 & avg_TCF7 > 0.20 ~ "ZNF683+TCF7hi CD8+",
      candidate_cd8_cluster ~ "Naive/Memory-like CD8+",
      candidate_cd4_cluster & score_cycling > 0.15 ~ "Cycling cell",
      candidate_cd4_cluster ~ "Naive/Memory-like CD4+",
      candidate_b_cluster & score_b_lineage > 2.5 ~ "Naïve B cell",
      candidate_myeloid_cluster ~ "Classical monocyte",
      TRUE ~ "Low-confidence mixed"
    ),
    keep_for_strict_cd8 = candidate_cd8_cluster
  )

write_csv(
  cluster_audit,
  file.path(table_dir, "GSE224445_highres1_cluster_annotation_audit_merge8_markers.csv")
)

cd8a_count <- safe_feature_vec(rna_counts, "CD8A")
cd8b_count <- safe_feature_vec(rna_counts, "CD8B")
trac_count <- safe_feature_vec(rna_counts, "TRAC")
cd3d_count <- safe_feature_vec(rna_counts, "CD3D")
cd4_count <- safe_feature_vec(rna_counts, "CD4")
ms4a1_count <- safe_feature_vec(rna_counts, "MS4A1")
cd14_count <- safe_feature_vec(rna_counts, "CD14")
fcgr3a_count <- safe_feature_vec(rna_counts, "FCGR3A")
cxcr6_count <- safe_feature_vec(rna_counts, "CXCR6")
cxcr6_log <- safe_feature_vec(rna_data, "CXCR6")

adt_cd8 <- if (!is.null(adt_data) && "CD8:RPA-T8" %in% rownames(adt_data)) as.numeric(adt_data["CD8:RPA-T8", ]) else rep(0, ncol(obj))
adt_cd4 <- if (!is.null(adt_data) && "CD4:SK3" %in% rownames(adt_data)) as.numeric(adt_data["CD4:SK3", ]) else rep(0, ncol(obj))
adt_cd3 <- if (!is.null(adt_data) && "CD3:SK7" %in% rownames(adt_data)) as.numeric(adt_data["CD3:SK7", ]) else rep(0, ncol(obj))

cd8_clusters <- cluster_audit$highres_cluster[cluster_audit$keep_for_strict_cd8]
cell_cluster_is_cd8 <- as.character(obj$highres_cluster) %in% cd8_clusters
cell_t_evidence <- trac_count > 0 | cd3d_count > 0 | adt_cd3 > 0.5
cell_cd8_evidence <- cd8a_count > 0 | cd8b_count > 0 | adt_cd8 > 0.5
cell_non_cd8_conflict <- (cd4_count > 0 & cd8a_count == 0 & cd8b_count == 0 & adt_cd4 > adt_cd8) |
  (ms4a1_count > 0 & !cell_cd8_evidence) |
  ((cd14_count > 0 | fcgr3a_count > 0) & !cell_t_evidence)

obj$strict_cd8_highres <- cell_cluster_is_cd8 & cell_t_evidence & cell_cd8_evidence & !cell_non_cd8_conflict
obj$CXCR6_count <- cxcr6_count
obj$CXCR6_positive <- cxcr6_count > 0
obj$CXCR6_log_norm <- cxcr6_log
obj$CD8A_count_positive <- cd8a_count > 0
obj$CD8B_count_positive <- cd8b_count > 0
obj$CD8A_or_CD8B_count_positive <- cd8a_count > 0 | cd8b_count > 0
obj$TRAC_or_CD3D_count_positive <- trac_count > 0 | cd3d_count > 0
obj$article_annotation_highres <- cluster_audit$article_annotation_highres[match(as.character(obj$highres_cluster), cluster_audit$highres_cluster)]

strict_summary <- obj@meta.data %>%
  mutate(
    highres_cluster = as.character(highres_cluster),
    old_cluster = as.character(seurat_clusters),
    strict_cd8_highres = as.logical(strict_cd8_highres)
  ) %>%
  group_by(highres_cluster, old_cluster, article_annotation_highres) %>%
  summarise(
    n_cells = n(),
    n_strict_cd8 = sum(strict_cd8_highres),
    pct_strict_cd8 = n_strict_cd8 / n_cells,
    pct_CD8A_count = mean(CD8A_count_positive),
    pct_CD8B_count = mean(CD8B_count_positive),
    pct_CD8A_or_CD8B_count = mean(CD8A_or_CD8B_count_positive),
    pct_TRAC_or_CD3D_count = mean(TRAC_or_CD3D_count_positive),
    .groups = "drop"
  )
write_csv(
  strict_summary,
  file.path(table_dir, "GSE224445_highres1_strict_cd8_by_cluster.csv")
)

cd8_cells <- colnames(obj)[obj$strict_cd8_highres]
cd8_refined <- subset(obj, cells = cd8_cells)
DefaultAssay(cd8_refined) <- "RNA"
cd8_refined <- NormalizeData(cd8_refined, verbose = FALSE)
cd8_refined <- FindVariableFeatures(cd8_refined, selection.method = "vst", nfeatures = min(1000, nrow(cd8_refined)), verbose = FALSE)
cd8_refined <- ScaleData(cd8_refined, verbose = FALSE)
cd8_refined <- RunPCA(cd8_refined, npcs = min(20, nrow(cd8_refined) - 1, ncol(cd8_refined) - 1), verbose = FALSE)
cd8_dims <- 1:min(15, ncol(Embeddings(cd8_refined, "pca")))
cd8_refined <- RunUMAP(cd8_refined, dims = cd8_dims, reduction = "pca", reduction.name = "cd8_umap_highres1_strict", verbose = FALSE)
cd8_refined <- FindNeighbors(cd8_refined, dims = cd8_dims, verbose = FALSE)
cd8_refined <- FindClusters(cd8_refined, resolution = 0.8, verbose = FALSE)
cd8_refined$CXCR6_count <- safe_feature_vec(get_assay_layer(cd8_refined, "RNA", "counts"), "CXCR6")
cd8_refined$CXCR6_positive <- cd8_refined$CXCR6_count > 0
cd8_refined$CXCR6_log_norm <- safe_feature_vec(get_assay_layer(cd8_refined, "RNA", "data"), "CXCR6")

saveRDS(
  obj,
  file.path(processed_dir, "GSE224445_doubletfinder_rpca_highres1_annotated_seurat.rds")
)
saveRDS(
  cd8_refined,
  file.path(processed_dir, "GSE224445_doubletfinder_rpca_highres1_strict_cd8_seurat.rds")
)

cd8_meta <- cd8_refined@meta.data %>%
  tibble::rownames_to_column("cell")
write_csv(
  cd8_meta,
  file.path(table_dir, "GSE224445_highres1_strict_cd8_cell_metadata.csv")
)

tag_summary <- cd8_meta %>%
  group_by(group, tag_id, bead, Sample_Tag) %>%
  summarise(
    n_cd8 = n(),
    n_cxcr6_positive = sum(CXCR6_positive),
    cxcr6_positive_frequency = mean(CXCR6_positive),
    mean_cxcr6_log_norm = mean(CXCR6_log_norm),
    median_cxcr6_log_norm = median(CXCR6_log_norm),
    .groups = "drop"
  ) %>%
  arrange(factor(group, levels = c("TOT", "STA", "BPAR")), tag_id)
write_csv(tag_summary, file.path(table_dir, "GSE224445_highres1_strict_cd8_cxcr6_by_tag.csv"))

group_summary <- tag_summary %>%
  group_by(group) %>%
  summarise(
    n_tags = n(),
    total_cd8 = sum(n_cd8),
    pooled_cxcr6_positive_frequency = sum(n_cxcr6_positive) / sum(n_cd8),
    mean_tag_cxcr6_positive_frequency = mean(cxcr6_positive_frequency),
    sd_tag_cxcr6_positive_frequency = sd(cxcr6_positive_frequency),
    mean_tag_cxcr6_log_norm = mean(mean_cxcr6_log_norm),
    sd_tag_cxcr6_log_norm = sd(mean_cxcr6_log_norm),
    .groups = "drop"
  )
write_csv(group_summary, file.path(table_dir, "GSE224445_highres1_strict_cd8_cxcr6_by_group.csv"))

theme_gse <- theme_classic(base_size = 8) +
  theme(
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0, size = 9),
    legend.title = element_text(size = 7),
    legend.text = element_text(size = 7)
  )

pdf(file.path(figure_dir, "01_highres1_strict_cd8_validation_umap.pdf"), width = 10, height = 8)
print(
  (DimPlot(obj, reduction = "umap", group.by = "highres_cluster", label = TRUE, repel = TRUE, pt.size = 0.1) |
     DimPlot(obj, reduction = "umap", group.by = "strict_cd8_highres", pt.size = 0.1)) /
    (FeaturePlot(obj, reduction = "umap", features = "CD8A", cols = c("lightgrey", "firebrick3"), pt.size = 0.1, order = TRUE) |
       FeaturePlot(obj, reduction = "umap", features = "CXCR6", cols = c("lightgrey", "firebrick3"), pt.size = 0.1, order = TRUE)) &
    theme_gse
)
dev.off()

pdf(file.path(figure_dir, "02_highres1_strict_cd8_subset_umap_cxcr6.pdf"), width = 8.8, height = 8.2)
print(
  (DimPlot(cd8_refined, reduction = "cd8_umap_highres1_strict", group.by = "article_annotation_highres", label = TRUE, repel = TRUE, pt.size = 0.1) |
     DimPlot(cd8_refined, reduction = "cd8_umap_highres1_strict", group.by = "group", pt.size = 0.1)) /
    FeaturePlot(cd8_refined, reduction = "cd8_umap_highres1_strict", features = c("CXCR6", "CD8A"), cols = c("lightgrey", "firebrick3"), ncol = 2, pt.size = 0.1, order = TRUE) &
    theme_gse
)
dev.off()

freq_plot <- ggplot(tag_summary, aes(x = factor(group, levels = c("TOT", "STA", "BPAR")), y = cxcr6_positive_frequency, colour = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(shape = bead), size = 2.3, position = position_jitter(width = 0.08, height = 0), stroke = 0.7) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "CXCR6+ frequency in high-resolution strict CD8 T cells", colour = NULL, shape = "Bead") +
  theme_classic(base_size = 8)
expr_plot <- ggplot(tag_summary, aes(x = factor(group, levels = c("TOT", "STA", "BPAR")), y = mean_cxcr6_log_norm, colour = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(shape = bead), size = 2.3, position = position_jitter(width = 0.08, height = 0), stroke = 0.7) +
  labs(x = NULL, y = "Mean CXCR6 log-normalized expression", colour = NULL, shape = "Bead") +
  theme_classic(base_size = 8)
pdf(file.path(figure_dir, "03_highres1_strict_cd8_cxcr6_summary.pdf"), width = 11, height = 5)
print(freq_plot | expr_plot)
dev.off()

message("Wrote high-resolution strict CD8 refinement outputs to: ", table_dir)
