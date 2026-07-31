#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
repo_lib <- Sys.getenv("R_LIBS_USER", unset = "")
if (nzchar(repo_lib)) .libPaths(c(repo_lib, .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(Seurat)
  library(DoubletFinder)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
})

set.seed(224445)
options(future.globals.maxSize = 8 * 1024^3)

project_root <- Sys.getenv("PROJECT_ROOT", unset = "")
root_dir <- if (nzchar(project_root)) normalizePath(project_root, mustWork = TRUE) else normalizePath(getwd(), mustWork = TRUE)
raw_dir <- file.path(root_dir, "data", "external_pbmc_raw")
processed_dir <- file.path(root_dir, "data", "processed", "external_pbmc_patient_group_analysis")
figure_dir <- file.path(root_dir, "results", "figures", "external_pbmc_patient_group_analysis")
table_dir <- file.path(root_dir, "results", "tables", "external_pbmc_patient_group_analysis")
log_dir <- file.path(root_dir, "logs", "external_pbmc_patient_group_analysis")
for (d in c(processed_dir, figure_dir, table_dir, log_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

sink(file.path(log_dir, "GSE224445_doubletfinder_rpca_pipeline.log"), split = TRUE)
on.exit({
  sink(type = "output")
}, add = TRUE)

message_stamp <- function(...) {
  message(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), " | ", paste0(..., collapse = ""))
}

message_stamp("R version: ", R.version.string)
message_stamp("Seurat version: ", as.character(packageVersion("Seurat")))
message_stamp("DoubletFinder version: ", as.character(packageVersion("DoubletFinder")))

tag_map <- data.frame(
  bead = c(rep("Bead0", 4), rep("Bead1", 3), rep("Bead2", 4), "Bead3"),
  sample_tag = c(
    "SampleTag05_hs", "SampleTag06_hs", "SampleTag07_hs", "SampleTag08_hs",
    "SampleTag01_hs", "SampleTag02_hs", "SampleTag04_hs",
    "SampleTag05_hs", "SampleTag06_hs", "SampleTag07_hs", "SampleTag08_hs",
    "SampleTag09_hs"
  ),
  group = c("BPAR", "BPAR", "STA", "STA", "TOT", "TOT", "TOT", "BPAR", "BPAR", "STA", "STA", "TOT"),
  tag_id = c(
    "Bead0_Tag05", "Bead0_Tag06", "Bead0_Tag07", "Bead0_Tag08",
    "Bead1_Tag01", "Bead1_Tag02", "Bead1_Tag04",
    "Bead2_Tag05", "Bead2_Tag06", "Bead2_Tag07", "Bead2_Tag08",
    "Bead3_Tag09"
  ),
  stringsAsFactors = FALSE
)
write_csv(tag_map, file.path(table_dir, "GSE224445_doubletfinder_rpca_tag_group_map.csv"))

find_header_skip <- function(path) {
  con <- gzfile(path, open = "rt")
  on.exit(close(con))
  i <- 0L
  repeat {
    line <- readLines(con, n = 1)
    if (!length(line)) stop("No Cell_Index header found in ", path)
    i <- i + 1L
    if (startsWith(line, "Cell_Index,")) return(i - 1L)
  }
}

read_bd_csv <- function(path) {
  read.csv(gzfile(path), skip = find_header_skip(path), check.names = FALSE)
}

parse_sample_id <- function(path) {
  base <- basename(path)
  list(
    gsm = sub("^(GSM[0-9]+).*", "\\1", base),
    bead = sub("^GSM[0-9]+_(Bead[0-9]+).*", "\\1", base)
  )
}

feature_info <- function(feature_names) {
  parts <- strsplit(feature_names, "\\|")
  data.frame(
    original_name = feature_names,
    symbol = vapply(parts, `[`, character(1), 1),
    molecule_type = ifelse(grepl("\\|pAbO$", feature_names), "AbSeq", "RNA"),
    stringsAsFactors = FALSE
  )
}

aggregate_feature_matrix <- function(cell_by_feature, info, cell_names) {
  if (ncol(cell_by_feature) == 0) return(NULL)
  mat <- as.matrix(cell_by_feature)
  storage.mode(mat) <- "numeric"
  feat_by_cell <- rowsum(t(mat), group = info$symbol, reorder = FALSE)
  rownames(feat_by_cell) <- make.unique(rownames(feat_by_cell))
  colnames(feat_by_cell) <- cell_names
  Matrix(feat_by_cell, sparse = TRUE)
}

make_bead_object <- function(mol_file) {
  ids <- parse_sample_id(mol_file)
  tag_file <- list.files(raw_dir, pattern = paste0(ids$gsm, ".*SCMC_Sample_Tag_Calls.*csv\\.gz$"), full.names = TRUE)
  if (length(tag_file) != 1) stop("Could not find one tag file for ", ids$gsm)

  message_stamp("Loading ", ids$bead, " from ", basename(mol_file))
  mol <- read_bd_csv(mol_file)
  tags <- read_bd_csv(tag_file)
  if (!identical(mol$Cell_Index, tags$Cell_Index)) {
    mol <- merge(mol, tags[, "Cell_Index"], by = "Cell_Index")
    tags <- tags[match(mol$Cell_Index, tags$Cell_Index), ]
  }

  meta <- as.data.frame(tags)
  meta$cell_index <- as.character(meta$Cell_Index)
  meta$gsm <- ids$gsm
  meta$bead <- ids$bead
  meta <- left_join(meta, tag_map, by = c("bead", "Sample_Tag" = "sample_tag"))
  meta$sample_tag_status <- ifelse(is.na(meta$group), as.character(meta$Sample_Tag), "Assigned")
  keep <- !is.na(meta$group)
  meta <- meta[keep, , drop = FALSE]
  mol <- mol[keep, , drop = FALSE]
  cell_names <- paste(ids$bead, ids$gsm, meta$cell_index, sep = "_")
  rownames(meta) <- cell_names

  features <- feature_info(setdiff(names(mol), "Cell_Index"))
  feature_cols <- features$original_name
  count_df <- as.data.frame(mol[, feature_cols, drop = FALSE])
  rna_idx <- features$molecule_type == "RNA"
  ab_idx <- features$molecule_type == "AbSeq"
  rna_counts <- aggregate_feature_matrix(count_df[, rna_idx, drop = FALSE], features[rna_idx, ], cell_names)
  ab_counts <- aggregate_feature_matrix(count_df[, ab_idx, drop = FALSE], features[ab_idx, ], cell_names)

  obj <- CreateSeuratObject(counts = rna_counts, project = ids$bead, meta.data = meta, min.cells = 0, min.features = 0)
  obj$orig.ident <- ids$bead
  if (!is.null(ab_counts)) obj[["ADT"]] <- CreateAssayObject(counts = ab_counts)
  obj
}

prepare_for_doubletfinder <- function(obj, dims = 1:20) {
  DefaultAssay(obj) <- "RNA"
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = min(2000, nrow(obj)), verbose = FALSE)
  obj <- ScaleData(obj, verbose = FALSE)
  npcs <- min(30, nrow(obj) - 1, ncol(obj) - 1)
  obj <- RunPCA(obj, npcs = npcs, verbose = FALSE)
  dims_use <- dims[dims <= ncol(Embeddings(obj, "pca"))]
  obj <- FindNeighbors(obj, dims = dims_use, verbose = FALSE)
  obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
  obj <- RunUMAP(obj, dims = dims_use[dims_use <= 10], verbose = FALSE)
  if ("ADT" %in% Assays(obj)) {
    DefaultAssay(obj) <- "ADT"
    obj <- NormalizeData(obj, normalization.method = "CLR", margin = 2, verbose = FALSE)
    DefaultAssay(obj) <- "RNA"
  }
  obj
}

pick_pk <- function(bcmvn) {
  if (is.null(bcmvn) || nrow(bcmvn) == 0 || !"BCmetric" %in% colnames(bcmvn)) {
    return(list(pK = 0.01, method = "fallback_no_bcmvn"))
  }
  bc_values <- suppressWarnings(as.numeric(bcmvn$BCmetric))
  pk_values <- suppressWarnings(as.numeric(as.character(bcmvn$pK)))
  ok <- is.finite(bc_values) & is.finite(pk_values)
  if (!any(ok)) return(list(pK = 0.01, method = "fallback_no_finite_metric"))
  max_idx <- which(ok)[which.max(bc_values[ok])]
  list(pK = pk_values[max_idx], method = "max_BCmetric")
}

run_doubletfinder_per_tag <- function(obj, sample_name) {
  message_stamp("DoubletFinder start: ", sample_name, " (", ncol(obj), " cells)")
  obj <- prepare_for_doubletfinder(obj)
  dims_use <- 1:min(10, ncol(Embeddings(obj, "pca")))

  sweep.res.list <- tryCatch(
    paramSweep(obj, PCs = dims_use, sct = FALSE, num.cores = 1),
    error = function(e) {
      message_stamp("paramSweep failed for ", sample_name, ": ", e$message)
      NULL
    }
  )
  sweep.stats <- if (!is.null(sweep.res.list)) summarizeSweep(sweep.res.list, GT = FALSE) else NULL
  bcmvn <- if (!is.null(sweep.stats)) find.pK(sweep.stats) else NULL
  pk_pick <- pick_pk(bcmvn)

  homotypic.prop <- modelHomotypic(obj@meta.data$seurat_clusters)
  doublet_rate <- nrow(obj@meta.data) / 500 * 0.4 / 100
  nExp_poi <- max(1, round(doublet_rate * nrow(obj@meta.data)))
  nExp_poi.adj <- max(1, round(nExp_poi * (1 - homotypic.prop)))

  obj <- doubletFinder(
    obj,
    PCs = dims_use,
    pN = 0.25,
    pK = pk_pick$pK,
    nExp = nExp_poi,
    sct = FALSE
  )
  pann_col <- grep("^pANN", colnames(obj@meta.data), value = TRUE)
  if (length(pann_col) < 1) stop("No pANN column found after first DoubletFinder pass for ", sample_name)
  pann_col <- tail(pann_col, 1)

  obj <- doubletFinder(
    obj,
    PCs = dims_use,
    pN = 0.25,
    pK = pk_pick$pK,
    nExp = nExp_poi.adj,
    reuse.pANN = pann_col,
    sct = FALSE
  )
  class_col <- tail(grep("^DF.classifications", colnames(obj@meta.data), value = TRUE), 1)
  if (!length(class_col)) stop("No DF.classifications column found for ", sample_name)
  obj$doubletfinder_class <- obj@meta.data[[class_col]]
  obj$doubletfinder_pANN <- obj@meta.data[[pann_col]]
  obj$doubletfinder_pK <- pk_pick$pK
  obj$doubletfinder_pK_method <- pk_pick$method
  obj$doubletfinder_expected_rate <- doublet_rate
  obj$doubletfinder_nExp <- nExp_poi
  obj$doubletfinder_nExp_adj <- nExp_poi.adj
  obj$doubletfinder_homotypic_prop <- homotypic.prop

  df_meta <- obj@meta.data
  df_meta$cell <- rownames(df_meta)
  write_tsv(df_meta, file.path(table_dir, paste0(sample_name, "_doubletfinder.tsv")))
  if (!is.null(bcmvn)) write_tsv(bcmvn, file.path(table_dir, paste0(sample_name, "_doubletfinder_pK_sweep.tsv")))

  n_doublet <- sum(obj$doubletfinder_class == "Doublet", na.rm = TRUE)
  n_singlet <- sum(obj$doubletfinder_class == "Singlet", na.rm = TRUE)
  message_stamp("DoubletFinder done: ", sample_name, " singlets=", n_singlet, " doublets=", n_doublet, " pK=", pk_pick$pK)

  summary <- data.frame(
    tag_id = sample_name,
    group = unique(obj$group),
    bead = unique(obj$bead),
    cells_input = ncol(obj),
    singlets = n_singlet,
    doublets = n_doublet,
    doublet_fraction = n_doublet / ncol(obj),
    doublet_rate_formula = doublet_rate,
    nExp_poi = nExp_poi,
    nExp_poi_adj = nExp_poi.adj,
    homotypic_prop = homotypic.prop,
    pK = pk_pick$pK,
    pK_method = pk_pick$method,
    stringsAsFactors = FALSE
  )

  list(object = obj, singlet = subset(obj, subset = doubletfinder_class == "Singlet"), summary = summary)
}

score_features <- function(obj, features, assay = "RNA", layer = "data") {
  present <- intersect(features, rownames(obj[[assay]]))
  if (length(present) == 0) return(rep(0, ncol(obj)))
  Matrix::colMeans(GetAssayData(obj, assay = assay, layer = layer)[present, , drop = FALSE])
}

cluster_marker_summary <- function(object, genes, assay = "RNA", layer = "data") {
  genes <- intersect(genes, rownames(object[[assay]]))
  if (length(genes) == 0) return(data.frame())
  mat <- GetAssayData(object, assay = assay, layer = layer)[genes, , drop = FALSE]
  cl <- as.character(object$seurat_clusters)
  out <- sapply(sort(unique(cl)), function(k) Matrix::rowMeans(mat[, cl == k, drop = FALSE]))
  as.data.frame(t(out)) %>%
    tibble::rownames_to_column("seurat_clusters") %>%
    pivot_longer(-seurat_clusters, names_to = "gene", values_to = "avg_log_norm")
}

annotate_clusters <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  obj$score_t <- score_features(obj, c("CD3D", "CD3E", "CD3G", "TRAC", "TRBC2", "CD247"))
  obj$score_cd8 <- score_features(obj, c("CD8A", "CD8B", "RUNX3", "GZMK", "GZMA", "GZMB"))
  obj$score_cd4 <- score_features(obj, c("CD4", "IL7R", "CCR7", "FOXP3", "IL2RA", "CTLA4"))
  obj$score_nk <- score_features(obj, c("NKG7", "GNLY", "KLRD1", "NCAM1", "KLRF1", "TRDC"))
  obj$score_b <- score_features(obj, c("MS4A1", "CD79A", "CD79B", "CD19", "PAX5", "JCHAIN"))
  obj$score_myeloid <- score_features(obj, c("LYZ", "CD14", "FCGR3A", "S100A9", "FCN1", "C1QA"))
  if ("ADT" %in% Assays(obj)) {
    DefaultAssay(obj) <- "ADT"
    obj$score_t <- obj$score_t + score_features(obj, c("CD3:SK7"), assay = "ADT")
    obj$score_cd8 <- obj$score_cd8 + score_features(obj, c("CD8:RPA-T8"), assay = "ADT")
    obj$score_cd4 <- obj$score_cd4 + score_features(obj, c("CD4:SK3"), assay = "ADT")
    DefaultAssay(obj) <- "RNA"
  }

  cluster_scores <- obj@meta.data %>%
    mutate(seurat_clusters = as.character(seurat_clusters)) %>%
    group_by(seurat_clusters) %>%
    summarise(
      n_cells = n(),
      across(starts_with("score_"), mean, .names = "mean_{.col}"),
      .groups = "drop"
    )

  marker_genes <- c(
    "CD3D", "CD3E", "TRAC", "CD4", "CD8A", "CD8B", "TCF7", "CCR7", "IL7R",
    "FOXP3", "IL2RA", "NKG7", "GNLY", "KLRD1", "TRDC", "KLRB1", "CXCR6",
    "GZMK", "GZMB", "IFNG", "MS4A1", "CD79A", "JCHAIN", "LYZ", "CD14",
    "FCGR3A", "C1QA", "MKI67"
  )
  marker_avg <- cluster_marker_summary(obj, marker_genes)
  marker_wide <- marker_avg %>%
    pivot_wider(names_from = gene, values_from = avg_log_norm, values_fill = 0)
  cluster_scores <- left_join(cluster_scores, marker_wide, by = "seurat_clusters")

  getv <- function(row, nm) if (nm %in% names(row)) as.numeric(row[[nm]]) else 0
  ann <- apply(cluster_scores, 1, function(row) {
    t_score <- getv(row, "mean_score_t")
    cd8_score <- getv(row, "mean_score_cd8")
    cd4_score <- getv(row, "mean_score_cd4")
    nk_score <- getv(row, "mean_score_nk")
    b_score <- getv(row, "mean_score_b")
    my_score <- getv(row, "mean_score_myeloid")
    cxcr6 <- getv(row, "CXCR6")
    gzmb <- getv(row, "GZMB")
    gNLY <- getv(row, "GNLY")
    tcf7 <- getv(row, "TCF7")
    mki67 <- getv(row, "MKI67")
    foxp3 <- getv(row, "FOXP3")

    if (my_score > max(t_score, b_score, nk_score, cd8_score, cd4_score)) return("Myeloid")
    if (b_score > max(t_score, my_score, nk_score, cd8_score, cd4_score)) return(ifelse(getv(row, "JCHAIN") > 0.5, "Plasma/B cell", "B cell"))
    if (t_score <= 0 && nk_score > max(b_score, my_score)) return("NK-like")
    if (cd4_score > cd8_score && cd4_score >= t_score * 0.35) {
      if (foxp3 > 0.1) return("FOXP3+ regulatory CD4 T")
      return("CD4 T")
    }
    if (cd8_score >= cd4_score && (t_score > 0 || cd8_score > 0)) {
      if (mki67 > 0.2) return("Cycling cytotoxic CD8 T")
      if (cxcr6 > 0.2 && gzmb > 0.2) return("CXCR6+ effector CD8 T")
      if (cxcr6 > 0.2) return("CXCR6+ CD8 T")
      if (tcf7 > 0.2) return("Naive/Memory-like CD8 T")
      if (gNLY > 0.5 && nk_score > t_score) return("NK/gamma-delta-like")
      return("Cytotoxic CD8 T")
    }
    if (nk_score > 0) return("NK/gamma-delta-like")
    "Other"
  })

  cluster_scores$cluster_annotation_repro <- ann
  cd8_labels <- c(
    "Cycling cytotoxic CD8 T",
    "CXCR6+ effector CD8 T",
    "CXCR6+ CD8 T",
    "Naive/Memory-like CD8 T",
    "Cytotoxic CD8 T"
  )
  cluster_scores$is_repro_cd8_cluster <- cluster_scores$cluster_annotation_repro %in% cd8_labels
  obj$cluster_annotation_repro <- cluster_scores$cluster_annotation_repro[match(as.character(obj$seurat_clusters), cluster_scores$seurat_clusters)]
  obj$is_repro_cd8_cluster <- cluster_scores$is_repro_cd8_cluster[match(as.character(obj$seurat_clusters), cluster_scores$seurat_clusters)]

  list(object = obj, cluster_scores = cluster_scores, marker_average = marker_avg)
}

mol_files <- sort(list.files(raw_dir, pattern = "Combined_SCMC_DBEC_MolsPerCell\\.csv\\.gz$", full.names = TRUE))
if (length(mol_files) != 4) stop("Expected 4 molecule files in ", raw_dir)

bead_objects <- lapply(mol_files, make_bead_object)
names(bead_objects) <- vapply(bead_objects, function(x) unique(x$bead), character(1))
bead_counts <- data.frame(bead = names(bead_objects), n_assigned_cells = vapply(bead_objects, ncol, numeric(1)))
write_csv(bead_counts, file.path(table_dir, "GSE224445_doubletfinder_rpca_cells_loaded_by_bead.csv"))

sample_objects <- list()
for (bead_name in names(bead_objects)) {
  obj <- bead_objects[[bead_name]]
  for (tag in sort(unique(obj$tag_id))) {
    sample_objects[[tag]] <- subset(obj, cells = colnames(obj)[obj$tag_id == tag])
  }
}
sample_counts <- data.frame(
  tag_id = names(sample_objects),
  bead = vapply(sample_objects, function(x) unique(x$bead), character(1)),
  group = vapply(sample_objects, function(x) unique(x$group), character(1)),
  cells_before_doubletfinder = vapply(sample_objects, ncol, numeric(1)),
  stringsAsFactors = FALSE
)
write_csv(sample_counts, file.path(table_dir, "GSE224445_doubletfinder_rpca_cells_loaded_by_tag.csv"))

df_results <- lapply(names(sample_objects), function(tag) run_doubletfinder_per_tag(sample_objects[[tag]], tag))
names(df_results) <- names(sample_objects)
doublet_summary <- bind_rows(lapply(df_results, `[[`, "summary")) %>%
  arrange(bead, tag_id)
write_csv(doublet_summary, file.path(table_dir, "GSE224445_doubletfinder_rpca_doublet_summary.csv"))

singlet_objects <- lapply(df_results, `[[`, "singlet")
names(singlet_objects) <- names(df_results)
for (nm in names(singlet_objects)) {
  DefaultAssay(singlet_objects[[nm]]) <- "RNA"
  singlet_objects[[nm]] <- NormalizeData(singlet_objects[[nm]], verbose = FALSE)
  singlet_objects[[nm]] <- FindVariableFeatures(singlet_objects[[nm]], selection.method = "vst", nfeatures = min(2000, nrow(singlet_objects[[nm]])), verbose = FALSE)
  singlet_objects[[nm]] <- ScaleData(singlet_objects[[nm]], verbose = FALSE)
  npcs <- min(30, nrow(singlet_objects[[nm]]) - 1, ncol(singlet_objects[[nm]]) - 1)
  singlet_objects[[nm]] <- RunPCA(singlet_objects[[nm]], npcs = npcs, verbose = FALSE)
  if ("ADT" %in% Assays(singlet_objects[[nm]])) {
    DefaultAssay(singlet_objects[[nm]]) <- "ADT"
    singlet_objects[[nm]] <- NormalizeData(singlet_objects[[nm]], normalization.method = "CLR", margin = 2, verbose = FALSE)
    DefaultAssay(singlet_objects[[nm]]) <- "RNA"
  }
}

message_stamp("Starting RPCA integration on ", length(singlet_objects), " singlet sample-tag objects")
features <- SelectIntegrationFeatures(object.list = singlet_objects, nfeatures = min(2000, nrow(singlet_objects[[1]])))
anchors <- FindIntegrationAnchors(
  object.list = singlet_objects,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:20,
  verbose = FALSE
)
merged <- IntegrateData(anchorset = anchors, dims = 1:20, verbose = FALSE)
DefaultAssay(merged) <- "integrated"
merged <- ScaleData(merged, verbose = FALSE)
merged <- RunPCA(merged, npcs = 30, verbose = FALSE)
merged <- RunUMAP(merged, dims = 1:20, reduction = "pca", verbose = FALSE)
merged <- FindNeighbors(merged, dims = 1:20, verbose = FALSE)
merged <- FindClusters(merged, resolution = 0.5, verbose = FALSE)
merged <- JoinLayers(merged, assay = "RNA")
if ("ADT" %in% Assays(merged) && inherits(merged[["ADT"]], "Assay5")) merged <- JoinLayers(merged, assay = "ADT")

cx_counts <- as.numeric(GetAssayData(merged, assay = "RNA", layer = "counts")["CXCR6", ])
cx_log <- as.numeric(GetAssayData(merged, assay = "RNA", layer = "data")["CXCR6", ])
merged$CXCR6_count <- cx_counts
merged$CXCR6_positive <- cx_counts > 0
merged$CXCR6_log_norm <- cx_log

ann <- annotate_clusters(merged)
merged <- ann$object
write_csv(ann$cluster_scores, file.path(table_dir, "GSE224445_doubletfinder_rpca_cluster_annotation_scores.csv"))
write_csv(ann$marker_average, file.path(table_dir, "GSE224445_doubletfinder_rpca_cluster_marker_average_expression.csv"))

merged_meta <- merged@meta.data %>%
  tibble::rownames_to_column("cell")
write_csv(merged_meta, file.path(table_dir, "GSE224445_doubletfinder_rpca_integrated_cell_metadata.csv"))
saveRDS(merged, file.path(processed_dir, "GSE224445_doubletfinder_rpca_integrated_seurat.rds"))

DefaultAssay(merged) <- "RNA"
all_markers <- FindAllMarkers(merged, assay = "RNA", only.pos = TRUE, logfc.threshold = 0.25, min.pct = 0.10)
write_csv(all_markers, file.path(table_dir, "GSE224445_doubletfinder_rpca_all_cluster_markers.csv"))

pdf(file.path(figure_dir, "01_gse224445_doubletfinder_rpca_integrated_umap_overview.pdf"), width = 14, height = 10)
print(
  (DimPlot(merged, reduction = "umap", group.by = "group", pt.size = 0.18) |
     DimPlot(merged, reduction = "umap", group.by = "seurat_clusters", label = TRUE, repel = TRUE, pt.size = 0.18)) /
    (DimPlot(merged, reduction = "umap", group.by = "cluster_annotation_repro", label = TRUE, repel = TRUE, pt.size = 0.18) |
       FeaturePlot(merged, features = "CXCR6", reduction = "umap", cols = c("lightgrey", "firebrick3"), order = TRUE, pt.size = 0.16))
)
dev.off()

marker_genes <- intersect(
  c("CD3D", "TRAC", "CD4", "CD8A", "CD8B", "TCF7", "CCR7", "IL7R", "NKG7", "GNLY", "KLRD1", "TRDC", "KLRB1", "CXCR6", "GZMK", "GZMB", "IFNG", "MS4A1", "CD79A", "JCHAIN", "LYZ", "CD14", "FCGR3A", "C1QA", "MKI67"),
  rownames(merged)
)
pdf(file.path(figure_dir, "02_gse224445_doubletfinder_rpca_marker_featureplots.pdf"), width = 12, height = 30)
print(
  FeaturePlot(
    merged,
    features = marker_genes,
    cols = c("lightgrey", "firebrick3"),
    ncol = 2,
    label = FALSE,
    pt.size = 0.12,
    order = TRUE
  )
)
dev.off()

pdf(file.path(figure_dir, "03_gse224445_doubletfinder_rpca_annotation_dotplot.pdf"), width = 13, height = 7)
print(
  DotPlot(merged, features = marker_genes, group.by = "cluster_annotation_repro", assay = "RNA", scale = TRUE) +
    RotatedAxis() +
    scale_colour_gradient(low = "lightgrey", high = "firebrick3") +
    ggtitle("Marker expression by reproducible cluster annotation")
)
dev.off()

heat_df <- ann$marker_average %>%
  filter(gene %in% marker_genes) %>%
  left_join(ann$cluster_scores[, c("seurat_clusters", "cluster_annotation_repro")], by = "seurat_clusters") %>%
  group_by(cluster_annotation_repro, gene) %>%
  summarise(avg_log_norm = mean(avg_log_norm), .groups = "drop") %>%
  pivot_wider(names_from = gene, values_from = avg_log_norm, values_fill = 0)
heat_mat <- as.matrix(heat_df[, -1, drop = FALSE])
rownames(heat_mat) <- heat_df$cluster_annotation_repro
heat_mat <- t(scale(t(heat_mat)))
heat_mat[is.na(heat_mat)] <- 0
pdf(file.path(figure_dir, "04_gse224445_doubletfinder_rpca_annotation_heatmap.pdf"), width = 10, height = 7)
pheatmap(
  heat_mat,
  color = colorRampPalette(c("#3B5A7A", "white", "firebrick3"))(101),
  border_color = NA,
  fontsize = 8,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  main = "Reproducible cluster annotation marker heatmap"
)
dev.off()

cd8_clusters <- ann$cluster_scores$seurat_clusters[ann$cluster_scores$is_repro_cd8_cluster]
if (!length(cd8_clusters)) stop("No reproducible CD8 clusters were identified.")
message_stamp("CD8 clusters retained: ", paste(cd8_clusters, collapse = ", "))

cd8 <- subset(merged, subset = is_repro_cd8_cluster)
DefaultAssay(cd8) <- "RNA"
cd8 <- NormalizeData(cd8, verbose = FALSE)
cd8 <- FindVariableFeatures(cd8, selection.method = "vst", nfeatures = min(1000, nrow(cd8)), verbose = FALSE)
cd8 <- ScaleData(cd8, verbose = FALSE)
cd8 <- RunPCA(cd8, npcs = min(20, nrow(cd8) - 1, ncol(cd8) - 1), verbose = FALSE)
cd8_dims <- 1:min(15, ncol(Embeddings(cd8, "pca")))
cd8 <- RunUMAP(cd8, dims = cd8_dims, reduction = "pca", reduction.name = "cd8_umap_doubletfinder_rpca", verbose = FALSE)
cd8 <- FindNeighbors(cd8, dims = cd8_dims, verbose = FALSE)
cd8 <- FindClusters(cd8, resolution = 0.4, verbose = FALSE)
cd8$CXCR6_count <- as.numeric(GetAssayData(cd8, assay = "RNA", layer = "counts")["CXCR6", ])
cd8$CXCR6_positive <- cd8$CXCR6_count > 0
cd8$CXCR6_log_norm <- as.numeric(GetAssayData(cd8, assay = "RNA", layer = "data")["CXCR6", ])
saveRDS(cd8, file.path(processed_dir, "GSE224445_doubletfinder_rpca_cd8_cluster_based_seurat.rds"))

cd8_meta <- cd8@meta.data %>%
  tibble::rownames_to_column("cell")
write_csv(cd8_meta, file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cluster_based_cell_metadata.csv"))

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
write_csv(tag_summary, file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cxcr6_by_tag.csv"))

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
write_csv(group_summary, file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cxcr6_by_group.csv"))

bead_group_summary <- tag_summary %>%
  group_by(bead, group) %>%
  summarise(
    n_tags = n(),
    total_cd8 = sum(n_cd8),
    pooled_cxcr6_positive_frequency = sum(n_cxcr6_positive) / sum(n_cd8),
    mean_tag_cxcr6_positive_frequency = mean(cxcr6_positive_frequency),
    mean_tag_cxcr6_log_norm = mean(mean_cxcr6_log_norm),
    .groups = "drop"
  )
write_csv(bead_group_summary, file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cxcr6_by_bead_group.csv"))

pairwise <- do.call(rbind, lapply(c("TOT", "STA"), function(ref) {
  data.frame(
    comparison = paste("BPAR_vs", ref, sep = "_"),
    wilcox_p_frequency = wilcox.test(
      tag_summary$cxcr6_positive_frequency[tag_summary$group == "BPAR"],
      tag_summary$cxcr6_positive_frequency[tag_summary$group == ref],
      exact = FALSE
    )$p.value,
    wilcox_p_mean_expression = wilcox.test(
      tag_summary$mean_cxcr6_log_norm[tag_summary$group == "BPAR"],
      tag_summary$mean_cxcr6_log_norm[tag_summary$group == ref],
      exact = FALSE
    )$p.value
  )
}))
write_csv(pairwise, file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cxcr6_tag_level_wilcox.csv"))

theme_gse224445 <- theme_classic(base_size = 9) +
  theme(
    axis.text = element_text(colour = "black"),
    axis.title = element_text(colour = "black"),
    plot.title = element_text(face = "bold", hjust = 0)
  )

pdf(file.path(figure_dir, "05_gse224445_doubletfinder_rpca_cd8_umap_annotation_group_cxcr6.pdf"), width = 14, height = 8)
print(
  (DimPlot(cd8, reduction = "cd8_umap_doubletfinder_rpca", group.by = "cluster_annotation_repro", label = TRUE, repel = TRUE, pt.size = 0.22) |
     DimPlot(cd8, reduction = "cd8_umap_doubletfinder_rpca", group.by = "group", pt.size = 0.22)) /
    FeaturePlot(cd8, reduction = "cd8_umap_doubletfinder_rpca", features = c("CXCR6", "CD8A"), cols = c("lightgrey", "firebrick3"), ncol = 2, pt.size = 0.22, order = TRUE)
)
dev.off()

freq_plot <- ggplot(tag_summary, aes(x = factor(group, levels = c("TOT", "STA", "BPAR")), y = cxcr6_positive_frequency, colour = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(shape = bead), size = 3, position = position_jitter(width = 0.08, height = 0)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
  labs(x = NULL, y = "CXCR6+ frequency in doublet-filtered CD8 T cells", colour = "Group", shape = "Bead") +
  theme_gse224445
expr_plot <- ggplot(tag_summary, aes(x = factor(group, levels = c("TOT", "STA", "BPAR")), y = mean_cxcr6_log_norm, colour = group)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.15) +
  geom_point(aes(shape = bead), size = 3, position = position_jitter(width = 0.08, height = 0)) +
  labs(x = NULL, y = "Mean CXCR6 log-normalized expression", colour = "Group", shape = "Bead") +
  theme_gse224445
pdf(file.path(figure_dir, "06_gse224445_doubletfinder_rpca_cd8_cxcr6_summary_plots.pdf"), width = 11, height = 5)
print(freq_plot | expr_plot)
dev.off()

sink(file.path(log_dir, "GSE224445_doubletfinder_rpca_sessionInfo.txt"))
print(sessionInfo())
sink()

message_stamp("Completed GSE224445 reproducible DoubletFinder + RPCA integration workflow.")
message_stamp("Key output: ", file.path(table_dir, "GSE224445_doubletfinder_rpca_cd8_cxcr6_by_group.csv"))
