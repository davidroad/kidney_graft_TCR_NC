#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
input_rdata <- if (length(args) >= 1) args[[1]] else "intermediate_objects/merged_4_with_vdj_9_24_2025.RData"
output_dir <- if (length(args) >= 2) args[[2]] else "intermediate_objects"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_rdata_env <- function(path) {
  if (!file.exists(path)) stop("Input RData not found: ", path, call. = FALSE)
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  env
}

get_first_seurat <- function(env, preferred) {
  for (name in preferred) {
    if (exists(name, envir = env) && inherits(get(name, envir = env), "Seurat")) {
      return(get(name, envir = env))
    }
  }
  obj_names <- ls(env)
  is_seurat <- vapply(obj_names, function(x) inherits(get(x, envir = env), "Seurat"), logical(1))
  if (!any(is_seurat)) stop("No Seurat object found in input RData.", call. = FALSE)
  get(obj_names[which(is_seurat)[[1]]], envir = env)
}

first_existing <- function(x, candidates) {
  hit <- candidates[candidates %in% x]
  if (length(hit) == 0) return(NULL)
  hit[[1]]
}

subset_by_regex <- function(obj, column, pattern) {
  values <- as.character(obj@meta.data[[column]])
  cells <- rownames(obj@meta.data)[grepl(pattern, values, ignore.case = TRUE)]
  if (length(cells) == 0) {
    stop("No cells matched ", column, " with pattern: ", pattern, call. = FALSE)
  }
  subset(obj, cells = cells)
}

reintegrate_by_sample <- function(obj) {
  sample_col <- first_existing(colnames(obj@meta.data), c("sample_id", "batch", "orig.ident"))
  if (is.null(sample_col)) return(obj)

  object_list <- SplitObject(obj, split.by = sample_col)
  if (length(object_list) < 2) return(obj)

  object_list <- lapply(object_list, function(x) {
    DefaultAssay(x) <- "RNA"
    x <- NormalizeData(x, verbose = FALSE)
    x <- FindVariableFeatures(x, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
    x
  })
  features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
  object_list <- lapply(object_list, function(x) {
    x <- ScaleData(x, features = features, verbose = FALSE)
    x <- RunPCA(x, features = features, verbose = FALSE)
    x
  })
  anchors <- FindIntegrationAnchors(
    object.list = object_list,
    anchor.features = features,
    reduction = "rpca",
    dims = 1:30
  )
  integrated <- IntegrateData(anchorset = anchors, dims = 1:30)
  DefaultAssay(integrated) <- "integrated"
  integrated <- ScaleData(integrated, verbose = FALSE)
  integrated <- RunPCA(integrated, npcs = 50, verbose = FALSE)
  integrated <- RunUMAP(integrated, dims = 1:30, reduction = "pca")
  integrated <- FindNeighbors(integrated, dims = 1:20, verbose = FALSE)
  integrated <- FindClusters(integrated, resolution = 0.5, verbose = FALSE)
  integrated
}

env <- load_rdata_env(input_rdata)
all_cell_obj <- get_first_seurat(env, c("merged.all.batchcorrected", "merged_combined_regressed"))

annotation_col <- Sys.getenv("ALL_CELL_ANNOTATION_COLUMN", "cell_type_annotated")
if (!annotation_col %in% colnames(all_cell_obj@meta.data)) {
  stop("Annotation column not found: ", annotation_col, call. = FALSE)
}

t_cell_regex <- Sys.getenv("T_CELL_KEEP_REGEX", "CD4|CD8|Treg|T cell")
message("Extracting T cells using ", annotation_col, " pattern: ", t_cell_regex)
t_cell_obj <- subset_by_regex(all_cell_obj, annotation_col, t_cell_regex)

run_reintegration <- tolower(Sys.getenv("REINTEGRATE_T_CELLS", "true")) %in% c("1", "true", "yes")
merged_all_T_batchcorrected <- if (run_reintegration) reintegrate_by_sample(t_cell_obj) else t_cell_obj

save(
  merged_all_T_batchcorrected,
  file = file.path(output_dir, "merged_4_with_vdj_T_cells_10_8_2025.RData")
)

filter_col <- Sys.getenv("T_CELL_FILTER_COLUMN", "cell_type_annotated_granularity")
filter_regex <- Sys.getenv("T_CELL_FILTER_REGEX", "CD4|CD8|Treg|T cell")
if (filter_col %in% colnames(merged_all_T_batchcorrected@meta.data)) {
  merged_all_T_batchcorrected_filter <- subset_by_regex(
    merged_all_T_batchcorrected,
    filter_col,
    filter_regex
  )
} else {
  merged_all_T_batchcorrected_filter <- merged_all_T_batchcorrected
  warning("Filter column not found; saving unfiltered T-cell object: ", filter_col)
}

meta <- merged_all_T_batchcorrected_filter@meta.data
clone_col <- first_existing(colnames(meta), c("t_cdr3s_aa", "cdr3s_aa", "t_clonotype_id", "clonotype_id"))
tissue_col <- first_existing(colnames(meta), c("tissue", "label", "sample_type", "batch"))
sample_col <- first_existing(colnames(meta), c("sample_id", "batch", "orig.ident"))

if (is.null(clone_col) || is.null(tissue_col)) {
  warning("Cannot build shared-clone subset because clone or tissue metadata are missing.")
  combined <- data.frame()
  save(combined, merged_all_T_batchcorrected_filter,
       file = file.path(output_dir, "merged_4_with_vdj_T_cells_filter_10_28_2025.RData"))
} else {
  clone_key <- as.character(meta[[clone_col]])
  tissue <- as.character(meta[[tissue_col]])
  sample <- if (!is.null(sample_col)) as.character(meta[[sample_col]]) else rownames(meta)
  valid <- !is.na(clone_key) & clone_key != "" & clone_key != "NA"

  clone_frame <- data.frame(
    cell = rownames(meta),
    clone_key = clone_key,
    tissue = tissue,
    sample = sample,
    stringsAsFactors = FALSE
  )
  clone_frame <- clone_frame[valid, , drop = FALSE]
  clone_ids <- unique(clone_frame$clone_key)
  shared_clone_summary <- do.call(rbind, lapply(clone_ids, function(id) {
    rows <- clone_frame[clone_frame$clone_key == id, , drop = FALSE]
    data.frame(
      clone_key = id,
      n_cells = nrow(rows),
      n_samples = length(unique(rows$sample)),
      has_graft = any(grepl("graft", rows$tissue, ignore.case = TRUE)),
      has_pbmc = any(grepl("pbmc|blood", rows$tissue, ignore.case = TRUE)),
      stringsAsFactors = FALSE
    )
  }))

  if (is.null(shared_clone_summary) || nrow(shared_clone_summary) == 0) {
    combined <- data.frame()
    warning("No valid clonotype metadata found for shared-clone subset.")
  } else {
    shared_ids <- shared_clone_summary$clone_key[
      shared_clone_summary$has_graft & shared_clone_summary$has_pbmc
    ]
    shared_cells <- clone_frame$cell[clone_frame$clone_key %in% shared_ids]
    combined <- shared_clone_summary
    if (length(shared_cells) > 0) {
      all_clone_test_filter <- subset(merged_all_T_batchcorrected_filter, cells = shared_cells)
      save(all_clone_test_filter, file = file.path(output_dir, "all_clone_filter_10_28.Rdata"))
    } else {
      warning("No shared graft/PBMC clonotypes were detected.")
    }
  }

  save(combined, merged_all_T_batchcorrected_filter,
       file = file.path(output_dir, "merged_4_with_vdj_T_cells_filter_10_28_2025.RData"))
}

writeLines(capture.output(sessionInfo()), file.path(output_dir, "02_extract_t_cells_sessionInfo.txt"))

