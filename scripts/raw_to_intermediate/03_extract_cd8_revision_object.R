#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
input_rdata <- if (length(args) >= 1) args[[1]] else "intermediate_objects/merged_4_with_vdj_T_cells_filter_10_28_2025.RData"
output_dir <- if (length(args) >= 2) args[[2]] else "intermediate_objects"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

env <- new.env(parent = emptyenv())
if (!file.exists(input_rdata)) stop("Input RData not found: ", input_rdata, call. = FALSE)
load(input_rdata, envir = env)

if (exists("merged_all_T_batchcorrected_filter", envir = env)) {
  t_obj <- get("merged_all_T_batchcorrected_filter", envir = env)
} else {
  obj_names <- ls(env)
  is_seurat <- vapply(obj_names, function(x) inherits(get(x, envir = env), "Seurat"), logical(1))
  if (!any(is_seurat)) stop("No Seurat object found in input RData.", call. = FALSE)
  t_obj <- get(obj_names[which(is_seurat)[[1]]], envir = env)
}

annotation_col <- Sys.getenv("CD8_ANNOTATION_COLUMN", "cell_type_annotated_granularity")
if (!annotation_col %in% colnames(t_obj@meta.data)) {
  annotation_col <- Sys.getenv("FALLBACK_CD8_ANNOTATION_COLUMN", "cell_type_annotated")
}
if (!annotation_col %in% colnames(t_obj@meta.data)) {
  stop("No CD8 annotation column found. Set CD8_ANNOTATION_COLUMN.", call. = FALSE)
}

cd8_regex <- Sys.getenv("CD8_KEEP_REGEX", "CD8")
cells <- rownames(t_obj@meta.data)[
  grepl(cd8_regex, as.character(t_obj@meta.data[[annotation_col]]), ignore.case = TRUE)
]
if (length(cells) == 0) stop("No CD8 cells matched pattern: ", cd8_regex, call. = FALSE)

merged_all_T_batchcorrected_filter_CD8 <- subset(t_obj, cells = cells)
save(
  merged_all_T_batchcorrected_filter_CD8,
  file = file.path(output_dir, "merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "03_extract_cd8_sessionInfo.txt"))

