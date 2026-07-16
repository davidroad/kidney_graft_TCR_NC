#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- if (length(args) >= 1) args[[1]] else "scripts/raw_to_intermediate/sample_manifest.csv"
output_dir <- if (length(args) >= 2) args[[2]] else "intermediate_objects"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

stop_if_missing <- function(path, label) {
  if (is.na(path) || path == "" || !file.exists(path)) {
    stop(label, " not found: ", path, call. = FALSE)
  }
}

load_seurat_from_rdata <- function(rdata_file, object_name = "") {
  env <- new.env(parent = emptyenv())
  load(rdata_file, envir = env)
  if (!is.na(object_name) && object_name != "" && exists(object_name, envir = env)) {
    return(get(object_name, envir = env))
  }
  obj_names <- ls(env)
  is_seurat <- vapply(obj_names, function(x) inherits(get(x, envir = env), "Seurat"), logical(1))
  if (!any(is_seurat)) {
    stop("No Seurat object found in ", rdata_file, call. = FALSE)
  }
  get(obj_names[which(is_seurat)[[1]]], envir = env)
}

add_tcr_metadata <- function(seurat_obj, vdj_dir, prefix = "t") {
  if (is.na(vdj_dir) || vdj_dir == "") return(seurat_obj)

  contig_file <- file.path(vdj_dir, "filtered_contig_annotations.csv")
  clonotype_file <- file.path(vdj_dir, "clonotypes.csv")
  if (!file.exists(contig_file) || !file.exists(clonotype_file)) {
    warning("Skipping TCR metadata because VDJ files are missing in: ", vdj_dir)
    return(seurat_obj)
  }

  tcr <- read.csv(contig_file, stringsAsFactors = FALSE)
  clono <- read.csv(clonotype_file, stringsAsFactors = FALSE)
  required_contig <- c("barcode", "raw_clonotype_id")
  required_clono <- c("clonotype_id", "cdr3s_aa")
  if (!all(required_contig %in% colnames(tcr)) || !all(required_clono %in% colnames(clono))) {
    warning("Skipping TCR metadata because expected VDJ columns are missing in: ", vdj_dir)
    return(seurat_obj)
  }

  tcr <- tcr[!duplicated(tcr$barcode), c("barcode", "raw_clonotype_id")]
  colnames(tcr)[colnames(tcr) == "raw_clonotype_id"] <- "clonotype_id"
  tcr <- merge(tcr, clono[, c("clonotype_id", "cdr3s_aa")], by = "clonotype_id", all.x = TRUE)

  meta <- tcr[, c("clonotype_id", "cdr3s_aa")]
  colnames(meta) <- paste(prefix, colnames(meta), sep = "_")

  raw_barcode <- tcr$barcode
  no_suffix <- sub("-1$", "", raw_barcode)
  candidates <- list(raw_barcode, no_suffix, paste0(no_suffix, "-1"))
  best <- candidates[[which.max(vapply(candidates, function(x) sum(x %in% colnames(seurat_obj)), integer(1)))]]
  rownames(meta) <- best
  meta <- meta[intersect(rownames(meta), colnames(seurat_obj)), , drop = FALSE]
  if (nrow(meta) == 0) {
    warning("No matching barcodes between VDJ metadata and Seurat object for: ", vdj_dir)
    return(seurat_obj)
  }

  AddMetaData(seurat_obj, metadata = meta)
}

load_one_sample <- function(row) {
  sample_id <- row[["sample_id"]]
  if (!is.na(row[["existing_seurat_rdata"]]) && row[["existing_seurat_rdata"]] != "") {
    stop_if_missing(row[["existing_seurat_rdata"]], "existing_seurat_rdata")
    obj <- load_seurat_from_rdata(row[["existing_seurat_rdata"]], row[["existing_seurat_object"]])
  } else {
    stop_if_missing(row[["matrix_dir"]], "matrix_dir")
    counts <- Read10X(row[["matrix_dir"]])
    if (is.list(counts)) {
      counts <- if ("Gene Expression" %in% names(counts)) {
        counts[["Gene Expression"]]
      } else {
        counts[[1]]
      }
    }
    obj <- CreateSeuratObject(counts = counts, project = sample_id, min.cells = 3, min.features = 200)
  }

  obj$sample_id <- sample_id
  obj$patient_id <- row[["patient_id"]]
  obj$tissue <- row[["tissue"]]
  obj$batch <- row[["batch_label"]]
  obj <- add_tcr_metadata(obj, row[["vdj_dir"]], prefix = "t")
  obj[["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-")
  obj
}

required_cols <- c(
  "sample_id", "patient_id", "tissue", "batch_label",
  "matrix_dir", "vdj_dir", "existing_seurat_rdata", "existing_seurat_object"
)
manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
missing_cols <- setdiff(required_cols, colnames(manifest))
if (length(missing_cols) > 0) {
  stop("Manifest is missing columns: ", paste(missing_cols, collapse = ", "), call. = FALSE)
}

object_list <- lapply(seq_len(nrow(manifest)), function(i) {
  message("Loading sample: ", manifest$sample_id[[i]])
  obj <- load_one_sample(manifest[i, , drop = FALSE])
  obj <- NormalizeData(obj, verbose = FALSE)
  obj <- FindVariableFeatures(obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  obj
})
names(object_list) <- manifest$sample_id

merged_combined_regressed <- Reduce(function(x, y) merge(x, y), object_list)

features <- SelectIntegrationFeatures(object.list = object_list, nfeatures = 3000)
object_list <- lapply(object_list, function(obj) {
  obj <- ScaleData(obj, features = features, verbose = FALSE)
  obj <- RunPCA(obj, features = features, verbose = FALSE)
  obj
})

anchors <- FindIntegrationAnchors(
  object.list = object_list,
  anchor.features = features,
  reduction = "rpca",
  dims = 1:30
)

merged.all.batchcorrected <- IntegrateData(anchorset = anchors, dims = 1:30)
DefaultAssay(merged.all.batchcorrected) <- "integrated"
merged.all.batchcorrected <- ScaleData(merged.all.batchcorrected, verbose = FALSE)
merged.all.batchcorrected <- RunPCA(merged.all.batchcorrected, npcs = 50, verbose = FALSE)
merged.all.batchcorrected <- RunUMAP(merged.all.batchcorrected, dims = 1:30, reduction = "pca")
merged.all.batchcorrected <- FindNeighbors(merged.all.batchcorrected, dims = 1:30, verbose = FALSE)
merged.all.batchcorrected <- FindClusters(merged.all.batchcorrected, resolution = 0.5, verbose = FALSE)

save(
  merged_combined_regressed,
  merged.all.batchcorrected,
  file = file.path(output_dir, "merged_4_with_vdj_9_24_2025.RData")
)

writeLines(capture.output(sessionInfo()), file.path(output_dir, "01_build_all_cell_sessionInfo.txt"))
