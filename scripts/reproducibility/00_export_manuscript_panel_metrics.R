#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  idx <- match(flag, args)
  if (is.na(idx) || idx == length(args)) return(default)
  args[[idx + 1L]]
}

all_path <- arg_value("--all-cd45")
tcell_path <- arg_value("--tcell")
shared_path <- arg_value("--shared")
out_dir <- arg_value("--out", "inputs/plotting_inputs")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load_seurat <- function(path, preferred) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  env <- new.env(parent = emptyenv())
  loaded <- load(path, envir = env)
  candidates <- unique(c(preferred, loaded))
  for (name in candidates) {
    if (exists(name, envir = env, inherits = FALSE)) {
      x <- get(name, envir = env, inherits = FALSE)
      if (inherits(x, "Seurat")) return(x)
    }
  }
  stop("No Seurat object found in ", path)
}

first_column <- function(meta, candidates, default = NA_character_) {
  hit <- candidates[candidates %in% colnames(meta)][1]
  if (!length(hit) || is.na(hit)) rep(default, nrow(meta)) else as.character(meta[[hit]])
}

anonymous_ids <- function(n, prefix) sprintf("%s_%06d", prefix, seq_len(n))

write_csv <- function(x, name) {
  path <- file.path(out_dir, name)
  con <- if (grepl("\\.gz$", path)) gzfile(path, open = "wt") else path
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  write.csv(x, con, row.names = FALSE, quote = TRUE)
}

anonymize_context <- function(meta) {
  tissue <- first_column(meta, c("label", "group"))
  raw_sample <- first_column(meta, c("batch", "orig.ident"))
  raw_patient <- first_column(meta, c("patient"))
  missing_patient <- is.na(raw_patient) | !nzchar(raw_patient)
  raw_patient[missing_patient] <- sub(".*?([0-9]{4}).*", "\\1", raw_sample[missing_patient])
  patient_levels <- sort(unique(raw_patient[!is.na(raw_patient) & nzchar(raw_patient)]))
  patient_map <- setNames(sprintf("P%d", seq_along(patient_levels)), patient_levels)
  patient <- unname(patient_map[raw_patient])
  sample <- ifelse(!is.na(patient) & !is.na(tissue), paste(patient, tissue, sep = "_"), raw_sample)
  list(patient = patient, sample = sample, tissue = tissue)
}

export_embedding <- function(object, prefix, include_expression = character()) {
  emb <- Embeddings(object, reduction = "umap")
  meta <- object@meta.data[rownames(emb), , drop = FALSE]
  context <- anonymize_context(meta)
  out <- data.frame(
    cell_id = anonymous_ids(nrow(emb), prefix),
    UMAP_1 = emb[, 1],
    UMAP_2 = emb[, 2],
    cluster = first_column(meta, c("cell_type_annotated_granularity", "cell_type_annotated", "integrated_snn_res.0.5", "seurat_clusters")),
    cluster_id = first_column(meta, c("integrated_snn_res.0.5", "seurat_clusters")),
    tissue = context$tissue,
    sample = context$sample,
    patient = context$patient,
    stringsAsFactors = FALSE
  )
  genes <- intersect(include_expression, rownames(object))
  if (length(genes)) {
    expr <- FetchData(object, vars = genes, layer = "data")
    expr <- expr[rownames(emb), genes, drop = FALSE]
    out <- cbind(out, as.data.frame(expr, check.names = FALSE))
  }
  out
}

export_proportions <- function(embedding, prefix) {
  keep <- !is.na(embedding$sample) & !is.na(embedding$cluster)
  counts <- aggregate(
    list(n_cells = rep(1L, sum(keep))),
    embedding[keep, c("patient", "sample", "tissue", "cluster", "cluster_id"), drop = FALSE],
    sum
  )
  totals <- aggregate(n_cells ~ sample, counts, sum)
  names(totals)[2] <- "sample_total_cells"
  counts <- merge(counts, totals, by = "sample", all.x = TRUE, sort = FALSE)
  counts$proportion <- counts$n_cells / counts$sample_total_cells
  write_csv(counts, paste0(prefix, "_sample_cluster_proportions.csv"))
}

all_object <- load_seurat(all_path, c("merged.all.batchcorrected", "merged_all_batchcorrected"))
if (!is.null(all_object)) {
  all_embedding <- export_embedding(all_object, "cd45")
  write_csv(all_embedding, "fig1_all_cd45_umap_annotations.csv.gz")
  export_proportions(all_embedding, "fig1_all_cd45")

  marker_genes <- c(
    "LYZ", "CD14", "TCF7", "LEF1", "CD3D", "CD8A", "TCL1A", "MS4A1",
    "GNLY", "NKG7", "LAG3", "GZMB", "FOXP3", "CTLA4", "TRDV1", "TRDC",
    "TSHZ2", "RALGPS2", "XCL1", "XCL2", "CDKN1C", "CX3CR1", "JCHAIN",
    "MZB1", "C1QA", "C1QB", "CSF3R", "NAMPT", "TYMS", "STMN1", "AFF3",
    "IGHM", "IGHD", "CD1C", "FCER1A", "IL3RA", "CLEC4C", "HBB", "ALAS2",
    "CPA3", "HDC"
  )
  marker_genes <- intersect(marker_genes, rownames(all_object))
  group <- factor(all_embedding$cluster)
  expr <- FetchData(all_object, vars = marker_genes, layer = "data")
  dot_rows <- lapply(levels(group), function(level) {
    idx <- group == level
    data.frame(
      cluster = level,
      gene = marker_genes,
      average_expression = colMeans(expr[idx, marker_genes, drop = FALSE]),
      percent_expressing = colMeans(expr[idx, marker_genes, drop = FALSE] > 0) * 100,
      stringsAsFactors = FALSE
    )
  })
  dot <- do.call(rbind, dot_rows)
  dot$scaled_average_expression <- ave(dot$average_expression, dot$gene, FUN = function(x) {
    z <- as.numeric(scale(x)); z[is.na(z)] <- 0; z
  })
  write_csv(dot, "fig1d_marker_dotplot_metrics.csv")
  rm(all_object, expr); gc()
}

tcell_object <- load_seurat(tcell_path, c("merged_all_T_batchcorrected_filter"))
if (!is.null(tcell_object)) {
  tcell_embedding <- export_embedding(tcell_object, "tcell")
  write_csv(tcell_embedding, "fig2_tcell_umap_annotations.csv.gz")
  export_proportions(tcell_embedding, "fig2_tcell")
  rm(tcell_embedding); gc()
}

shared_genes <- c("TCF7", "CXCR6", "GZMB", "IFNG", "PRF1", "TOX", "PDCD1", "HAVCR2", "KLF2", "CISH", "JUN", "MKI67", "CDKN2A")
shared_object <- load_seurat(shared_path, c("all_clone_test_filter"))
if (!is.null(shared_object)) {
  shared_embedding <- export_embedding(shared_object, "shared", shared_genes)
  write_csv(shared_embedding, "fig2_shared_clone_umap_expression.csv.gz")

  meta <- shared_object@meta.data
  context <- anonymize_context(meta)
  clone_col <- c("CTaa", "t_cdr3s_aa")[c("CTaa", "t_cdr3s_aa") %in% colnames(meta)][1]
  if (!is.na(clone_col)) {
    clone <- data.frame(
      patient = context$patient,
      sample = context$sample,
      tissue = context$tissue,
      clonotype = as.character(meta[[clone_col]]),
      stringsAsFactors = FALSE
    )
    clone <- clone[!is.na(clone$clonotype) & nzchar(clone$clonotype), , drop = FALSE]
    abundance <- aggregate(list(n_cells = rep(1L, nrow(clone))), clone, sum)
    totals <- aggregate(n_cells ~ sample, abundance, sum)
    names(totals)[2] <- "sample_total_shared_clone_cells"
    abundance <- merge(abundance, totals, by = "sample", all.x = TRUE, sort = FALSE)
    abundance$relative_frequency <- abundance$n_cells / abundance$sample_total_shared_clone_cells
    write_csv(abundance, "fig2c_shared_clonotype_abundance.csv")
  }

  Idents(shared_object) <- factor(first_column(meta, c("label", "group")), levels = c("PBMC", "Graft"))
  if (all(c("PBMC", "Graft") %in% levels(Idents(shared_object)))) {
    de <- FindMarkers(shared_object, ident.1 = "Graft", ident.2 = "PBMC", min.pct = 0.25, logfc.threshold = 0)
    de$gene <- rownames(de)
    write_csv(de, "fig2f_shared_clone_graft_vs_pbmc_de.csv")
  }
}

message("Metrics export complete: ", normalizePath(out_dir))
