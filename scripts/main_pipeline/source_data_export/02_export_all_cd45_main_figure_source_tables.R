args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript 02_export_all_cd45_main_figure_source_tables.R <all_cd45_seurat.RData> <output_dir>",
    call. = FALSE
  )
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_path <- tryCatch(normalizePath(script_path, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)

if (!is.na(script_path)) {
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE)
  local_lib <- normalizePath(file.path(repo_root, "..", "Rlib"), winslash = "/", mustWork = FALSE)
  if (dir.exists(local_lib)) .libPaths(c(local_lib, .libPaths()))
}

suppressPackageStartupMessages(library(SeuratObject))

input_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading all-CD45 object: ", input_file)
env <- new.env(parent = emptyenv())
loaded <- load(input_file, envir = env)

is_seurat <- vapply(loaded, function(nm) "Seurat" %in% class(get(nm, envir = env)), logical(1))
seurat_names <- loaded[is_seurat]
if (!length(seurat_names)) stop("No Seurat object found in input RData.", call. = FALSE)

object_name <- if ("merged.all.batchcorrected" %in% seurat_names) {
  "merged.all.batchcorrected"
} else {
  seurat_names[[1]]
}
seu <- get(object_name, envir = env)
message("Using Seurat object: ", object_name)

write_tsv <- function(x, path, compress = FALSE) {
  con <- if (compress) gzfile(path, open = "wt") else file(path, open = "wt")
  on.exit(close(con), add = TRUE)
  write.table(x, file = con, sep = "\t", quote = FALSE, row.names = FALSE, na = "")
}

first_existing <- function(df, candidates) {
  hit <- candidates[candidates %in% colnames(df)]
  if (length(hit)) hit[[1]] else NA_character_
}

get_expr_matrix <- function(seu, genes, assay = "RNA") {
  assay_obj <- seu@assays[[assay]]
  genes_present <- genes[genes %in% rownames(assay_obj)]
  out <- matrix(NA_real_, nrow = length(genes), ncol = ncol(seu))
  rownames(out) <- genes
  colnames(out) <- colnames(seu)

  if (!length(genes_present)) return(out)

  mat <- tryCatch(
    GetAssayData(seu, assay = assay, layer = "data"),
    error = function(e) tryCatch(GetAssayData(seu, assay = assay, slot = "data"), error = function(e2) NULL)
  )
  if (is.null(mat)) {
    mat <- tryCatch(
      GetAssayData(seu, assay = assay, layer = "counts"),
      error = function(e) tryCatch(GetAssayData(seu, assay = assay, slot = "counts"), error = function(e2) NULL)
    )
  }
  if (is.null(mat)) stop("Could not retrieve RNA data or counts matrix.", call. = FALSE)

  out[genes_present, ] <- as.matrix(mat[genes_present, colnames(seu), drop = FALSE])
  out
}

md <- seu@meta.data
umap <- as.data.frame(seu@reductions$umap@cell.embeddings)
if (ncol(umap) < 2) stop("UMAP reduction has fewer than two dimensions.", call. = FALSE)
colnames(umap)[1:2] <- c("UMAP_1", "UMAP_2")

cell_ids <- rownames(md)
public_cell_id <- sprintf("CD45cell_%06d", seq_along(cell_ids))
names(public_cell_id) <- cell_ids

batch_col <- first_existing(md, c("batch", "group", "sample", "orig.ident"))
tissue_col <- first_existing(md, c("label", "tissue", "Tissue", "origin"))
cluster_col <- first_existing(md, c("seurat_clusters", "integrated_snn_res.0.5", "RNA_snn_res.0.5"))
annotation_col <- first_existing(md, c("cell_type_annotated", "cell_type_annotated_granularity", "cell_type_annotated_previous"))

batch_raw <- if (!is.na(batch_col)) as.character(md[[batch_col]]) else rep("sample", nrow(md))
tissue <- if (!is.na(tissue_col)) as.character(md[[tissue_col]]) else ifelse(grepl("graft", batch_raw, ignore.case = TRUE), "Graft", "PBMC")
cluster_id <- if (!is.na(cluster_col)) as.character(md[[cluster_col]]) else NA_character_
cell_type <- if (!is.na(annotation_col)) as.character(md[[annotation_col]]) else NA_character_

subject_raw <- sub("^.*_", "", batch_raw)
subject_raw[subject_raw == batch_raw] <- batch_raw[subject_raw == batch_raw]
subject_levels <- sort(unique(subject_raw[!is.na(subject_raw)]))
subject_map <- setNames(sprintf("P%02d", seq_along(subject_levels)), subject_levels)
subject_id <- unname(subject_map[subject_raw])
subject_id[is.na(subject_id)] <- "P00"
sample_id <- paste0(tissue, "_", subject_id)

all_cell_umap <- data.frame(
  cell_id = public_cell_id[cell_ids],
  subject_id = subject_id,
  sample_id = sample_id,
  tissue = tissue,
  cluster_id = cluster_id,
  cell_type = cell_type,
  UMAP_1 = umap[cell_ids, "UMAP_1"],
  UMAP_2 = umap[cell_ids, "UMAP_2"],
  stringsAsFactors = FALSE
)

write_tsv(all_cell_umap, file.path(output_dir, "Fig1BCE_umap.tsv.gz"), compress = TRUE)

freq_input <- all_cell_umap
freq_input$cell_type[is.na(freq_input$cell_type) | freq_input$cell_type == ""] <- "unannotated"
freq_counts <- aggregate(
  cell_id ~ sample_id + subject_id + tissue + cluster_id + cell_type,
  data = freq_input,
  FUN = length
)
colnames(freq_counts)[ncol(freq_counts)] <- "n_cells"
sample_totals <- aggregate(
  cell_id ~ sample_id + subject_id + tissue,
  data = freq_input,
  FUN = length
)
colnames(sample_totals)[ncol(sample_totals)] <- "total_cd45_cells"
freq_counts <- merge(freq_counts, sample_totals, by = c("sample_id", "subject_id", "tissue"), all.x = TRUE)
freq_counts$proportion <- freq_counts$n_cells / freq_counts$total_cd45_cells
freq_counts <- freq_counts[order(freq_counts$sample_id, freq_counts$cluster_id), ]
write_tsv(freq_counts, file.path(output_dir, "Fig1FG_celltype_freq.tsv"))

marker_genes <- c(
  "HDC", "CPA3", "ALAS2", "HBB", "CLEC4C", "IL3RA", "FCER1A", "CD1C",
  "IGHD", "IGHM", "AFF3", "STMN1", "TYMS", "NAMPT", "CSF3R", "C1QB",
  "C1QA", "MZB1", "JCHAIN", "CX3CR1", "CDKN1C", "XCL2", "XCL1",
  "RALGPS2", "TSHZ2", "TRDC", "TRDV1", "CTLA4", "FOXP3", "GZMB",
  "LAG3", "NKG7", "GNLY", "MS4A1", "TCL1A", "CD8A", "CD3D", "LEF1",
  "TCF7", "CD14", "LYZ"
)

expr <- get_expr_matrix(seu, marker_genes, assay = "RNA")
group <- ifelse(is.na(cell_type) | cell_type == "", paste0("cluster_", cluster_id), cell_type)
group_levels <- unique(group)

dotplot <- do.call(
  rbind,
  lapply(marker_genes, function(gene) {
    vals <- expr[gene, ]
    do.call(
      rbind,
      lapply(group_levels, function(grp) {
        idx <- group == grp
        data.frame(
          gene = gene,
          cell_type = grp,
          n_cells = sum(idx),
          avg_expression = if (all(is.na(vals[idx]))) NA_real_ else mean(vals[idx], na.rm = TRUE),
          pct_expressing = if (all(is.na(vals[idx]))) NA_real_ else 100 * mean(vals[idx] > 0, na.rm = TRUE),
          stringsAsFactors = FALSE
        )
      })
    )
  })
)
write_tsv(dotplot, file.path(output_dir, "Fig1D_dotplot.tsv"))

readme <- c(
  "# All-CD45 main-figure source data",
  "",
  "These tables were exported from the all-CD45 Seurat RData object using `scripts/source_data_export/02_export_all_cd45_main_figure_source_tables.R`.",
  "",
  "Identifiers were neutralized before export:",
  "",
  "- original sample/patient IDs were converted to `P01`, `P02`, ...",
  "- sample IDs were converted to `PBMC_Pxx` or `Graft_Pxx`",
  "- original cell barcodes were converted to `CD45cell_000001`, `CD45cell_000002`, ...",
  "",
  "Files:",
  "",
  "- `Fig1BCE_umap.tsv.gz`: cell-level UMAP coordinates, neutralized sample/tissue labels, cluster IDs, and annotated cell types for Fig. 1B, Fig. 1C, and Fig. 1E.",
  "- `Fig1D_dotplot.tsv`: average expression and percent expressing values for the selected marker genes shown in the all-CD45 dot plot.",
  "- `Fig1FG_celltype_freq.tsv`: sample-level cell counts and proportions by cluster and annotated cell type for Fig. 1F and Fig. 1G.",
  "",
  "The marker gene list was taken from the first-round manuscript figure text for the all-CD45 marker dot plot."
)
writeLines(readme, file.path(output_dir, "README.md"), useBytes = TRUE)

message("Wrote source tables to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
