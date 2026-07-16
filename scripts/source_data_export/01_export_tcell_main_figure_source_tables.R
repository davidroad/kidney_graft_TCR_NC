args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    "Usage: Rscript 01_export_tcell_main_figure_source_tables.R <tcell_seurat.RData> <output_dir>",
    call. = FALSE
  )
}

file_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_path <- if (length(file_arg)) sub("^--file=", "", file_arg[[1]]) else NA_character_
script_path <- tryCatch(normalizePath(script_path, winslash = "/", mustWork = TRUE), error = function(e) NA_character_)

if (!is.na(script_path)) {
  repo_root <- normalizePath(file.path(dirname(script_path), "..", ".."), winslash = "/", mustWork = FALSE)
  local_lib <- file.path(dirname(repo_root), "Rlib")
  prep_lib <- normalizePath(file.path(repo_root, "..", "Rlib"), winslash = "/", mustWork = FALSE)
  for (lib in unique(c(local_lib, prep_lib))) {
    if (dir.exists(lib)) .libPaths(c(lib, .libPaths()))
  }
}

suppressPackageStartupMessages(library(SeuratObject))

input_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Loading T-cell object: ", input_file)
env <- new.env(parent = emptyenv())
loaded <- load(input_file, envir = env)

is_seurat <- vapply(loaded, function(nm) inherits(get(nm, envir = env), "Seurat"), logical(1))
if (!any(is_seurat)) stop("No Seurat object found in input RData.", call. = FALSE)

object_name <- loaded[which(is_seurat)[1]]
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
  genes_present <- genes[genes %in% rownames(seu[[assay]])]
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
public_cell_id <- sprintf("Tcell_%06d", seq_along(cell_ids))
names(public_cell_id) <- cell_ids

patient_col <- first_existing(md, c("patient", "Patient", "donor", "subject"))
tissue_col <- first_existing(md, c("label", "tissue", "Tissue", "origin"))
batch_col <- first_existing(md, c("batch", "group", "sample", "orig.ident"))
cluster_col <- first_existing(md, c("seurat_clusters", "integrated_snn_res.0.5", "RNA_snn_res.0.5"))
state_col <- first_existing(md, c("cell_type_annotated_granularity", "cell_type_annotated_previous"))
state_previous_col <- first_existing(md, c("cell_type_annotated_previous"))
clonotype_col <- first_existing(md, c("t_clonotype_id", "clonotype_id", "raw_clonotype_id"))
cdr3_col <- first_existing(md, c("t_cdr3s_aa", "CTaa"))
clone_size_col <- first_existing(md, c("cloneSize"))
clonal_freq_col <- first_existing(md, c("clonalFrequency"))
strict_clone_col <- first_existing(md, c("CTstrict", "t_cdr3s_aa", "CTaa", "t_clonotype_id"))

patient_raw <- if (!is.na(patient_col)) as.character(md[[patient_col]]) else rep("subject", nrow(md))
tissue <- if (!is.na(tissue_col)) as.character(md[[tissue_col]]) else rep(NA_character_, nrow(md))
batch_raw <- if (!is.na(batch_col)) as.character(md[[batch_col]]) else rep(NA_character_, nrow(md))

patient_levels <- sort(unique(patient_raw[!is.na(patient_raw)]))
patient_map <- setNames(sprintf("P%02d", seq_along(patient_levels)), patient_levels)
subject_id <- unname(patient_map[patient_raw])
subject_id[is.na(subject_id)] <- "P00"

sample_id <- ifelse(
  is.na(tissue) | tissue == "",
  paste0("Sample_", subject_id),
  paste0(tissue, "_", subject_id)
)

cluster_id <- if (!is.na(cluster_col)) as.character(md[[cluster_col]]) else NA_character_
tcell_state <- if (!is.na(state_col)) as.character(md[[state_col]]) else NA_character_
tcell_state_previous <- if (!is.na(state_previous_col)) as.character(md[[state_previous_col]]) else NA_character_
clonotype_raw <- if (!is.na(clonotype_col)) as.character(md[[clonotype_col]]) else NA_character_
cdr3_aa <- if (!is.na(cdr3_col)) as.character(md[[cdr3_col]]) else NA_character_
clone_size <- if (!is.na(clone_size_col)) as.character(md[[clone_size_col]]) else NA_character_
clonal_frequency <- if (!is.na(clonal_freq_col)) as.integer(md[[clonal_freq_col]]) else NA_integer_
strict_clone <- if (!is.na(strict_clone_col)) as.character(md[[strict_clone_col]]) else clonotype_raw

valid_clone <- !is.na(strict_clone) & strict_clone != "" & strict_clone != "NA"
clone_pair_key <- ifelse(valid_clone, paste(subject_id, strict_clone, sep = "||"), NA_character_)
unique_clone_keys <- sort(unique(clone_pair_key[!is.na(clone_pair_key)]))
public_clone_map <- setNames(sprintf("clonotype_%06d", seq_along(unique_clone_keys)), unique_clone_keys)
public_clonotype_id <- unname(public_clone_map[clone_pair_key])

clone_tissue_list <- tapply(
  tissue[!is.na(clone_pair_key)],
  clone_pair_key[!is.na(clone_pair_key)],
  function(x) sort(unique(x[!is.na(x) & x != ""]))
)
shared_clone_keys <- names(clone_tissue_list)[
  vapply(clone_tissue_list, function(x) all(c("PBMC", "Graft") %in% x), logical(1))
]
shared_between_pbmc_graft <- clone_pair_key %in% shared_clone_keys
shared_between_pbmc_graft[is.na(shared_between_pbmc_graft)] <- FALSE

genes <- c("TCF7", "CXCR6", "GZMB")
expr <- get_expr_matrix(seu, genes, assay = "RNA")
expr_df <- as.data.frame(t(expr))
expr_df <- expr_df[cell_ids, genes, drop = FALSE]

base_cell_table <- data.frame(
  cell_id = public_cell_id[cell_ids],
  subject_id = subject_id,
  sample_id = sample_id,
  tissue = tissue,
  cluster_id = cluster_id,
  tcell_state = tcell_state,
  tcell_state_previous = tcell_state_previous,
  clonotype_id = public_clonotype_id,
  cdr3_aa = cdr3_aa,
  clone_size = clone_size,
  clonal_frequency = clonal_frequency,
  shared_between_pbmc_graft = shared_between_pbmc_graft,
  UMAP_1 = umap[cell_ids, "UMAP_1"],
  UMAP_2 = umap[cell_ids, "UMAP_2"],
  TCF7 = expr_df$TCF7,
  CXCR6 = expr_df$CXCR6,
  GZMB = expr_df$GZMB,
  stringsAsFactors = FALSE
)

write_tsv(
  base_cell_table,
  file.path(output_dir, "Fig2A_G_umap_features.tsv.gz"),
  compress = TRUE
)

cluster_input <- base_cell_table
cluster_input$tcell_state[is.na(cluster_input$tcell_state) | cluster_input$tcell_state == ""] <- "unannotated"
cluster_counts <- aggregate(
  cell_id ~ sample_id + subject_id + tissue + cluster_id + tcell_state,
  data = cluster_input,
  FUN = length
)
colnames(cluster_counts)[ncol(cluster_counts)] <- "n_cells"
sample_totals <- aggregate(
  cell_id ~ sample_id + subject_id + tissue,
  data = cluster_input,
  FUN = length
)
colnames(sample_totals)[ncol(sample_totals)] <- "total_t_cells"
cluster_counts <- merge(cluster_counts, sample_totals, by = c("sample_id", "subject_id", "tissue"), all.x = TRUE)
cluster_counts$proportion <- cluster_counts$n_cells / cluster_counts$total_t_cells
cluster_counts <- cluster_counts[order(cluster_counts$sample_id, cluster_counts$cluster_id), ]
write_tsv(cluster_counts, file.path(output_dir, "Fig2B_tcell_frequency.tsv"))

shared_cell_table <- base_cell_table[base_cell_table$shared_between_pbmc_graft, ]
write_tsv(
  shared_cell_table,
  file.path(output_dir, "Fig2DEH_shared_umap.tsv.gz"),
  compress = TRUE
)

clone_input <- base_cell_table[!is.na(base_cell_table$clonotype_id) & base_cell_table$clonotype_id != "", ]
clone_counts <- aggregate(
  cell_id ~ clonotype_id + subject_id + sample_id + tissue + cdr3_aa + shared_between_pbmc_graft,
  data = clone_input,
  FUN = length
)
colnames(clone_counts)[ncol(clone_counts)] <- "n_cells"
clone_counts <- clone_counts[order(clone_counts$subject_id, clone_counts$clonotype_id, clone_counts$tissue), ]
write_tsv(clone_counts, file.path(output_dir, "Fig2CDE_clonotype_counts.tsv"))

feature_table <- base_cell_table[, c(
  "cell_id", "subject_id", "sample_id", "tissue", "cluster_id", "tcell_state",
  "UMAP_1", "UMAP_2", genes
)]
write_tsv(
  feature_table,
  file.path(output_dir, "Fig2G_feature_expr.tsv.gz"),
  compress = TRUE
)

readme <- c(
  "# T-cell main-figure source data",
  "",
  "These tables were exported from the T-cell Seurat RData object using `scripts/source_data_export/01_export_tcell_main_figure_source_tables.R`.",
  "",
  "Identifiers were neutralized before export:",
  "",
  "- original patient IDs were converted to `P01`, `P02`, ...",
  "- sample IDs were converted to `PBMC_Pxx` or `Graft_Pxx`",
  "- original cell barcodes were converted to `Tcell_000001`, `Tcell_000002`, ...",
  "- clonotypes were converted to public `clonotype_000001`, `clonotype_000002`, ... IDs within subject-clonotype pairs",
  "",
  "Files:",
  "",
  "- `Fig2A_G_umap_features.tsv.gz`: cell-level UMAP coordinates, T-cell annotations, TCR summary fields, and TCF7/CXCR6/GZMB expression.",
  "- `Fig2B_tcell_frequency.tsv`: sample-level T-cell counts and proportions by cluster/state.",
  "- `Fig2CDE_clonotype_counts.tsv`: clonotype counts by sample/tissue.",
  "- `Fig2DEH_shared_umap.tsv.gz`: cell-level source data for shared PBMC/graft clonotype UMAP panels.",
  "- `Fig2G_feature_expr.tsv.gz`: compact cell-level source data for TCF7, CXCR6, and GZMB feature plots.",
  "",
  "The `shared_between_pbmc_graft` flag is defined within each subject as a clonotype observed in both PBMC and graft compartments."
)
writeLines(readme, file.path(output_dir, "README.md"), useBytes = TRUE)

message("Wrote source tables to: ", normalizePath(output_dir, winslash = "/", mustWork = FALSE))
