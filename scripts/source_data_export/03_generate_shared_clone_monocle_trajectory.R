args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop(
    paste(
      "Usage:",
      "Rscript scripts/source_data_export/03_generate_shared_clone_monocle_trajectory.R",
      "<all_clone_filter_10_28.Rdata> <output_dir>"
    ),
    call. = FALSE
  )
}

suppressPackageStartupMessages({
  library(Seurat)
  library(monocle3)
  library(ggplot2)
})

input_file <- normalizePath(args[[1]], winslash = "/", mustWork = TRUE)
output_dir <- args[[2]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

env <- new.env(parent = emptyenv())
loaded <- load(input_file, envir = env)

if (!"all_clone_test_filter" %in% loaded) {
  stop("Expected object `all_clone_test_filter` was not found in the input RData.", call. = FALSE)
}

shared_clone_seurat <- get("all_clone_test_filter", envir = env)

get_counts <- function(seu) {
  tryCatch(
    GetAssayData(seu, assay = "RNA", layer = "counts"),
    error = function(e) GetAssayData(seu, assay = "RNA", slot = "counts")
  )
}

counts <- get_counts(shared_clone_seurat)
cell_metadata <- shared_clone_seurat@meta.data
gene_metadata <- data.frame(gene_short_name = rownames(counts), row.names = rownames(counts))
stopifnot(identical(colnames(counts), rownames(cell_metadata)))

cds <- new_cell_data_set(
  counts,
  cell_metadata = cell_metadata,
  gene_metadata = gene_metadata
)

cds <- preprocess_cds(cds, num_dim = 50)
cds <- reduce_dimension(cds, preprocess_method = "PCA")

seurat_umap <- Embeddings(shared_clone_seurat, reduction = "umap")
cds_cells <- colnames(cds)
seurat_umap_aligned <- seurat_umap[cds_cells, , drop = FALSE]
stopifnot(!any(is.na(seurat_umap_aligned)))

cds@int_colData$reducedDims$UMAP <- seurat_umap_aligned

# These parameters document the Monocle 3 workflow used for the shared-clone
# pseudotime panel. The original run ordered cells interactively in Monocle 3;
# the start/root point was selected from the TCF7-defined stem-like T-cell
# state. Users regenerating the exact panel should root the trajectory in this
# state before ordering progression toward differentiated effector states.
cds <- cluster_cells(cds, resolution = 0.00105)
cds <- learn_graph(cds)
cds <- order_cells(cds)

cell_out <- as.data.frame(colData(cds))
umap_out <- as.data.frame(reducedDims(cds)$UMAP)
colnames(umap_out)[1:2] <- c("UMAP_1", "UMAP_2")

source_data <- data.frame(
  cell_barcode = rownames(cell_out),
  patient = cell_out$patient,
  tissue = cell_out$label,
  batch = cell_out$batch,
  cell_type_annotated_granularity = cell_out$cell_type_annotated_granularity,
  shared_clone_list = cell_out$shared_clone_list,
  UMAP_1 = umap_out[rownames(cell_out), "UMAP_1"],
  UMAP_2 = umap_out[rownames(cell_out), "UMAP_2"],
  pseudotime = pseudotime(cds),
  monocle_partition = partitions(cds),
  monocle_cluster = clusters(cds),
  stringsAsFactors = FALSE
)

write.table(
  source_data,
  file = file.path(output_dir, "Fig2I_monocle_pseudotime_source_data.tsv"),
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  na = ""
)

pdf(file.path(output_dir, "Fig2I_monocle_by_pseudotime.pdf"), width = 8, height = 6)
print(
  plot_cells(
    cds,
    color_cells_by = "pseudotime",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.9,
    graph_label_size = 1.5,
    label_roots = FALSE
  ) +
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
)
dev.off()

cds$clusters <- cds$cell_type_annotated_granularity
pdf(file.path(output_dir, "Fig2I_monocle_by_clusters.pdf"), width = 8, height = 6)
print(
  plot_cells(
    cds,
    color_cells_by = "clusters",
    label_cell_groups = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.9,
    graph_label_size = 1.5,
    label_roots = FALSE
  ) +
    theme(
      axis.line.x = element_blank(),
      axis.line.y = element_blank(),
      axis.text = element_blank(),
      axis.ticks = element_blank(),
      axis.title = element_blank()
    )
)
dev.off()
