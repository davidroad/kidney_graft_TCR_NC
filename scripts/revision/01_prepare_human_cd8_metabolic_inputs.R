#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  repo_lib <- Sys.getenv("R_LIBS_USER", unset = "")
  if (nzchar(repo_lib)) .libPaths(c(repo_lib, .libPaths()))
  library(Seurat)
  library(Matrix)
  library(ggplot2)
  library(edgeR)
})

set.seed(1)

input_rdata <- "merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData"
out_dir <- file.path("results", "human_cd8_metabolic_analysis")
input_dir <- file.path(out_dir, "input")
fig_dir <- file.path(out_dir, "figures")
edgeR_dir <- file.path(out_dir, "edgeR")
humess_dir <- file.path(out_dir, "humess")
troppo_dir <- file.path(out_dir, "troppo")
for (d in c(out_dir, input_dir, fig_dir, edgeR_dir, humess_dir, troppo_dir)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

message("Loading Seurat object: ", input_rdata)
env <- new.env(parent = emptyenv())
load(input_rdata, envir = env)
obj_name <- ls(env)[1]
obj <- env[[obj_name]]
stopifnot(inherits(obj, "Seurat"))
DefaultAssay(obj) <- "RNA"

get_assay_layer <- function(object, assay = "RNA", layer = "counts") {
  out <- tryCatch(
    GetAssayData(object, assay = assay, layer = layer),
    error = function(e) GetAssayData(object, assay = assay, slot = layer)
  )
  out
}

write_matrix_tsv <- function(mat, file) {
  df <- data.frame(gene = rownames(mat), as.matrix(mat), check.names = FALSE)
  con <- if (grepl("\\.gz$", file)) gzfile(file, "wt") else file
  on.exit(if (inherits(con, "connection")) close(con), add = TRUE)
  write.table(df, con, sep = "\t", quote = FALSE, row.names = FALSE)
}

sanitize_id <- function(x) {
  x <- as.character(x)
  x <- gsub("[^A-Za-z0-9_]+", "_", x)
  x <- gsub("_+", "_", x)
  x <- gsub("^_|_$", "", x)
  x
}

aggregate_sparse_by_group <- function(mat, group) {
  ok <- !is.na(group) & nzchar(as.character(group))
  group <- sanitize_id(group[ok])
  lv <- sort(unique(group))
  j <- match(group, lv)
  design <- sparseMatrix(
    i = seq_along(j),
    j = j,
    x = 1,
    dims = c(length(j), length(lv)),
    dimnames = list(colnames(mat)[ok], lv)
  )
  mat[, ok, drop = FALSE] %*% design
}

counts <- get_assay_layer(obj, layer = "counts")
data_layer <- get_assay_layer(obj, layer = "data")

meta <- obj@meta.data
meta$barcode <- rownames(meta)
if (!"patient" %in% colnames(meta)) {
  stop("metadata column `patient` is required")
}
if (!"label" %in% colnames(meta)) {
  stop("metadata column `label` is required")
}
meta$condition <- as.character(meta$label)
meta$patient <- as.character(meta$patient)
meta$group_clean <- as.character(meta$group)
missing_group <- is.na(meta$group_clean) | !nzchar(meta$group_clean)
meta$group_clean[missing_group] <- paste0(meta$condition[missing_group], "_", meta$patient[missing_group])
meta$sample_id <- paste(meta$patient, meta$condition, sep = "_")
meta$sample_id <- sanitize_id(meta$sample_id)
meta$celltype_granular <- as.character(meta$cell_type_annotated_granularity)
meta$celltype_granular[is.na(meta$celltype_granular)] <- "Unknown_CD8"
meta$celltype_granular_clean <- sanitize_id(meta$celltype_granular)
meta$sample_celltype_id <- paste(meta$sample_id, meta$celltype_granular_clean, sep = "__")

umap <- Embeddings(obj, "umap")
meta$UMAP_1 <- umap[rownames(meta), 1]
meta$UMAP_2 <- umap[rownames(meta), 2]

glycolysis_genes <- unique(c(
  "HK1", "HK2", "GPI", "PFKL", "PFKP", "ALDOA", "GAPDH", "PGK1",
  "PGAM1", "ENO1", "ENO2", "PKM", "LDHA", "LDHB", "SLC2A1",
  "G6PC3", "AKR1A1"
))
stem_markers <- c("TCF7", "CXCR6")
all_score_genes <- unique(c(glycolysis_genes, stem_markers))
present_genes <- intersect(all_score_genes, rownames(data_layer))
missing_genes <- setdiff(all_score_genes, rownames(data_layer))
if (length(missing_genes) > 0) {
  message("Missing score genes: ", paste(missing_genes, collapse = ", "))
}

score_from_genes <- function(mat, genes) {
  genes <- intersect(genes, rownames(mat))
  if (length(genes) == 0) {
    return(rep(NA_real_, ncol(mat)))
  }
  as.numeric(Matrix::colMeans(mat[genes, , drop = FALSE]))
}

meta$glycolysis_score <- score_from_genes(data_layer, glycolysis_genes)
for (g in stem_markers) {
  meta[[paste0(g, "_expr")]] <- if (g %in% rownames(data_layer)) {
    as.numeric(data_layer[g, ])
  } else {
    NA_real_
  }
}
meta$TCF7_positive <- !is.na(meta$TCF7_expr) & meta$TCF7_expr > 0
meta$CXCR6_positive <- !is.na(meta$CXCR6_expr) & meta$CXCR6_expr > 0
meta$TCF7_CXCR6_double_positive <- meta$TCF7_positive & meta$CXCR6_positive
meta$stem_like_annotation <- meta$celltype_granular == "ZNF683+TCF7hi CD8+"

metadata_file <- file.path(input_dir, "cd8_cell_metadata_with_umap.tsv.gz")
scores_file <- file.path(input_dir, "cd8_metabolic_scores.tsv.gz")
write.table(meta, gzfile(metadata_file, "wt"), sep = "\t", quote = FALSE, row.names = FALSE)
score_cols <- c(
  "barcode", "patient", "condition", "sample_id", "group_clean",
  "celltype_granular", "UMAP_1", "UMAP_2", "glycolysis_score",
  "TCF7_expr", "CXCR6_expr", "TCF7_positive", "CXCR6_positive",
  "TCF7_CXCR6_double_positive", "stem_like_annotation"
)
write.table(meta[, score_cols], gzfile(scores_file, "wt"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Building pseudobulk matrices")
pb_sample_counts <- aggregate_sparse_by_group(counts, meta$sample_id)
sample_meta <- aggregate(
  barcode ~ sample_id + patient + condition,
  data = meta,
  FUN = length
)
colnames(sample_meta)[colnames(sample_meta) == "barcode"] <- "n_cells"
sample_meta <- sample_meta[match(colnames(pb_sample_counts), sample_meta$sample_id), ]

pb_sample_celltype_counts <- aggregate_sparse_by_group(counts, meta$sample_celltype_id)
celltype_meta <- aggregate(
  barcode ~ sample_celltype_id + sample_id + patient + condition + celltype_granular,
  data = meta,
  FUN = length
)
colnames(celltype_meta)[colnames(celltype_meta) == "barcode"] <- "n_cells"
celltype_meta <- celltype_meta[match(colnames(pb_sample_celltype_counts), celltype_meta$sample_celltype_id), ]

write.table(sample_meta, file.path(input_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(celltype_meta, file.path(input_dir, "sample_celltype_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write_matrix_tsv(pb_sample_counts, file.path(input_dir, "pseudobulk_counts_by_sample.tsv.gz"))
write_matrix_tsv(pb_sample_celltype_counts, file.path(input_dir, "pseudobulk_counts_by_sample_celltype.tsv.gz"))

dge <- DGEList(counts = pb_sample_counts)
dge <- calcNormFactors(dge)
pb_cpm <- cpm(dge, log = FALSE)
pb_logcpm <- cpm(dge, log = TRUE, prior.count = 1)
write_matrix_tsv(pb_cpm, file.path(input_dir, "pseudobulk_cpm_by_sample.tsv.gz"))
write_matrix_tsv(pb_logcpm, file.path(input_dir, "pseudobulk_logcpm_by_sample.tsv.gz"))

message("Running paired edgeR: Graft vs PBMC")
sample_meta$condition <- factor(sample_meta$condition, levels = c("PBMC", "Graft"))
sample_meta$patient <- factor(sample_meta$patient)
keep <- filterByExpr(dge, group = sample_meta$condition)
dge2 <- dge[keep, , keep.lib.sizes = FALSE]
dge2 <- calcNormFactors(dge2)
design <- model.matrix(~ patient + condition, data = sample_meta)
dge2 <- estimateDisp(dge2, design)
fit <- glmQLFit(dge2, design)
qlf <- glmQLFTest(fit, coef = "conditionGraft")
de <- topTags(qlf, n = Inf)$table
de$gene <- rownames(de)
de$is_glycolysis_gene <- de$gene %in% glycolysis_genes
de <- de[, c("gene", setdiff(colnames(de), "gene"))]
write.table(de, file.path(edgeR_dir, "paired_graft_vs_pbmc_edgeR.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

message("Writing HUMESS input files")
humess_counts <- file.path(humess_dir, "umi_counts_by_sample.tsv.gz")
write_matrix_tsv(pb_sample_counts, humess_counts)
write.table(sample_meta, file.path(humess_dir, "sample_sheet.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
humess_samplesheet_native <- file.path(humess_dir, "samplesheet.native.tsv")
write.table(
  sample_meta[, c("sample_id", "condition")],
  humess_samplesheet_native,
  sep = "\t",
  quote = FALSE,
  row.names = FALSE,
  col.names = FALSE
)
humess_abundance_native <- file.path(humess_dir, "abundance.native.tsv")
humess_df <- data.frame(gene = rownames(pb_sample_counts), as.matrix(pb_sample_counts), check.names = FALSE)
colnames(humess_df)[1] <- ""
write.table(humess_df, humess_abundance_native, sep = "\t", quote = FALSE, row.names = FALSE)
expressed_dir <- file.path(humess_dir, "expressed_genes_by_sample")
dir.create(expressed_dir, showWarnings = FALSE)
global_keep <- filterByExpr(DGEList(counts = pb_sample_counts), group = sample_meta$condition)
global_keep_genes <- rownames(pb_sample_counts)[global_keep]
writeLines(global_keep_genes, file.path(humess_dir, "expressed_genes_edgeR_global.txt"))
for (sid in colnames(pb_cpm)) {
  expressed <- rownames(pb_cpm)[pb_cpm[, sid] >= 1 & rownames(pb_cpm) %in% global_keep_genes]
  writeLines(expressed, file.path(expressed_dir, paste0(sid, ".expressed_genes.txt")))
}
humess_de <- de[, c("gene", "logFC", "FDR")]
write.table(humess_de, file.path(humess_dir, "graft_vs_pbmc_deg_for_humess.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
writeLines(
  paste(
    "PBMC",
    "Graft",
    normalizePath(file.path(humess_dir, "graft_vs_pbmc_deg_for_humess.tsv"), mustWork = FALSE),
    sep = "\t"
  ),
  file.path(humess_dir, "comparisons.native.tsv")
)
humess_config <- c(
  "# Template config for external HUMESS workflow.",
  "# Adjust keys after cloning HUMESS because exact config names can change by release.",
  paste0("counts_table: ", normalizePath(humess_counts, mustWork = FALSE)),
  paste0("sample_sheet: ", normalizePath(file.path(humess_dir, "sample_sheet.tsv"), mustWork = FALSE)),
  paste0("expressed_gene_dir: ", normalizePath(expressed_dir, mustWork = FALSE)),
  "species: human",
  "model_template: Recon3D",
  "filtering: edgeR_filterByExpr_plus_sample_CPM_ge_1",
  "comparison:",
  "  reference: PBMC",
  "  case: Graft",
  "paired_by: patient"
)
writeLines(humess_config, file.path(humess_dir, "config_template.yaml"))

message("Writing Troppo/human_ts_models input files")
write_matrix_tsv(pb_sample_counts, file.path(troppo_dir, "counts_by_sample.tsv.gz"))
write_matrix_tsv(pb_logcpm, file.path(troppo_dir, "logcpm_by_sample.tsv.gz"))
write.table(sample_meta, file.path(troppo_dir, "sample_metadata.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
tas <- pb_logcpm
for (i in seq_len(nrow(tas))) {
  x <- tas[i, ]
  center <- median(x, na.rm = TRUE)
  spread <- mad(x, constant = 1, na.rm = TRUE)
  if (is.na(spread) || spread == 0) spread <- sd(x, na.rm = TRUE)
  if (is.na(spread) || spread == 0) spread <- 1
  tas[i, ] <- (x - center) / spread
}
write_matrix_tsv(tas, file.path(troppo_dir, "transcript_activity_scores_z.tsv.gz"))
human_ts_expr <- data.frame(sample_id = colnames(tas), t(as.matrix(tas)), check.names = FALSE)
write.csv(human_ts_expr, file.path(troppo_dir, "human_ts_models_expression_scores.csv"), quote = FALSE, row.names = FALSE)
human_ts_config <- c(
  paste0("ROOT_FOLDER = ", normalizePath(troppo_dir, mustWork = FALSE)),
  paste0("MODEL_PATH = ", normalizePath(file.path("external", "human_ts_models", "shared", "models", "Human-GEM_latest_consistent.xml"), mustWork = FALSE)),
  paste0("DATA_PATH = ", normalizePath(file.path(troppo_dir, "human_ts_models_expression_scores.csv"), mustWork = FALSE)),
  paste0("SAMPLE_INFO = ", normalizePath(file.path(troppo_dir, "sample_metadata.tsv"), mustWork = FALSE)),
  paste0("CS_MODEL_DF_FOLDER = ", normalizePath(file.path(troppo_dir, "reconstructions"), mustWork = FALSE)),
  "CS_MODEL_NAMES = cd8_paired_human_ts_models",
  "INDEX_COLUMNS = 0",
  "NTHREADS = 4",
  "PROTECTED_REACTION_LIST = biomass_human,HMR_10023,HMR_10024",
  paste0(
    "SOURCES_TO_ADD = ",
    paste(
      normalizePath(file.path("external", "troppo"), mustWork = FALSE),
      normalizePath(file.path("external", "troppo", "src"), mustWork = FALSE),
      normalizePath(file.path("external", "human_ts_models"), mustWork = FALSE),
      normalizePath(file.path("external", "human_ts_models", "shared", "src"), mustWork = FALSE),
      sep = ":"
    )
  )
)
writeLines(human_ts_config, file.path(troppo_dir, "human_ts_models_config.txt"))
troppo_notes <- c(
  "Troppo/human_ts_models input notes",
  "",
  "- `counts_by_sample.tsv.gz`: raw paired pseudobulk counts.",
  "- `logcpm_by_sample.tsv.gz`: edgeR normalized logCPM expression.",
  "- `transcript_activity_scores_z.tsv.gz`: robust z-scored logCPM per gene, for adaptation to TAS/RAS preprocessing.",
  "- `human_ts_models_expression_scores.csv`: samples x genes input for `shared/scripts/model_reconstruction.py`.",
  "- `human_ts_models_config.txt`: config file for `model_reconstruction.py`.",
  "- `sample_metadata.tsv`: patient and condition metadata for paired PBMC/Graft samples.",
  "",
  "Use these files after installing the BioSystemsUM `troppo` and `human_ts_models` repositories."
)
writeLines(troppo_notes, file.path(troppo_dir, "README.txt"))

message("Generating figures")
fig_meta <- meta
fig_meta$condition <- factor(fig_meta$condition, levels = c("PBMC", "Graft"))
point_size <- 0.12
theme_umap <- theme_classic(base_size = 10) +
  theme(axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank())

p1 <- ggplot(fig_meta, aes(UMAP_1, UMAP_2, color = glycolysis_score)) +
  geom_point(size = point_size, alpha = 0.75) +
  scale_color_viridis_c(option = "magma") +
  theme_umap +
  labs(title = "CD8 glycolysis score", color = "Score", x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(fig_dir, "umap_glycolysis_score.pdf"), p1, width = 5.2, height = 4.4)

p2 <- ggplot(fig_meta, aes(UMAP_1, UMAP_2, color = condition)) +
  geom_point(size = point_size, alpha = 0.65) +
  theme_umap +
  labs(title = "Paired CD8 samples", color = "Condition", x = "UMAP 1", y = "UMAP 2")
ggsave(file.path(fig_dir, "umap_condition.pdf"), p2, width = 5.2, height = 4.4)

p3 <- ggplot(fig_meta, aes(UMAP_1, UMAP_2, color = celltype_granular)) +
  geom_point(size = point_size, alpha = 0.65) +
  theme_umap +
  labs(title = "CD8 subtype annotation", color = "Subtype", x = "UMAP 1", y = "UMAP 2") +
  guides(color = guide_legend(override.aes = list(size = 2, alpha = 1)))
ggsave(file.path(fig_dir, "umap_cd8_subtypes.pdf"), p3, width = 7.2, height = 4.8)

for (g in c("TCF7_expr", "CXCR6_expr")) {
  p <- ggplot(fig_meta, aes(UMAP_1, UMAP_2, color = .data[[g]])) +
    geom_point(size = point_size, alpha = 0.75) +
    scale_color_viridis_c(option = "plasma") +
    theme_umap +
    labs(title = sub("_expr", "", g), color = "Expr", x = "UMAP 1", y = "UMAP 2")
  ggsave(file.path(fig_dir, paste0("umap_", g, ".pdf")), p, width = 5.2, height = 4.4)
}

sample_scores <- aggregate(
  glycolysis_score ~ patient + condition + sample_id,
  data = fig_meta,
  FUN = mean
)
p4 <- ggplot(sample_scores, aes(condition, glycolysis_score, group = patient, color = patient)) +
  geom_line(linewidth = 0.5, alpha = 0.8) +
  geom_point(size = 2) +
  theme_classic(base_size = 11) +
  labs(title = "Paired patient-level glycolysis score", x = NULL, y = "Mean single-cell score")
ggsave(file.path(fig_dir, "paired_glycolysis_score_by_patient.pdf"), p4, width = 4.8, height = 3.8)

p5 <- ggplot(fig_meta, aes(stem_like_annotation, glycolysis_score, fill = stem_like_annotation)) +
  geom_boxplot(outlier.size = 0.1) +
  facet_wrap(~ condition, scales = "free_y") +
  theme_classic(base_size = 11) +
  labs(title = "Glycolysis in stem-like annotated CD8 cells", x = "ZNF683+TCF7hi annotation", y = "Glycolysis score") +
  theme(legend.position = "none")
ggsave(file.path(fig_dir, "glycolysis_stem_like_boxplot.pdf"), p5, width = 6.0, height = 3.8)

volcano <- de
volcano$neg_log10_fdr <- -log10(pmax(volcano$FDR, .Machine$double.xmin))
volcano$highlight <- ifelse(volcano$is_glycolysis_gene, "Glycolysis", "Other")
p6 <- ggplot(volcano, aes(logFC, neg_log10_fdr, color = highlight)) +
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_manual(values = c(Glycolysis = "#D55E00", Other = "grey70")) +
  theme_classic(base_size = 11) +
  labs(title = "Paired Graft vs PBMC CD8 transcriptome", x = "logFC Graft/PBMC", y = "-log10 FDR", color = NULL)
ggsave(file.path(fig_dir, "edgeR_paired_graft_vs_pbmc_volcano.pdf"), p6, width = 5.2, height = 4.2)

summary_lines <- c(
  paste0("object_name\t", obj_name),
  paste0("n_cells\t", ncol(obj)),
  paste0("n_genes\t", nrow(obj)),
  paste0("samples\t", paste(colnames(pb_sample_counts), collapse = ",")),
  paste0("glycolysis_genes_present\t", paste(intersect(glycolysis_genes, rownames(obj)), collapse = ",")),
  paste0("glycolysis_genes_missing\t", paste(setdiff(glycolysis_genes, rownames(obj)), collapse = ",")),
  paste0("edgeR_genes_tested\t", nrow(de)),
  paste0("output_dir\t", normalizePath(out_dir, mustWork = FALSE))
)
writeLines(summary_lines, file.path(out_dir, "run_summary.tsv"))
message("Done. Outputs written to: ", out_dir)
