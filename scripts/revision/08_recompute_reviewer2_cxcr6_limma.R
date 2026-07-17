#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript 08_recompute_reviewer2_cxcr6_limma.R ",
    "<primary_cd8.RData> <gse224445_strict_cd8.rds> <output.csv>"
  )
}

primary_file <- normalizePath(args[[1]], mustWork = TRUE)
gse_file <- normalizePath(args[[2]], mustWork = TRUE)
output_file <- args[[3]]
dir.create(dirname(output_file), recursive = TRUE, showWarnings = FALSE)

load_single_seurat <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  seurat_names <- object_names[
    vapply(object_names, function(name) inherits(env[[name]], "Seurat"), logical(1))
  ]
  if (length(seurat_names) != 1L) {
    stop("Expected one Seurat object in ", path)
  }
  env[[seurat_names[[1]]]]
}

sample_matrix <- function(object, sample_id, genes) {
  DefaultAssay(object) <- "RNA"
  counts <- LayerData(
    object,
    assay = "RNA",
    layer = "counts",
    features = genes
  )[genes, colnames(object), drop = FALSE]
  sample_id <- sample_id[colnames(object)]
  ids <- unique(sample_id)
  y <- matrix(
    NA_real_,
    nrow = length(genes),
    ncol = length(ids),
    dimnames = list(genes, ids)
  )
  n_cells <- setNames(integer(length(ids)), ids)

  for (id in ids) {
    cells <- which(sample_id == id)
    n_cells[[id]] <- length(cells)
    positive <- Matrix::rowSums(counts[, cells, drop = FALSE] > 0)
    y[, id] <- qlogis((positive + 0.5) / (length(cells) + 1))
  }
  list(y = y, n_cells = n_cells)
}

extract_cxcr6 <- function(
  fit,
  coefficient,
  model,
  comparison,
  n_samples,
  n_pairs = NA_integer_,
  caution = ""
) {
  result <- topTable(
    fit,
    coef = coefficient,
    number = Inf,
    sort.by = "none"
  )
  result$gene <- rownames(result)
  result <- result[result$gene == "CXCR6", , drop = FALSE]
  if (nrow(result) != 1L) stop("CXCR6 was not recovered from the fitted model.")

  data.frame(
    comparison = comparison,
    model = model,
    response = paste(
      "sample-level gene-positive frequency with",
      "0.5-count empirical-logit"
    ),
    coefficient_log_odds = result$logFC,
    p_value = result$P.Value,
    FDR_across_common_genes = result$adj.P.Val,
    n_feature_rows = nrow(fit$coefficients),
    n_samples = n_samples,
    n_pairs = n_pairs,
    eBayes = "trend=TRUE; robust=FALSE",
    caution = caution,
    stringsAsFactors = FALSE
  )
}

gse <- readRDS(gse_file)
primary <- load_single_seurat(primary_file)
common_genes <- sort(intersect(rownames(gse), rownames(primary)))
if (!"CXCR6" %in% common_genes || length(common_genes) < 50L) {
  stop("Insufficient shared gene rows for the joint empirical-Bayes model.")
}

gse_sample <- paste(
  as.character(gse$bead),
  as.character(gse$Sample_Tag),
  sep = "_"
)
names(gse_sample) <- colnames(gse)
gse_matrix <- sample_matrix(gse, gse_sample, common_genes)
gse_meta <- unique(data.frame(
  sample = gse_sample,
  condition = as.character(gse$group),
  stringsAsFactors = FALSE
))
rownames(gse_meta) <- gse_meta$sample
gse_meta <- gse_meta[colnames(gse_matrix$y), , drop = FALSE]

primary_group <- as.character(primary$group)
primary_group[
  is.na(primary_group) & as.character(primary$patient) == "0701"
] <- "PBMC_0701"
if (any(is.na(primary_group))) stop("Unresolved primary-study sample labels.")
names(primary_group) <- colnames(primary)
primary_matrix <- sample_matrix(primary, primary_group, common_genes)
primary_meta <- unique(data.frame(
  sample = primary_group,
  patient = as.character(primary$patient),
  condition = ifelse(grepl("^PBMC_", primary_group), "PBMC", "Graft"),
  stringsAsFactors = FALSE
))
rownames(primary_meta) <- primary_meta$sample
primary_meta <- primary_meta[colnames(primary_matrix$y), , drop = FALSE]

primary_pbmc <- rownames(primary_meta)[primary_meta$condition == "PBMC"]
cross_y <- cbind(
  gse_matrix$y,
  primary_matrix$y[, primary_pbmc, drop = FALSE]
)
cross_meta <- rbind(
  data.frame(
    sample = colnames(gse_matrix$y),
    condition = gse_meta[colnames(gse_matrix$y), "condition"],
    dataset = "GSE224445"
  ),
  data.frame(
    sample = primary_pbmc,
    condition = "OurPBMC",
    dataset = "ThisStudy"
  )
)
rownames(cross_meta) <- cross_meta$sample
cross_meta <- cross_meta[colnames(cross_y), , drop = FALSE]
cross_meta$condition <- factor(
  cross_meta$condition,
  levels = c("TOT", "STA", "BPAR", "OurPBMC")
)

cross_design <- model.matrix(~ 0 + condition, cross_meta)
colnames(cross_design) <- levels(cross_meta$condition)
cross_contrasts <- makeContrasts(
  OurPBMC_vs_TOT = OurPBMC - TOT,
  OurPBMC_vs_STA = OurPBMC - STA,
  OurPBMC_vs_BPAR = OurPBMC - BPAR,
  levels = cross_design
)
cross_fit <- eBayes(
  contrasts.fit(lmFit(cross_y, cross_design), cross_contrasts),
  trend = TRUE,
  robust = FALSE
)
cross_caution <- paste(
  "Exploratory cross-cohort contrast: dataset/platform is confounded",
  "with clinical group and cannot be adjusted separately."
)
results <- rbind(
  extract_cxcr6(
    cross_fit,
    "OurPBMC_vs_TOT",
    "unpaired cross-cohort limma",
    "PBMC:Rejection vs PBMC:TOT",
    ncol(cross_y),
    caution = cross_caution
  ),
  extract_cxcr6(
    cross_fit,
    "OurPBMC_vs_STA",
    "unpaired cross-cohort limma",
    "PBMC:Rejection vs PBMC:STA",
    ncol(cross_y),
    caution = cross_caution
  ),
  extract_cxcr6(
    cross_fit,
    "OurPBMC_vs_BPAR",
    "unpaired cross-cohort limma",
    "PBMC:Rejection vs PBMC:BPAR",
    ncol(cross_y),
    caution = cross_caution
  )
)

paired_meta <- primary_meta[colnames(primary_matrix$y), , drop = FALSE]
paired_meta$patient <- factor(paired_meta$patient)
paired_meta$condition <- factor(
  paired_meta$condition,
  levels = c("PBMC", "Graft")
)
paired_design <- model.matrix(~ patient + condition, paired_meta)
paired_fit <- eBayes(
  lmFit(primary_matrix$y, paired_design),
  trend = TRUE,
  robust = FALSE
)
results <- rbind(
  results,
  extract_cxcr6(
    paired_fit,
    "conditionGraft",
    "paired within-cohort limma ~ patient + condition",
    "Graft:Rejection vs PBMC:Rejection",
    ncol(primary_matrix$y),
    n_pairs = 4L,
    caution = "Within-cohort paired contrast; patient is the blocking term."
  )
)

results$BH_across_4_requested_CXCR6_contrasts <- p.adjust(
  results$p_value,
  method = "BH"
)
results$p_value_4dp <- formatC(results$p_value, format = "f", digits = 4)
results$FDR_across_common_genes_4dp <- formatC(
  results$FDR_across_common_genes,
  format = "f",
  digits = 4
)
write.csv(results, output_file, row.names = FALSE, quote = TRUE)
message("Wrote: ", output_file)
