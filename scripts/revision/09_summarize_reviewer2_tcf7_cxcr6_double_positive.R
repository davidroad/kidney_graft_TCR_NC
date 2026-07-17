#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript 09_summarize_reviewer2_tcf7_cxcr6_double_positive.R ",
    "<primary_cd8.RData> <gse224445_strict_cd8.rds> <output_directory>"
  )
}
primary_file <- normalizePath(args[[1]], mustWork = TRUE)
gse_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

load_single_seurat <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  seurat_names <- object_names[
    vapply(object_names, function(name) inherits(env[[name]], "Seurat"), logical(1))
  ]
  if (length(seurat_names) != 1L) stop("Expected one Seurat object in ", path)
  env[[seurat_names[[1]]]]
}

summarize_object <- function(object, sample_id, condition, dataset) {
  DefaultAssay(object) <- "RNA"
  counts <- LayerData(
    object,
    assay = "RNA",
    layer = "counts",
    features = c("TCF7", "CXCR6")
  )[, colnames(object), drop = FALSE]
  cells <- data.frame(
    sample_id = unname(sample_id[colnames(object)]),
    condition = unname(condition[colnames(object)]),
    TCF7_positive = as.numeric(counts["TCF7", ] > 0),
    CXCR6_positive = as.numeric(counts["CXCR6", ] > 0),
    stringsAsFactors = FALSE
  )
  cells$double_positive <- cells$TCF7_positive * cells$CXCR6_positive
  groups <- split(
    seq_len(nrow(cells)),
    interaction(cells$sample_id, cells$condition, drop = TRUE)
  )
  do.call(rbind, lapply(groups, function(index) {
    data.frame(
      dataset = dataset,
      sample_id = cells$sample_id[index[[1]]],
      condition = cells$condition[index[[1]]],
      n_cd8_cells = length(index),
      n_TCF7_positive = sum(cells$TCF7_positive[index]),
      n_CXCR6_positive = sum(cells$CXCR6_positive[index]),
      n_TCF7_CXCR6_double_positive = sum(cells$double_positive[index]),
      TCF7_CXCR6_double_positive_percent = 100 * mean(cells$double_positive[index]),
      stringsAsFactors = FALSE
    )
  }))
}

primary <- load_single_seurat(primary_file)
primary_group <- as.character(primary$group)
primary_group[
  is.na(primary_group) & as.character(primary$patient) == "0701"
] <- "PBMC_0701"
if (any(is.na(primary_group))) stop("Unresolved primary-study sample labels.")
names(primary_group) <- colnames(primary)
primary_condition <- ifelse(
  grepl("^PBMC_", primary_group),
  "PBMC:Rejection",
  "Graft:Rejection"
)
names(primary_condition) <- colnames(primary)
primary_summary <- summarize_object(
  primary,
  primary_group,
  primary_condition,
  "This study"
)

gse <- readRDS(gse_file)
gse_sample <- paste(
  as.character(gse$bead),
  as.character(gse$Sample_Tag),
  sep = "_"
)
names(gse_sample) <- colnames(gse)
gse_condition <- paste0("PBMC:", as.character(gse$group))
names(gse_condition) <- colnames(gse)
gse_summary <- summarize_object(
  gse,
  gse_sample,
  gse_condition,
  "GSE224445"
)

condition_levels <- c(
  "PBMC:TOT",
  "PBMC:STA",
  "PBMC:BPAR",
  "PBMC:Rejection",
  "Graft:Rejection"
)
result <- rbind(gse_summary, primary_summary)
result$condition <- factor(result$condition, levels = condition_levels)
result <- result[order(result$condition, result$sample_id), ]
write.csv(
  result,
  file.path(
    output_dir,
    "reviewer2_tcf7_cxcr6_double_positive_frequency_five_groups.csv"
  ),
  row.names = FALSE,
  quote = TRUE
)

colours <- c(
  "PBMC:TOT" = "#4E79A7",
  "PBMC:STA" = "#59A14F",
  "PBMC:BPAR" = "#F28E2B",
  "PBMC:Rejection" = "#B07AA1",
  "Graft:Rejection" = "#E15759"
)
labels <- c(
  "PBMC:TOT" = "PBMC\nTOT",
  "PBMC:STA" = "PBMC\nSTA",
  "PBMC:BPAR" = "PBMC\nBPAR",
  "PBMC:Rejection" = "PBMC\nRejection",
  "Graft:Rejection" = "Graft\nRejection"
)
set.seed(20260717)
plot <- ggplot(
  result,
  aes(condition, TCF7_CXCR6_double_positive_percent, colour = condition)
) +
  geom_boxplot(
    width = 0.48,
    outlier.shape = NA,
    colour = "#4B4B4B",
    fill = "white",
    linewidth = 0.42
  ) +
  geom_point(
    size = 2.5,
    position = position_jitter(width = 0.07, height = 0),
    alpha = 0.95
  ) +
  stat_summary(
    fun = mean,
    geom = "crossbar",
    width = 0.28,
    colour = "black",
    linewidth = 0.4
  ) +
  scale_colour_manual(values = colours, guide = "none") +
  scale_x_discrete(labels = labels) +
  scale_y_continuous(expand = expansion(mult = c(0.03, 0.12))) +
  labs(
    x = NULL,
    y = "TCF7+CXCR6+ among CD8+ T cells (%)"
  ) +
  theme_classic(base_size = 12, base_family = "DejaVu Sans") +
  theme(
    axis.text.x = element_text(size = 10.5, lineheight = 0.9, hjust = 0.5),
    plot.margin = margin(8, 10, 8, 8)
  )

ggsave(
  file.path(
    output_dir,
    "reviewer2_tcf7_cxcr6_double_positive_frequency_five_groups.pdf"
  ),
  plot,
  width = 7.8,
  height = 5.0,
  device = cairo_pdf
)
ggsave(
  file.path(
    output_dir,
    "reviewer2_tcf7_cxcr6_double_positive_frequency_five_groups.png"
  ),
  plot,
  width = 7.8,
  height = 5.0,
  dpi = 300
)
message("Wrote double-positive source table and figure to: ", output_dir)
