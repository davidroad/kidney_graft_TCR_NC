#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3L) {
  stop(
    "Usage: Rscript 05_generate_suppfig6a_cxcr6_umaps.R ",
    "<primary_cd8.RData> <gse224445_strict_cd8.rds> <output_directory>"
  )
}

primary_file <- normalizePath(args[[1]], mustWork = TRUE)
gse_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dir <- normalizePath(args[[3]], mustWork = FALSE)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

conditions <- c(
  "PBMC:TOT",
  "PBMC:STA",
  "PBMC:BPAR",
  "PBMC:Rejection",
  "Graft:Rejection"
)
condition_titles <- c(
  "PBMC:TOT" = "PBMC;TOT",
  "PBMC:STA" = "PBMC;STA",
  "PBMC:BPAR" = "PBMC;BPAR",
  "PBMC:Rejection" = "PBMC",
  "Graft:Rejection" = "Graft"
)

load_single_seurat <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  seurat_names <- object_names[
    vapply(object_names, function(name) inherits(env[[name]], "Seurat"), logical(1))
  ]
  if (length(seurat_names) != 1L) {
    stop("Expected one Seurat object in ", path, "; found: ", paste(seurat_names, collapse = ", "))
  }
  env[[seurat_names[[1]]]]
}

extract_plot_data <- function(object, reduction, group, dataset) {
  DefaultAssay(object) <- "RNA"
  coordinates <- Embeddings(object, reduction = reduction)[, 1:2, drop = FALSE]
  colnames(coordinates) <- c("UMAP_1", "UMAP_2")
  if (!"CXCR6" %in% rownames(object)) stop("CXCR6 is absent from ", dataset)
  expression <- as.numeric(LayerData(
    object,
    assay = "RNA",
    layer = "data",
    features = "CXCR6"
  )[1, rownames(coordinates)])
  data.frame(
    cell = rownames(coordinates),
    UMAP_1 = coordinates[, 1],
    UMAP_2 = coordinates[, 2],
    condition = unname(group[rownames(coordinates)]),
    dataset = dataset,
    CXCR6 = expression,
    stringsAsFactors = FALSE
  )
}

primary <- load_single_seurat(primary_file)
primary_group <- if ("group" %in% colnames(primary@meta.data)) as.character(primary$group) else rep(NA_character_, ncol(primary))
primary_tissue <- rep(NA_character_, ncol(primary))
for (field in c("tissue", "label", "condition")) {
  if (field %in% colnames(primary@meta.data)) {
    values <- as.character(primary[[field, drop = TRUE]])
    primary_tissue[is.na(primary_tissue) & grepl("PBMC", values, ignore.case = TRUE)] <- "PBMC"
    primary_tissue[is.na(primary_tissue) & grepl("Graft|Kidney", values, ignore.case = TRUE)] <- "Graft"
  }
}
primary_tissue[grepl("^PBMC_", primary_group)] <- "PBMC"
primary_tissue[grepl("^Graft_", primary_group)] <- "Graft"
if (anyNA(primary_tissue)) stop("Unable to resolve PBMC or graft metadata for every primary-study CD8 cell.")
primary_condition <- ifelse(primary_tissue == "PBMC", "PBMC:Rejection", "Graft:Rejection")
names(primary_condition) <- colnames(primary)
primary_data <- extract_plot_data(primary, "umap", primary_condition, "GSE319007")
rm(primary)
invisible(gc())

gse <- readRDS(gse_file)
gse_condition <- paste0("PBMC:", as.character(gse$group))
names(gse_condition) <- colnames(gse)
gse_data <- extract_plot_data(gse, "cd8_umap_highres1_strict", gse_condition, "GSE224445")
rm(gse)
invisible(gc())

plot_data <- rbind(gse_data, primary_data)
plot_data$condition <- factor(plot_data$condition, levels = conditions)

scale_within_dataset <- function(values) {
  cap <- unname(quantile(values[is.finite(values)], 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(values, na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) return(rep(0, length(values)))
  pmin(pmax(values / cap, 0), 1)
}

plot_data$CXCR6_relative <- NA_real_
for (dataset_name in unique(plot_data$dataset)) {
  selected <- plot_data$dataset == dataset_name
  plot_data$CXCR6_relative[selected] <- scale_within_dataset(plot_data$CXCR6[selected])
}

dataset_limits <- lapply(split(plot_data, plot_data$dataset), function(data) {
  x_range <- range(data$UMAP_1, finite = TRUE)
  y_range <- range(data$UMAP_2, finite = TRUE)
  span <- max(diff(x_range), diff(y_range))
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  list(
    x = x_mid + c(-0.53, 0.53) * span,
    y = y_mid + c(-0.53, 0.53) * span
  )
})

panel_theme <- theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5),
    plot.margin = margin(2, 2, 2, 2),
    panel.border = element_rect(colour = "#9A9A9A", fill = NA, size = 0.28)
  )

make_panel <- function(condition) {
  dataset_name <- if (grepl("^PBMC:(TOT|STA|BPAR)$", condition)) {
    "GSE224445"
  } else {
    "GSE319007"
  }
  data <- plot_data[
    plot_data$condition == condition & plot_data$dataset == dataset_name,
    ,
    drop = FALSE
  ]
  data <- data[order(data$CXCR6_relative), ]
  limits <- dataset_limits[[dataset_name]]
  ggplot(data, aes(UMAP_1, UMAP_2, colour = CXCR6_relative)) +
    geom_point(size = 0.6, alpha = 0.9, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#D9D9D9", "#F5B39D", "#D7553A", "#7F0000"),
      values = c(0, 0.08, 0.42, 1),
      limits = c(0, 1),
      breaks = c(0, 0.5, 1),
      oob = scales::squish,
      guide = if (condition == "Graft:Rejection") {
        guide_colorbar(
          title = "Relative log-normalized CXCR6\n(99th percentile = 1)",
          title.position = "top"
        )
      } else {
        "none"
      }
    ) +
    coord_fixed(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = unname(condition_titles[[condition]])) +
    panel_theme
}

panels <- lapply(conditions, make_panel)
figure <- wrap_plots(panels, nrow = 1, guides = "collect") +
  plot_annotation(
    title = "CXCR6 expression in CD8+ T cells",
    subtitle = "GSE224445 scale: TOT, STA and BPAR; GSE319007 scale: PBMC and graft"
  ) &
  theme(legend.position = "right")

ggsave(
  file.path(output_dir, "suppfig6a_cxcr6_umaps.pdf"),
  figure,
  width = 12,
  height = 3.2,
  units = "in",
  device = cairo_pdf
)

metrics_connection <- gzfile(
  file.path(output_dir, "suppfig6a_cxcr6_umap_plotting_data.csv.gz"),
  open = "wt"
)
write.csv(
  plot_data[c("UMAP_1", "UMAP_2", "condition", "dataset", "CXCR6_relative")],
  metrics_connection,
  row.names = FALSE
)
close(metrics_connection)

message("Supplementary Fig. 6A CXCR6 UMAP output written to ", normalizePath(output_dir))
