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
    "Usage: Rscript 07_generate_reviewer2_five_group_featureplots.R ",
    "<primary_cd8.RData> <gse224445_strict_cd8.rds> <output_directory>"
  )
}
primary_file <- normalizePath(args[[1]], mustWork = TRUE)
gse_file <- normalizePath(args[[2]], mustWork = TRUE)
output_dirs <- normalizePath(args[[3]], mustWork = FALSE)
invisible(lapply(output_dirs, dir.create, recursive = TRUE, showWarnings = FALSE))

genes <- c("CXCR6", "TCF7", "IFNG", "GZMB", "TOX", "PDCD1")
conditions <- c("PBMC:TOT", "PBMC:STA", "PBMC:BPAR", "PBMC:Rejection", "Graft:Rejection")
condition_titles <- c("PBMC:TOT" = "PBMC\nTOT", "PBMC:STA" = "PBMC\nSTA", "PBMC:BPAR" = "PBMC\nBPAR", "PBMC:Rejection" = "PBMC\nRejection", "Graft:Rejection" = "Graft\nRejection")

load_single_seurat <- function(path) {
  env <- new.env(parent = emptyenv())
  object_names <- load(path, envir = env)
  seurat_names <- object_names[vapply(object_names, function(nm) inherits(env[[nm]], "Seurat"), logical(1))]
  if (length(seurat_names) != 1L) {
    stop("Expected exactly one Seurat object in ", path, "; found: ", paste(seurat_names, collapse = ", "))
  }
  env[[seurat_names[[1]]]]
}

extract_plot_data <- function(object, reduction, group, dataset, genes) {
  DefaultAssay(object) <- "RNA"
  coords <- Embeddings(object, reduction = reduction)[, 1:2, drop = FALSE]
  colnames(coords) <- c("umap_1", "umap_2")
  available <- intersect(genes, rownames(object))
  expression <- as.matrix(LayerData(object, assay = "RNA", layer = "data", features = available))
  expression <- t(expression[, rownames(coords), drop = FALSE])
  out <- data.frame(
    cell = rownames(coords),
    umap_1 = coords[, 1],
    umap_2 = coords[, 2],
    condition = unname(group[rownames(coords)]),
    dataset = dataset,
    stringsAsFactors = FALSE
  )
  for (gene in genes) {
    out[[gene]] <- if (gene %in% colnames(expression)) expression[, gene] else NA_real_
  }
  out
}

primary <- load_single_seurat(primary_file)
primary_group <- as.character(primary$group)
primary_group[is.na(primary_group) & as.character(primary$patient) == "0701"] <- "PBMC_0701"
if (any(is.na(primary_group))) {
  stop("Unresolved primary-study group labels remain after restoring PBMC_0701 from patient metadata.")
}
primary_condition <- ifelse(grepl("^PBMC_", primary_group), "PBMC:Rejection", "Graft:Rejection")
names(primary_condition) <- colnames(primary)
primary_df <- extract_plot_data(primary, "umap", primary_condition, "This study", genes)
rm(primary)
invisible(gc())

gse <- readRDS(gse_file)
gse_condition <- paste0("PBMC:", as.character(gse$group))
names(gse_condition) <- colnames(gse)
gse_df <- extract_plot_data(gse, "cd8_umap_highres1_strict", gse_condition, "GSE224445", genes)
rm(gse)
invisible(gc())

all_df <- rbind(gse_df, primary_df)
all_df$condition <- factor(all_df$condition, levels = conditions)

dataset_limits <- lapply(split(all_df, all_df$dataset), function(dat) {
  x_range <- range(dat$umap_1, finite = TRUE)
  y_range <- range(dat$umap_2, finite = TRUE)
  span <- max(diff(x_range), diff(y_range))
  x_mid <- mean(x_range)
  y_mid <- mean(y_range)
  list(
    x = x_mid + c(-0.53, 0.53) * span,
    y = y_mid + c(-0.53, 0.53) * span
  )
})

relative_expression <- function(dat, gene) {
  values <- dat[[gene]]
  if (all(is.na(values))) return(rep(NA_real_, length(values)))
  cap <- unname(quantile(values[is.finite(values)], 0.99, na.rm = TRUE))
  if (!is.finite(cap) || cap <= 0) cap <- max(values, na.rm = TRUE)
  if (!is.finite(cap) || cap <= 0) return(rep(0, length(values)))
  pmin(pmax(values / cap, 0), 1)
}

for (dataset in unique(all_df$dataset)) {
  rows <- all_df$dataset == dataset
  for (gene in genes) {
    all_df[rows, paste0(gene, "_relative")] <- relative_expression(all_df[rows, , drop = FALSE], gene)
  }
}

panel_theme <- theme_void(base_family = "DejaVu Sans") +
  theme(
    plot.title = element_text(size = 10.5, face = "bold", hjust = 0.5, margin = margin(b = 2)),
    plot.margin = margin(2, 2, 2, 2),
    panel.border = element_rect(colour = "#9A9A9A", fill = NA, linewidth = 0.28)
  )

make_panel <- function(gene, condition, show_title = FALSE) {
  dataset <- if (grepl("^PBMC:(TOT|STA|BPAR)$", condition)) "GSE224445" else "This study"
  dat <- all_df[all_df$condition == condition & all_df$dataset == dataset, , drop = FALSE]
  limits <- dataset_limits[[dataset]]
  rel_col <- paste0(gene, "_relative")
  title <- if (show_title) unname(condition_titles[[condition]]) else NULL

  if (!nrow(dat) || all(is.na(dat[[rel_col]]))) {
    return(
      ggplot() +
        annotate("text", x = mean(limits$x), y = mean(limits$y), label = if (gene == "TOX" && condition == "PBMC:STA") "Not measured in\nGSE224445 targeted panel" else "", size = 3.0, colour = "#555555") +
        coord_fixed(xlim = limits$x, ylim = limits$y, expand = FALSE) +
        labs(title = title) +
        panel_theme
    )
  }

  dat <- dat[order(dat[[rel_col]], na.last = TRUE), , drop = FALSE]
  point_size <- if (dataset == "GSE224445") 0.80 else 0.50
  ggplot(dat, aes(x = umap_1, y = umap_2, colour = .data[[rel_col]])) +
    geom_point(size = point_size, alpha = 0.9, stroke = 0) +
    scale_colour_gradientn(
      colours = c("#D9D9D9", "#F5B39D", "#D7553A", "#7F0000"),
      values = c(0, 0.08, 0.42, 1),
      limits = c(0, 1),
      oob = scales::squish,
      guide = "none"
    ) +
    coord_fixed(xlim = limits$x, ylim = limits$y, expand = FALSE) +
    labs(title = title) +
    panel_theme
}

header_external <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "External PBMC (GSE224445)", fontface = "bold", size = 4.0, family = "DejaVu Sans") +
  theme_void()
header_primary <- ggplot() +
  annotate("text", x = 0.5, y = 0.5, label = "This study", fontface = "bold", size = 4.0, family = "DejaVu Sans") +
  theme_void()
header_spacer <- plot_spacer()

header <- header_spacer + header_external + header_primary +
  plot_layout(widths = c(0.34, 3, 2))

rows <- lapply(seq_along(genes), function(i) {
  gene <- genes[[i]]
  label_plot <- ggplot() +
    annotate("text", x = 0.5, y = 0.5, label = gene, fontface = "bold", size = 4.0, family = "DejaVu Sans") +
    theme_void()
  panels <- lapply(conditions, function(condition) make_panel(gene, condition, show_title = i == 1L))
  wrap_plots(c(list(label_plot), panels), nrow = 1, widths = c(0.34, rep(1, 5)))
})

figure <- wrap_plots(c(list(header), rows), ncol = 1, heights = c(0.18, rep(1, length(rows)))) +
  plot_annotation(
    title = "CD8+ T-cell marker expression by clinical group",
    theme = theme(
      plot.title = element_text(family = "DejaVu Sans", face = "bold", size = 14, hjust = 0.5, margin = margin(b = 6))
    )
  )

basename_pdf <- "reviewer2_requested_marker_featureplots_five_groups.pdf"
basename_png <- "reviewer2_requested_marker_featureplots_five_groups.png"
for (output_dir in output_dirs) {
  ggsave(file.path(output_dir, basename_pdf), figure, width = 15.8, height = 18.5, units = "in", device = cairo_pdf)
  ggsave(file.path(output_dir, basename_png), figure, width = 15.8, height = 18.5, units = "in", dpi = 220)
}

message("Wrote five-group reviewer feature plots to: ", paste(output_dirs, collapse = "; "))
