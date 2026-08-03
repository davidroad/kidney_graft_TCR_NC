#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2L) {
  stop("Usage: Rscript 05_generate_suppfig3_4_from_tcr_object.R <T-cell.RData> <output-dir>")
}
input_file <- normalizePath(args[[1]], mustWork = TRUE)
out_dir <- args[[2]]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
load(input_file)
if (!exists("combined") || !is.list(combined)) stop("Input must contain the scRepertoire 'combined' list.")
if (!exists("merged_all_T_batchcorrected_filter") || !inherits(merged_all_T_batchcorrected_filter, "Seurat")) {
  stop("Input must contain 'merged_all_T_batchcorrected_filter'.")
}

annotation_order <- c(
  "Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+", "CXCR6+ effector CD8+",
  "CCR2+ CD4+", "Naive/Memory-like CD8+", "Naive CD4+", "Treg",
  "CX3CR1+ effector CD8+", "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+"
)
cluster_labels <- paste0(0:10, ". ", annotation_order)
meta <- merged_all_T_batchcorrected_filter@meta.data
meta$cluster <- factor(
  match(as.character(meta$cell_type_annotated_granularity), annotation_order) - 1L,
  levels = 0:10, labels = cluster_labels
)
clone_levels <- c(
  "Single (0 < X <= 1)", "Small (1 < X <= 5)", "Medium (5 < X <= 20)",
  "Large (20 < X <= 100)", "Hyperexpanded (100 < X <= 500)"
)
clone_labels <- c("Single (<=1)", "Small (1-5)", "Medium (5-20)", "Large (20-100)", "Hyperexpanded (100-500)")
occ <- as.data.frame(table(meta$cluster, factor(as.character(meta$cloneSize), levels = clone_levels)))
colnames(occ) <- c("cluster", "clone_size", "cell_number")
occ$clone_size <- factor(occ$clone_size, levels = clone_levels, labels = clone_labels)

p3a <- ggplot(occ, aes(cluster, cell_number, fill = clone_size)) +
  geom_col(width = 0.92, colour = "black", linewidth = 0.12) +
  geom_text(aes(label = ifelse(cell_number >= 10, cell_number, "")),
            position = position_stack(vjust = 0.5), size = 2.2) +
  scale_fill_manual(values = c(
    "Single (<=1)" = "#111111", "Small (1-5)" = "#542057", "Medium (5-20)" = "#D41473",
    "Large (20-100)" = "#F08A12", "Hyperexpanded (100-500)" = "#FFF37A"
  ), drop = FALSE) +
  labs(x = NULL, y = "Cell number", fill = "Clone size") +
  theme_classic(base_size = 9) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, vjust = 1), legend.position = "right")
ggsave(file.path(out_dir, "suppfig3a_clonal_occupancy_by_tcell_cluster.pdf"), p3a,
       width = 12, height = 7.5, units = "in", device = cairo_pdf)
write.csv(occ, file.path(out_dir, "suppfig3a_clonal_occupancy_counts.csv"), row.names = FALSE)

sample_order <- c(
  "PBMC_0825_1", "PBMC_0712_2", "PBMC_0701_3", "PBMC_0708_4",
  "graft_0825_5", "graft_0712_6", "graft_0701_7", "graft_0708_8"
)
sample_labels <- c(
  "PBMC 08252021", "PBMC 07122021", "PBMC 07012022", "PBMC 07082022",
  "Graft 08252021", "Graft 07122021", "Graft 07012022", "Graft 07082022"
)
clone_sets <- lapply(combined[sample_order], function(z) unique(z$CTaa[!is.na(z$CTaa) & nzchar(z$CTaa)]))
overlap <- expand.grid(x = sample_order, y = sample_order, stringsAsFactors = FALSE)
overlap$x_index <- match(overlap$x, sample_order)
overlap$y_index <- match(overlap$y, sample_order)
overlap$raw <- mapply(function(x, y) length(intersect(clone_sets[[x]], clone_sets[[y]])), overlap$x, overlap$y)
overlap$jaccard <- mapply(function(x, y) {
  den <- length(union(clone_sets[[x]], clone_sets[[y]]))
  if (den == 0L) 0 else length(intersect(clone_sets[[x]], clone_sets[[y]])) / den
}, overlap$x, overlap$y)
overlap <- overlap[overlap$y_index > overlap$x_index, , drop = FALSE]
overlap$x_label <- factor(overlap$x, levels = sample_order, labels = sample_labels)
overlap$y_label <- factor(overlap$y, levels = rev(sample_order), labels = rev(sample_labels))

base_heatmap <- list(
  geom_tile(colour = "white", linewidth = 0.25),
  scale_x_discrete(drop = FALSE), scale_y_discrete(drop = FALSE),
  coord_equal(), theme_minimal(base_size = 7),
  theme(panel.grid = element_blank(), axis.title = element_blank(),
        axis.text.x = element_text(angle = 35, hjust = 1), legend.position = "right")
)
p_j <- ggplot(overlap, aes(x_label, y_label, fill = jaccard)) + base_heatmap +
  geom_text(aes(label = ifelse(jaccard >= 0.001, sprintf("%.3f", jaccard), "0")), size = 2) +
  scale_fill_viridis_c(option = "B", name = "Jaccard")
p_r <- ggplot(overlap, aes(x_label, y_label, fill = raw)) + base_heatmap +
  geom_text(aes(label = raw), size = 2) +
  scale_fill_viridis_c(option = "B", name = "Raw")
p3b <- p_j + p_r + plot_layout(ncol = 2)
ggsave(file.path(out_dir, "suppfig3b_tcr_clonotype_overlap.pdf"), p3b,
       width = 13.5, height = 6.5, units = "in", device = cairo_pdf)

pair_defs <- list(
  list(patient = "0701", pbmc = "PBMC_0701_3", graft = "graft_0701_7", pbmc_lab = "PBMC 07012022", graft_lab = "Graft 07012022"),
  list(patient = "0708", pbmc = "PBMC_0708_4", graft = "graft_0708_8", pbmc_lab = "PBMC 07082022", graft_lab = "Graft 07082022"),
  list(patient = "0712", pbmc = "PBMC_0712_2", graft = "graft_0712_6", pbmc_lab = "PBMC 07122021", graft_lab = "Graft 07122021"),
  list(patient = "0825", pbmc = "PBMC_0825_1", graft = "graft_0825_5", pbmc_lab = "PBMC 08252021", graft_lab = "Graft 08252021")
)

make_scatter <- function(def) {
  count_one <- function(sample_name) {
    z <- combined[[sample_name]]
    tab <- as.data.frame(table(z$CTstrict), stringsAsFactors = FALSE)
    colnames(tab) <- c("clonotype", "n")
    tab[tab$clonotype != "" & !is.na(tab$clonotype), , drop = FALSE]
  }
  x <- count_one(def$pbmc); colnames(x)[2] <- "n_pbmc"
  y <- count_one(def$graft); colnames(y)[2] <- "n_graft"
  d <- merge(x, y, by = "clonotype", all = TRUE)
  d$n_pbmc[is.na(d$n_pbmc)] <- 0
  d$n_graft[is.na(d$n_graft)] <- 0
  d$prop_pbmc <- d$n_pbmc / sum(d$n_pbmc)
  d$prop_graft <- d$n_graft / sum(d$n_graft)
  d$total_n <- d$n_pbmc + d$n_graft
  d$class <- ifelse(d$n_pbmc > 0 & d$n_graft > 0, "Dual Expanded",
    ifelse(d$n_graft > 0 & d$n_graft > 1, paste(def$graft_lab, "Expanded"),
    ifelse(d$n_graft > 0, paste(def$graft_lab, "Singlet"),
    ifelse(d$n_pbmc > 1, paste(def$pbmc_lab, "Expanded"), paste(def$pbmc_lab, "Singlet")))))
  d$patient <- def$patient
  d$clonotype_id <- sprintf("%s_clonotype_%06d", def$patient, seq_len(nrow(d)))
  p <- ggplot(d, aes(prop_pbmc, prop_graft, size = total_n, colour = class)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", linewidth = 0.3, colour = "grey55") +
    geom_point(alpha = 0.82) +
    scale_colour_manual(values = c(
      "Dual Expanded" = "black",
      setNames(c("#542057", "#C58BBB", "#E57800", "#F4CF65"),
               c(paste(def$graft_lab, "Expanded"), paste(def$graft_lab, "Singlet"),
                 paste(def$pbmc_lab, "Expanded"), paste(def$pbmc_lab, "Singlet")))
    )) +
    scale_size_area(max_size = 5) +
    labs(x = "Proportion in PBMC", y = "Proportion in Graft", colour = "Class", size = "Total n") +
    theme_classic(base_size = 8) +
    theme(legend.position = "right")
  list(plot = p, data = d[, c("patient", "clonotype_id", "n_pbmc", "n_graft", "prop_pbmc", "prop_graft", "total_n", "class")])
}

scatter_results <- lapply(pair_defs, make_scatter)
p4 <- wrap_plots(lapply(scatter_results, `[[`, "plot"), ncol = 2)
ggsave(file.path(out_dir, "suppfig4_clonotype_expansion_graft_vs_pbmc.pdf"), p4,
       width = 13.5, height = 11, units = "in", device = cairo_pdf)
write.csv(do.call(rbind, lapply(scatter_results, `[[`, "data")),
          file.path(out_dir, "suppfig4_clonotype_expansion_plot_data.csv"), row.names = FALSE)

write.csv(overlap[, c("x", "y", "raw", "jaccard")],
          file.path(out_dir, "suppfig3b_clonotype_overlap_values.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(out_dir, "suppfig3_4_sessionInfo.txt"))
