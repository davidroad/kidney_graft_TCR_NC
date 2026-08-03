#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(scRepertoire)
})

args <- commandArgs(trailingOnly = TRUE)
manifest_file <- if (length(args) >= 1) args[[1]] else "scripts/reproducibility/fig2c_filtered_contig_manifest_template.csv"
output_dir <- if (length(args) >= 2) args[[2]] else "results/figure_2/fig2c_patient_panels"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

manifest <- read.csv(manifest_file, stringsAsFactors = FALSE, check.names = FALSE)
required_columns <- c("sample", "filtered_contig_annotations")
if (!all(required_columns %in% colnames(manifest))) {
  stop("Manifest must contain: ", paste(required_columns, collapse = ", "))
}

sample_order <- c(
  "PBMC_0825", "PBMC_0712", "PBMC_0701", "PBMC_0708",
  "graft_0825", "graft_0712", "graft_0701", "graft_0708"
)
idx <- match(sample_order, manifest$sample)
if (anyNA(idx)) {
  stop("Manifest is missing: ", paste(sample_order[is.na(idx)], collapse = ", "))
}
manifest <- manifest[idx, , drop = FALSE]

missing_files <- manifest$filtered_contig_annotations[
  !file.exists(manifest$filtered_contig_annotations)
]
if (length(missing_files)) {
  stop("Missing filtered contig annotation file(s): ", paste(missing_files, collapse = ", "))
}

contig_list <- lapply(manifest$filtered_contig_annotations, read.csv)
combined <- combineTCR(
  contig_list,
  samples = sample_order,
  ID = as.character(seq_along(sample_order))
)

comparisons <- list(
  "0701" = c("PBMC_0701_3", "graft_0701_7"),
  "0708" = c("PBMC_0708_4", "graft_0708_8"),
  "0712" = c("PBMC_0712_2", "graft_0712_6"),
  "0825" = c("PBMC_0825_1", "graft_0825_5")
)

for (patient in names(comparisons)) {
  out_file <- file.path(output_dir, paste0("fig2c_alluvial_", patient, ".pdf"))
  grDevices::pdf(out_file)
  print(clonalCompare(
    input.data = combined,
    cloneCall = "aa",
    samples = comparisons[[patient]],
    top.clones = 10,
    graph = "alluvial"
  ))
  grDevices::dev.off()
}

message("Fig. 2C patient alluvial panels written to ", normalizePath(output_dir))
