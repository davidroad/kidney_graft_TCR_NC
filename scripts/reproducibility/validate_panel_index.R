#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
root <- normalizePath(if (length(args)) args[[1]] else ".", mustWork = TRUE)
index_path <- file.path(root, "provenance", "figure_panel_reproducibility_index.tsv")
index <- read.delim(index_path, check.names = FALSE, stringsAsFactors = FALSE, na.strings = c("NA", ""))

required_columns <- c(
  "figure_panel", "analysis_object", "unit_of_analysis", "upstream_accession",
  "private_upstream_object", "released_metrics", "reproduce_script", "released_output",
  "statistical_method", "status", "notes"
)
stopifnot(all(required_columns %in% colnames(index)))
stopifnot(!anyDuplicated(index$figure_panel))

expected_primary <- c(paste0("Fig1", LETTERS[1:7]), paste0("Fig2", LETTERS[1:9]))
missing_primary <- setdiff(expected_primary, index$figure_panel)
if (length(missing_primary)) stop("Missing primary panel rows: ", paste(missing_primary, collapse = ", "))

split_paths <- function(value) {
  if (is.na(value) || !nzchar(value)) return(character())
  trimws(strsplit(value, ";", fixed = TRUE)[[1]])
}

for (i in which(index$status == "reproducible")) {
  for (field in c("released_metrics", "reproduce_script", "released_output")) {
    paths <- split_paths(index[[field]][i])
    if (!length(paths)) stop(index$figure_panel[i], ": missing ", field)
    missing <- paths[!file.exists(file.path(root, paths))]
    if (length(missing)) stop(index$figure_panel[i], ": missing files in ", field, ": ", paste(missing, collapse = ", "))
  }
}

large_objects <- list.files(root, recursive = TRUE, full.names = TRUE, pattern = "\\.(RData|rdata|rds|h5ad)$")
if (length(large_objects)) stop("Large processed objects are present in the public deposit: ", paste(large_objects, collapse = ", "))

message(
  "Panel index validation passed: ", nrow(index), " indexed modules; ",
  sum(index$status == "reproducible"), " directly reproducible; no large processed objects."
)
