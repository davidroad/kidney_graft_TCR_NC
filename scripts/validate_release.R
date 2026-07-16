#!/usr/bin/env Rscript

repo <- normalizePath(getwd(), mustWork = TRUE)
required <- c(
  "README.md",
  "DATA_AVAILABILITY.md",
  "environment.yml",
  "scripts/recompute_fig2b_from_latest_tcell_object.R",
  "scripts/06_reviewer3_shared_clone_t_cell_cycle.R",
  "scripts/03_make_article_humess_visuals.R",
  "tables/fig2b_latest_11cluster_sample_level_proportions.csv",
  "tables/fig2b_limma_ebayes_logit_proportions.csv"
)
missing <- required[!file.exists(file.path(repo, required))]
if (length(missing)) stop("Missing release files: ", paste(missing, collapse = ", "))

stats <- read.csv(file.path(repo, "tables", "fig2b_limma_ebayes_logit_proportions.csv"))
source <- read.csv(file.path(repo, "tables", "fig2b_latest_11cluster_sample_level_proportions.csv"))
stopifnot(
  nrow(stats) == 11L,
  !anyNA(stats[c("p_value", "FDR")]),
  all(stats$FDR + 1e-12 >= stats$p_value),
  !all(abs(stats$FDR - stats$p_value) < 1e-12),
  identical(sort(unique(source$patient)), paste0("P", 1:4)),
  all(c("PBMC", "Graft") %in% source$condition)
)

text_files <- list.files(repo, recursive = TRUE, full.names = TRUE, pattern = "\\.(R|md|csv|tsv|txt|yml)$")
text_files <- setdiff(text_files, file.path(repo, "scripts", "validate_release.R"))
text <- unlist(lapply(text_files, readLines, warn = FALSE), use.names = FALSE)
forbidden <- c("/home/0.collaboration", "/home/26_immune_NC", "/data_sys/collab2", "Patient 701", "Patient 708", "Patient 712", "Patient 825")
hits <- forbidden[vapply(forbidden, function(x) any(grepl(x, text, fixed = TRUE)), logical(1))]
if (length(hits)) stop("Internal identifiers remain in the release: ", paste(hits, collapse = ", "))

message("Release validation passed: 11-cluster joint limma model, P1-P4 source data, and no internal paths.")
