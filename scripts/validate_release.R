#!/usr/bin/env Rscript

repo <- normalizePath(getwd(), mustWork = TRUE)
required <- c(
  "README.md",
  "DATA_AVAILABILITY.md",
  "environment.yml",
  "scripts/revision/recompute_fig1g_limma.R",
  "scripts/revision/plot_fig1g_limma_paired_barplot.py",
  "scripts/revision/recompute_fig2b_from_latest_tcell_object.R",
  "scripts/revision/06_reviewer3_shared_clone_t_cell_cycle.R",
  "scripts/revision/03_make_article_humess_visuals.R",
  "tables/fig1g_current_cluster_proportions_by_patient.csv",
  "tables/fig1g_limma_ebayes_logit_proportions.csv",
  "tables/fig2b_latest_11cluster_sample_level_proportions.csv",
  "tables/fig2b_limma_ebayes_logit_proportions.csv",
  "tables/reviewer3_shared_clone_t_cell_cycle_patient_summary.csv"
)
missing <- required[!file.exists(file.path(repo, required))]
if (length(missing)) stop("Missing release files: ", paste(missing, collapse = ", "))

stats1 <- read.csv(file.path(repo, "tables", "fig1g_limma_ebayes_logit_proportions.csv"))
source1 <- read.csv(file.path(repo, "tables", "fig1g_current_cluster_proportions_by_patient.csv"))
stats <- read.csv(file.path(repo, "tables", "fig2b_limma_ebayes_logit_proportions.csv"))
source <- read.csv(file.path(repo, "tables", "fig2b_latest_11cluster_sample_level_proportions.csv"))
shared_t <- read.csv(file.path(repo, "tables", "reviewer3_shared_clone_t_cell_cycle_patient_summary.csv"))
stopifnot(
  nrow(stats1) == 28L,
  !anyNA(stats1[c("p_value", "FDR")]),
  all(stats1$FDR + 1e-12 >= stats1$p_value),
  !all(abs(stats1$FDR - stats1$p_value) < 1e-12),
  identical(sort(unique(as.character(source1$patient))), paste0("P", 1:4)),
  nrow(stats) == 11L,
  !anyNA(stats[c("p_value", "FDR")]),
  all(stats$FDR + 1e-12 >= stats$p_value),
  !all(abs(stats$FDR - stats$p_value) < 1e-12),
  identical(sort(unique(as.character(source$patient))), paste0("P", 1:4)),
  all(c("PBMC", "Graft") %in% source$condition),
  nrow(shared_t) == 16L,
  sum(shared_t$n_cells) == 4110L,
  identical(sort(unique(as.character(shared_t$patient))), paste0("P", 1:4)),
  setequal(as.character(shared_t$label), c("PBMC", "Graft")),
  setequal(as.character(shared_t$CXCR6_status), c("CXCR6-", "CXCR6+"))
)

text_files <- list.files(repo, recursive = TRUE, full.names = TRUE, pattern = "\\.(R|md|csv|tsv|txt|yml)$")
text_files <- setdiff(text_files, file.path(repo, "scripts", "validate_release.R"))
text <- unlist(lapply(text_files, readLines, warn = FALSE), use.names = FALSE)
forbidden <- c("/home/0.collaboration", "/home/26_immune_NC", "/data_sys/collab2", "Patient 701", "Patient 708", "Patient 712", "Patient 825")
hits <- forbidden[vapply(forbidden, function(x) any(grepl(x, text, fixed = TRUE)), logical(1))]
if (length(hits)) stop("Internal identifiers remain in the release: ", paste(hits, collapse = ", "))

message("Release validation passed: joint Fig. 1G/Fig. 2B limma models, 4,110 shared-clonotype T cells, P1-P4 source data, and no internal paths.")
