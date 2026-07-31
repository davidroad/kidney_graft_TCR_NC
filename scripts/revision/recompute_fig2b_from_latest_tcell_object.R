#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) default else args[[pos + 1L]]
}

repo_dir <- normalizePath(arg_value("--repo", getwd()), mustWork = TRUE)
table_input <- arg_value("--table-input", NULL)
obj_file <- arg_value(
  "--input",
  Sys.getenv(
    "FIG2B_INPUT_RDATA",
    file.path(repo_dir, "data", "processed", "merged_4_with_vdj_T_cells_filter_10_28_2025.RData")
  )
)
out <- arg_value("--out", file.path(repo_dir, "results", "fig2b"))

if (is.null(table_input) && !file.exists(obj_file)) {
  stop(
    "Fig. 2B input object not found: ", obj_file, "\n",
    "Provide it with --input /path/to/file.RData or FIG2B_INPUT_RDATA, ",
    "or use --table-input with the released sample-level count table."
  )
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# Original Fig. 2B definition: 11 annotated T-cell clusters within each
# matched graft/PBMC T-cell sample.
cluster_levels <- c(
  "Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+",
  "CXCR6+ effector CD8+", "CCR2+ CD4+", "Naive/Memory-like CD8+",
  "Naive CD4+", "Treg", "CX3CR1+ effector CD8+",
  "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+"
)

if (!is.null(table_input)) {
  if (!file.exists(table_input)) stop("Fig. 2B table input not found: ", table_input)
  tab <- read.csv(table_input, check.names = FALSE)
  required <- c("patient", "condition", "tcluster", "n", "denom")
  if (!all(required %in% names(tab))) {
    stop("Fig. 2B table input is missing required columns.")
  }
  source_label <- basename(table_input)
  object_label <- "released sample-level count table"
} else {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("The Seurat package is required when an RData object is used.")
  }
  e <- new.env(parent = emptyenv())
  load(obj_file, envir = e)
  obj_name <- intersect(ls(e), "merged_all_T_batchcorrected_filter")[1]
  if (is.na(obj_name)) {
    stop("Object merged_all_T_batchcorrected_filter was not found in the input RData.")
  }
  obj <- get(obj_name, e)
  md <- obj@meta.data
  if (!all(c("batch", "label") %in% names(md))) {
    stop("Input metadata must contain batch and label.")
  }
  md$cluster <- factor(as.character(Seurat::Idents(obj)), levels = cluster_levels)
  md$patient_internal <- sub("^(graft|PBMC)_", "", as.character(md$batch))
  md$condition <- ifelse(grepl("^graft", md$batch, ignore.case = TRUE), "Graft", "PBMC")
  md <- md[!is.na(md$cluster) & md$condition %in% c("Graft", "PBMC"), , drop = FALSE]

  patient_levels <- sort(unique(md$patient_internal))
  patient_map <- setNames(sprintf("P%d", seq_along(patient_levels)), patient_levels)
  md$patient <- unname(patient_map[md$patient_internal])

  counts <- as.data.frame(table(md$patient, md$condition, md$cluster), stringsAsFactors = FALSE)
  names(counts) <- c("patient", "condition", "tcluster", "n")
  denom <- aggregate(n ~ patient + condition, counts, sum)
  names(denom)[3] <- "denom"
  tab <- merge(counts, denom, by = c("patient", "condition"), all.x = TRUE)
  source_label <- basename(obj_file)
  object_label <- obj_name
}

tab$patient <- as.character(tab$patient)
condition_key <- tolower(trimws(as.character(tab$condition)))
tab$condition <- ifelse(condition_key == "pbmc", "PBMC",
  ifelse(condition_key == "graft", "Graft", NA_character_))
tab$tcluster <- as.character(tab$tcluster)
tab$n <- as.numeric(tab$n)
tab$denom <- as.numeric(tab$denom)
if (anyNA(tab$condition) || anyNA(tab$n) || anyNA(tab$denom) ||
    any(tab$n < 0) || any(tab$denom <= 0) || any(tab$n > tab$denom)) {
  stop("Fig. 2B table contains invalid condition labels or cell counts.")
}
if (!setequal(unique(tab$tcluster), cluster_levels)) {
  stop("Fig. 2B table does not contain exactly the expected 11 T-cell clusters.")
}
if (anyDuplicated(tab[c("patient", "condition", "tcluster")])) {
  stop("Fig. 2B table contains duplicated patient-condition-cluster rows.")
}
sample_totals <- aggregate(n ~ patient + condition, tab, sum)
sample_denoms <- unique(tab[c("patient", "condition", "denom")])
if (anyDuplicated(sample_denoms[c("patient", "condition")])) {
  stop("More than one denominator was found for at least one Fig. 2B sample.")
}
denom_check <- merge(sample_totals, sample_denoms,
  by = c("patient", "condition"), all = TRUE)
if (anyNA(denom_check) || any(abs(denom_check$n - denom_check$denom) > 1e-8)) {
  stop("Fig. 2B cluster counts do not sum to the reported sample denominators.")
}
tab$proportion_pct <- 100 * tab$n / tab$denom
tab <- tab[order(tab$patient, match(tab$condition, c("PBMC", "Graft")), match(tab$tcluster, cluster_levels)), ]
write.csv(tab, file.path(out, "fig2b_latest_11cluster_sample_level_proportions.csv"), row.names = FALSE)

# Fit all 11 cluster responses together so eBayes estimates a common variance
# prior across clusters. A 0.5-count empirical-logit correction handles zeros
# without collapsing distinct small non-zero counts to a fixed proportion.
sample_meta <- unique(tab[c("patient", "condition")])
sample_meta <- sample_meta[order(sample_meta$patient, match(sample_meta$condition, c("PBMC", "Graft"))), ]
sample_meta$patient <- factor(sample_meta$patient)
sample_meta$condition <- factor(sample_meta$condition, levels = c("PBMC", "Graft"))
sample_meta$sample_id <- paste(sample_meta$patient, sample_meta$condition, sep = "__")

response <- matrix(
  NA_real_,
  nrow = length(cluster_levels),
  ncol = nrow(sample_meta),
  dimnames = list(cluster_levels, sample_meta$sample_id)
)
for (i in seq_len(nrow(sample_meta))) {
  z <- tab[
    tab$patient == as.character(sample_meta$patient[i]) &
      tab$condition == as.character(sample_meta$condition[i]),
  ]
  z <- z[match(cluster_levels, z$tcluster), ]
  response[, i] <- log((z$n + 0.5) / (z$denom - z$n + 0.5))
}
if (anyNA(response)) stop("The cluster-by-sample response matrix contains missing values.")

design <- model.matrix(~ patient + condition, data = sample_meta)
fit <- lmFit(response, design)
fit <- eBayes(fit, trend = TRUE, robust = FALSE)
tt <- topTable(
  fit,
  coef = "conditionGraft",
  number = Inf,
  sort.by = "none",
  adjust.method = "BH"
)

stats <- data.frame(
  analysis = "Fig2B_latest_11cluster_frequency",
  feature = rownames(tt),
  n_pairs = length(unique(sample_meta$patient)),
  logit_effect_graft_vs_pbmc = tt$logFC,
  moderated_t = tt$t,
  p_value = tt$P.Value,
  FDR = tt$adj.P.Val,
  model = paste(
    "limma lmFit + eBayes(trend=TRUE, robust=FALSE) across 11 clusters;",
    "design ~ patient + condition; count-based empirical logit with",
    "0.5 pseudocount; BH across clusters"
  ),
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(out, "fig2b_limma_ebayes_logit_proportions.csv"), row.names = FALSE)

writeLines(
  c(
    paste("Source:", source_label),
    paste("Object:", object_label),
    "Unit of inference: four matched patients (P1-P4).",
    "Features: 11 annotated T-cell clusters used in Fig. 2B.",
    "Model: one 11 x 8 cluster-by-sample matrix fitted with limma lmFit/eBayes.",
    "Design: ~ patient + condition; contrast: Graft versus PBMC.",
    paste(
      "Transformation: count-based empirical logit",
      "log((cluster count + 0.5)/(total count - cluster count + 0.5))."
    ),
    "Multiplicity: Benjamini-Hochberg across 11 clusters."
  ),
  file.path(out, "fig2b_provenance.txt")
)
writeLines(capture.output(sessionInfo()), file.path(out, "fig2b_sessionInfo.txt"))
message("Wrote corrected 11-cluster Fig. 2B source and limma/eBayes tables to ", out)
