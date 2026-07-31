#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages(library(limma))

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) default else args[[pos + 1L]]
}

repo <- normalizePath(arg_value("--repo", getwd()), mustWork = TRUE)
input <- arg_value("--input", file.path(repo, "tables", "fig1g_current_cluster_proportions_by_patient.csv"))
out <- arg_value("--out", file.path(repo, "results", "fig1g"))
dir.create(out, recursive = TRUE, showWarnings = FALSE)

d <- read.csv(input, check.names = FALSE)
feature_col <- "integrated_snn_res.0.5"
required <- c("patient", "condition", feature_col, "n", "denom", "proportion_pct")
if (!all(required %in% names(d))) stop("Fig. 1G source table is missing required columns.")

condition_key <- tolower(trimws(as.character(d$condition)))
d$condition <- ifelse(condition_key == "pbmc", "PBMC",
  ifelse(condition_key == "graft", "Graft", NA_character_))
if (anyNA(d$condition)) stop("Unexpected condition labels in Fig. 1G source table.")
d[[feature_col]] <- as.character(d[[feature_col]])
d$n <- as.numeric(d$n)
d$denom <- as.numeric(d$denom)
if (anyNA(d$n) || anyNA(d$denom) || any(d$n < 0) || any(d$denom <= 0) ||
    any(d$n > d$denom)) {
  stop("Fig. 1G n and denom columns must contain valid non-negative cell counts.")
}

patients <- sort(unique(as.character(d$patient)))
features <- sort(unique(as.character(d[[feature_col]])), method = "radix")
sample_denom <- unique(d[c("patient", "condition", "denom")])
sample_denom$key <- paste(sample_denom$patient, sample_denom$condition, sep = "__")
if (anyDuplicated(sample_denom$key)) {
  stop("More than one denominator was found for at least one Fig. 1G sample.")
}
denom_lookup <- setNames(sample_denom$denom, sample_denom$key)

grid <- expand.grid(patient = patients, condition = c("PBMC", "Graft"),
  feature = features, stringsAsFactors = FALSE)
names(grid)[3] <- feature_col
d <- merge(grid, d, by = c("patient", "condition", feature_col), all.x = TRUE)
missing_n <- is.na(d$n)
d$n[missing_n] <- 0
d$denom[missing_n] <- unname(denom_lookup[
  paste(d$patient[missing_n], d$condition[missing_n], sep = "__")
])
if (anyNA(d$denom)) stop("A Fig. 1G sample denominator could not be restored.")
d$proportion_pct <- 100 * d$n / d$denom
d$patient <- factor(d$patient, levels = patients)
d$condition <- factor(d$condition, levels = c("PBMC", "Graft"))
d$feature <- as.character(d[[feature_col]])
d$sample <- paste(d$patient, d$condition, sep = "__")
d$y <- log((d$n + 0.5) / (d$denom - d$n + 0.5))

sample_meta <- unique(d[c("patient", "condition", "sample")])
sample_meta <- sample_meta[order(sample_meta$patient, sample_meta$condition), ]
wide <- reshape(d[c("feature", "sample", "y")], idvar = "feature",
  timevar = "sample", direction = "wide")
rownames(wide) <- wide$feature
wide$feature <- NULL
wide <- wide[features, paste0("y.", sample_meta$sample), drop = FALSE]
colnames(wide) <- sample_meta$sample
if (anyNA(wide)) stop("Incomplete Fig. 1G feature-by-sample matrix.")

design <- model.matrix(~ patient + condition, data = sample_meta)
fit <- eBayes(lmFit(as.matrix(wide), design), trend = TRUE, robust = FALSE)
tt <- topTable(fit, coef = "conditionGraft", number = Inf,
  sort.by = "none", adjust.method = "BH")
stats <- data.frame(
  analysis = "Fig1G",
  feature = rownames(tt),
  n_pairs = length(patients),
  logit_effect_graft_vs_pbmc = tt$logFC,
  moderated_t = tt$t,
  p_value = tt$P.Value,
  FDR = tt$adj.P.Val,
  model = paste("limma lmFit + eBayes(trend=TRUE, robust=FALSE);",
    "one 28-cluster matrix; design ~ patient + condition;",
    "count-based empirical logit with 0.5 pseudocount"),
  row.names = NULL
)
write.csv(stats, file.path(out, "fig1g_limma_ebayes_logit_proportions.csv"), row.names = FALSE)
message("Wrote corrected Fig. 1G limma/eBayes table to ", out)
