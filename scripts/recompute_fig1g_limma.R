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
required <- c("patient", "condition", feature_col, "proportion_pct")
if (!all(required %in% names(d))) stop("Fig. 1G source table is missing required columns.")

condition_key <- tolower(trimws(as.character(d$condition)))
d$condition <- ifelse(condition_key == "pbmc", "PBMC",
  ifelse(condition_key == "graft", "Graft", NA_character_))
if (anyNA(d$condition)) stop("Unexpected condition labels in Fig. 1G source table.")

patients <- sort(unique(as.character(d$patient)))
features <- sort(unique(as.character(d[[feature_col]])), method = "radix")
grid <- expand.grid(patient = patients, condition = c("PBMC", "Graft"),
  feature = features, stringsAsFactors = FALSE)
names(grid)[3] <- feature_col
d <- merge(grid, d, by = c("patient", "condition", feature_col), all.x = TRUE)
d$proportion_pct[is.na(d$proportion_pct)] <- 0
d$patient <- factor(d$patient, levels = patients)
d$condition <- factor(d$condition, levels = c("PBMC", "Graft"))
d$feature <- as.character(d[[feature_col]])
d$sample <- paste(d$patient, d$condition, sep = "__")
d$y <- qlogis(pmin(pmax(d$proportion_pct / 100, 0.005), 0.995))

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
    "one 28-cluster matrix; design ~ patient + condition; bounded logit proportion"),
  row.names = NULL
)
write.csv(stats, file.path(out, "fig1g_limma_ebayes_logit_proportions.csv"), row.names = FALSE)
writeLines(c(
  "Unit of inference: four matched patients (P1-P4).",
  "Model: one 28 x 8 cluster-by-sample matrix fitted with limma lmFit/eBayes.",
  "Design: ~ patient + condition; contrast: Graft versus PBMC.",
  "Transformation: logit of sample-level proportions bounded at 0.005/0.995.",
  "Multiplicity: Benjamini-Hochberg across 28 clusters."
), file.path(out, "fig1g_provenance.txt"))
writeLines(capture.output(sessionInfo()), file.path(out, "fig1g_sessionInfo.txt"))
message("Wrote corrected Fig. 1G limma/eBayes table to ", out)
