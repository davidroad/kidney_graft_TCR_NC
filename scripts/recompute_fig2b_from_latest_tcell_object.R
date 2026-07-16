#!/usr/bin/env Rscript

options(stringsAsFactors = FALSE)
suppressPackageStartupMessages({
  library(Seurat)
  library(limma)
})

args <- commandArgs(trailingOnly = TRUE)
arg_value <- function(flag, default = NULL) {
  pos <- match(flag, args)
  if (is.na(pos) || pos == length(args)) default else args[[pos + 1L]]
}

repo_dir <- normalizePath(arg_value("--repo", getwd()), mustWork = TRUE)
obj_file <- arg_value(
  "--input",
  Sys.getenv(
    "FIG2B_INPUT_RDATA",
    file.path(repo_dir, "data", "processed", "merged_4_with_vdj_T_cells_filter_10_28_2025.RData")
  )
)
out <- arg_value("--out", file.path(repo_dir, "results", "fig2b"))

if (!file.exists(obj_file)) {
  stop(
    "Fig. 2B input object not found: ", obj_file, "\n",
    "Provide it with --input /path/to/file.RData or FIG2B_INPUT_RDATA."
  )
}
dir.create(out, recursive = TRUE, showWarnings = FALSE)

e <- new.env(parent = emptyenv())
load(obj_file, envir = e)
obj_name <- intersect(ls(e), "merged_all_T_batchcorrected_filter")[1]
if (is.na(obj_name)) stop("Object merged_all_T_batchcorrected_filter was not found in the input RData.")
obj <- get(obj_name, e)

# Original Fig. 2B definition: 11 annotated T-cell clusters within each
# matched graft/PBMC T-cell sample.
cluster_levels <- c(
  "Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+",
  "CXCR6+ effector CD8+", "CCR2+ CD4+", "Naive/Memory-like CD8+",
  "Naive CD4+", "Treg", "CX3CR1+ effector CD8+",
  "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+"
)

md <- obj@meta.data
if (!all(c("batch", "label") %in% names(md))) stop("Input metadata must contain batch and label.")
md$cluster <- factor(as.character(Idents(obj)), levels = cluster_levels)
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
tab$proportion_pct <- 100 * tab$n / tab$denom
tab <- tab[order(tab$patient, match(tab$condition, c("PBMC", "Graft")), match(tab$tcluster, cluster_levels)), ]
write.csv(tab, file.path(out, "fig2b_latest_11cluster_sample_level_proportions.csv"), row.names = FALSE)

# Fit all 11 cluster responses together so eBayes estimates a common variance
# prior across clusters. The 0.5-count empirical-logit correction handles zeros.
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
fit <- eBayes(fit, trend = FALSE, robust = FALSE)
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
    "limma lmFit + eBayes across 11 clusters; design ~ patient + condition;",
    "empirical logit log((n+0.5)/(denom-n+0.5)); BH across clusters"
  ),
  stringsAsFactors = FALSE
)
write.csv(stats, file.path(out, "fig2b_limma_ebayes_logit_proportions.csv"), row.names = FALSE)

writeLines(
  c(
    paste("Source:", basename(obj_file)),
    paste("Object:", obj_name),
    "Unit of inference: four matched patients (P1-P4).",
    "Features: 11 annotated T-cell clusters used in Fig. 2B.",
    "Model: one 11 x 8 cluster-by-sample matrix fitted with limma lmFit/eBayes.",
    "Design: ~ patient + condition; contrast: Graft versus PBMC.",
    "Transformation: empirical logit log((n+0.5)/(denom-n+0.5)).",
    "Multiplicity: Benjamini-Hochberg across 11 clusters."
  ),
  file.path(out, "fig2b_provenance.txt")
)
writeLines(capture.output(sessionInfo()), file.path(out, "fig2b_sessionInfo.txt"))
message("Wrote corrected 11-cluster Fig. 2B source and limma/eBayes tables to ", out)
