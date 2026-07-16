# R 4.5 Seurat + limma runbook

Verified analysis versions: R 4.5.3, Seurat 5.5.1, limma 3.66.0, and statmod
1.5.2. Create the public environment with `environment.yml`.

## Fig. 1G

```bash
Rscript scripts/recompute_fig1g_limma.R --repo . --out results/fig1g
python scripts/plot_fig1g_limma_paired_barplot.py \
  --source tables/fig1g_current_cluster_proportions_by_patient.csv \
  --stats tables/fig1g_limma_ebayes_logit_proportions.csv \
  --out results/fig1g
```

The model fits one 28-cluster by 8-sample matrix using bounded logit sample
proportions, design `~ patient + condition`, one
`eBayes(trend=TRUE, robust=FALSE)` call, and BH correction across 28 clusters.

## Fig. 2B

```bash
Rscript scripts/recompute_fig2b_from_latest_tcell_object.R \
  --repo . --input data/processed/merged_4_with_vdj_T_cells_filter_10_28_2025.RData \
  --out results/fig2b
python scripts/plot_fig2b_limma_barplot.py \
  --source tables/fig2b_latest_11cluster_sample_level_proportions.csv \
  --stats tables/fig2b_limma_ebayes_logit_proportions.csv \
  --out results/fig2b
```

The model fits one 11-cluster by 8-sample matrix using the same transformation,
paired design, one empirical-Bayes step, and BH correction across 11 clusters.

Always preserve the source table, statistical result, provenance, session
information, and script version together.
