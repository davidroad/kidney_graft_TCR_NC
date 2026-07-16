# Kidney Graft TCR Analysis

This repository contains reproducible analysis code and lightweight source data for a paired human kidney graft/PBMC single-cell RNA-seq and TCR-seq study. The study examines how TCR clonotypes shared between blood and graft connect circulating stem-like and graft-infiltrating T-cell states, with focused analyses of CXCR6-associated differentiation, T-cell cluster composition, transcriptional metabolic programs, and in situ cell-cycle activity.

## Repository contents

- `scripts/`: R scripts for Fig. 2B cluster-frequency modeling, shared-clonotype T-cell cycle scoring, and HUMESS-derived figure generation.
- `tables/`: de-identified sample-level source data and statistical results.
- `provenance/`: analysis notes and R session information.
- `DATA_AVAILABILITY.md`: public accessions and required processed inputs.
- `environment.yml`: reproducible software environment.

## Main analysis definitions

- Fig. 2B contains 11 annotated T-cell clusters measured in four matched graft/PBMC pairs. The corrected model fits one `11 clusters x 8 samples` empirical-logit matrix with `limma::lmFit` and `eBayes`, using `~ patient + condition`; Benjamini-Hochberg correction is applied across all 11 clusters.
- The Fig. 2D shared-clonotype object contains 4,110 T cells, including CD8, CD4, and Treg annotations. It is therefore described as a shared-clonotype T-cell analysis, not a CD8-restricted analysis.
- Patient/sample is the unit of inference. Cell-level tests are retained only as explicitly labeled exploratory source data.

## Running the analyses

Create the software environment:

```bash
conda env create -f environment.yml
conda activate kidney-graft-tcr-nc
```

Recompute Fig. 2B:

```bash
Rscript scripts/recompute_fig2b_from_latest_tcell_object.R \
  --input data/processed/merged_4_with_vdj_T_cells_filter_10_28_2025.RData \
  --out results/fig2b
```

Render the corrected Fig. 2B panel:

```bash
python scripts/plot_fig2b_limma_barplot.py --out results/fig2b
```

Recompute shared-clonotype T-cell cycle scores and figures:

```bash
Rscript scripts/06_reviewer3_shared_clone_t_cell_cycle.R \
  --input data/processed/all_clone_filter_10_28.Rdata \
  --out results/shared_clone_t_cell_cycle
```

Generate HUMESS-derived revision figures after staging the inputs described in `DATA_AVAILABILITY.md`:

```bash
Rscript scripts/03_make_article_humess_visuals.R --repo .
```

Validate the lightweight public release:

```bash
Rscript scripts/validate_release.R
```

## Corrected Fig. 2B result

At FDR < 0.05, graft samples showed higher proportions of `CXCL13+CXCR6+ effector CD8+` and `CXCR6+ effector CD8+` clusters, and lower proportions of `CCR2+ CD4+` and `Naive CD4+` clusters. These estimates are based on four matched patients and should be interpreted as paired cohort evidence rather than population-level prevalence estimates.

## Citation

Please cite the associated Nature Communications article and Figshare record when available. Repository metadata will be updated with the final DOI before publication.
