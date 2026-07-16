# Kidney Graft TCR Analysis

Long-term alloimmune responses can persist for years after transplantation, but the cellular organization that sustains them remains incompletely understood. The accompanying study combines paired single-cell RNA sequencing and TCR sequencing from rejected human kidney allografts and peripheral blood collected approximately a decade after transplantation. TCR clonotypes shared between blood and graft connect circulating TCF1+ stem-like CD8 T cells with graft-infiltrating CXCR6+ cytotoxic effector states, supporting a model in which a persistent stem-like reservoir repeatedly replenishes short-lived effectors during ongoing rejection.

This repository provides the reproducible human computational analyses and lightweight source data supporting the single-cell/TCR and revision components of the study. It includes paired cluster-composition modeling, shared-clonotype and cell-cycle analyses, CXCR6-focused comparisons, and transcriptional metabolic analyses. Large processed objects and raw sequencing data are provided through the linked Figshare and GEO records rather than stored in GitHub.

## Repository contents

- `scripts/`: reproducible scripts for Fig. 1G/Fig. 2B paired limma modeling and plotting, shared-clonotype T-cell cycle scoring, and HUMESS-derived figure generation.
- `scripts/raw_to_intermediate/`: functional programs for rebuilding the primary all-cell, T-cell/shared-clonotype, and CD8 intermediate objects.
- `scripts/source_data_export/`: functional programs for exporting primary figure source tables and the shared-clone trajectory.
- `scripts/primary_scrna_tcr_analysis.md`: main manuscript scRNA/TCR workflow and data contract.
- `tables/`: de-identified sample-level source data and statistical results.
- `provenance/`: analysis notes and R session information.
- `DATA_AVAILABILITY.md`: public accessions and required processed inputs.
- `environment.yml`: reproducible software environment.

## Main analysis definitions

- Fig. 1G contains 28 CD45+ clusters measured in four matched graft/PBMC pairs. One `28 clusters x 8 samples` bounded-logit matrix is fitted with `limma::lmFit` and one `eBayes` call using `~ patient + condition`; BH correction is applied across 28 clusters.
- Fig. 2B contains 11 annotated T-cell clusters measured in the same four matched pairs. One `11 clusters x 8 samples` bounded-logit matrix is fitted using the same donor-blocked model; BH correction is applied across 11 clusters.
- The Fig. 2D shared-clonotype object contains 4,110 T cells, including CD8, CD4, and Treg annotations. Population-level cell-cycle metrics use this full shared-clonotype population, while the annotated CXCR6+ effector CD8 states retain the manuscript's CD8-focused biological context.
- Patient/sample is the unit of inference. Cell-level tests are retained only as explicitly labeled exploratory source data.

The full manuscript pipeline is documented in [`scripts/primary_scrna_tcr_analysis.md`](scripts/primary_scrna_tcr_analysis.md). The scripts in `scripts/revision/` are downstream revision and figure-support analyses; they are not substitutes for the main eight-sample integration and shared-TCR workflow.

## Running the analyses

Create the software environment:

```bash
conda env create -f environment.yml
conda activate kidney-graft-tcr-nc
```

Recompute and plot Fig. 1G:

```bash
Rscript scripts/revision/recompute_fig1g_limma.R --repo . --out results/fig1g
python scripts/revision/plot_fig1g_limma_paired_barplot.py --out results/fig1g
```

Recompute Fig. 2B:

```bash
Rscript scripts/revision/recompute_fig2b_from_latest_tcell_object.R \
  --input data/processed/merged_4_with_vdj_T_cells_filter_10_28_2025.RData \
  --out results/fig2b
```

Render the corrected Fig. 2B panel:

```bash
python scripts/revision/plot_fig2b_limma_barplot.py --out results/fig2b
```

Recompute shared-clonotype T-cell cycle scores and figures:

```bash
Rscript scripts/revision/06_reviewer3_shared_clone_t_cell_cycle.R \
  --input data/processed/all_clone_filter_10_28.Rdata \
  --out results/shared_clone_t_cell_cycle
```

Generate HUMESS-derived revision figures after staging the inputs described in `DATA_AVAILABILITY.md`:

```bash
Rscript scripts/revision/03_make_article_humess_visuals.R --repo .
```

Validate the lightweight public release:

```bash
Rscript scripts/validate_release.R
```

## Corrected paired results

For Fig. 1G, clusters 0, 2, 4, 5, 6, 8, 9, and 11 have FDR < 0.05. For Fig. 2B, `CXCL13+CXCR6+ effector CD8+` is higher in graft, whereas `CCR2+ CD4+` and `Naive CD4+` are lower in graft at FDR < 0.05. `CXCR6+ effector CD8+` has raw p = 0.0331 but FDR = 0.0910. These estimates are based on four matched patients and should be interpreted as paired cohort evidence rather than population-level prevalence estimates.

## Citation

Please cite the associated Nature Communications article and the Figshare data and reproducibility record:

> Dai Y. *Kidney graft TCR clonotype analysis: revision figures, source data, and reproducibility files*. Figshare. 2026. https://doi.org/10.6084/m9.figshare.33001412.v1
