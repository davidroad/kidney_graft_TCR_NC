# Kidney Graft TCR Analysis

Long-term alloimmune responses can persist for years after transplantation, but the cellular organization that sustains them remains incompletely understood. The accompanying study combines paired single-cell RNA sequencing and TCR sequencing from rejected human kidney allografts and peripheral blood collected approximately a decade after transplantation. TCR clonotypes shared between blood and graft connect circulating TCF1+ stem-like CD8 T cells with graft-infiltrating CXCR6+ cytotoxic effector states, supporting a model in which a persistent stem-like reservoir repeatedly replenishes short-lived effectors during ongoing rejection.

This repository provides the reproducible human computational workflow for the primary single-cell/TCR analyses in the study, together with downstream revision modules. The primary route rebuilds the integrated all-CD45 and T-cell objects, TCR clonotypes, shared graft/PBMC clonotypes, and main-figure source tables. Revision modules then add paired limma models, cell-cycle, CXCR6, metabolic, HUMESS, and external-PBMC analyses. Raw sequencing data are available from the documented GEO accessions, while lightweight plotting metrics and source tables are provided through Figshare; large processed objects are not stored in either public code package.

## Repository contents

- `scripts/revision/`: downstream revision and figure-support analyses, including paired limma modeling, cell-cycle, CXCR6, metabolic, HUMESS, and external-PBMC analyses.
- `scripts/raw_to_intermediate/`: functional programs for rebuilding the primary all-cell, T-cell/shared-clonotype, and CD8 intermediate objects.
- `scripts/source_data_export/`: functional programs for exporting primary figure source tables and the shared-clone trajectory.
- `scripts/reproducibility/`: metrics exporters and plotting wrappers for rebuilding manuscript computational panels from the lightweight Figshare metrics.
- `scripts/primary_scrna_tcr_analysis.md`: main manuscript scRNA/TCR workflow and data contract.
- `provenance/`: analysis notes and R session information.
- `DATA_AVAILABILITY.md`: public accessions and required processed inputs.
- `environment.yml`: reproducible software environment.

## Main analysis definitions

- Fig. 1G contains 28 CD45+ clusters measured in four matched graft/PBMC pairs. One `28 clusters x 8 samples` bounded-logit matrix is fitted with `limma::lmFit` and one `eBayes` call using `~ patient + condition`; BH correction is applied across 28 clusters.
- Fig. 2B contains 11 annotated T-cell clusters measured in the same four matched pairs. One `11 clusters x 8 samples` bounded-logit matrix is fitted using the same donor-blocked model; BH correction is applied across 11 clusters.
- The Fig. 2D shared-clonotype object contains 4,110 T cells, including CD8, CD4, and Treg annotations. Population-level cell-cycle metrics use this full shared-clonotype population, while the annotated CXCR6+ effector CD8 states retain the manuscript's CD8-focused biological context.
- The Fig. 2I Monocle 3 trajectory was rooted in the TCF7-defined stem-like T-cell state, with pseudotime progressing toward differentiated effector states.
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

Generate the reviewer-requested five-group marker plots and TCF7+CXCR6+ frequency summaries after reconstructing both CD8 objects:

```bash
Rscript scripts/revision/07_generate_reviewer2_five_group_featureplots.R \
  data/processed/merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData \
  data/processed/GSE224445_strict_cd8_seurat.rds results/reviewer2
Rscript scripts/revision/09_summarize_reviewer2_tcf7_cxcr6_double_positive.R \
  data/processed/merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData \
  data/processed/GSE224445_strict_cd8_seurat.rds results/reviewer2
Rscript scripts/revision/08_recompute_reviewer2_cxcr6_limma.R \
  data/processed/merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData \
  data/processed/GSE224445_strict_cd8_seurat.rds \
  results/reviewer2/requested_cxcr6_contrasts_limma_ebayes.csv
```

The cross-cohort limma contrasts use all genes shared by the two assays for variance moderation, but remain exploratory because platform and cohort are confounded. `TOX` is absent from the GSE224445 targeted panel and is not imputed.

GSE224445 strict CD8 cells required T-lineage and CD8 RNA/ADT evidence without conflicting CD4, B-cell, or myeloid signal; our cohort used the six marker-annotated CD8 clusters. Feature-plot RNA log-normalized expression was scaled per gene within each dataset to the 99th percentile (0-1), so colors represent relative within-dataset intensity.

After downloading the lightweight plotting metrics from Figshare, regenerate primary computational panels without a Seurat object:

```bash
Rscript scripts/reproducibility/01_plot_figure1_from_metrics.R \
  figshare/inputs/plotting_inputs results/figure_1
Rscript scripts/reproducibility/02_plot_figure2_from_metrics.R \
  figshare/inputs/plotting_inputs results/figure_2
```

The metrics exporter documents how those public inputs are derived from locally reconstructed Seurat objects. Large `.RData`/`.rds` objects and raw matrices are not stored in this repository.

## Corrected paired results

For Fig. 1G, clusters 0, 2, 4, 5, 6, 8, 9, and 11 have FDR < 0.05. For Fig. 2B, `CXCL13+CXCR6+ effector CD8+` is higher in graft, whereas `CCR2+ CD4+` and `Naive CD4+` are lower in graft at FDR < 0.05. `CXCR6+ effector CD8+` has raw p = 0.0331 but FDR = 0.0910. These estimates are based on four matched patients and should be interpreted as paired cohort evidence rather than population-level prevalence estimates.

## Citation

Please cite the associated Nature Communications article and the Figshare data and reproducibility record:

> Dai Y. *Kidney graft TCR clonotype analysis: primary scRNA/TCR workflow, source data, and revision analyses*. Figshare. 2026. https://doi.org/10.6084/m9.figshare.33001412
