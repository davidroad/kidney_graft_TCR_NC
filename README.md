# Kidney Graft TCR Analysis

Finalized for Nature Communications: 2026-07-31

Repository: https://github.com/davidroad/kidney_graft_TCR_NC

This repository contains the analysis code for the human scRNA-seq and TCR-seq results. The scripts build processed objects, export source data, run panel-level analyses, and reproduce figures from released metrics.

Figure source data and standalone panel outputs are provided in the linked Figshare package.

## Repository contents

- `scripts/raw_to_intermediate/`: loads the scRNA-seq and TCR-seq inputs and builds the integrated all-CD45, T-cell, shared-clonotype, and CD8 objects.
- `scripts/source_data_export/`: exports primary figure source tables and the Fig. 2I Monocle trajectory.
- `scripts/reproducibility/`: exports lightweight metrics and plots Fig. 1, Fig. 2, and Supplementary Fig. 13 panels from Figshare data.
- `scripts/analysis/`: runs the Fig. 1G/Fig. 2B statistics and the GSE224445 PBMC analysis used in Supplementary Fig. 6A.
- `DATA_AVAILABILITY.md`: accessions and expected processed inputs.
- `SHA256SUMS.txt`: SHA-256 checksums for file-integrity verification.

## Figure and script map

| Figure | Main analysis | Script location |
|---|---|---|
| Fig. 1 | All-CD45 atlas, sample and tissue distributions, marker summaries, and cluster proportions | `scripts/source_data_export/02_export_all_cd45_main_figure_source_tables.R`; `scripts/reproducibility/01_plot_figure1_from_metrics.R`; `scripts/analysis/recompute_fig1g_limma.R` |
| Fig. 2 | T-cell atlas, shared graft-PBMC clonotypes, expression patterns, differential expression, and trajectory analysis | `scripts/source_data_export/01_export_tcell_main_figure_source_tables.R`; `scripts/source_data_export/03_generate_shared_clone_monocle_trajectory.R`; `scripts/reproducibility/02_plot_figure2_from_metrics.R`; `scripts/analysis/recompute_fig2b_from_latest_tcell_object.R` |
| Supplementary Fig. 6A | GSE224445 PBMC processing, CD8 T-cell selection, and CXCR6 visualization | `scripts/analysis/03_process_external_pbmc_dataset.R`; `scripts/analysis/04_refine_external_cd8_annotation.R`; `scripts/analysis/05_generate_suppfig6a_cxcr6_umaps.R` |
| Supplementary Fig. 13 | CD8 T-cell glycolysis-score visualization and gene-level paired results | `scripts/reproducibility/03_plot_suppfig13_from_metrics.R` |

## Custom analyses

- **Object construction and integration:** `scripts/raw_to_intermediate/` loads the single-cell and TCR inputs and builds the all-CD45, T-cell, shared-clonotype, and CD8 objects.
- **Matched cluster-frequency analysis:** `scripts/analysis/recompute_fig1g_limma.R` and `scripts/analysis/recompute_fig2b_from_latest_tcell_object.R` generate the Fig. 1G and Fig. 2B statistical tables.
- **External PBMC analysis:** `scripts/analysis/03_process_external_pbmc_dataset.R` and `04_refine_external_cd8_annotation.R` process and annotate GSE224445; `05_generate_suppfig6a_cxcr6_umaps.R` generates the CXCR6 visualization inputs.
- **Shared-clonotype trajectory:** `scripts/source_data_export/03_generate_shared_clone_monocle_trajectory.R` generates the Fig. 2I Monocle trajectory.
- **Panel reconstruction from released data:** `scripts/reproducibility/` redraws the supported panels using the lightweight Figshare metrics.

## Reproduce in article order

Download and extract the Figshare source-data package (version 3) next to this repository as `figshare/`:

https://doi.org/10.6084/m9.figshare.33131777.v3

### Software environment

Create the software environment:

```bash
conda env create -f environment.yml
conda activate kidney-graft-tcr-nc
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
```

### Fig. 1

```bash
Rscript scripts/reproducibility/01_plot_figure1_from_metrics.R \
  figshare/inputs/plotting_inputs results/figure_1
Rscript scripts/analysis/recompute_fig1g_limma.R \
  --input figshare/tables/figure_1/fig1g_current_cluster_proportions_by_patient.csv \
  --out results/fig1g
```

### Fig. 2

```bash
Rscript scripts/reproducibility/02_plot_figure2_from_metrics.R \
  figshare/inputs/plotting_inputs results/figure_2
Rscript scripts/analysis/recompute_fig2b_from_latest_tcell_object.R \
  --table-input figshare/tables/figure_2/fig2b_latest_11cluster_sample_level_proportions.csv \
  --out results/fig2b
```

Fig. 2I uses the shared-clonotype intermediate object and the documented Monocle 3 workflow:

```bash
Rscript scripts/source_data_export/03_generate_shared_clone_monocle_trajectory.R \
  data/processed/all_clone_filter_10_28.Rdata results/figure_2i
```

### Supplementary Fig. 6A

Generate the CXCR6 UMAP portion after reconstructing or supplying the two processed CD8 objects:

```bash
Rscript scripts/analysis/05_generate_suppfig6a_cxcr6_umaps.R \
  data/processed/merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData \
  data/processed/GSE224445_strict_cd8_seurat.rds results/suppfig6a
```

### Supplementary Fig. 13A-C

```bash
Rscript scripts/reproducibility/03_plot_suppfig13_from_metrics.R \
  figshare/inputs/plotting_inputs figshare/tables results/supplementary_figure_13
```

## Data access

- GSE319007: matched scRNA-seq and TCR-seq.
- GSE319298: CUT&Tag.
- GSE224445: external PBMC scRNA-seq.
- Figshare source-data record: https://doi.org/10.6084/m9.figshare.33131777.v3
