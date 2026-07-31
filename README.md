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

## Required inputs

- **Fig. 1:** the Figshare Fig. 1 plotting metrics; Fig. 1G additionally uses the released sample-level cluster table.
- **Fig. 2A-H:** the Figshare Fig. 2 plotting metrics; Fig. 2B additionally uses the released sample-level T-cell cluster table.
- **Fig. 2I:** the processed shared-clonotype object `all_clone_filter_10_28.Rdata`.
- **Supplementary Fig. 6A:** the processed GSE319007 CD8 object and the processed GSE224445 CD8 object.
- **Supplementary Fig. 13A-C:** the released glycolysis-score metrics, paired sample values, and gene-level result table from Figshare.
- **Upstream object construction:** staged GSE319007 single-cell matrices, matching TCR files, and a sample manifest are used by `scripts/raw_to_intermediate/`; downloaded GSE224445 files are used by `scripts/analysis/03_process_external_pbmc_dataset.R`.

Software dependencies are listed in `environment.yml`. The released metrics and tables are available at https://doi.org/10.6084/m9.figshare.33131777.v3.

## Data access

- GSE319007: matched scRNA-seq and TCR-seq.
- GSE319298: CUT&Tag.
- GSE224445: external PBMC scRNA-seq.
- Figshare source-data record: https://doi.org/10.6084/m9.figshare.33131777.v3
