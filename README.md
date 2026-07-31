# Kidney Graft TCR Analysis

Finalized for Nature Communications: 2026-07-30

Repository: https://github.com/davidroad/kidney_graft_TCR_NC

This code package supports the human computational panels present in `Figures 07 30 2026.pdf`. It reconstructs the primary scRNA-seq and TCR-seq analysis route, exports lightweight source data, reproduces selected panels from released metrics, and records the paired cluster-frequency statistics.

The repository is code-only. Figure files and source-data tables are supplied in the linked Figshare package. Reviewer-only analyses and scripts whose outputs are absent from the 07-30 full figure set are not included.

## Repository contents

- `scripts/raw_to_intermediate/`: rebuilds the integrated all-CD45, T-cell, shared-clonotype, and CD8 objects.
- `scripts/source_data_export/`: exports primary figure source tables and the Fig. 2I Monocle trajectory.
- `scripts/reproducibility/`: exports lightweight metrics and plots Fig. 1, Fig. 2, and Supplementary Fig. 13 panels from Figshare data.
- `scripts/revision/`: final paired Fig. 1G/Fig. 2B statistics and the GSE224445 external PBMC route used in Supplementary Fig. 6A.
- `provenance/`: panel mapping, analysis notes, session information, and checksums.
- `DATA_AVAILABILITY.md`: accessions and expected processed inputs.

## Core analysis definitions

Fig. 1G contains 28 clusters across four matched graft-PBMC pairs. Fig. 2B contains 11 T-cell clusters across the same pairs. Zero-count cluster-sample combinations are retained. Each proportion is transformed as

`log((x + 0.5) / (N - x + 0.5))`,

where `x` is the cluster count and `N` is the total count in that sample. Each panel is fitted as one complete cluster-by-sample matrix using `limma`, design `~ patient + condition`, one empirical-Bayes moderation step, and Benjamini-Hochberg correction across all clusters in the panel.

The values printed in Fig. 1G and Fig. 2B are moderated two-sided P values. The released source-data tables contain both the P values and Benjamini-Hochberg FDR values.

The primary integration route uses Seurat reciprocal PCA anchors to align the eight samples while retaining patient and tissue metadata. Full parameters are documented in `scripts/raw_to_intermediate/`.

Supplementary Fig. 6A uses GSE224445 and GSE319007 CD8 T cells. Log-normalized CXCR6 expression is divided by the 99th percentile separately within each dataset for UMAP color display and limited to 0-1. The frequency graphs are not percentile-scaled and use detectable CXCR6 transcript expression.

Supplementary Fig. 13A-B uses the 17-gene glycolysis score. Supplementary Fig. 13C contains gene-level paired pseudobulk edgeR results; it is not a HUMESS activity panel.

## Running the code

Create the software environment:

```bash
conda env create -f environment.yml
conda activate kidney-graft-tcr-nc
Rscript -e 'remotes::install_github("chris-mcginnis-ucsf/DoubletFinder")'
```

Recompute final Fig. 1G and Fig. 2B statistical tables:

```bash
Rscript scripts/revision/recompute_fig1g_limma.R \
  --input figshare/tables/figure_1/fig1g_current_cluster_proportions_by_patient.csv \
  --out results/fig1g
Rscript scripts/revision/recompute_fig2b_from_latest_tcell_object.R \
  --table-input figshare/tables/figure_2/fig2b_latest_11cluster_sample_level_proportions.csv \
  --out results/fig2b
```

Generate only the CXCR6 UMAP portion of Supplementary Fig. 6A after reconstructing the two CD8 objects:

```bash
Rscript scripts/revision/05_generate_suppfig6a_cxcr6_umaps.R \
  data/processed/merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData \
  data/processed/GSE224445_strict_cd8_seurat.rds results/suppfig6a
```

After downloading the Figshare plotting inputs:

```bash
Rscript scripts/reproducibility/01_plot_figure1_from_metrics.R \
  figshare/inputs/plotting_inputs results/figure_1
Rscript scripts/reproducibility/02_plot_figure2_from_metrics.R \
  figshare/inputs/plotting_inputs results/figure_2
Rscript scripts/reproducibility/03_plot_suppfig13_from_metrics.R \
  figshare/inputs/plotting_inputs figshare/tables results/supplementary_figure_13
```

The complete panel map is `provenance/figure_panel_reproducibility_index.tsv`.

## Data access

- GSE319007: matched scRNA-seq and TCR-seq.
- GSE319298: CUT&Tag.
- GSE224445: external PBMC scRNA-seq.
- Figshare source-data record: https://doi.org/10.6084/m9.figshare.33001412
