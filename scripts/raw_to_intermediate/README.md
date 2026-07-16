# Raw-to-intermediate workflow

This folder gives a cleaned route from deposited scRNA-seq/TCR-seq data to the principal intermediate objects used by the manuscript analyses.

## Inputs

Create a local `sample_manifest.csv` that is not committed to GitHub. Each row should provide these fields:

- `sample_id`: sample identifier
- `patient_id`: anonymized patient identifier
- `tissue`: `PBMC` or `Graft`
- `batch_label`: label used during integration
- `matrix_dir`: Cell Ranger filtered feature matrix directory
- `vdj_dir`: matching TCR `vdj_t` directory containing `filtered_contig_annotations.csv` and `clonotypes.csv`
- `existing_seurat_rdata`: optional prebuilt per-sample Seurat RData path
- `existing_seurat_object`: optional object name in that RData file

The manifest is a local input configuration and must remain outside version control.

## Outputs

- `merged_4_with_vdj_9_24_2025.RData`: merged all-cell/all-immune object with TCR metadata.
- `merged_4_with_vdj_T_cells_10_8_2025.RData`: extracted T-cell object.
- `merged_4_with_vdj_T_cells_filter_10_28_2025.RData`: filtered/annotated T-cell object.
- `all_clone_filter_10_28.Rdata`: shared graft/PBMC clonotype subset.
- `merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData`: CD8 subset used for revision metabolic analysis.

The date-stamped filenames match the analysis provenance. These large objects are not part of the default figshare source-data deposit.
