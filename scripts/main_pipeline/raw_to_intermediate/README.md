# Raw-to-intermediate workflow

This folder gives a cleaned route from deposited scRNA-seq/TCR-seq data to the principal intermediate objects used by the manuscript analyses.

## Inputs

Prepare `sample_manifest.csv` from `sample_manifest_template.csv`. Each row should provide:

- sample identifier
- patient identifier
- tissue/source label
- Cell Ranger filtered feature matrix directory
- matching TCR `vdj_t` directory containing `filtered_contig_annotations.csv` and `clonotypes.csv`

If starting from a prebuilt per-sample Seurat object, fill `existing_seurat_rdata` and `existing_seurat_object`.

## Outputs

- `merged_4_with_vdj_9_24_2025.RData`: merged all-cell/all-immune object with TCR metadata.
- `merged_4_with_vdj_T_cells_10_8_2025.RData`: extracted T-cell object.
- `merged_4_with_vdj_T_cells_filter_10_28_2025.RData`: filtered/annotated T-cell object.
- `all_clone_filter_10_28.Rdata`: shared graft/PBMC clonotype subset.
- `merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData`: CD8 subset used for revision metabolic analysis.

The date-stamped filenames match the analysis provenance. These large objects are not part of the default figshare source-data deposit.

