# Main manuscript pipeline

The main human scRNA-seq/TCR-seq analysis is rooted in the canonical
`merge_8_all_R2.R` program. It covers the eight paired PBMC/graft samples,
Cell Ranger-derived Seurat inputs, doublet filtering, integration, all-CD45
clustering, T-cell re-clustering, TCR clonotype annotation, shared-clonotype
selection, and the downstream objects used for the primary figures.

The monolithic program is retained for provenance and for reproducing the
historical analysis sequence. Its server-specific library loading and absolute
paths have been removed; configure `SCVDJ_RAW_ROOT` and `SCVDJ_WORK_ROOT` for a
local environment. It still expects the date-stamped per-tissue intermediate
workspaces used by the original sequence. The split programs below provide the
cleaner route when rebuilding directly from a sample manifest and GEO-derived
inputs.

## Split workflow

1. `raw_to_intermediate/01_build_all_cell_object_from_raw_10x.R` reads the
   eight-sample manifest, loads filtered 10x matrices or existing per-sample
   Seurat objects, attaches VDJ metadata, integrates the samples, and writes
   `merged_4_with_vdj_9_24_2025.RData`.
2. `raw_to_intermediate/02_extract_t_cells_and_shared_clone_objects.R`
   extracts the T-cell object, applies the documented cluster filter, and
   writes the filtered T-cell and shared-clonotype objects.
3. `raw_to_intermediate/03_extract_cd8_object_for_revision.R` derives the CD8
   object used by the revision metabolic/HUMESS analyses.
4. `source_data_export/01_export_tcell_main_figure_source_tables.R` and
   `source_data_export/02_export_all_cd45_main_figure_source_tables.R` export
   neutralized source tables for the primary figures.
5. `source_data_export/03_generate_shared_clone_monocle_trajectory.R` rebuilds
   the shared-clonotype Monocle 3 trajectory used for the main text.

The shared-clonotype object contains 4,110 T cells across CD8, CD4, and Treg
annotations. CD8-focused figures use the appropriate annotated CD8 states;
the population definition is not silently relabeled as CD8-only.

## Inputs and outputs

Use `raw_to_intermediate/sample_manifest_template.csv` to create a local
manifest. Required raw inputs are documented in `DATA_AVAILABILITY.md` and include
the filtered feature-barcode matrix and `vdj_t` directory for each of the four
matched graft/PBMC pairs. Raw FASTQ/10x data, RData objects, and generated
figures are intentionally excluded from GitHub.

The revision scripts in the parent `scripts/` directory are downstream of this
main pipeline. They are not substitutes for the all-CD45/T-cell integration
and shared-TCR workflow documented here.
