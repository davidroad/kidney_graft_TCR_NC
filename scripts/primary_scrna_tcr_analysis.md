# Primary scRNA-seq and TCR-seq workflow

The workflow is organized into four script modules:

- `scripts/raw_to_intermediate/`: loads the single-cell matrices, TCR annotations, and sample manifest, then builds the integrated all-CD45, T-cell, shared-clonotype, and CD8 objects.
- `scripts/source_data_export/`: exports the source tables used by the main figures and generates the Fig. 2I Monocle trajectory.
- `scripts/analysis/`: runs the paired cluster-frequency analyses and processes GSE224445 for the Supplementary Fig. 6A CXCR6 analysis.
- `scripts/reproducibility/`: redraws article panels from the lightweight metrics released through Figshare.

The main object-building scripts are:

1. `01_build_primary_all_cell_object.R`
2. `02_extract_t_cells_and_shared_clonotypes.R`
3. `03_extract_cd8_object.R`

Raw scRNA-seq and TCR-seq data are available under `GSE319007`. The panel-level source data are available through the linked Figshare record.
