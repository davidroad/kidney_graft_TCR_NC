# Primary scRNA/TCR Analysis

The original exploratory analysis was developed as a local monolithic script. The public repository does not include that path-specific exploratory script because it contains local server paths and many unused plotting branches.

For reproducibility, the public code package instead provides:

- `scripts/raw_to_intermediate/01_build_primary_all_cell_object.R`
- `scripts/raw_to_intermediate/02_extract_t_cells_and_shared_clonotypes.R`
- `scripts/raw_to_intermediate/03_extract_cd8_revision_object.R`
- `scripts/source_data_export/03_generate_shared_clone_monocle_trajectory.R`

These scripts document the route from deposited scRNA-seq/TCR-seq data to the principal intermediate objects used by the primary human scRNA/TCR analysis:

- all-cell/all-immune Seurat object with TCR metadata
- extracted T-cell object
- filtered T-cell object
- shared graft/PBMC clonotype subset
- CD8 subset for downstream revision analyses

The Monocle 3 pseudotime panel uses the shared-clone intermediate object
`all_clone_filter_10_28.Rdata`, which contains the Seurat object
`all_clone_test_filter` with 4,110 shared-clone T cells. The Monocle 3 code
converts the RNA count matrix and metadata to a `cell_data_set`, inserts the
existing Seurat UMAP coordinates for visualization consistency, runs
`cluster_cells`, `learn_graph`, and `order_cells`, and plots cells colored by
pseudotime and annotated T-cell state. The pseudotime start/root point was
selected from the TCF7-defined stem-like T-cell state through the interactive
Monocle 3 `order_cells()` interface. Pseudotime then progressed toward
differentiated effector states. The intermediate Monocle `cell_data_set`
object was not retained as a deposited object because it can be regenerated
from this intermediate and code.

Primary raw single-cell and TCR-seq data are deposited in GEO under `GSE319007`. Panel-level source data for the primary scRNA/TCR figures are provided in the figshare source-data package under `source_data_by_analysis/primary_scrna_tcr_analysis/`.
