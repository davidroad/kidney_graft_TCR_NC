# Shared-clonotype T-cell cycle analysis

- Input object: `all_clone_test_filter` from `all_clone_filter_10_28.Rdata`.
- Population: 4,110 Fig. 2D shared-clonotype T cells, including CD8, CD4, and Treg annotations.
- Biological unit: patient; four matched patients are reported as `P1-P4`.
- Stratification: tissue compartment and detectable CXCR6 RNA expression.
- Scores: Seurat S and G2/M scores plus specified proliferation, arrest, and survival/apoptosis modules.
- Inference: paired patient-level Wilcoxon tests; cell-level tests are exploratory only.
- Feature plots: MKI67 and CDKN2A shown on the saved original UMAP and split by PBMC versus graft.

No patient-level cell-cycle metric differs significantly between CXCR6-positive and CXCR6-negative shared-clonotype T cells after BH correction. The result supports low basal in situ cycling in the sampled shared-clonotype population and is not presented as a CD8-restricted analysis or as a replacement for the ex vivo proliferation assay.
