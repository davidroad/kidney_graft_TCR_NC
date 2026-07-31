# Data availability

Raw scRNA-seq and TCR-seq data are available under `GSE319007`. CUT&Tag data are available under `GSE319298`. The external PBMC analysis uses `GSE224445`.

Large processed objects are not committed to the code repository. The main scripts expect reconstructed objects such as:

```text
data/processed/
  merged_4_with_vdj_8_24_2025.RData
  merged_4_with_vdj_T_cells_filter_10_28_2025.RData
  all_clone_filter_10_28.Rdata
  merged_4_with_vdj_CD8_T_cells_filter_10_28_2025.RData
  GSE224445_strict_cd8_seurat.rds
```

Lightweight plotting metrics and statistical tables are supplied through Figshare:

https://doi.org/10.6084/m9.figshare.33135038.v1
