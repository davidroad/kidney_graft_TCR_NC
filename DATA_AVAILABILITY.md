# Data availability

Raw scRNA-seq and TCR-seq data are associated with GEO accession `GSE319007`. CUT&Tag data are associated with `GSE319298`. The external PBMC analysis uses `GSE224445`.

Large processed Seurat objects and HUMESS model outputs are not committed to GitHub. They should be obtained from the associated Figshare record and placed in the following layout:

```text
data/processed/
  merged_4_with_vdj_T_cells_filter_10_28_2025.RData
  all_clone_filter_10_28.Rdata

results/reviewer2_metabolic_profile/
  input/cd8_metabolic_scores.tsv.gz
  input/pseudobulk_logcpm_by_sample.tsv.gz
  input/pseudobulk_counts_by_sample.tsv.gz
  edgeR/paired_graft_vs_pbmc_edgeR.tsv
  humess/humess_run/comparisons/PBMC__vs__Graft/ReporterMetabolites/
    reporter_metabolites_Graft_vs_PBMC_all.tsv
  humess/humess_run/models/
```

Expected RData object names:

- `merged_4_with_vdj_T_cells_filter_10_28_2025.RData`: `merged_all_T_batchcorrected_filter`
- `all_clone_filter_10_28.Rdata`: `all_clone_test_filter`

The Figshare DOI and article DOI must be inserted here before the repository is cited in the final publication.
