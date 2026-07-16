# PPT Bioinformatics/Biostatistics Issue Map

Updated: 2026-07-16

This file maps the revision requests in the PPT to the reproducible artifact that answers each request. Original manuscript figures remain in their original directories; revision panels are listed separately.

## Cell-cycle request: G2/M scoring

- Question: perform G2/M cell-cycle scoring using all CD8+ T cells or shared TCR clonotypes, optionally stratified by CXCR6 or TCF1 state.
- Current analysis choice: all 4,110 shared-TCR-clone CD8+ cells from `all_clone_filter_10_28.Rdata`, stratified by tissue and CXCR6 status.
- Method: Seurat `CellCycleScoring` using the updated Seurat S-phase and G2/M gene sets intersected with genes present in the RNA assay.
- Figure: `Reviewer3_Major2_cell_cycle_in_situ/figures/12_reviewer3_shared_clone_cd8_cell_cycle_summary.pdf`, panels H-I.
- Source tables: `reviewer3_shared_clone_cd8_cell_cycle_patient_summary.csv` and `reviewer3_shared_clone_cd8_cell_cycle_patient_level_stats.csv`.
- Reproducible script: `Reviewer3_Major2_cell_cycle_in_situ/scripts/06_reviewer3_shared_clone_cell_cycle.R`.

## Shared-clone feature plots / Fig. 2G geometry

- Question: generate marker feature plots on the shared TCR clonotypes with the original Fig. 2G UMAP geometry and PBMC/Graft split.
- Current source: `all_clone_filter_10_28.Rdata`, object `all_clone_test_filter` (4,110 cells).
- UMAP rule: use the saved `umap` reduction directly; do not rerun UMAP or reconstruct a 3,545-cell subset.
- Original reference plot: `/home/0.collaboration/20_sc_VDJ/merge_4_with_vdj/shared_clone_T_cells_feature_plot_gene_list0_split_by_tissue_final.pdf`.
- Revision output: `figures/13_reviewer3_shared_clone_cd8_cell_cycle_marker_featureplots.pdf`.
- Reproducible script: `scripts/rebuild_fig2g_original_umap.R`.

## Fig. 2B cluster-frequency analysis

- Definition: frequency of the 11 annotated T-cell clusters among total T cells, by matched graft/PBMC sample.
- Statistical method: limma linear model with empirical-Bayes variance moderation (`lmFit` + `eBayes`) on logit-transformed sample-level proportions, with patient and condition terms.
- Figure: `figures/fig2b_limma_paired_barplot.pdf`.
- Source tables: `fig2b_latest_11cluster_sample_level_proportions.csv` and `fig2b_limma_ebayes_logit_proportions.csv`.
- Script: `scripts/recompute_fig2b_from_latest_tcell_object.R`.

## GSE224445 external PBMC analysis

- Question: provide CXCR6 feature plots, CXCR6+ CD8+ percentages, exact datapoints, and requested p-values for TOT, STA, BPAR, PBMC rejection, and graft rejection.
- Outputs: `figures/gse224445_cxcr6_boxplot_reviewer2_summary.pdf`, `figures/gse224445_cxcr6_summary_plots.pdf`, and `figures/r2_cxcr6_frequency_revised.pdf`.
- Source tables: `prism_ready_cxcr6_frequency_datapoints.csv`, `requested_cxcr6_contrasts_exact.csv`, and `requested_feature_gene_availability.csv`.
- Integration method: sample-tag-aware preprocessing, per-tag DoubletFinder, Seurat RPCA integration, marker-based strict CD8 selection.

## HUMESS revision

- Question: split PBMC/Graft plots, exact patient-level values, and explanation of the metric.
- Figure: `figures/nm_humess_support_overview_v2.pdf`.
- Source metric table: `tables/panel_e_exact_pseudobulk_glycolysis_metric.tsv`.
- Interpretation: the plotted HUMESS metric is calculated from the specified pseudobulk glycolysis-feature comparison; exact patient values are retained in the source table and the calculation is described in the response text/caption.

## Integration and batch correction response

- Main eight-sample dataset: Seurat reciprocal-PCA integration across the matched graft/PBMC objects, followed by scaling, PCA, UMAP, and SNN clustering.
- External GSE224445 dataset: per-sample-tag preprocessing followed by 2,000 integration features, `FindIntegrationAnchors(reduction = "rpca", dims = 1:20)`, and `IntegrateData(dims = 1:20)`.
- Response text: `response_text/reviewer_single_cell_comments_Dai_response_final.md` and `response_text/Dai_new_analysis_methods_addendum.md`.

## Packaging rule

The GitHub package contains scripts, source tables, methods, response text, and provenance. The Figshare package contains the release figures, source tables, methods, response text, scripts, and checksums. Historical/reference figures are not overwritten.
