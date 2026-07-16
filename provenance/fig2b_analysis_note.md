# Corrected Fig. 2B analysis

## Analysis contract

- Biological unit: one matched patient; four patients contribute one PBMC and one graft sample each.
- Features: 11 annotated T-cell clusters used in Fig. 2B.
- Response: cluster count divided by the total annotated T-cell count within each sample.
- Transformation: empirical logit `log((n + 0.5) / (denominator - n + 0.5))`.
- Model: one 11-by-8 feature-by-sample matrix fitted with `limma::lmFit`; design `~ patient + condition`.
- Variance moderation: one `eBayes` call across all 11 clusters.
- Multiplicity: Benjamini-Hochberg correction across 11 clusters.

## Validation

The output contains 11 rows, four matched pairs per row, no missing statistics, and FDR values that differ from raw p values. Public sample identifiers are `P1-P4`.

## Result

Four clusters meet FDR < 0.05: `CXCL13+CXCR6+ effector CD8+` and `CXCR6+ effector CD8+` are higher in graft, whereas `CCR2+ CD4+` and `Naive CD4+` are lower in graft. Given the four-patient cohort and compositional outcome, these are paired cohort estimates and should not be generalized as population prevalence estimates.
