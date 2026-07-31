# Fig. 2B analysis

- Biological unit: four matched patients, each with PBMC and graft samples.
- Features: 11 annotated T-cell clusters.
- Source matrix: 88 sample-cluster rows (8 samples × 11 clusters), including five explicit zero-count rows.
- Response: count-based empirical logit, `log[(x + 0.5)/(N - x + 0.5)]`, where `x` is the cluster count and `N` is the sample total.
- Model: one 11-by-8 matrix fitted with `limma::lmFit`; design `~ patient + condition`.
- Variance moderation: one `eBayes(trend=TRUE, robust=FALSE)` call across 11 rows.
- Multiplicity: Benjamini-Hochberg correction across 11 contrasts.

Four clusters meet FDR < 0.05: `CXCL13+CXCR6+ effector CD8+` and `CXCR6+ effector CD8+` are higher in graft, whereas `CCR2+ CD4+` and `Naive CD4+` are lower in graft. Given the four-pair cohort and compositional outcome, these are paired cohort estimates.
