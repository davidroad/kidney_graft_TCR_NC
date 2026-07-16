# Corrected Fig. 2B analysis

- Biological unit: four matched patients, each with PBMC and graft samples.
- Features: 11 annotated T-cell clusters.
- Response: sample-level cluster proportion, bounded at 0.005/0.995 and logit transformed.
- Model: one 11-by-8 matrix fitted with `limma::lmFit`; design `~ patient + condition`.
- Variance moderation: one `eBayes(trend=TRUE, robust=FALSE)` call across 11 rows.
- Multiplicity: Benjamini-Hochberg correction across 11 contrasts.

Three clusters meet FDR < 0.05: `CXCL13+CXCR6+ effector CD8+` is higher in
graft, whereas `CCR2+ CD4+` and `Naive CD4+` are lower in graft.
`CXCR6+ effector CD8+` has raw p = 0.0331 and FDR = 0.0910. Given the four-pair
cohort and compositional outcome, these are paired cohort estimates.
