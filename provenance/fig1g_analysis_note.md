# Fig. 1G analysis

- Biological unit: four matched patients, each with PBMC and graft samples.
- Features: 28 CD45+ clusters.
- Source matrix: 224 sample-cluster rows (8 samples × 28 clusters), including 20 explicit zero-count rows.
- Response: count-based empirical logit, `log[(x + 0.5)/(N - x + 0.5)]`, where `x` is the cluster count and `N` is the sample total.
- Model: one 28-by-8 matrix fitted with `limma::lmFit`; design `~ patient + condition`.
- Variance moderation: one `eBayes(trend=TRUE, robust=FALSE)` call across 28 rows.
- Multiplicity: Benjamini-Hochberg correction across 28 contrasts.

Clusters 0, 2, 5, 6, 7, 8, 9, 11, 14, 16, and 22 meet FDR < 0.05.
