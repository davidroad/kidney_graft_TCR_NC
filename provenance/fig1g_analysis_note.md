# Corrected Fig. 1G analysis

- Biological unit: four matched patients, each with PBMC and graft samples.
- Features: 28 CD45+ clusters.
- Response: sample-level cluster proportion, bounded at 0.005/0.995 and logit transformed.
- Model: one 28-by-8 matrix fitted with `limma::lmFit`; design `~ patient + condition`.
- Variance moderation: one `eBayes(trend=TRUE, robust=FALSE)` call across 28 rows.
- Multiplicity: Benjamini-Hochberg correction across 28 contrasts.

Clusters 0, 2, 4, 5, 6, 8, 9, and 11 meet FDR < 0.05. The revised figure
shows all four patient-level values and the matched PBMC/graft connections.
