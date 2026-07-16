#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot paired Fig. 1G cluster frequencies and joint-limma FDR labels.")
    parser.add_argument("--source", default="tables/fig1g_current_cluster_proportions_by_patient.csv")
    parser.add_argument("--stats", default="tables/fig1g_limma_ebayes_logit_proportions.csv")
    parser.add_argument("--out", default="results/fig1g")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    d = pd.read_csv(args.source)
    stats = pd.read_csv(args.stats)
    feature_col = "integrated_snn_res.0.5"
    d["condition"] = d["condition"].astype(str).str.strip().str.lower().map({"pbmc": "PBMC", "graft": "Graft"})
    if d["condition"].isna().any():
        raise ValueError("Unexpected condition label in Fig. 1G source table")
    d[feature_col] = d[feature_col].astype(str)
    d["patient"] = d["patient"].astype(str)
    patients = sorted(d["patient"].unique())
    features = sorted(d[feature_col].unique(), key=lambda value: int(float(value)))

    grid = pd.MultiIndex.from_product(
        [patients, ["PBMC", "Graft"], features],
        names=["patient", "condition", feature_col],
    ).to_frame(index=False)
    d = grid.merge(d[["patient", "condition", feature_col, "proportion_pct"]],
                   on=["patient", "condition", feature_col], how="left")
    d["proportion_pct"] = d["proportion_pct"].fillna(0.0)
    fdr = stats.assign(feature=stats["feature"].astype(str)).set_index("feature")["FDR"].to_dict()

    colors = {"PBMC": "#3B8EA5", "Graft": "#D95F59"}
    x = np.arange(len(features), dtype=float)
    offset = {"PBMC": -0.19, "Graft": 0.19}
    width = 0.33
    fig, ax = plt.subplots(figsize=(16, 6.8))
    means = d.groupby([feature_col, "condition"])["proportion_pct"].mean()
    for condition in ["PBMC", "Graft"]:
        values = [means.loc[(feature, condition)] for feature in features]
        ax.bar(x + offset[condition], values, width=width, color=colors[condition],
               alpha=0.72, edgecolor="white", linewidth=0.5, label=condition, zorder=1)

    for i, feature in enumerate(features):
        sub = d[d[feature_col] == feature]
        for patient in patients:
            pair = sub[sub["patient"] == patient].set_index("condition")["proportion_pct"]
            values = [pair["PBMC"], pair["Graft"]]
            ax.plot([i + offset["PBMC"], i + offset["Graft"]], values,
                    color="#555555", alpha=0.68, linewidth=0.75, zorder=2)
            ax.scatter([i + offset["PBMC"], i + offset["Graft"]], values,
                       c=[colors["PBMC"], colors["Graft"]], s=14,
                       edgecolors="white", linewidths=0.35, zorder=3)

    ymax = max(1.0, d["proportion_pct"].max())
    for i, feature in enumerate(features):
        if fdr.get(feature, 1.0) < 0.05:
            local_max = d.loc[d[feature_col] == feature, "proportion_pct"].max()
            mark = "**" if fdr[feature] < 0.01 else "*"
            ax.text(i, local_max + 0.025 * ymax, mark, ha="center", va="bottom", fontsize=10)

    ax.set_xlim(-0.7, len(features) - 0.3)
    ax.set_ylim(0, ymax * 1.12)
    ax.set_xticks(x)
    ax.set_xticklabels(features, fontsize=8)
    ax.set_xlabel("CD45+ cluster", labelpad=11)
    ax.set_ylabel("Cells in cluster per sample (%)")
    ax.set_title("Fig. 1G  Paired cluster proportions", pad=12)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#E6E6E6", linewidth=0.65, zorder=0)
    fig.text(0.08, 0.015,
             "Bars: mean; points and lines: four matched patients; * limma BH FDR < 0.05; ** FDR < 0.01",
             fontsize=8, color="#444444")
    fig.tight_layout(rect=(0, 0.07, 1, 1))
    fig.savefig(out / "fig1g_limma_paired_barplot.pdf", bbox_inches="tight")
    fig.savefig(out / "fig1g_limma_paired_barplot.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
