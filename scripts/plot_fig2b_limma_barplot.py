#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description="Plot paired Fig. 2B cluster frequencies and joint-limma FDR labels.")
    parser.add_argument("--source", default="tables/fig2b_latest_11cluster_sample_level_proportions.csv")
    parser.add_argument("--stats", default="tables/fig2b_limma_ebayes_logit_proportions.csv")
    parser.add_argument("--out", default="results/fig2b")
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    data = pd.read_csv(args.source)
    stats = pd.read_csv(args.stats).set_index("feature")

    cluster_order = [
        "Memory-like CD4+", "CXCL13+CXCR6+ effector CD8+", "CXCR6+ effector CD8+",
        "CCR2+ CD4+", "Naive/Memory-like CD8+", "Naive CD4+", "Treg",
        "CX3CR1+ effector CD8+", "ZNF683+TCF7hi CD8+", "CXCL13+ CD4+", "ISG-high CD8+",
    ]
    clusters = [x for x in cluster_order if x in set(data["tcluster"])]
    data["condition"] = pd.Categorical(data["condition"], ["PBMC", "Graft"])
    means = (
        data.groupby(["tcluster", "condition"], observed=False)["proportion_pct"]
        .mean().unstack().reindex(clusters)
    )

    fig, ax = plt.subplots(figsize=(14.5, 7.2), dpi=180)
    x = np.arange(len(clusters))
    width = 0.34
    colors = {"PBMC": "#222222", "Graft": "#c93434"}
    for j, condition in enumerate(["PBMC", "Graft"]):
        ax.bar(
            x + (j - 0.5) * width, means[condition].to_numpy(), width,
            color=colors[condition], alpha=0.82, edgecolor="black", linewidth=0.55,
            label=condition, zorder=2,
        )

    for i, cluster in enumerate(clusters):
        pairs = []
        for patient in sorted(data["patient"].unique()):
            values = {}
            for condition in ["PBMC", "Graft"]:
                q = data[
                    (data.tcluster == cluster)
                    & (data.patient == patient)
                    & (data.condition == condition)
                ]["proportion_pct"]
                values[condition] = float(q.iloc[0]) if len(q) else 0.0
            pairs.append((values["PBMC"], values["Graft"]))
            ax.plot(
                [x[i] - width / 2, x[i] + width / 2], pairs[-1],
                color="#777777", lw=0.75, alpha=0.7, zorder=3,
            )
            ax.scatter(x[i] - width / 2, pairs[-1][0], s=22, facecolor="white", edgecolor=colors["PBMC"], linewidth=0.8, zorder=4)
            ax.scatter(x[i] + width / 2, pairs[-1][1], s=22, facecolor="white", edgecolor=colors["Graft"], linewidth=0.8, zorder=4)

        fdr = float(stats.loc[cluster, "FDR"])
        if np.isfinite(fdr) and fdr < 0.05:
            star = "***" if fdr < 0.001 else "**" if fdr < 0.01 else "*"
            ymax = max(max(z) for z in pairs)
            ysig = ymax + max(0.65, ymax * 0.12)
            ax.plot(
                [x[i] - width / 2, x[i] - width / 2, x[i] + width / 2, x[i] + width / 2],
                [ysig - 0.18, ysig, ysig, ysig - 0.18], color="black", lw=0.8,
                clip_on=False, zorder=5,
            )
            ax.text(x[i], ysig + 0.04, star, ha="center", va="bottom", fontsize=10, fontweight="bold", zorder=6)

    wrapped = [
        "Memory-like\nCD4+", "CXCL13+CXCR6+\neffector CD8+", "CXCR6+\neffector CD8+",
        "CCR2+\nCD4+", "Naive/Memory-like\nCD8+", "Naive\nCD4+", "Treg",
        "CX3CR1+\neffector CD8+", "ZNF683+TCF7hi\nCD8+", "CXCL13+\nCD4+", "ISG-high\nCD8+",
    ]
    label_map = dict(zip(cluster_order, wrapped))
    ax.set_xticks(x)
    ax.set_xticklabels([label_map[c] for c in clusters], ha="center", fontsize=7.8, linespacing=1.35)
    ax.set_xlabel("T-cell cluster", labelpad=17)
    ax.set_ylabel("Frequency among total T cells (%)")
    ax.set_title("Fig. 2B  T-cell cluster frequencies", pad=18)
    ax.legend(frameon=False, ncol=2, loc="upper right")
    ax.spines[["top", "right"]].set_visible(False)
    ax.grid(axis="y", color="#dddddd", lw=0.5, zorder=0)
    fig.text(
        0.09, 0.025,
        "Bars: mean; points: matched patients; lines: paired values; stars: joint limma/eBayes BH FDR",
        fontsize=8, color="#444444",
    )
    fig.subplots_adjust(left=0.075, right=0.985, top=0.88, bottom=0.25)
    fig.savefig(out / "fig2b_limma_paired_barplot.pdf", bbox_inches="tight")
    fig.savefig(out / "fig2b_limma_paired_barplot.png", dpi=300, bbox_inches="tight")


if __name__ == "__main__":
    main()
