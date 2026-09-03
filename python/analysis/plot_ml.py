#!/usr/bin/env python3
"""M23 ML figures: metric comparison + advisory agreement."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "legend.fontsize": 9,
    "figure.titlesize": 12,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(data / "m23_results.csv")

    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M23 — 10 m GNSS Bias: Classical Baselines vs ML Challenger", fontweight="bold")
    x = range(len(df))
    axs[0].bar(x, df["precision"], width=0.35, label="precision", color="royalblue")
    axs[0].bar([i + 0.35 for i in x], df["recall"], width=0.35, label="recall", color="crimson")
    axs[0].set_xticks([i + 0.175 for i in x])
    axs[0].set_xticklabels(df["model"], rotation=12)
    axs[0].set_title("Precision / Recall on Held-Out Seeds")
    axs[0].grid(True)
    axs[0].legend()
    axs[1].bar(x, df["false_alarm_rate"], width=0.35, label="false-alarm", color="darkorange")
    axs[1].bar([i + 0.35 for i in x], df["miss_rate"], width=0.35, label="miss", color="purple")
    axs[1].set_xticks([i + 0.175 for i in x])
    axs[1].set_xticklabels(df["model"], rotation=12)
    axs[1].set_title("Error Rates (accuracy NOT reported)")
    axs[1].grid(True)
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(figs / "m23_metrics.png")
    plt.close(fig)
    print("Wrote M23 figures")


if __name__ == "__main__":
    main()
