#!/usr/bin/env python3
"""M21 mission framework plots: scorecard + dispersion per mission."""

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
    "lines.linewidth": 1.5,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})


def plot_scorecard(scorecard: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M21 — Regression Mission Scorecard (20 seeded runs each)", fontweight="bold")
    axs[0].bar(scorecard["mission_id"], scorecard["success_rate"], color="forestgreen")
    axs[0].axhline(0.9, color="k", ls=":", label=" loosest bound (0.90)")
    axs[0].set_title("Success Rate")
    axs[0].set_ylabel("Fraction converged")
    axs[0].tick_params(axis="x", rotation=15)
    axs[0].grid(True)
    axs[0].legend()
    axs[1].bar(scorecard["mission_id"], scorecard["p95_range_m"], color="royalblue")
    axs[1].set_title("P95 Final Range")
    axs[1].set_ylabel("Range (m)")
    axs[1].tick_params(axis="x", rotation=15)
    axs[1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_dispersion(data_dir: Path, out: Path):
    fig, ax = plt.subplots(figsize=(10, 5), dpi=300)
    fig.suptitle("M21 — Final-Range Dispersion per Mission", fontweight="bold")
    frames = []
    for csv in sorted(data_dir.glob("m21_rendezvous_*.csv")) + sorted(data_dir.glob("m21_docking_*.csv")):
        df = pd.read_csv(csv, comment="#")
        frames.append((csv.stem, df["final_range_m"]))
    for i, (name, series) in enumerate(frames):
        ax.scatter([i] * len(series), series, s=12, alpha=0.7)
    ax.set_xticks(list(range(len(frames))))
    ax.set_xticklabels([name.replace("m21_", "") for name, _ in frames], rotation=15)
    ax.set_title("Per-run Final Range")
    ax.set_ylabel("Final range (m)")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_scorecard(pd.read_csv(data / "m21_scorecard.csv"), figs / "m21_scorecard.png")
    plot_dispersion(data, figs / "m21_dispersion.png")
    print("Wrote M21 figures")


if __name__ == "__main__":
    main()
