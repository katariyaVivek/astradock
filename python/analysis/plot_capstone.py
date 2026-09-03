#!/usr/bin/env python3
"""M25 capstone plots: mission timeline + Monte Carlo histogram."""

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
    cap = pd.read_csv(data / "m25_capstone.csv")
    mc = pd.read_csv(data / "m25_monte_carlo.csv")

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle("M25 — Capstone Mission Timeline (seed 42)", fontweight="bold")
    axs[0, 0].semilogy(cap["time_s"], cap["range_m"], color="royalblue")
    axs[0, 0].set_title("Range to Target")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Range (m)")
    axs[0, 0].grid(True, which="both")
    axs[0, 1].semilogy(cap["time_s"], cap["att_err_deg"], color="crimson")
    axs[0, 1].set_title("Attitude Error")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Error (deg)")
    axs[0, 1].grid(True, which="both")
    axs[1, 0].semilogy(cap["time_s"], cap["nis"], color="purple", alpha=0.8)
    axs[1, 0].axhline(14.4494, color="k", ls=":", label="95% gate")
    axs[1, 0].set_title("GNSS NIS (dropout 900-960 s reads 0)")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].grid(True, which="both")
    axs[1, 0].legend()
    axs[1, 1].hist(mc["final_range_m"], bins=20, color="forestgreen", edgecolor="white")
    axs[1, 1].set_title("Monte Carlo Final Range (100 runs)")
    axs[1, 1].set_xlabel("Final range (m)")
    axs[1, 1].set_ylabel("Runs")
    axs[1, 1].grid(True)
    fig.tight_layout()
    fig.savefig(figs / "m25_capstone.png")
    plt.close(fig)
    print("Wrote M25 figures")


if __name__ == "__main__":
    main()
