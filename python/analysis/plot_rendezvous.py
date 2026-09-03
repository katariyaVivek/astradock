#!/usr/bin/env python3
"""M18 rendezvous plots: approach profiles + ECI truth run."""

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


def plot_sweep(df: pd.DataFrame, out: Path, title: str):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle(title, fontweight="bold")
    axs[0, 0].plot(df["time_s"] / 3600.0, df["range_m"], color="royalblue")
    axs[0, 0].set_title("Range to Target")
    axs[0, 0].set_xlabel("Time (h)")
    axs[0, 0].set_ylabel("Range (m)")
    axs[0, 0].grid(True)
    axs[0, 1].plot(df["time_s"] / 3600.0, df["closing_mps"], color="crimson")
    axs[0, 1].set_title("Closing Speed")
    axs[0, 1].set_xlabel("Time (h)")
    axs[0, 1].set_ylabel("m/s")
    axs[0, 1].grid(True)
    axs[1, 0].plot(df["time_s"] / 3600.0, df["accel_mps2"] * 1000.0, color="forestgreen")
    axs[1, 0].set_title("Commanded Acceleration")
    axs[1, 0].set_xlabel("Time (h)")
    axs[1, 0].set_ylabel("mm/s^2")
    axs[1, 0].grid(True)
    axs[1, 1].plot(df["time_s"] / 3600.0, df["leg"], color="darkorange")
    axs[1, 1].set_title("Active Leg")
    axs[1, 1].set_xlabel("Time (h)")
    axs[1, 1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_eci_comparison(nominal: pd.DataFrame, dropout: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M18 — ECI Truth Rendezvous: Dropout Robustness", fontweight="bold")
    axs[0].plot(nominal["time_s"] / 3600.0, nominal["range_m"], label="no dropout",
                color="royalblue")
    axs[0].plot(dropout["time_s"] / 3600.0, dropout["range_m"], label="60 s dropout",
                color="crimson", ls="--", alpha=0.8)
    axs[0].set_title("Range to Target")
    axs[0].set_xlabel("Time (h)")
    axs[0].set_ylabel("Range (m)")
    axs[0].grid(True)
    axs[0].legend()
    axs[1].semilogy(nominal["time_s"] / 3600.0,
                    (nominal["est_range_m"] - nominal["range_m"]).abs(), label="no dropout",
                    color="royalblue", alpha=0.7)
    axs[1].semilogy(dropout["time_s"] / 3600.0,
                    (dropout["est_range_m"] - dropout["range_m"]).abs(), label="dropout",
                    color="crimson", alpha=0.7)
    axs[1].set_title("Estimate-vs-Scored Range Gap")
    axs[1].set_xlabel("Time (h)")
    axs[1].set_ylabel("|est - scored| (m)")
    axs[1].grid(True, which="both")
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_sweep(pd.read_csv(data / "m18_nominal.csv"), figs / "m18_nominal_approach.png",
               "M18 — Nominal 6 km to 50 m Approach (CW Plant)")
    plot_eci_comparison(pd.read_csv(data / "m18_eci_truth.csv"),
                        pd.read_csv(data / "m18_dropout.csv"),
                        figs / "m18_eci_dropout.png")
    print("Wrote M18 figures")


if __name__ == "__main__":
    main()
