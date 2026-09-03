#!/usr/bin/env python3
"""M22 vision plots: approach errors + noise sensitivity."""

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


def plot_approach(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M22 — Closing Approach 15 m to 5 m (0.25 px noise)", fontweight="bold")
    axs[0].plot(df["range_m"], df["pos_err_m"], "o-", color="royalblue", ms=4)
    axs[0].set_title("Position Error vs Range")
    axs[0].set_xlabel("Range (m)")
    axs[0].set_ylabel("Error (m)")
    axs[0].grid(True)
    axs[1].plot(df["range_m"], df["rms_px"], "o-", color="forestgreen", ms=4)
    axs[1].set_title("RMS Reprojection")
    axs[1].set_xlabel("Range (m)")
    axs[1].set_ylabel("RMS (px)")
    axs[1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_sensitivity(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M22 — Pixel-Noise Sensitivity at 10 m", fontweight="bold")
    axs[0].plot(df["range_m"], df["pos_err_m"], "o-", color="crimson", ms=5)
    axs[0].set_title("Position Error vs Pixel Noise")
    axs[0].set_xlabel("Noise std (px)")
    axs[0].set_ylabel("Error (m)")
    axs[0].grid(True)
    axs[1].plot(df["range_m"], df["rms_px"], "o-", color="darkorange", ms=5)
    axs[1].set_title("RMS Reprojection vs Noise")
    axs[1].set_xlabel("Noise std (px)")
    axs[1].set_ylabel("RMS (px)")
    axs[1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_approach(pd.read_csv(data / "m22_approach.csv"), figs / "m22_approach.png")
    plot_sensitivity(pd.read_csv(data / "m22_sensitivity.csv"), figs / "m22_sensitivity.png")
    print("Wrote M22 figures")


if __name__ == "__main__":
    main()
