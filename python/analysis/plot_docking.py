#!/usr/bin/env python3
"""M19 docking plots: nominal approach + abort-case summary."""

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


def plot_approach(df: pd.DataFrame, out: Path, title: str):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle(title, fontweight="bold")
    axs[0, 0].plot(df["time_s"], df["axial_m"], color="royalblue")
    axs[0, 0].set_title("Axial Separation (positive = apart)")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Axial (m)")
    axs[0, 0].grid(True)
    axs[0, 1].plot(df["time_s"], df["lateral_m"], color="crimson")
    axs[0, 1].axhline(0.25, color="k", ls=":", label="capture limit")
    axs[0, 1].set_title("Lateral Error")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Lateral (m)")
    axs[0, 1].grid(True)
    axs[0, 1].legend()
    axs[1, 0].plot(df["time_s"], df["closing_mps"], color="forestgreen")
    axs[1, 0].axhline(0.1, color="k", ls=":", label="capture limit")
    axs[1, 0].set_title("Closing Speed")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].set_ylabel("m/s")
    axs[1, 0].grid(True)
    axs[1, 0].legend()
    axs[1, 1].plot(df["time_s"], df["contact_force_N"], color="darkorange")
    axs[1, 1].axhline(500.0, color="k", ls=":", label="crush limit")
    axs[1, 1].set_title("Contact Force")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].set_ylabel("Force (N)")
    axs[1, 1].grid(True)
    axs[1, 1].legend()
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_approach(pd.read_csv(data / "m19_nominal.csv"), figs / "m19_nominal_approach.png",
                  "M19 — Nominal Final Approach to Capture (50 m to latch)")
    plot_approach(pd.read_csv(data / "m19_dropout.csv"), figs / "m19_dropout.png",
                  "M19 — Dropout Ride-Through (60 s GNSS outage)")
    print("Wrote M19 figures")


if __name__ == "__main__":
    main()
