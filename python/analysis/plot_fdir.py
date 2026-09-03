#!/usr/bin/env python3
"""M20 FDIR plots: bias-jump detection + dropout coast."""

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


def plot_case(df: pd.DataFrame, out: Path, title: str, fault_t: float):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle(title, fontweight="bold")
    axs[0, 0].semilogy(df["time_s"], df["gnss_nis"], color="royalblue")
    axs[0, 0].axhline(14.4494, color="k", ls=":", label="95% gate (df=6)")
    axs[0, 0].axvline(fault_t, color="crimson", ls="--", label="fault onset")
    axs[0, 0].set_title("GNSS NIS")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].grid(True, which="both")
    axs[0, 0].legend()
    axs[0, 1].plot(df["time_s"], df["pos_err_m"], color="forestgreen")
    axs[0, 1].axvline(fault_t, color="crimson", ls="--", label="fault onset")
    axs[0, 1].set_title("Position Error")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Error (m)")
    axs[0, 1].grid(True)
    axs[0, 1].legend()
    axs[1, 0].plot(df["time_s"], df["triggered"], color="darkorange")
    axs[1, 0].set_title("Monitor Trigger")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].grid(True)
    axs[1, 1].plot(df["time_s"], df["excluded"], color="purple")
    axs[1, 1].set_title("Channel Excluded (recovery)")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_case(pd.read_csv(data / "m20_bias_jump.csv"), figs / "m20_bias_jump.png",
              "M20 — 50 m GNSS Bias Jump: Detect, Isolate, Exclude", 60.0)
    plot_case(pd.read_csv(data / "m20_dropout.csv"), figs / "m20_dropout.png",
              "M20 — 60 s GNSS Dropout: Coast + Reacquire", 60.0)
    print("Wrote M20 figures")


if __name__ == "__main__":
    main()
