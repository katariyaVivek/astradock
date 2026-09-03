#!/usr/bin/env python3
"""M15 relative-dynamics plots: CW breakdown scaling + 5 km trajectory overlay."""

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


def plot_breakdown(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M15 — CW Linearization Breakdown vs Separation Scale", fontweight="bold")
    axs[0].loglog(df["separation_m"], df["cw_error_1orbit_m"], "o-", color="royalblue", label="1 orbit")
    axs[0].loglog(df["separation_m"], df["cw_error_3orbits_m"], "s-", color="crimson", label="3 orbits")
    axs[0].set_title("Absolute CW Position Error")
    axs[0].set_xlabel("Initial separation (m)")
    axs[0].set_ylabel("Error (m)")
    axs[0].grid(True, which="both")
    axs[0].legend()
    rel1 = df["cw_error_1orbit_m"] / df["separation_m"] * 100.0
    rel3 = df["cw_error_3orbits_m"] / df["separation_m"] * 100.0
    axs[1].loglog(df["separation_m"], rel1, "o-", color="royalblue", label="1 orbit")
    axs[1].loglog(df["separation_m"], rel3, "s-", color="crimson", label="3 orbits")
    axs[1].axhline(1.0, color="k", ls=":", label="1% threshold")
    axs[1].set_title("Relative CW Position Error")
    axs[1].set_xlabel("Initial separation (m)")
    axs[1].set_ylabel("Error / separation (%)")
    axs[1].grid(True, which="both")
    axs[1].legend()
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_trajectory(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 2, figsize=(12, 5), dpi=300)
    fig.suptitle("M15 — 5 km Bounded Ellipse: Nonlinear Truth vs CW Prediction", fontweight="bold")
    axs[0].plot(df["truth_y_m"], df["truth_x_m"], color="royalblue", label="ECI truth (LVLH)")
    axs[0].plot(df["pred_y_m"], df["pred_x_m"], "--", color="crimson", label="CW prediction")
    axs[0].set_title("In-Plane Relative Motion (x radial vs y along-track)")
    axs[0].set_xlabel("y along-track (m)")
    axs[0].set_ylabel("x radial (m)")
    axs[0].grid(True)
    axs[0].legend()
    axs[0].set_aspect("equal", adjustable="datalim")
    axs[1].semilogy(df["time_s"] / 3600.0, df["err_norm_m"], color="darkorange")
    axs[1].set_title("CW Position Error Growth")
    axs[1].set_xlabel("Time (hours)")
    axs[1].set_ylabel("Error norm (m)")
    axs[1].grid(True, which="both")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_breakdown(pd.read_csv(data / "m15_cw_breakdown.csv"), figs / "m15_cw_breakdown.png")
    plot_trajectory(
        pd.read_csv(data / "m15_cw_trajectory_5km.csv"), figs / "m15_cw_trajectory_5km.png"
    )
    print("Wrote M15 figures")


if __name__ == "__main__":
    main()
