#!/usr/bin/env python3
"""M17 closed-loop telemetry plots: hold, maneuver, Monte Carlo."""

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


def plot_run(df: pd.DataFrame, out: Path, title: str):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle(title, fontweight="bold")
    axs[0, 0].semilogy(df["time_s"], df["scored_error_rad"] * 180.0 / 3.141592653589793,
                       label="simulated vs reference", color="royalblue")
    axs[0, 0].semilogy(df["time_s"], df["estimate_error_rad"] * 180.0 / 3.141592653589793,
                       label="simulated vs estimate", color="forestgreen", alpha=0.8)
    axs[0, 0].set_title("Attitude Errors")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Error (deg)")
    axs[0, 0].grid(True, which="both")
    axs[0, 0].legend()
    axs[0, 1].plot(df["time_s"], df["desired_torque_Nm"], label="desired", color="gray",
                   ls="--", alpha=0.7)
    axs[0, 1].plot(df["time_s"], df["achieved_torque_Nm"], label="achieved", color="crimson")
    axs[0, 1].set_title("Control Torque Norm")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Torque (N m)")
    axs[0, 1].grid(True)
    axs[0, 1].legend()
    axs[1, 0].plot(df["time_s"], df["wheel_momentum_total_Nms"], color="darkorange")
    axs[1, 0].set_title("Total Wheel Momentum")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].set_ylabel("H (N m s)")
    axs[1, 0].grid(True)
    axs[1, 1].semilogy(df["time_s"], df["nis"], color="purple", alpha=0.8)
    axs[1, 1].set_title("Star Tracker NIS (df=3)")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].grid(True, which="both")
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_monte_carlo(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5), dpi=300)
    fig.suptitle("M17 — 100-Run Monte Carlo: Closed-Loop Attitude Hold", fontweight="bold")
    axs[0].hist(df["settle_time_s"], bins=20, color="royalblue", edgecolor="white")
    axs[0].set_title("Settle Time Distribution")
    axs[0].set_xlabel("Settle time (s)")
    axs[0].set_ylabel("Runs")
    axs[0].grid(True)
    axs[1].semilogy(df["seed"], df["final_error_rad"] * 180.0 / 3.141592653589793,
                    "o", color="forestgreen", ms=3)
    axs[1].set_title("Final Error per Seed")
    axs[1].set_xlabel("Seed")
    axs[1].set_ylabel("Final error (deg)")
    axs[1].grid(True, which="both")
    axs[2].plot(df["seed"], df["saturated_fraction"], "o", color="darkorange", ms=3)
    axs[2].set_title("Saturated Fraction per Seed")
    axs[2].set_xlabel("Seed")
    axs[2].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_run(pd.read_csv(data / "m17_attitude_hold.csv"), figs / "m17_attitude_hold.png",
             "M17 — Closed-Loop Attitude Hold (estimates-only control)")
    plot_run(pd.read_csv(data / "m17_maneuver.csv"), figs / "m17_maneuver.png",
             "M17 — 60-deg Maneuver at t = 40 s (estimates-only control)")
    plot_monte_carlo(pd.read_csv(data / "m17_monte_carlo.csv"), figs / "m17_monte_carlo.png")
    print("Wrote M17 figures")


if __name__ == "__main__":
    main()
