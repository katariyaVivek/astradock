#!/usr/bin/env python3
"""M16 control telemetry plots: attitude hold/slew + station keeping."""

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


def plot_attitude(df: pd.DataFrame, out: Path, title: str):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle(title, fontweight="bold")
    axs[0, 0].semilogy(df["time_s"], df["attitude_error_deg"], color="royalblue")
    axs[0, 0].set_title("Attitude Error")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Error (deg)")
    axs[0, 0].grid(True, which="both")
    axs[0, 1].semilogy(df["time_s"], df["rate_norm_deg_s"], color="crimson")
    axs[0, 1].set_title("Body Rate Norm")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Rate (deg/s)")
    axs[0, 1].grid(True, which="both")
    axs[1, 0].plot(df["time_s"], df["torque_norm_Nm"], color="forestgreen")
    axs[1, 0].set_title("Commanded Torque Norm (0.5 Nm cap)")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].set_ylabel("Torque (N m)")
    axs[1, 0].grid(True)
    axs[1, 1].plot(df["time_s"], df["saturated"], color="darkorange")
    axs[1, 1].set_title("Saturation Flag")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_station_keeping(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(1, 3, figsize=(13, 4.5), dpi=300)
    fig.suptitle("M16 — 100 m Radial Station Keeping (CW Plant)", fontweight="bold")
    axs[0].semilogy(df["time_s"], df["position_norm_m"], color="royalblue")
    axs[0].set_title("Relative Position Norm")
    axs[0].set_xlabel("Time (s)")
    axs[0].set_ylabel("Offset (m)")
    axs[0].grid(True, which="both")
    axs[1].semilogy(df["time_s"], df["velocity_norm_mps"], color="crimson")
    axs[1].set_title("Relative Velocity Norm")
    axs[1].set_xlabel("Time (s)")
    axs[1].set_ylabel("Speed (m/s)")
    axs[1].grid(True, which="both")
    axs[2].plot(df["time_s"], df["accel_norm_mps2"] * 1000.0, color="forestgreen")
    axs[2].set_title("Commanded Acceleration")
    axs[2].set_xlabel("Time (s)")
    axs[2].set_ylabel("Accel (mm/s^2)")
    axs[2].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_attitude(
        pd.read_csv(data / "m16_attitude_hold.csv"),
        figs / "m16_attitude_hold.png",
        "M16 — Detumble + Attitude Hold (tumble [0.2, -0.1, 0.15] rad/s)",
    )
    plot_attitude(
        pd.read_csv(data / "m16_attitude_slew.csv"),
        figs / "m16_attitude_slew.png",
        "M16 — 90-degree Yaw Slew",
    )
    plot_station_keeping(pd.read_csv(data / "m16_station_keeping.csv"), figs / "m16_station_keeping.png")
    print("Wrote M16 figures")


if __name__ == "__main__":
    main()
