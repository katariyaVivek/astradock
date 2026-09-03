#!/usr/bin/env python3
"""M14 actuator telemetry plots: momentum buildup/saturation + thruster burn."""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 10, "axes.titlesize": 11,
    "axes.labelsize": 10, "legend.fontsize": 9, "figure.titlesize": 12,
    "lines.linewidth": 1.5, "grid.alpha": 0.3, "grid.linestyle": "--",
})


def plot_wheel_saturation(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle("M14 — Reaction Wheel Momentum Buildup & Saturation", fontweight="bold")
    axs[0, 0].plot(df["time_s"], df["commanded_torque_Nm"], label="commanded", color="gray", ls="--")
    axs[0, 0].plot(df["time_s"], df["achieved_torque_Nm"], label="achieved", color="royalblue")
    axs[0, 0].set_title("Motor Torque (command vs achieved)")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Torque (N m)")
    axs[0, 0].grid(True)
    axs[0, 0].legend()
    axs[0, 1].plot(df["time_s"], df["wheel_speed_rad_s"], color="crimson")
    axs[0, 1].axhline(100.0, color="k", ls=":", label="speed limit")
    axs[0, 1].set_title("Wheel Speed")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Speed (rad/s)")
    axs[0, 1].grid(True)
    axs[0, 1].legend()
    axs[1, 0].plot(df["time_s"], df["stored_momentum_Nms"], color="forestgreen")
    axs[1, 0].axhline(0.1, color="k", ls=":", label="capacity I*w_max")
    axs[1, 0].set_title("Stored Momentum")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].set_ylabel("H (N m s)")
    axs[1, 0].grid(True)
    axs[1, 0].legend()
    axs[1, 1].plot(df["time_s"], df["authority_lost"], label="authority lost", color="darkorange")
    axs[1, 1].plot(df["time_s"], df["speed_saturated"], label="speed saturated", color="purple", alpha=0.7)
    axs[1, 1].set_title("Saturation Flags")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].grid(True)
    axs[1, 1].legend()
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def plot_thruster_burn(df: pd.DataFrame, out: Path):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle("M14 — Offset Thruster Burn & Propellant Depletion", fontweight="bold")
    axs[0, 0].plot(df["time_s"], df["thrust_achieved_N"], color="royalblue")
    axs[0, 0].set_title("Achieved Thrust")
    axs[0, 0].set_xlabel("Time (s)")
    axs[0, 0].set_ylabel("Thrust (N)")
    axs[0, 0].grid(True)
    axs[0, 1].plot(df["time_s"], df["torque_body_z_Nm"], color="crimson")
    axs[0, 1].set_title("Body Torque (moment arm)")
    axs[0, 1].set_xlabel("Time (s)")
    axs[0, 1].set_ylabel("Torque z (N m)")
    axs[0, 1].grid(True)
    axs[1, 0].plot(df["time_s"], df["spacecraft_mass_kg"], color="forestgreen")
    axs[1, 0].set_title("Wet Mass")
    axs[1, 0].set_xlabel("Time (s)")
    axs[1, 0].set_ylabel("Mass (kg)")
    axs[1, 0].grid(True)
    axs[1, 1].plot(df["time_s"], df["propellant_used_kg"], color="darkorange")
    axs[1, 1].set_title("Cumulative Propellant")
    axs[1, 1].set_xlabel("Time (s)")
    axs[1, 1].set_ylabel("Used (kg)")
    axs[1, 1].grid(True)
    fig.tight_layout()
    fig.savefig(out)
    plt.close(fig)


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    figs = root / "artifacts" / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    plot_wheel_saturation(pd.read_csv(data / "m14_wheel_saturation.csv"), figs / "m14_wheel_saturation.png")
    plot_thruster_burn(pd.read_csv(data / "m14_thruster_burn.csv"), figs / "m14_thruster_burn.png")
    print("Wrote m14_wheel_saturation.png and m14_thruster_burn.png")


if __name__ == "__main__":
    main()
