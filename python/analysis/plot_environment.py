#!/usr/bin/env python3
"""Visualization and analysis scripts for AstraDock Milestone M11.

Generates publication-grade figures demonstrating environmental perturbation dynamics:
  1. Plot A: J2 Orbital Element Precession (m11_j2_precession.png)
  2. Plot B: Atmospheric Drag Orbital Decay & Energy Dissipation (m11_drag_decay.png)
  3. Plot C: Third-Body Lunar Gravitational Acceleration (m11_third_body.png)
  4. Plot D: Gravity-Gradient Torque & Pitch Libration (m11_gravity_gradient.png)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

# Styling configuration
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.size": 10,
    "axes.titlesize": 11,
    "axes.labelsize": 10,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "legend.fontsize": 9,
    "figure.titlesize": 12,
    "lines.linewidth": 1.5,
    "grid.alpha": 0.3,
    "grid.linestyle": "--",
})


def plot_j2_precession(df: pd.DataFrame, output_path: Path):
    """Plot A: J2 Orbital Precession showing RAAN regression, ArgP advancement, and element stability."""
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)
    fig.suptitle("M11 — Earth J2 Oblateness Perturbation & Orbital Element Precession", fontweight="bold")

    t_hr = df["time_s"] / 3600.0

    # 1. RAAN Regression
    axs[0, 0].plot(t_hr, df["raan_deg"], color="royalblue", label="Simulated RAAN")
    axs[0, 0].set_title("Right Ascension of Ascending Node (RAAN)")
    axs[0, 0].set_xlabel("Time (hours)")
    axs[0, 0].set_ylabel(r"$\Omega$ (deg)")
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    # 2. Argument of Periapsis Advancement
    axs[0, 1].plot(t_hr, df["argument_of_periapsis_deg"], color="forestgreen", label="Simulated ArgP")
    axs[0, 1].set_title("Argument of Periapsis Precession")
    axs[0, 1].set_xlabel("Time (hours)")
    axs[0, 1].set_ylabel(r"$\omega$ (deg)")
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    # 3. Semi-major Axis Oscillations (No secular drift)
    axs[1, 0].plot(t_hr, (df["semi_major_axis_m"] - df["semi_major_axis_m"].iloc[0]) / 1000.0, color="crimson")
    axs[1, 0].set_title("Semi-Major Axis Variations")
    axs[1, 0].set_xlabel("Time (hours)")
    axs[1, 0].set_ylabel(r"$\Delta a$ (km)")
    axs[1, 0].grid(True)

    # 4. Inclination Variations (No secular drift)
    axs[1, 1].plot(t_hr, df["inclination_deg"] - df["inclination_deg"].iloc[0], color="purple")
    axs[1, 1].set_title("Inclination Variations")
    axs[1, 1].set_xlabel("Time (hours)")
    axs[1, 1].set_ylabel(r"$\Delta i$ (deg)")
    axs[1, 1].grid(True)

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_drag_decay(df: pd.DataFrame, output_path: Path):
    """Plot B: Atmospheric drag altitude loss, semi-major axis shrinkage, and energy dissipation."""
    fig, axs = plt.subplots(3, 1, figsize=(10, 9), dpi=300, sharex=True)
    fig.suptitle("M11 — Atmospheric Drag Orbit Decay & Energy Dissipation (300 km LEO)", fontweight="bold")

    t_hr = df["time_s"] / 3600.0

    # 1. Altitude Decay
    axs[0].plot(t_hr, df["altitude_km"], color="darkorange", label="Geometric Altitude")
    axs[0].set_ylabel("Altitude (km)")
    axs[0].set_title("Spacecraft Altitude Decay")
    axs[0].grid(True)
    axs[0].legend()

    # 2. Semi-Major Axis Decay
    axs[1].plot(t_hr, df["semi_major_axis_km"], color="teal", label="Semi-Major Axis")
    axs[1].set_ylabel("SMA (km)")
    axs[1].set_title("Semi-Major Axis Shrinkage")
    axs[1].grid(True)
    axs[1].legend()

    # 3. Specific Orbital Energy Dissipation
    delta_energy = df["specific_energy_m2_s2"] - df["specific_energy_m2_s2"].iloc[0]
    axs[2].plot(t_hr, delta_energy, color="firebrick", label=r"$\Delta \varepsilon = \varepsilon(t) - \varepsilon_0$")
    axs[2].set_ylabel(r"$\Delta \varepsilon$ ($\mathrm{m}^2/\mathrm{s}^2$)")
    axs[2].set_xlabel("Time (hours)")
    axs[2].set_title("Mechanical Energy Dissipation")
    axs[2].grid(True)
    axs[2].legend()

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_third_body(df: pd.DataFrame, output_path: Path):
    """Plot C: Lunar third-body gravitational tidal acceleration along geostationary orbit."""
    fig, axs = plt.subplots(2, 1, figsize=(10, 7), dpi=300, sharex=True)
    fig.suptitle("M11 — Lunar Third-Body Gravitational Perturbation (GEO Orbit)", fontweight="bold")

    t_days = df["time_s"] / 86400.0

    # 1. Acceleration Components
    axs[0].plot(t_days, df["a_3b_x_mps2"] * 1e6, color="royalblue", label=r"$a_{3B,x}$ (ECI X)")
    axs[0].plot(t_days, df["a_3b_y_mps2"] * 1e6, color="forestgreen", label=r"$a_{3B,y}$ (ECI Y)")
    axs[0].plot(t_days, df["a_3b_z_mps2"] * 1e6, color="purple", label=r"$a_{3B,z}$ (ECI Z)")
    axs[0].set_ylabel(r"Acceleration ($\mu\mathrm{m}/\mathrm{s}^2$)")
    axs[0].set_title("Third-Body Tidal Acceleration Components")
    axs[0].grid(True)
    axs[0].legend(loc="upper right")

    # 2. Total Magnitude
    axs[1].plot(t_days, df["a_3b_mag_mps2"] * 1e6, color="crimson", label=r"$\|\mathbf{a}_{3B}\|$")
    axs[1].set_ylabel(r"Magnitude ($\mu\mathrm{m}/\mathrm{s}^2$)")
    axs[1].set_xlabel("Time (days)")
    axs[1].set_title("Total Lunar Gravitational Perturbation Magnitude")
    axs[1].grid(True)
    axs[1].legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {output_path}")


def plot_gravity_gradient(df: pd.DataFrame, output_path: Path):
    """Plot D: Gravity-gradient restoring torque and resulting attitude pitch libration."""
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), dpi=300, sharex=True)
    fig.suptitle("M11 — Gravity-Gradient Torque & Pitch Libration (Dumbbell Satellite)", fontweight="bold")

    t_min = df["time_s"] / 60.0

    # 1. Pitch Angle (Libration)
    axs[0].plot(t_min, df["pitch_deg"], color="mediumblue", label=r"Pitch Angle $\theta$")
    axs[0].set_ylabel(r"Pitch $\theta$ (deg)")
    axs[0].set_title("Attitude Pitch Angle Oscillation (Libration)")
    axs[0].grid(True)
    axs[0].legend(loc="upper right")

    # 2. Pitch Angular Velocity
    axs[1].plot(t_min, df["omega_y_rad_s"] * 1e3, color="darkcyan", label=r"$\omega_y$")
    axs[1].set_ylabel(r"$\omega_y$ ($\mathrm{mrad}/\mathrm{s}$)")
    axs[1].set_title("Pitch Angular Velocity")
    axs[1].grid(True)
    axs[1].legend(loc="upper right")

    # 3. Gravity-Gradient Restoring Torque
    axs[2].plot(t_min, df["torque_y_Nm"] * 1e6, color="crimson", label=r"$\tau_{gg,y}$")
    axs[2].set_ylabel(r"Torque ($\mu\mathrm{N}\cdot\mathrm{m}$)")
    axs[2].set_xlabel("Time (minutes)")
    axs[2].set_title("Restoring Gravity-Gradient Torque")
    axs[2].grid(True)
    axs[2].legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(output_path, dpi=300)
    plt.close(fig)
    print(f"Saved: {output_path}")


def main():
    root = Path(__file__).resolve().parent.parent.parent
    data_dir = root / "data"
    fig_dir = root / "artifacts" / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    j2_csv = data_dir / "environment_j2_orbit.csv"
    drag_csv = data_dir / "environment_drag_decay.csv"
    tb_csv = data_dir / "environment_third_body.csv"
    gg_csv = data_dir / "environment_gravity_gradient.csv"

    if j2_csv.exists():
        df_j2 = pd.read_csv(j2_csv)
        plot_j2_precession(df_j2, fig_dir / "m11_j2_precession.png")

    if drag_csv.exists():
        df_drag = pd.read_csv(drag_csv)
        plot_drag_decay(df_drag, fig_dir / "m11_drag_decay.png")

    if tb_csv.exists():
        df_tb = pd.read_csv(tb_csv)
        plot_third_body(df_tb, fig_dir / "m11_third_body.png")

    if gg_csv.exists():
        df_gg = pd.read_csv(gg_csv)
        plot_gravity_gradient(df_gg, fig_dir / "m11_gravity_gradient.png")

    print("\nAll M11 environmental figures generated successfully in artifacts/figures/")


if __name__ == "__main__":
    main()
