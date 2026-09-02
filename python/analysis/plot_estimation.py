"""AstraDock — Milestone M13A: State Estimation Visualizations.

Generates publication-grade figures from the EKF telemetry CSV:
  A. Position: truth vs GNSS vs EKF estimate (ECI components)
  B. Position error: raw GNSS vs EKF
  C. Velocity: truth vs GNSS vs EKF
  D. Covariance: 1-sigma position and velocity uncertainty
  E. Innovations: GNSS position innovation + NIS time series
"""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def setup_style() -> None:
    plt.rcParams.update({
        "font.family": "serif",
        "font.size": 10,
        "axes.labelsize": 11,
        "axes.titlesize": 12,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "legend.fontsize": 9,
        "figure.titlesize": 13,
        "lines.linewidth": 1.2,
        "grid.alpha": 0.35,
        "grid.linestyle": "--",
    })


def plot_position(data, output_dir):
    print("Generating Figure: m13_position.png...")
    t = data["time_s"]
    fig, axs = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    comps = ["x", "y", "z"]
    tkeys = ["truth_pos_x_m", "truth_pos_y_m", "truth_pos_z_m"]
    ekeys = ["est_pos_x_m", "est_pos_y_m", "est_pos_z_m"]
    mkeys = ["meas_pos_x_m", "meas_pos_y_m", "meas_pos_z_m"]
    valid = data["meas_valid"] > 0.5
    for i, (ax, c, tk, ek, mk) in enumerate(zip(axs, comps, tkeys, ekeys, mkeys)):
        ax.plot(t, data[tk], label="Truth " + c, color="#1f77b4", lw=1.0)
        ax.plot(t[valid], data[mk][valid], label="GNSS " + c, color="#ff7f0e",
                alpha=0.3, lw=0.5, marker=".", markersize=2)
        ax.plot(t, data[ek], label="EKF " + c, color="#2ca02c", lw=0.9)
        ax.set_ylabel(c + " [m]")
        ax.grid(True)
        if i == 0:
            ax.legend(loc="upper right", ncol=3)
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("M13A - Position: Truth vs GNSS vs EKF Estimate")
    fig.tight_layout()
    fig.savefig(output_dir / "m13_position.png", dpi=150)
    plt.close(fig)


def plot_position_error(data, output_dir):
    print("Generating Figure: m13_position_error.png...")
    t = data["time_s"]
    fig, ax = plt.subplots(figsize=(13, 5))
    ax.plot(t, data["position_error_m"], label="EKF position error", color="#2ca02c")
    ax.axhline(np.sqrt(3) * 10.0, color="#ff7f0e", ls="--", lw=1,
               label="sqrt(3)*sigma_r (raw GNSS)")
    ax.axvspan(3000, 3900, color="gray", alpha=0.25, label="GNSS dropout")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Position error [m]")
    ax.set_title("M13A - Position Error: EKF vs Raw GNSS Noise Floor")
    ax.set_yscale("log")
    ax.set_ylim(bottom=0.1)
    ax.legend(loc="upper right")
    ax.grid(True)
    fig.tight_layout()
    fig.savefig(output_dir / "m13_position_error.png", dpi=150)
    plt.close(fig)


def plot_velocity(data, output_dir):
    print("Generating Figure: m13_velocity.png...")
    t = data["time_s"]
    fig, axs = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    comps = ["x", "y", "z"]
    tkeys = ["truth_vel_x_mps", "truth_vel_y_mps", "truth_vel_z_mps"]
    ekeys = ["est_vel_x_mps", "est_vel_y_mps", "est_vel_z_mps"]
    for i, (ax, c, tk, ek) in enumerate(zip(axs, comps, tkeys, ekeys)):
        ax.plot(t, data[tk], label="Truth " + c, color="#1f77b4", lw=1.0)
        ax.plot(t, data[ek], label="EKF " + c, color="#2ca02c", lw=0.9)
        ax.set_ylabel(c + " [m/s]")
        ax.grid(True)
        if i == 0:
            ax.legend(loc="upper right")
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("M13A - Velocity: Truth vs EKF Estimate")
    fig.tight_layout()
    fig.savefig(output_dir / "m13_velocity.png", dpi=150)
    plt.close(fig)


def plot_covariance(data, output_dir):
    print("Generating Figure: m13_covariance.png...")
    t = data["time_s"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6), sharex=True)
    ax1.plot(t, data["sigma_pos_max_m"], color="#d62728", label="Max sigma_pos")
    ax1.axhline(10.0, color="#ff7f0e", ls=":", alpha=0.7, label="GNSS sigma_r = 10 m")
    ax1.axvspan(3000, 3900, color="gray", alpha=0.25, label="GNSS dropout")
    ax1.set_ylabel("Position 1-sigma [m]")
    ax1.set_title("M13A - Position Uncertainty")
    ax1.legend(loc="upper left")
    ax1.grid(True)

    ax2.plot(t, data["sigma_vel_max_mps"], color="#9467bd", label="Max sigma_vel")
    ax2.axhline(0.05, color="#ff7f0e", ls=":", alpha=0.7, label="GNSS sigma_v = 0.05 m/s")
    ax2.axvspan(3000, 3900, color="gray", alpha=0.25)
    ax2.set_ylabel("Velocity 1-sigma [m/s]")
    ax2.set_xlabel("Time [s]")
    ax2.set_title("M13A - Velocity Uncertainty")
    ax2.legend(loc="upper left")
    ax2.grid(True)

    fig.tight_layout()
    fig.savefig(output_dir / "m13_covariance.png", dpi=150)
    plt.close(fig)


def plot_innovations(data, output_dir):
    print("Generating Figure: m13_innovations.png...")
    t = data["time_s"]
    valid = data["meas_valid"] > 0.5
    nis = data["nis"]
    nis_valid = nis[valid]
    t_valid = t[valid]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6), sharex=True)

    ikeys = ["innov_pos_x_m", "innov_pos_y_m", "innov_pos_z_m"]
    labels = ["x", "y", "z"]
    colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]
    for key, label, color in zip(ikeys, labels, colors):
        ax1.plot(t_valid, data[key][valid], label=label, color=color, alpha=0.6, lw=0.5)
    ax1.axhline(0, color="black", lw=0.5)
    ax1.axvspan(3000, 3900, color="gray", alpha=0.25)
    ax1.set_ylabel("Innovation [m]")
    ax1.set_title("M13A - GNSS Position Innovation")
    ax1.legend(loc="upper right", ncol=3)
    ax1.grid(True)

    ax2.plot(t_valid, nis_valid, color="#d62728", lw=0.7)
    ax2.axhline(6.0, color="black", ls="--", lw=1, label="chi2_6 mean = 6")
    ax2.axvspan(3000, 3900, color="gray", alpha=0.25, label="dropout")
    ax2.set_ylabel("NIS")
    ax2.set_xlabel("Time [s]")
    ax2.set_title("M13A - Normalized Innovation Squared (Consistency Diagnostic)")
    ax2.legend(loc="upper right")
    ax2.grid(True)

    fig.tight_layout()
    fig.savefig(output_dir / "m13_innovations.png", dpi=150)
    plt.close(fig)


def main():
    setup_style()
    data_path = Path("data/m13_ekf_telemetry.csv")
    output_dir = Path("artifacts/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        print("ERROR:", data_path, "not found. Run astradock_ekf_demo first.")
        return

    data = np.genfromtxt(data_path, delimiter=",", names=True)
    print("Loaded", len(data), "rows from", data_path.name)

    plot_position(data, output_dir)
    plot_position_error(data, output_dir)
    plot_velocity(data, output_dir)
    plot_covariance(data, output_dir)
    plot_innovations(data, output_dir)

    print("All estimation figures written to", output_dir)


if __name__ == "__main__":
    main()
