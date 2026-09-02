"""AstraDock -- Milestone M13B: IMU-Aided Navigation Visualizations.

Generates publication-grade figures from the M13B telemetry CSV
(data/m13b_imu_ekf_telemetry.csv, produced by astradock_imu_ekf_demo):
  A. Position: truth vs GNSS vs IMU+GNSS EKF estimate
  B. Velocity: truth vs GNSS vs IMU+GNSS EKF estimate
  C. GNSS dropout: prediction-only interval highlighted with covariance growth
  D. IMU specific force: measured body-frame components over time
  E. Estimation error vs +-1-sigma covariance envelopes
  F. GNSS innovation and NIS consistency diagnostic
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


DROPOUT_START = 600.0
DROPOUT_END = 700.0


def mark_dropout(ax) -> None:
    ax.axvspan(DROPOUT_START, DROPOUT_END, color="gray", alpha=0.25)


def plot_position(data, output_dir: Path) -> None:
    print("Generating Figure: m13b_position.png ...")
    t = data["time_s"]
    valid = data["gnss_valid"] > 0.5
    fig, axs = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    for i, (ax, comp) in enumerate(zip(axs, ["x", "y", "z"])):
        tkey = f"truth_position_eci_{comp}_m"
        ekey = f"estimated_position_eci_{comp}_m"
        mkey = f"gnss_fix_{comp}_m" if f"gnss_fix_{comp}_m" in data.dtype.names else None
        ax.plot(t, data[tkey], label="Truth", color="#1f77b4", lw=1.0)
        if mkey:
            ax.plot(t[valid], data[mkey][valid], ".", color="#ff7f0e",
                    ms=2, alpha=0.35, label="GNSS fix")
        ax.plot(t, data[ekey], label="IMU+GNSS EKF", color="#2ca02c", lw=0.9)
        mark_dropout(ax)
        ax.set_ylabel(f"{comp} [m]")
        ax.grid(True)
        if i == 0:
            ax.legend(loc="upper left", ncol=3)
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("M13B - Position: Truth vs GNSS vs IMU-Aided EKF")
    fig.tight_layout()
    fig.savefig(output_dir / "m13b_position.png", dpi=150)
    plt.close(fig)


def plot_velocity(data, output_dir: Path) -> None:
    print("Generating Figure: m13b_velocity.png ...")
    t = data["time_s"]
    fig, axs = plt.subplots(3, 1, figsize=(13, 8), sharex=True)
    for i, (ax, comp) in enumerate(zip(axs, ["x", "y", "z"])):
        ax.plot(t, data[f"truth_velocity_eci_{comp}_mps"], label="Truth",
                color="#1f77b4", lw=1.0)
        ax.plot(t, data[f"estimated_velocity_eci_{comp}_mps"],
                label="IMU+GNSS EKF", color="#2ca02c", lw=0.9)
        mark_dropout(ax)
        ax.set_ylabel(f"v{comp} [m/s]")
        ax.grid(True)
        if i == 0:
            ax.legend(loc="upper left", ncol=2)
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("M13B - Velocity: Truth vs IMU-Aided EKF Estimate")
    fig.tight_layout()
    fig.savefig(output_dir / "m13b_velocity.png", dpi=150)
    plt.close(fig)


def plot_dropout(data, output_dir: Path) -> None:
    print("Generating Figure: m13b_dropout_covariance.png ...")
    t = data["time_s"]
    sigma_norm = np.sqrt(
        data["sigma_position_x_m"] ** 2
        + data["sigma_position_y_m"] ** 2
        + data["sigma_position_z_m"] ** 2)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True)

    ax1.plot(t, data["position_error_m"], color="#2ca02c", lw=0.8,
             label="Position error (dead reckoning through dropout)")
    ax1.axhline(np.sqrt(3) * 10.0, color="#ff7f0e", ls="--", lw=1,
                label=r"$\sqrt{3}\,\sigma_r$ raw GNSS floor")
    mark_dropout(ax1)
    ax1.set_ylabel("Position error [m]")
    ax1.set_title("M13B - GNSS Dropout [600 s, 700 s]: IMU Dead Reckoning")
    ax1.legend(loc="upper right")
    ax1.grid(True)

    ax2.semilogy(t, sigma_norm, color="#d62728", lw=1.0,
                 label=r"$|\boldsymbol{\sigma}_{pos}|$ filter uncertainty")
    ax2.axvspan(DROPOUT_START, DROPOUT_END, color="gray", alpha=0.25,
                label="GNSS unavailable (prediction only)")
    ax2.set_ylabel(r"$|\sigma_{pos}|$ [m]")
    ax2.set_xlabel("Time [s]")
    ax2.set_title("Covariance growth during outage and contraction after recovery")
    ax2.legend(loc="lower right")
    ax2.grid(True)

    fig.tight_layout()
    fig.savefig(output_dir / "m13b_dropout_covariance.png", dpi=150)
    plt.close(fig)


def plot_specific_force(data, output_dir: Path) -> None:
    print("Generating Figure: m13b_specific_force.png ...")
    t = data["time_s"]
    fig, axs = plt.subplots(3, 1, figsize=(13, 7), sharex=True)
    for ax, comp in zip(axs, ["x", "y", "z"]):
        key = f"imu_specific_force_body_{comp}_mps2"
        ax.plot(t, data[key], color="#9467bd", lw=0.5,
                label=f"measured f_{comp} (BODY)")
        ax.axhline(0.0, color="black", lw=0.5)
        mark_dropout(ax)
        ax.set_ylabel(f"f_{comp} [m/s$^2$]")
        ax.grid(True)
        ax.legend(loc="upper right")
    axs[-1].set_xlabel("Time [s]")
    fig.suptitle(
        "M13B - Measured Specific Force (free fall: |f| $\\approx$ noise; "
        "near-zero mean confirms f = a - g $\\approx$ 0 in orbit)")
    fig.tight_layout()
    fig.savefig(output_dir / "m13b_specific_force.png", dpi=150)
    plt.close(fig)


def plot_error_envelopes(data, output_dir: Path) -> None:
    print("Generating Figure: m13b_error_envelopes.png ...")
    t = data["time_s"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 7), sharex=True)

    pos_err = data["position_error_m"]
    sig_pos = np.sqrt(
        data["sigma_position_x_m"] ** 2
        + data["sigma_position_y_m"] ** 2
        + data["sigma_position_z_m"] ** 2)
    vel_err = data["velocity_error_mps"]
    sig_vel = np.sqrt(
        data["sigma_velocity_x_mps"] ** 2
        + data["sigma_velocity_y_mps"] ** 2
        + data["sigma_velocity_z_mps"] ** 2)

    ax1.plot(t, pos_err, color="#2ca02c", lw=0.7, label="|position error|")
    ax1.plot(t, sig_pos, color="#d62728", lw=0.9, label=r"$|\sigma_{pos}|$ ($\pm1\sigma$ envelope norm)")
    mark_dropout(ax1)
    ax1.set_ylabel("Position [m]")
    ax1.set_title("M13B - Estimation Error vs 1-Sigma Uncertainty (position)")
    ax1.legend(loc="upper right")
    ax1.grid(True)

    ax2.plot(t, vel_err, color="#2ca02c", lw=0.7, label="|velocity error|")
    ax2.plot(t, sig_vel, color="#d62728", lw=0.9, label=r"$|\sigma_{vel}|$")
    mark_dropout(ax2)
    ax2.set_ylabel("Velocity [m/s]")
    ax2.set_xlabel("Time [s]")
    ax2.set_title("Estimation Error vs 1-Sigma Uncertainty (velocity)")
    ax2.legend(loc="upper right")
    ax2.grid(True)

    fig.tight_layout()
    fig.savefig(output_dir / "m13b_error_envelopes.png", dpi=150)
    plt.close(fig)


def plot_innovations(data, output_dir: Path) -> None:
    print("Generating Figure: m13b_innovations.png ...")
    t = data["time_s"]
    valid = (data["gnss_valid"] > 0.5) & ~np.isnan(data["nis"])

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 6), sharex=True)

    ax1.plot(t[valid], data["innovation_norm_m"][valid], ".-", color="#1f77b4",
             lw=0.7, ms=3, label="||position innovation||")
    ax1.axhline(0.0, color="black", lw=0.5)
    mark_dropout(ax1)
    ax1.set_ylabel(r"$\|y_r\|$ [m]")
    ax1.set_title("M13B - GNSS Position Innovation Norm")
    ax1.legend(loc="upper right")
    ax1.grid(True)

    ax2.plot(t[valid], data["nis"][valid], ".-", color="#d62728", lw=0.7, ms=3)
    ax2.axhline(6.0, color="black", ls="--", lw=1, label=r"$\chi^2_6$ mean = 6")
    ax2.axhspan(1.64, 14.45, color="gray", alpha=0.15, label="95% band")
    mark_dropout(ax2)
    ax2.set_ylabel("NIS")
    ax2.set_xlabel("Time [s]")
    ax2.set_title("Normalized Innovation Squared (consistency diagnostic)")
    ax2.legend(loc="upper right")
    ax2.grid(True)

    fig.tight_layout()
    fig.savefig(output_dir / "m13b_innovations.png", dpi=150)
    plt.close(fig)


def main() -> None:
    setup_style()
    data_path = Path("data/m13b_imu_ekf_telemetry.csv")
    output_dir = Path("artifacts/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    if not data_path.exists():
        print("ERROR:", data_path, "not found. Run astradock_imu_ekf_demo first.")
        return

    data = np.genfromtxt(data_path, delimiter=",", names=True)
    print("Loaded", len(data), "rows from", data_path.name)

    plot_position(data, output_dir)
    plot_velocity(data, output_dir)
    plot_dropout(data, output_dir)
    plot_specific_force(data, output_dir)
    plot_error_envelopes(data, output_dir)
    plot_innovations(data, output_dir)

    print("All M13B figures written to", output_dir)


if __name__ == "__main__":
    main()
