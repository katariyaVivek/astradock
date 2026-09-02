"""AstraDock -- Milestone M13C: Attitude Error-State & Range Update Visualizations.

Generates publication-quality figures from M13C telemetry datasets:
  Plot A: Attitude error: raw gyro vs star-tracker vs EKF estimate
  Plot B: Gyro bias: true vs estimated bias across X, Y, Z
  Plot C: Attitude covariance envelopes (+-1-sigma, +-3-sigma) with dropout highlighted
  Plot D: Attitude NIS consistency diagnostic vs chi^2(3) theoretical bounds
  Plot E: Range measurement update: truth, predicted, measured, and estimation error
  Plot F: Range NIS consistency diagnostic vs chi^2(1) theoretical bounds
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


def mark_attitude_dropout(ax) -> None:
    ax.axvspan(30.0, 45.0, color="#d62728", alpha=0.15, label="Star Tracker Outage (15 s)")


def mark_range_dropout(ax) -> None:
    ax.axvspan(20.0, 35.0, color="#d62728", alpha=0.15, label="Range Sensor Outage (15 s)")


def plot_attitude_error(data, output_dir: Path) -> None:
    t = data["time_s"]
    raw_err_deg = np.rad2deg(data["truth_attitude_error_rad"])
    ekf_err_deg = np.rad2deg(data["estimated_attitude_error_rad"])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, raw_err_deg, label="Open-Loop Gyro Integration", color="#d62728", lw=1.5, ls="--")
    ax.plot(t, ekf_err_deg, label="Attitude Error-State EKF (MEKF)", color="#2ca02c", lw=1.5)
    mark_attitude_dropout(ax)

    ax.set_yscale("log")
    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Attitude Error Norm [deg] (log scale)")
    ax.set_title("Plot A: Attitude Estimation Error — MEKF vs Open-Loop Gyro Integration")
    ax.grid(True, which="both")
    ax.legend(loc="upper left")

    fig.tight_layout()
    fig.savefig(output_dir / "m13c_attitude_error.png", dpi=300)
    plt.close(fig)
    print("Saved m13c_attitude_error.png")


def plot_gyro_bias(data, output_dir: Path) -> None:
    t = data["time_s"]
    true_bias_mrad_s = [5.0, -3.0, 2.0]
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    for i, (ax, comp, tb) in enumerate(zip(axs, ["x", "y", "z"], true_bias_mrad_s)):
        est_b_mrad = data[f"estimated_bias_{comp}"] * 1000.0
        sig_b_mrad = data[f"sigma_bias_{comp}"] * 1000.0

        ax.axhline(tb, color="#1f77b4", ls="--", lw=1.5, label=f"True Bias {comp.upper()} ({tb:.1f} mrad/s)")
        ax.plot(t, est_b_mrad, color="#ff7f0e", lw=1.2, label=f"Estimated Bias {comp.upper()}")
        ax.fill_between(t, est_b_mrad - sig_b_mrad, est_b_mrad + sig_b_mrad, color="#ff7f0e", alpha=0.2, label="+-1-sigma")
        mark_attitude_dropout(ax)

        ax.set_ylabel(f"b_{comp} [mrad/s]")
        ax.grid(True)
        if i == 0:
            ax.legend(loc="upper right", ncol=2)

    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("Plot B: Dynamic Gyroscope Bias Estimation Across Principal Axes")
    fig.tight_layout()
    fig.savefig(output_dir / "m13c_gyro_bias.png", dpi=300)
    plt.close(fig)
    print("Saved m13c_gyro_bias.png")


def plot_attitude_covariance(data, output_dir: Path) -> None:
    t = data["time_s"]
    fig, axs = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    for i, (ax, comp) in enumerate(zip(axs, ["x", "y", "z"])):
        # Approximate component error via 1/sqrt(3) scaling of total error for envelope comparison
        sig_1 = np.rad2deg(data[f"sigma_attitude_{comp}"])
        sig_3 = 3.0 * sig_1

        ax.fill_between(t, -sig_3, sig_3, color="#aec7e8", alpha=0.4, label="+-3-sigma Envelope")
        ax.fill_between(t, -sig_1, sig_1, color="#1f77b4", alpha=0.3, label="+-1-sigma Envelope")
        mark_attitude_dropout(ax)

        ax.set_ylabel(f"sigma_{comp} [deg]")
        ax.grid(True)
        if i == 0:
            ax.legend(loc="upper left")

    axs[-1].set_xlabel("Time [s]")
    fig.suptitle("Plot C: Attitude Error-State Covariance Envelopes & Sensor Outage Growth")
    fig.tight_layout()
    fig.savefig(output_dir / "m13c_attitude_covariance.png", dpi=300)
    plt.close(fig)
    print("Saved m13c_attitude_covariance.png")


def plot_attitude_nis(data, output_dir: Path) -> None:
    t = data["time_s"]
    nis = data["attitude_NIS"]
    valid = nis > 0.0

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t[valid], nis[valid], "o", color="#2ca02c", ms=3, alpha=0.6, label="Star Tracker NIS")

    # Chi-square(3) bounds: mean = 3.0, 95% = 7.815, 99% = 11.345
    ax.axhline(3.0, color="black", ls="-", lw=1.2, label="Theoretical Mean E[NIS] = 3.0")
    ax.axhline(7.815, color="#ff7f0e", ls="--", lw=1.2, label="95% Confidence Bound (7.815)")
    ax.axhline(11.345, color="#d62728", ls=":", lw=1.2, label="99% Confidence Bound (11.345)")
    mark_attitude_dropout(ax)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Normalized Innovation Squared (NIS)")
    ax.set_ylim(0, 16)
    ax.set_title("Plot D: Attitude Star Tracker Innovation Consistency Diagnostic (df = 3)")
    ax.grid(True)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output_dir / "m13c_attitude_nis.png", dpi=300)
    plt.close(fig)
    print("Saved m13c_attitude_nis.png")


def plot_range_update(data, output_dir: Path) -> None:
    t = data["time_s"]
    true_r = data["truth_range_m"]
    pred_r = data["predicted_range_m"]
    meas_r = data["measured_range_m"]
    pos_err = data["position_error_m"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    valid_meas = meas_r > 0.0
    ax1.plot(t, true_r, color="#1f77b4", lw=1.5, label="Truth Range")
    ax1.plot(t, pred_r, color="#2ca02c", lw=1.2, ls="--", label="Filter Predicted Range")
    ax1.plot(t[valid_meas], meas_r[valid_meas], ".", color="#ff7f0e", ms=3, alpha=0.5, label="Range Sensor Fix")
    mark_range_dropout(ax1)
    ax1.set_ylabel("Range [m]")
    ax1.set_title("Plot E: Nonlinear Range Measurement Tracking & Relative Approach")
    ax1.grid(True)
    ax1.legend(loc="upper right")

    ax2.plot(t, pos_err, color="#9467bd", lw=1.5, label="Chaser 3D Position Error Norm")
    mark_range_dropout(ax2)
    ax2.set_xlabel("Time [s]")
    ax2.set_ylabel("Position Error [m]")
    ax2.grid(True)
    ax2.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output_dir / "m13c_range_update.png", dpi=300)
    plt.close(fig)
    print("Saved m13c_range_update.png")


def plot_range_nis(data, output_dir: Path) -> None:
    t = data["time_s"]
    nis = data["range_NIS"]
    valid = nis > 0.0

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t[valid], nis[valid], "o", color="#1f77b4", ms=3, alpha=0.6, label="Scalar Range NIS")

    # Chi-square(1) bounds: mean = 1.0, 95% = 3.841, 99% = 6.635
    ax.axhline(1.0, color="black", ls="-", lw=1.2, label="Theoretical Mean E[NIS] = 1.0")
    ax.axhline(3.841, color="#ff7f0e", ls="--", lw=1.2, label="95% Confidence Bound (3.841)")
    ax.axhline(6.635, color="#d62728", ls=":", lw=1.2, label="99% Confidence Bound (6.635)")
    mark_range_dropout(ax)

    ax.set_xlabel("Time [s]")
    ax.set_ylabel("Range NIS")
    ax.set_ylim(0, 10)
    ax.set_title("Plot F: Range Measurement Innovation Consistency Diagnostic (df = 1)")
    ax.grid(True)
    ax.legend(loc="upper right")

    fig.tight_layout()
    fig.savefig(output_dir / "m13c_range_nis.png", dpi=300)
    plt.close(fig)
    print("Saved m13c_range_nis.png")


def main() -> None:
    setup_style()
    output_dir = Path("artifacts/figures")
    output_dir.mkdir(parents=True, exist_ok=True)

    att_csv = Path("data/m13c_attitude_telemetry.csv")
    range_csv = Path("data/m13c_range_telemetry.csv")

    if att_csv.exists():
        att_data = np.genfromtxt(att_csv, delimiter=",", names=True)
        plot_attitude_error(att_data, output_dir)
        plot_gyro_bias(att_data, output_dir)
        plot_attitude_covariance(att_data, output_dir)
        plot_attitude_nis(att_data, output_dir)
    else:
        print(f"Warning: {att_csv} not found.")

    if range_csv.exists():
        range_data = np.genfromtxt(range_csv, delimiter=",", names=True)
        plot_range_update(range_data, output_dir)
        plot_range_nis(range_data, output_dir)
    else:
        print(f"Warning: {range_csv} not found.")


if __name__ == "__main__":
    main()
