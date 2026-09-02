"""AstraDock — Milestone M12: Sensor Simulation Visualizations.

Generates publication-grade figures illustrating simulated spacecraft sensor telemetry:
  1. m12_imu_telemetry.png: Gyroscope angular rate & Accelerometer specific force
  2. m12_gnss_telemetry.png: ECI position and velocity truth vs GNSS measurements
  3. m12_star_tracker_telemetry.png: Attitude orientation error angle evolution
  4. m12_range_telemetry.png: Truth vs measured target range and residuals
  5. m12_sensor_sampling_rates.png: Multi-rate sampling timeline & dropout window
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


def plot_imu(data_path: Path, output_dir: Path) -> None:
    print("Generating Figure: m12_imu_telemetry.png...")
    data = np.genfromtxt(data_path, delimiter=",", names=True)
    # Downsample for crisp plotting (every 50th sample)
    ds = slice(None, None, 50)
    t = data["time_s"][ds]
    valid = data["valid"][ds] == 1

    fig, axs = plt.subplots(2, 2, figsize=(13, 8), sharex=True)

    # 1. Truth vs Measured Angular Velocity
    axs[0, 0].plot(t[valid], data["truth_omega_x_rad_s"][ds][valid], label=r"Truth $\omega_x$", color="#1f77b4")
    axs[0, 0].plot(t[valid], data["meas_omega_x_rad_s"][ds][valid], label=r"Measured $\omega_x$", color="#ff7f0e", alpha=0.6, lw=0.8)
    axs[0, 0].axvspan(2000, 2100, color="gray", alpha=0.25, label="Dropout Window")
    axs[0, 0].set_ylabel(r"Angular Rate $\omega_x$ [rad/s]")
    axs[0, 0].set_title("Gyroscope Body Rate (X-Axis)")
    axs[0, 0].grid(True)
    axs[0, 0].legend(loc="upper right")

    axs[1, 0].plot(t[valid], data["truth_omega_z_rad_s"][ds][valid], label=r"Truth $\omega_z$", color="#2ca02c")
    axs[1, 0].plot(t[valid], data["meas_omega_z_rad_s"][ds][valid], label=r"Measured $\omega_z$", color="#d62728", alpha=0.6, lw=0.8)
    axs[1, 0].axvspan(2000, 2100, color="gray", alpha=0.25)
    axs[1, 0].set_xlabel("Time [s]")
    axs[1, 0].set_ylabel(r"Angular Rate $\omega_z$ [rad/s]")
    axs[1, 0].set_title("Gyroscope Body Rate (Z-Axis)")
    axs[1, 0].grid(True)
    axs[1, 0].legend(loc="upper right")

    # 2. Accelerometer Specific Force: Free Fall vs Measurement
    truth_spec_norm = np.linalg.norm(np.column_stack([
        data["truth_spec_force_x_mps2"][ds],
        data["truth_spec_force_y_mps2"][ds],
        data["truth_spec_force_z_mps2"][ds]
    ]), axis=1)

    meas_spec_norm = np.linalg.norm(np.column_stack([
        data["meas_spec_force_x_mps2"][ds],
        data["meas_spec_force_y_mps2"][ds],
        data["meas_spec_force_z_mps2"][ds]
    ]), axis=1)

    axs[0, 1].plot(t[valid], truth_spec_norm[valid] * 1e6, label=r"Truth Specific Force $f_{\text{true}}$ ($\mu$m/s$^2$)", color="#9467bd", lw=1.5)
    axs[0, 1].axvspan(2000, 2100, color="gray", alpha=0.25, label="Dropout Window")
    axs[0, 1].set_ylabel(r"Specific Force [$\mu$m/s$^2$]")
    axs[0, 1].set_title(r"Orbital Free-Fall Truth: $f = a - g \approx 0$ (Drag $\sim 0.9\,\mu$m/s$^2$)")
    axs[0, 1].grid(True)
    axs[0, 1].legend(loc="upper right")

    axs[1, 1].plot(t[valid], meas_spec_norm[valid] * 1e3, label=r"Measured Specific Force $f_{\text{meas}}$", color="#8c564b", alpha=0.7, lw=0.8)
    axs[1, 1].axvspan(2000, 2100, color="gray", alpha=0.25)
    axs[1, 1].set_xlabel("Time [s]")
    axs[1, 1].set_ylabel(r"Measured Specific Force [mm/s$^2$]")
    axs[1, 1].set_title(r"Accelerometer Output with Bias ($2.69$ mm/s$^2$) & Noise ($\sigma=5$ mm/s$^2$)")
    axs[1, 1].grid(True)
    axs[1, 1].legend(loc="upper right")

    plt.suptitle("M12A — 6-Axis IMU Simulation Telemetry (100 Hz, BODY Frame)", y=0.98)
    plt.tight_layout()
    fig.savefig(output_dir / "m12_imu_telemetry.png", dpi=300)
    plt.close(fig)


def plot_gnss(data_path: Path, output_dir: Path) -> None:
    print("Generating Figure: m12_gnss_telemetry.png...")
    data = np.genfromtxt(data_path, delimiter=",", names=True)
    t = data["time_s"]
    valid = data["valid"] == 1

    fig, axs = plt.subplots(2, 2, figsize=(13, 8))

    # Position Error
    pos_err_x = data["meas_r_x_m"][valid] - data["truth_r_x_m"][valid]
    pos_err_y = data["meas_r_y_m"][valid] - data["truth_r_y_m"][valid]
    pos_err_z = data["meas_r_z_m"][valid] - data["truth_r_z_m"][valid]

    axs[0, 0].plot(t[valid], pos_err_x, label="Residual X (Bias=2.5 m)", color="#1f77b4", alpha=0.7)
    axs[0, 0].plot(t[valid], pos_err_y, label="Residual Y (Bias=-1.5 m)", color="#ff7f0e", alpha=0.7)
    axs[0, 0].plot(t[valid], pos_err_z, label="Residual Z (Bias=3.0 m)", color="#2ca02c", alpha=0.7)
    axs[0, 0].axvspan(2000, 2100, color="gray", alpha=0.25, label="Dropout Window")
    axs[0, 0].set_xlabel("Time [s]")
    axs[0, 0].set_ylabel("Position Error [m]")
    axs[0, 0].set_title("GNSS Position Measurement Residuals (ECI Frame)")
    axs[0, 0].grid(True)
    axs[0, 0].legend(loc="upper right")

    # Position Error Histogram
    axs[0, 1].hist(pos_err_x, bins=35, density=True, alpha=0.5, color="#1f77b4", label="X Residuals")
    axs[0, 1].hist(pos_err_y, bins=35, density=True, alpha=0.5, color="#ff7f0e", label="Y Residuals")
    axs[0, 1].hist(pos_err_z, bins=35, density=True, alpha=0.5, color="#2ca02c", label="Z Residuals")
    axs[0, 1].set_xlabel("Position Error [m]")
    axs[0, 1].set_ylabel("Probability Density")
    axs[0, 1].set_title(r"Position Error Distribution ($\sigma = 3.0$ m)")
    axs[0, 1].grid(True)
    axs[0, 1].legend(loc="upper right")

    # Velocity Error
    vel_err_x = data["meas_v_x_mps"][valid] - data["truth_v_x_mps"][valid]
    vel_err_y = data["meas_v_y_mps"][valid] - data["truth_v_y_mps"][valid]
    vel_err_z = data["meas_v_z_mps"][valid] - data["truth_v_z_mps"][valid]

    axs[1, 0].plot(t[valid], vel_err_x * 1e2, label="Residual X (Bias=2.0 cm/s)", color="#1f77b4", alpha=0.7)
    axs[1, 0].plot(t[valid], vel_err_y * 1e2, label="Residual Y (Bias=-1.0 cm/s)", color="#ff7f0e", alpha=0.7)
    axs[1, 0].plot(t[valid], vel_err_z * 1e2, label="Residual Z (Bias=1.5 cm/s)", color="#2ca02c", alpha=0.7)
    axs[1, 0].axvspan(2000, 2100, color="gray", alpha=0.25)
    axs[1, 0].set_xlabel("Time [s]")
    axs[1, 0].set_ylabel("Velocity Error [cm/s]")
    axs[1, 0].set_title("GNSS Velocity Measurement Residuals (ECI Frame)")
    axs[1, 0].grid(True)
    axs[1, 0].legend(loc="upper right")

    # Velocity Error Histogram
    axs[1, 1].hist(vel_err_x * 1e2, bins=35, density=True, alpha=0.5, color="#1f77b4", label="X Residuals")
    axs[1, 1].hist(vel_err_y * 1e2, bins=35, density=True, alpha=0.5, color="#ff7f0e", label="Y Residuals")
    axs[1, 1].hist(vel_err_z * 1e2, bins=35, density=True, alpha=0.5, color="#2ca02c", label="Z Residuals")
    axs[1, 1].set_xlabel("Velocity Error [cm/s]")
    axs[1, 1].set_ylabel("Probability Density")
    axs[1, 1].set_title(r"Velocity Error Distribution ($\sigma = 3.0$ cm/s)")
    axs[1, 1].grid(True)
    axs[1, 1].legend(loc="upper right")

    plt.suptitle("M12B — Spaceborne GNSS Receiver Telemetry (1 Hz, ECI Frame)", y=0.98)
    plt.tight_layout()
    fig.savefig(output_dir / "m12_gnss_telemetry.png", dpi=300)
    plt.close(fig)


def plot_star_tracker(data_path: Path, output_dir: Path) -> None:
    print("Generating Figure: m12_star_tracker_telemetry.png...")
    data = np.genfromtxt(data_path, delimiter=",", names=True)
    ds = slice(None, None, 5)
    t = data["time_s"][ds]
    valid = data["valid"][ds] == 1

    err_rad = data["orientation_error_rad"][ds][valid]
    err_arcsec = np.degrees(err_rad) * 3600.0

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))

    axs[0].plot(t[valid], err_arcsec, color="#1f77b4", lw=1.0, alpha=0.8, label="Attitude Error Angle")
    axs[0].axhline(np.mean(err_arcsec), color="red", linestyle="--", label=f"Mean Error = {np.mean(err_arcsec):.1f}\"")
    axs[0].axvspan(2000, 2100, color="gray", alpha=0.25, label="Dropout Window")
    axs[0].set_xlabel("Time [s]")
    axs[0].set_ylabel("Orientation Error Angle [arcseconds]")
    axs[0].set_title("Star Tracker Attitude Geodesic Error vs Time")
    axs[0].grid(True)
    axs[0].legend(loc="upper right")

    axs[1].hist(err_arcsec, bins=40, density=True, color="#2ca02c", alpha=0.7, edgecolor="black")
    axs[1].axvline(np.mean(err_arcsec), color="red", linestyle="--", label=f"Mean = {np.mean(err_arcsec):.1f}\"")
    axs[1].set_xlabel("Orientation Error Angle [arcseconds]")
    axs[1].set_ylabel("Probability Density")
    axs[1].set_title(r"Attitude Error Distribution (Bias $\sim 103''$, $\sigma \sim 41''$)")
    axs[1].grid(True)
    axs[1].legend(loc="upper right")

    plt.suptitle("M12C — Optical Star Tracker Attitude Measurement (10 Hz, ECI←BODY)", y=0.98)
    plt.tight_layout()
    fig.savefig(output_dir / "m12_star_tracker_telemetry.png", dpi=300)
    plt.close(fig)


def plot_range(data_path: Path, output_dir: Path) -> None:
    print("Generating Figure: m12_range_telemetry.png...")
    data = np.genfromtxt(data_path, delimiter=",", names=True)
    ds = slice(None, None, 5)
    t = data["time_s"][ds]
    valid = data["valid"][ds] == 1

    truth_range = data["truth_range_m"][ds][valid]
    meas_range = data["meas_range_m"][ds][valid]
    range_err = meas_range - truth_range

    fig, axs = plt.subplots(1, 2, figsize=(13, 5))

    axs[0].plot(t[valid], truth_range, label="Truth Range", color="#1f77b4", lw=1.5)
    axs[0].plot(t[valid], meas_range, label="Measured Range (with Noise/Bias)", color="#ff7f0e", lw=0.8, alpha=0.7)
    axs[0].axvspan(2000, 2100, color="gray", alpha=0.25, label="Dropout Window")
    axs[0].set_xlabel("Time [s]")
    axs[0].set_ylabel("Range [m]")
    axs[0].set_title("Relative Range to Target Spacecraft")
    axs[0].grid(True)
    axs[0].legend(loc="upper right")

    axs[1].plot(t[valid], range_err, color="#2ca02c", lw=0.8, alpha=0.7, label="Residual (Meas - Truth)")
    axs[1].axhline(0.5, color="red", linestyle="--", label="Configured Bias = +0.50 m")
    axs[1].axvspan(2000, 2100, color="gray", alpha=0.25)
    axs[1].set_xlabel("Time [s]")
    axs[1].set_ylabel("Range Residual [m]")
    axs[1].set_title(r"Range Measurement Residuals ($\sigma = 0.10$ m)")
    axs[1].grid(True)
    axs[1].legend(loc="upper right")

    plt.suptitle("M12D — Line-of-Sight Relative Range Sensor Telemetry (10 Hz)", y=0.98)
    plt.tight_layout()
    fig.savefig(output_dir / "m12_range_telemetry.png", dpi=300)
    plt.close(fig)


def plot_sampling_rates(output_dir: Path) -> None:
    print("Generating Figure: m12_sensor_sampling_rates.png...")
    # Generate multi-rate timeline over a 5-second zoom window [10.0, 15.0]
    t_span = (10.0, 15.0)

    t_imu = np.arange(t_span[0], t_span[1] + 1e-6, 0.01)       # 100 Hz
    t_st = np.arange(t_span[0], t_span[1] + 1e-6, 0.1)         # 10 Hz
    t_range = np.arange(t_span[0], t_span[1] + 1e-6, 0.1)      # 10 Hz
    t_gnss = np.arange(t_span[0], t_span[1] + 1e-6, 1.0)       # 1 Hz

    fig, ax = plt.subplots(figsize=(12, 4.5))

    ax.scatter(t_imu, [4]*len(t_imu), s=15, color="#1f77b4", marker="|", label="IMU (100 Hz, BODY)")
    ax.scatter(t_st, [3]*len(t_st), s=40, color="#ff7f0e", marker="o", label="Star Tracker (10 Hz, SO(3))")
    ax.scatter(t_range, [2]*len(t_range), s=40, color="#2ca02c", marker="^", label="Range Sensor (10 Hz, Scalar)")
    ax.scatter(t_gnss, [1]*len(t_gnss), s=80, color="#d62728", marker="s", label="GNSS Receiver (1 Hz, ECI)")

    ax.set_yticks([1, 2, 3, 4])
    ax.set_yticklabels(["GNSS (1 Hz)", "Range (10 Hz)", "Star Tracker (10 Hz)", "IMU (100 Hz)"])
    ax.set_xlabel("Time [s]")
    ax.set_title("Multi-Rate Asynchronous Sensor Telemetry Timeline (5-Second Snapshot)")
    ax.set_xlim(9.8, 15.2)
    ax.grid(True, axis="x")
    ax.legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(output_dir / "m12_sensor_sampling_rates.png", dpi=300)
    plt.close(fig)


def main() -> None:
    setup_style()
    data_dir = Path("data")
    fig_dir = Path("artifacts/figures")
    fig_dir.mkdir(parents=True, exist_ok=True)

    plot_imu(data_dir / "m12_imu.csv", fig_dir)
    plot_gnss(data_dir / "m12_gnss.csv", fig_dir)
    plot_star_tracker(data_dir / "m12_star_tracker.csv", fig_dir)
    plot_range(data_dir / "m12_range.csv", fig_dir)
    plot_sampling_rates(fig_dir)

    print("All 5 sensor figures generated successfully in artifacts/figures/.")


if __name__ == "__main__":
    main()
