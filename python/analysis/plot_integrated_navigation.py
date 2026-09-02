"""AstraDock M13D — Integrated 15-State Navigation Visualization Suite.

Generates Plots A through H illustrating:
  - Plot A: Position tracking (Truth vs EKF vs GNSS)
  - Plot B: Velocity tracking (Truth vs EKF vs GNSS)
  - Plot C: Attitude error and covariance envelope
  - Plot D: Accelerometer and gyro bias estimation
  - Plot E: Covariance coupling across translation and attitude
  - Plot F: GNSS outage covariance growth and rapid recovery
  - Plot G: Star tracker outage and downstream translational effect
  - Plot H: Sensor combinations RMSE comparison (incremental value of sensors)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def load_telemetry(filepath: str) -> np.ndarray:
    return np.genfromtxt(filepath, delimiter=",", names=True)


def plot_position(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    coords = ["x", "y", "z"]
    for i, ax in enumerate(axes):
        c = coords[i]
        err = data[f"truth_position_eci_{c}_m"] - data[f"estimated_position_eci_{c}_m"]

        ax.plot(t, err, color="#1f77b4", label=f"EKF Error {c.upper()}")
        ax.axhline(0, color="gray", linestyle="--", alpha=0.6)
        ax.set_ylabel(f"{c.upper()} Error [m]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right")

    axes[0].set_title("Plot A: Integrated 15-State EKF Position Estimation Errors in ECI")
    axes[-1].set_xlabel("Simulation Time [s]")
    plt.tight_layout()
    fig.savefig(out_dir / "m13d_position.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_position.png")


def plot_velocity(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 8), sharex=True)

    coords = ["x", "y", "z"]
    for i, ax in enumerate(axes):
        c = coords[i]
        err = data[f"truth_velocity_eci_{c}_mps"] - data[f"estimated_velocity_eci_{c}_mps"]

        ax.plot(t, err, color="#2ca02c", label=f"EKF Vel Error {c.upper()}")
        ax.axhline(0, color="gray", linestyle="--", alpha=0.6)
        ax.set_ylabel(f"v_{c.upper()} Error [m/s]")
        ax.grid(True, alpha=0.3)
        ax.legend(loc="upper right")

    axes[0].set_title("Plot B: Integrated 15-State EKF Velocity Estimation Errors in ECI")
    axes[-1].set_xlabel("Simulation Time [s]")
    plt.tight_layout()
    fig.savefig(out_dir / "m13d_velocity.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_velocity.png")


def plot_attitude(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    att_err_deg = np.degrees(data["attitude_error_rad"])
    sig_att_deg = np.degrees(data["sigma_attitude"])

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, att_err_deg, color="#d62728", label="Actual Attitude Error (MEKF)")
    ax.plot(t, sig_att_deg, color="black", linestyle="--", label="Estimated 1σ Uncertainty")
    ax.plot(t, 3.0 * sig_att_deg, color="gray", linestyle=":", label="Estimated 3σ Uncertainty")

    ax.axvspan(30.0, 45.0, color="orange", alpha=0.2, label="Star Tracker Outage (30-45 s)")
    ax.set_xlabel("Simulation Time [s]")
    ax.set_ylabel("Attitude Error [deg]")
    ax.set_title("Plot C: Integrated MEKF Attitude Error and Covariance Growth During Outage")
    ax.grid(True, alpha=0.3)
    ax.legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(out_dir / "m13d_attitude.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_attitude.png")


def plot_biases(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8), sharex=True)

    # Accelerometer bias
    for c, col in [("x", "#1f77b4"), ("y", "#2ca02c"), ("z", "#d62728")]:
        truth = data[f"truth_ba_{c}_mps2"] * 1e3
        est = data[f"estimated_ba_{c}_mps2"] * 1e3
        ax1.plot(t, est, color=col, label=f"Est b_a_{c}")
        ax1.plot(t, truth, color=col, linestyle="--", alpha=0.5, label=f"Truth b_a_{c}")

    ax1.set_ylabel("Accel Bias [mm/s²]")
    ax1.set_title("Plot D: Accelerometer and Gyroscope Bias Estimation")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper right", ncol=3)

    # Gyro bias
    for c, col in [("x", "#1f77b4"), ("y", "#2ca02c"), ("z", "#d62728")]:
        truth = data[f"truth_bg_{c}_rad_s"] * 1e3
        est = data[f"estimated_bg_{c}_rad_s"] * 1e3
        ax2.plot(t, est, color=col, label=f"Est b_g_{c}")
        ax2.plot(t, truth, color=col, linestyle="--", alpha=0.5, label=f"Truth b_g_{c}")

    ax2.set_xlabel("Simulation Time [s]")
    ax2.set_ylabel("Gyro Bias [mrad/s]")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right", ncol=3)

    plt.tight_layout()
    fig.savefig(out_dir / "m13d_biases.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_biases.png")


def plot_covariance_coupling(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)

    # Position
    pos_err = np.sqrt(
        (data["truth_position_eci_x_m"] - data["estimated_position_eci_x_m"]) ** 2 +
        (data["truth_position_eci_y_m"] - data["estimated_position_eci_y_m"]) ** 2 +
        (data["truth_position_eci_z_m"] - data["estimated_position_eci_z_m"]) ** 2
    )
    axes[0].plot(t, pos_err, color="#1f77b4", label="Actual Position Error")
    axes[0].plot(t, data["sigma_position"], color="black", linestyle="--", label="1σ Covariance")
    axes[0].set_ylabel("Position [m]")
    axes[0].set_title("Plot E: Cross-Coupled Estimation Errors vs 1σ Covariance Envelopes")
    axes[0].grid(True, alpha=0.3)
    axes[0].legend(loc="upper right")

    # Velocity
    vel_err = np.sqrt(
        (data["truth_velocity_eci_x_mps"] - data["estimated_velocity_eci_x_mps"]) ** 2 +
        (data["truth_velocity_eci_y_mps"] - data["estimated_velocity_eci_y_mps"]) ** 2 +
        (data["truth_velocity_eci_z_mps"] - data["estimated_velocity_eci_z_mps"]) ** 2
    )
    axes[1].plot(t, vel_err, color="#2ca02c", label="Actual Velocity Error")
    axes[1].plot(t, data["sigma_velocity"], color="black", linestyle="--", label="1σ Covariance")
    axes[1].set_ylabel("Velocity [m/s]")
    axes[1].grid(True, alpha=0.3)
    axes[1].legend(loc="upper right")

    # Attitude
    att_err_deg = np.degrees(data["attitude_error_rad"])
    sig_att_deg = np.degrees(data["sigma_attitude"])
    axes[2].plot(t, att_err_deg, color="#d62728", label="Actual Attitude Error")
    axes[2].plot(t, sig_att_deg, color="black", linestyle="--", label="1σ Covariance")
    axes[2].set_ylabel("Attitude [deg]")
    axes[2].set_xlabel("Simulation Time [s]")
    axes[2].grid(True, alpha=0.3)
    axes[2].legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(out_dir / "m13d_covariance_coupling.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_covariance_coupling.png")


def plot_gnss_outage(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    pos_err = np.sqrt(
        (data["truth_position_eci_x_m"] - data["estimated_position_eci_x_m"]) ** 2 +
        (data["truth_position_eci_y_m"] - data["estimated_position_eci_y_m"]) ** 2 +
        (data["truth_position_eci_z_m"] - data["estimated_position_eci_z_m"]) ** 2
    )
    ax1.plot(t, pos_err, color="#1f77b4", label="Position Error")
    ax1.plot(t, data["sigma_position"], color="black", linestyle="--", label="Position 1σ")
    ax1.axvspan(60.0, 90.0, color="red", alpha=0.2, label="GNSS Outage (60-90 s)")
    ax1.set_ylabel("Position [m]")
    ax1.set_title("Plot F: 30-Second GNSS Outage Dead-Reckoning Drift and Reacquisition")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper left")

    ax2.plot(t, data["gnss_valid"], color="green", drawstyle="steps-post", label="GNSS Status")
    ax2.plot(t, data["gnss_NIS"], color="#9467bd", alpha=0.7, label="GNSS NIS (df=6)")
    ax2.axhline(6.0, color="gray", linestyle="--", label="Expected NIS (6.0)")
    ax2.set_xlabel("Simulation Time [s]")
    ax2.set_ylabel("GNSS Metrics")
    ax2.set_ylim(-0.5, 15.0)
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(out_dir / "m13d_gnss_outage.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_gnss_outage.png")


def plot_star_tracker_outage(data: np.ndarray, out_dir: Path):
    t = data["time_s"]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True)

    att_err_deg = np.degrees(data["attitude_error_rad"])
    sig_att_deg = np.degrees(data["sigma_attitude"])

    ax1.plot(t, att_err_deg, color="#d62728", label="Attitude Error")
    ax1.plot(t, sig_att_deg, color="black", linestyle="--", label="Attitude 1σ")
    ax1.axvspan(30.0, 45.0, color="orange", alpha=0.2, label="Star Tracker Outage (30-45 s)")
    ax1.set_ylabel("Attitude [deg]")
    ax1.set_title("Plot G: 15-Second Star Tracker Outage and Downstream Uncertainty Inflation")
    ax1.grid(True, alpha=0.3)
    ax1.legend(loc="upper left")

    # Downstream translational acceleration uncertainty inflation
    ax2.plot(t, data["sigma_velocity"], color="#2ca02c", label="Velocity 1σ Uncertainty [m/s]")
    ax2.set_xlabel("Simulation Time [s]")
    ax2.set_ylabel("Velocity 1σ [m/s]")
    ax2.grid(True, alpha=0.3)
    ax2.legend(loc="upper right")

    plt.tight_layout()
    fig.savefig(out_dir / "m13d_star_tracker_outage.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_star_tracker_outage.png")


def plot_sensor_combinations(comb_path: Path, out_dir: Path):
    data = np.genfromtxt(comb_path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    configs = [d["configuration"] for d in data]
    pos_rmse = [d["position_rmse_m"] for d in data]
    att_rmse_deg = [np.degrees(d["attitude_rmse_rad"]) for d in data]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

    x = np.arange(len(configs))
    ax1.bar(x, pos_rmse, color="#1f77b4", width=0.5)
    ax1.set_xticks(x)
    ax1.set_xticklabels(configs, rotation=20, ha="right")
    ax1.set_ylabel("Position RMSE [m]")
    ax1.set_title("Position Accuracy Across Configurations")
    ax1.grid(True, alpha=0.3, axis="y")

    ax2.bar(x, att_rmse_deg, color="#d62728", width=0.5)
    ax2.set_xticks(x)
    ax2.set_xticklabels(configs, rotation=20, ha="right")
    ax2.set_ylabel("Attitude RMSE [deg]")
    ax2.set_title("Attitude Accuracy Across Configurations")
    ax2.set_yscale("log")
    ax2.grid(True, alpha=0.3, axis="y")

    plt.suptitle("Plot H: Incremental Navigation Performance Across Sensor Configurations", y=1.02)
    plt.tight_layout()
    fig.savefig(out_dir / "m13d_sensor_combinations.png", dpi=300)
    plt.close(fig)
    print("  Saved m13d_sensor_combinations.png")


def main():
    print("Generating M13D Analysis Figures...")
    out_dir = Path("artifacts/figures")
    out_dir.mkdir(parents=True, exist_ok=True)

    tel_path = Path("data/m13d_integrated_telemetry.csv")
    comb_path = Path("data/m13d_sensor_combinations.csv")

    if not tel_path.exists():
        raise FileNotFoundError(f"Telemetry not found at {tel_path}")

    data = load_telemetry(str(tel_path))
    plot_position(data, out_dir)
    plot_velocity(data, out_dir)
    plot_attitude(data, out_dir)
    plot_biases(data, out_dir)
    plot_covariance_coupling(data, out_dir)
    plot_gnss_outage(data, out_dir)
    plot_star_tracker_outage(data, out_dir)

    if comb_path.exists():
        plot_sensor_combinations(comb_path, out_dir)

    print("All M13D figures generated successfully.")


if __name__ == "__main__":
    main()
