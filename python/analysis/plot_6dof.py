#!/usr/bin/env python3
"""Visualization and analysis scripts for AstraDock Milestone M10.

Generates publication-grade figures demonstrating integrated 6-DOF spacecraft simulation:
  1. Plot A: 3D ECI Orbit with Body Attitude Triads (m10_3d_orbit_attitude_triads.png)
  2. Plot B: Translational Telemetry & Orbital Conservation (m10_translational_telemetry.png)
  3. Plot C: Rotational Telemetry & Invariants (m10_rotational_telemetry.png)
  4. Plot D: Regression, Decoupling & Torque Verification (m10_regression_and_decoupling.png)
"""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
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


def quat_to_dcm(w: float, x: float, y: float, z: float) -> np.ndarray:
    """Converts a scalar-first unit quaternion into a 3x3 DCM (C_ECI_B)."""
    xx, yy, zz = x * x, y * y, z * z
    wx, wy, wz = w * x, w * y, w * z
    xy, xz, yz = x * y, x * z, y * z

    return np.array([
        [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz),       2.0 * (xz + wy)],
        [2.0 * (xy + wz),       1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
        [2.0 * (xz - wy),       2.0 * (yz + wx),       1.0 - 2.0 * (xx + yy)],
    ], dtype=np.float64)


def plot_3d_orbit_attitude_triads(df: pd.DataFrame, output_path: Path):
    """Plot A: 3D ECI Orbit Trajectory with Body Attitude Triads at discrete epochs."""
    fig = plt.figure(figsize=(10, 8), dpi=300)
    ax = fig.add_subplot(111, projection="3d")

    # Earth reference sphere
    r_earth_km = 6378.137
    u = np.linspace(0, 2 * np.pi, 50)
    v = np.linspace(0, np.pi, 50)
    xs = r_earth_km * np.outer(np.cos(u), np.sin(v))
    ys = r_earth_km * np.outer(np.sin(u), np.sin(v))
    zs = r_earth_km * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_wireframe(xs, ys, zs, color="steelblue", alpha=0.25, rstride=4, cstride=4)

    # Spacecraft trajectory
    pos_x = df["position_eci_x_m"].to_numpy() / 1000.0
    pos_y = df["position_eci_y_m"].to_numpy() / 1000.0
    pos_z = df["position_eci_z_m"].to_numpy() / 1000.0

    ax.plot(pos_x, pos_y, pos_z, color="#1f77b4", lw=2, label="500 km Orbit Trajectory")

    # Body frame triads at 5 key epochs
    total_time = df["time_s"].iloc[-1]
    target_fractions = [0.0, 0.25, 0.50, 0.75, 1.0]
    triad_len_km = 1200.0

    for frac in target_fractions:
        target_t = frac * total_time
        idx = int(np.argmin(np.abs(df["time_s"].to_numpy() - target_t)))

        rx, ry, rz = pos_x[idx], pos_y[idx], pos_z[idx]
        qw = df["quaternion_w"].iloc[idx]
        qx = df["quaternion_x"].iloc[idx]
        qy = df["quaternion_y"].iloc[idx]
        qz = df["quaternion_z"].iloc[idx]

        # DCM columns represent body unit vectors expressed in ECI
        C = quat_to_dcm(qw, qx, qy, qz)
        u_b = C[:, 0] * triad_len_km
        v_b = C[:, 1] * triad_len_km
        w_b = C[:, 2] * triad_len_km

        # Draw triad vectors: +X_B (red), +Y_B (green), +Z_B (blue)
        ax.quiver(rx, ry, rz, u_b[0], u_b[1], u_b[2], color="#d62728", lw=1.8, arrow_length_ratio=0.15)
        ax.quiver(rx, ry, rz, v_b[0], v_b[1], v_b[2], color="#2ca02c", lw=1.8, arrow_length_ratio=0.15)
        ax.quiver(rx, ry, rz, w_b[0], w_b[1], w_b[2], color="#1f77b4", lw=1.8, arrow_length_ratio=0.15)

        ax.scatter([rx], [ry], [rz], color="black", s=25, zorder=5)
        ax.text(rx, ry, rz + 300, f"t = {frac:.2f} T", fontsize=8, weight="bold")

    # Dummy lines for legend
    ax.plot([], [], color="#d62728", label=r"Body $+X_B$ (Principal Axis 1)")
    ax.plot([], [], color="#2ca02c", label=r"Body $+Y_B$ (Principal Axis 2)")
    ax.plot([], [], color="#1f77b4", label=r"Body $+Z_B$ (Principal Axis 3)")

    ax.set_xlabel("ECI X (km)")
    ax.set_ylabel("ECI Y (km)")
    ax.set_zlabel("ECI Z (km)")
    ax.set_title("AstraDock M10: 3D Integrated Orbit & Spacecraft Attitude Triads", pad=20, weight="bold")
    ax.legend(loc="upper right", framealpha=0.9)

    ax.set_box_aspect([1, 1, 0.8])
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Generated Plot A: {output_path}")


def plot_translational_telemetry(df: pd.DataFrame, output_path: Path):
    """Plot B: Translational position, velocity, and orbital invariants."""
    t_min = df["time_s"].to_numpy() / 60.0

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)

    # Top-Left: ECI Position
    axs[0, 0].plot(t_min, df["position_eci_x_m"] / 1000.0, label=r"$r_x$", color="#1f77b4")
    axs[0, 0].plot(t_min, df["position_eci_y_m"] / 1000.0, label=r"$r_y$", color="#ff7f0e")
    axs[0, 0].plot(t_min, df["position_eci_z_m"] / 1000.0, label=r"$r_z$", color="#2ca02c")
    axs[0, 0].set_ylabel("Position (km)")
    axs[0, 0].set_title("ECI Position Components")
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    # Top-Right: ECI Velocity
    axs[0, 1].plot(t_min, df["velocity_eci_x_mps"], label=r"$v_x$", color="#1f77b4")
    axs[0, 1].plot(t_min, df["velocity_eci_y_mps"], label=r"$v_y$", color="#ff7f0e")
    axs[0, 1].plot(t_min, df["velocity_eci_z_mps"], label=r"$v_z$", color="#2ca02c")
    axs[0, 1].set_ylabel("Velocity (m/s)")
    axs[0, 1].set_title("ECI Velocity Components")
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    # Bottom-Left: Specific Orbital Energy Drift
    init_energy = df["specific_orbital_energy_m2_s2"].iloc[0]
    energy_drift = np.maximum(np.abs(df["specific_orbital_energy_m2_s2"].to_numpy() - init_energy) / np.abs(init_energy), 1.0e-16)
    axs[1, 0].plot(t_min, energy_drift, color="#d62728", label=r"$|\Delta \mathcal{E}| / |\mathcal{E}_0|$")
    axs[1, 0].set_xlabel("Simulation Time (minutes)")
    axs[1, 0].set_ylabel("Relative Energy Drift")
    axs[1, 0].set_title("Orbital Specific Energy Conservation")
    axs[1, 0].set_yscale("log")
    axs[1, 0].set_ylim(bottom=1.0e-17, top=1.0e-12)
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    # Bottom-Right: Specific Angular Momentum Drift
    init_h = df["specific_angular_momentum_m2_s"].iloc[0]
    h_drift = np.maximum(np.abs(df["specific_angular_momentum_m2_s"].to_numpy() - init_h) / np.abs(init_h), 1.0e-16)
    axs[1, 1].plot(t_min, h_drift, color="#9467bd", label=r"$|\Delta h| / |h_0|$")
    axs[1, 1].set_xlabel("Simulation Time (minutes)")
    axs[1, 1].set_ylabel("Relative Momentum Drift")
    axs[1, 1].set_title("Orbital Angular Momentum Conservation")
    axs[1, 1].set_yscale("log")
    axs[1, 1].set_ylim(bottom=1.0e-17, top=1.0e-12)
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    fig.suptitle("AstraDock M10: Spacecraft Translational Dynamics & Orbital Invariants", fontsize=13, weight="bold")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Generated Plot B: {output_path}")


def plot_rotational_telemetry(df: pd.DataFrame, output_path: Path):
    """Plot C: Rotational quaternion, angular rates, kinetic energy, and quaternion norm."""
    t_s = df["time_s"].to_numpy()

    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)

    # Top-Left: Attitude Quaternion Components
    axs[0, 0].plot(t_s, df["quaternion_w"], label=r"$q_w$ (scalar)", color="#1f77b4")
    axs[0, 0].plot(t_s, df["quaternion_x"], label=r"$q_x$", color="#ff7f0e")
    axs[0, 0].plot(t_s, df["quaternion_y"], label=r"$q_y$", color="#2ca02c")
    axs[0, 0].plot(t_s, df["quaternion_z"], label=r"$q_z$", color="#d62728")
    axs[0, 0].set_ylabel("Quaternion Parameter")
    axs[0, 0].set_title("Attitude Quaternion Evolution")
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    # Top-Right: Body Angular Velocity Components
    axs[0, 1].plot(t_s, df["angular_velocity_body_x_rad_s"], label=r"$\omega_x$", color="#1f77b4")
    axs[0, 1].plot(t_s, df["angular_velocity_body_y_rad_s"], label=r"$\omega_y$", color="#ff7f0e")
    axs[0, 1].plot(t_s, df["angular_velocity_body_z_rad_s"], label=r"$\omega_z$", color="#2ca02c")
    axs[0, 1].set_ylabel("Angular Velocity (rad/s)")
    axs[0, 1].set_title("Body Angular Velocity Components (Polhode Motion)")
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    # Bottom-Left: Rotational Kinetic Energy Drift
    init_e_rot = df["rotational_kinetic_energy_J"].iloc[0]
    rot_drift = np.abs(df["rotational_kinetic_energy_J"].to_numpy() - init_e_rot) / init_e_rot
    axs[1, 0].plot(t_s, rot_drift, color="#e377c2", label=r"$|\Delta T_{rot}| / T_{rot,0}$")
    axs[1, 0].set_xlabel("Simulation Time (seconds)")
    axs[1, 0].set_ylabel("Relative Energy Drift")
    axs[1, 0].set_title("Rotational Kinetic Energy Conservation")
    axs[1, 0].set_yscale("log")
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    # Bottom-Right: Quaternion Norm Error
    norm_err = np.abs(df["quaternion_norm"].to_numpy() - 1.0)
    axs[1, 1].plot(t_s, norm_err, color="#8c564b", label=r"$|\|q\| - 1|$")
    axs[1, 1].set_xlabel("Simulation Time (seconds)")
    axs[1, 1].set_ylabel("Norm Error")
    axs[1, 1].set_title("Step-Boundary Unit Quaternion Normalization")
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    fig.suptitle("AstraDock M10: Spacecraft Rotational Dynamics & Invariants", fontsize=13, weight="bold")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Generated Plot C: {output_path}")


def plot_regression_and_decoupling(canon_df: pd.DataFrame, torque_df: pd.DataFrame, output_path: Path):
    """Plot D: Regression equivalence, decoupling, and external torque response."""
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), dpi=300)

    # Top-Left: Decoupled Orbital Path Comparison (Torque-Free vs Torque-Driven)
    canon_r = np.linalg.norm(canon_df[["position_eci_x_m", "position_eci_y_m", "position_eci_z_m"]].to_numpy(), axis=1) / 1000.0
    torque_r = np.linalg.norm(torque_df[["position_eci_x_m", "position_eci_y_m", "position_eci_z_m"]].to_numpy(), axis=1) / 1000.0

    axs[0, 0].plot(canon_df["time_s"].iloc[:len(torque_df)], canon_r[:len(torque_df)], label="Torque-Free (Scenario 1)", color="#1f77b4")
    axs[0, 0].plot(torque_df["time_s"], torque_r, "--", label=r"$\tau_z = 0.2$ N$\cdot$m (Scenario 2)", color="#d62728")
    axs[0, 0].set_ylabel("Orbital Radius (km)")
    axs[0, 0].set_title("Orbit Independence under External Body Torque")
    axs[0, 0].grid(True)
    axs[0, 0].legend()

    # Top-Right: External Torque Angular Velocity Ramp
    t_torque = torque_df["time_s"].to_numpy()
    analytical_wz = (0.2 / 30.0) * t_torque
    axs[0, 1].plot(t_torque, torque_df["angular_velocity_body_z_rad_s"], label=r"M10 Simulated $\omega_z$", color="#2ca02c", lw=2)
    axs[0, 1].plot(t_torque, analytical_wz, ":", label=r"Analytical $\frac{\tau_z}{I_{zz}} t$", color="black", lw=2)
    axs[0, 1].set_ylabel("Angular Velocity (rad/s)")
    axs[0, 1].set_title(r"Torque Response: Constant $\tau_z = 0.2$ N$\cdot$m on $I_{zz} = 30$ kg$\cdot$m$^2$")
    axs[0, 1].grid(True)
    axs[0, 1].legend()

    # Bottom-Left: Difference in orbital radius between torque-free and torque-driven
    delta_r_m = np.abs(
        canon_df["position_eci_x_m"].iloc[:len(torque_df)].to_numpy() - torque_df["position_eci_x_m"].to_numpy()
    )
    axs[1, 0].plot(t_torque, delta_r_m, color="#9467bd", label=r"$|\Delta r_{ECI, x}|$ Difference")
    axs[1, 0].set_xlabel("Simulation Time (seconds)")
    axs[1, 0].set_ylabel("Position Difference (m)")
    axs[1, 0].set_title("Cross-Coupling Absence: Machine-Precision Agreement")
    axs[1, 0].grid(True)
    axs[1, 0].legend()

    # Bottom-Right: Torque-Driven Attitude Evolution
    axs[1, 1].plot(t_torque, torque_df["quaternion_w"], label=r"$q_w$", color="#1f77b4")
    axs[1, 1].plot(t_torque, torque_df["quaternion_z"], label=r"$q_z$", color="#d62728")
    axs[1, 1].set_xlabel("Simulation Time (seconds)")
    axs[1, 1].set_ylabel("Quaternion Components")
    axs[1, 1].set_title("Attitude Spinning under Constant Body Torque")
    axs[1, 1].grid(True)
    axs[1, 1].legend()

    fig.suptitle("AstraDock M10: Subsystem Decoupling & Torque Verification", fontsize=13, weight="bold")
    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=300)
    plt.close()
    print(f"Generated Plot D: {output_path}")


def main():
    root = Path(__file__).resolve().parent.parent.parent
    data_dir = root / "data"
    figures_dir = root / "artifacts" / "figures"

    canon_csv = data_dir / "six_dof_canonical_orbit.csv"
    torque_csv = data_dir / "six_dof_constant_torque.csv"

    if not (canon_csv.exists() and torque_csv.exists()):
        print("Telemetry CSVs missing. Running 6-DOF demo executable first...")
        import subprocess
        demo_exe = root / "build" / "Release" / "astradock_6dof_demo.exe"
        subprocess.run([str(demo_exe)], check=True)

    canon_df = pd.read_csv(canon_csv)
    torque_df = pd.read_csv(torque_csv)

    plot_3d_orbit_attitude_triads(canon_df, figures_dir / "m10_3d_orbit_attitude_triads.png")
    plot_translational_telemetry(canon_df, figures_dir / "m10_translational_telemetry.png")
    plot_rotational_telemetry(canon_df, figures_dir / "m10_rotational_telemetry.png")
    plot_regression_and_decoupling(canon_df, torque_df, figures_dir / "m10_regression_and_decoupling.png")
    print("All M10 publication-grade figures generated successfully.")


if __name__ == "__main__":
    main()
