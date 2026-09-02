"""Plotting and visualization tools for M09 Rigid-Body Attitude Dynamics & Quaternion Kinematics."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib.pyplot as plt


def load_csv(filepath: Path) -> dict[str, list[float]]:
    data: dict[str, list[float]] = {}
    with open(filepath, encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for key, val in row.items():
                if key not in data:
                    data[key] = []
                data[key].append(float(val))
    return data


def plot_angular_velocity_evolution(asym_csv: Path, output_path: Path) -> None:
    """Plot A: Angular velocity components over time in asymmetric torque-free tumbling."""
    data = load_csv(asym_csv)
    t = data["time_s"]
    wx = data["wx_rad_s"]
    wy = data["wy_rad_s"]
    wz = data["wz_rad_s"]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, wx, label=r"$\omega_x$ (rad/s)", color="crimson", linewidth=1.8)
    ax.plot(t, wy, label=r"$\omega_y$ (rad/s)", color="forestgreen", linewidth=1.8)
    ax.plot(t, wz, label=r"$\omega_z$ (rad/s)", color="royalblue", linewidth=1.8)

    ax.set_xlabel("Time (s)", fontsize=11)
    ax.set_ylabel("Angular Velocity (rad/s)", fontsize=11)
    ax.set_title("AstraDock M09 — Asymmetric Rigid-Body Tumbling Angular Velocity Evolution\n(Gyroscopic Coupling Term: $\\boldsymbol{\\omega} \\times \\mathbf{I}\\boldsymbol{\\omega}$)", fontsize=12, fontweight="bold")
    ax.grid(True, linestyle="--", alpha=0.6)
    ax.legend(loc="upper right", fontsize=10)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot A: {output_path}")


def plot_quaternion_trajectory(asym_csv: Path, output_path: Path) -> None:
    """Plot B: Quaternion trajectory components and unit norm stability."""
    data = load_csv(asym_csv)
    t = data["time_s"]
    qw = data["qw"]
    qx = data["qx"]
    qy = data["qy"]
    qz = data["qz"]
    norm = data["quaternion_norm"]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 7), sharex=True, gridspec_kw={"height_ratios": [3, 1]})

    ax1.plot(t, qw, label=r"$q_w$ (scalar)", color="darkorange", linewidth=1.6)
    ax1.plot(t, qx, label=r"$q_x$", color="crimson", linewidth=1.6)
    ax1.plot(t, qy, label=r"$q_y$", color="forestgreen", linewidth=1.6)
    ax1.plot(t, qz, label=r"$q_z$", color="royalblue", linewidth=1.6)
    ax1.set_ylabel("Quaternion Parameter", fontsize=11)
    ax1.set_title("AstraDock M09 — Attitude Quaternion Kinematic Propagation ($dq/dt = \\frac{1}{2} q \\otimes \\boldsymbol{\\omega}$)", fontsize=12, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=10, ncol=4)

    norm_residuals = [max(abs(n - 1.0), 1e-16) for n in norm]
    ax2.plot(t, norm_residuals, color="purple", linewidth=1.5)
    ax2.set_xlabel("Time (s)", fontsize=11)
    ax2.set_ylabel(r"$|\|q\| - 1|$", fontsize=11)
    ax2.set_yscale("log")
    ax2.set_ylim([1e-17, 1e-13])
    ax2.set_title("Unit Norm Error (Post-Step Reprojection)", fontsize=10)
    ax2.grid(True, which="both", linestyle="--", alpha=0.6)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot B: {output_path}")


def plot_conservation_invariants(asym_csv: Path, output_path: Path) -> None:
    """Plot C: Conservation of rotational kinetic energy and inertial angular momentum."""
    data = load_csv(asym_csv)
    t = data["time_s"]
    energy = data["rotational_energy_J"]
    e0 = energy[0]
    e_rel_err = [max(abs(e - e0) / e0, 1e-16) for e in energy]

    h_ix = data["Hx_inertial_Nms"]
    h_iy = data["Hy_inertial_Nms"]
    h_iz = data["Hz_inertial_Nms"]
    h0_norm = (h_ix[0] ** 2 + h_iy[0] ** 2 + h_iz[0] ** 2) ** 0.5
    h_vec_err = [
        max((((hx - h_ix[0]) ** 2 + (hy - h_iy[0]) ** 2 + (hz - h_iz[0]) ** 2) ** 0.5) / h0_norm, 1e-16)
        for hx, hy, hz in zip(h_ix, h_iy, h_iz, strict=True)
    ]

    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(t, e_rel_err, label=r"Rotational Kinetic Energy Relative Error $|E(t) - E_0|/E_0$", color="forestgreen", linewidth=1.8)
    ax.plot(t, h_vec_err, label=r"Inertial Angular Momentum Relative Error $\|\mathbf{H}_I(t) - \mathbf{H}_I(0)\|/\|\mathbf{H}_I(0)\|$", color="darkblue", linewidth=1.8)

    ax.set_xlabel("Time (s)", fontsize=11)
    ax.set_ylabel("Relative Error", fontsize=11)
    ax.set_yscale("log")
    ax.set_title("AstraDock M09 — Physical Conservation Invariants in Torque-Free Tumbling", fontsize=12, fontweight="bold")
    ax.grid(True, which="both", linestyle="--", alpha=0.6)
    ax.legend(loc="upper left", fontsize=10)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot C: {output_path}")


def plot_torque_analytical_comparison(torque_csv: Path, output_path: Path) -> None:
    """Plot D: Constant torque acceleration vs exact analytical reference."""
    data = load_csv(torque_csv)
    t = data["time_s"]
    wx_num = data["wx_rad_s"]

    # Analytical: tau_x = 0.5, Ixx = 10.0 => alpha_x = 0.05 rad/s^2
    alpha_x = 0.5 / 10.0
    wx_ana = [alpha_x * ti for ti in t]
    diff = [max(abs(num - ana), 1e-16) for num, ana in zip(wx_num, wx_ana, strict=True)]

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 6), sharex=True, gridspec_kw={"height_ratios": [2.5, 1]})

    ax1.plot(t, wx_num, label="Numerical RK4 $\\omega_x(t)$", color="crimson", linewidth=2.0)
    ax1.plot(t, wx_ana, label="Analytical Reference $\\omega_x(t) = \\frac{\\tau_x}{I_{xx}} t$", color="black", linestyle="--", linewidth=1.5)
    ax1.set_ylabel(r"$\omega_x$ (rad/s)", fontsize=11)
    ax1.set_title(r"AstraDock M09 — Constant Principal-Axis Torque Response ($\tau_x = 0.5\,\mathrm{N\cdot m}, I_{xx} = 10\,\mathrm{kg\cdot m^2}$)", fontsize=12, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(loc="upper left", fontsize=10)

    ax2.plot(t, diff, color="purple", linewidth=1.5)
    ax2.set_xlabel("Time (s)", fontsize=11)
    ax2.set_ylabel(r"$|\omega_{num} - \omega_{ana}|$ (rad/s)", fontsize=10)
    ax2.set_yscale("log")
    ax2.set_title("Numerical Residual Error", fontsize=10)
    ax2.grid(True, which="both", linestyle="--", alpha=0.6)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot D: {output_path}")


def plot_body_vs_inertial_momentum(asym_csv: Path, output_path: Path) -> None:
    """Plot E: Body frame vs Inertial frame angular momentum vectors."""
    data = load_csv(asym_csv)
    t = data["time_s"]

    h_bx = data["Hx_body_Nms"]
    h_by = data["Hy_body_Nms"]
    h_bz = data["Hz_body_Nms"]

    h_ix = data["Hx_inertial_Nms"]
    h_iy = data["Hy_inertial_Nms"]
    h_iz = data["Hz_inertial_Nms"]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), sharey=True)

    # Subplot 1: Body Components
    ax1.plot(t, h_bx, label=r"$H_{B,x}$", color="crimson", linewidth=1.8)
    ax1.plot(t, h_by, label=r"$H_{B,y}$", color="forestgreen", linewidth=1.8)
    ax1.plot(t, h_bz, label=r"$H_{B,z}$", color="royalblue", linewidth=1.8)
    ax1.set_xlabel("Time (s)", fontsize=11)
    ax1.set_ylabel(r"Angular Momentum Components ($\mathrm{N\cdot m\cdot s}$)", fontsize=11)
    ax1.set_title(r"Body Frame $\mathbf{H}_B = \mathbf{I}\boldsymbol{\omega}$ (Varying)", fontsize=11, fontweight="bold")
    ax1.grid(True, linestyle="--", alpha=0.6)
    ax1.legend(loc="upper right", fontsize=10)

    # Subplot 2: Inertial Components
    ax2.plot(t, h_ix, label=r"$H_{I,x}$", color="crimson", linestyle="--", linewidth=1.8)
    ax2.plot(t, h_iy, label=r"$H_{I,y}$", color="forestgreen", linestyle="--", linewidth=1.8)
    ax2.plot(t, h_iz, label=r"$H_{I,z}$", color="royalblue", linestyle="--", linewidth=1.8)
    ax2.set_xlabel("Time (s)", fontsize=11)
    ax2.set_title(r"Inertial Frame $\mathbf{H}_I = C_{I\_B}\mathbf{H}_B$ (Conserved)", fontsize=11, fontweight="bold")
    ax2.grid(True, linestyle="--", alpha=0.6)
    ax2.legend(loc="upper right", fontsize=10)

    fig.suptitle("AstraDock M09 — Body-Frame vs Inertial-Frame Angular Momentum Distinction", fontsize=12, fontweight="bold")
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot E: {output_path}")


def main() -> None:
    figures_dir = Path("artifacts/figures")
    data_dir = Path("data")

    asym_csv = data_dir / "attitude_dynamics_asymmetric_tumble.csv"
    torque_csv = data_dir / "attitude_dynamics_constant_torque.csv"

    if asym_csv.exists():
        plot_angular_velocity_evolution(asym_csv, figures_dir / "m09_angular_velocity_evolution.png")
        plot_quaternion_trajectory(asym_csv, figures_dir / "m09_attitude_quaternion_trajectory.png")
        plot_conservation_invariants(asym_csv, figures_dir / "m09_conservation_invariants.png")
        plot_body_vs_inertial_momentum(asym_csv, figures_dir / "m09_body_vs_inertial_momentum.png")

    if torque_csv.exists():
        plot_torque_analytical_comparison(torque_csv, figures_dir / "m09_torque_analytical_comparison.png")

    print("All M09 attitude dynamics figures generated successfully.")


if __name__ == "__main__":
    main()
