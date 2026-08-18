"""Orbital elements visualization script for AstraDock M07.

Generates:
  1. artifacts/figures/m07_orbital_geometry.png (Plot A: 3D orbit geometry, angles, vectors)
  2. artifacts/figures/m07_elements_time_evolution.png (Plot B: Invariance of a, e, i, Omega, omega and evolution of nu)
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D


def load_csv(csv_path: Path) -> dict[str, np.ndarray]:
    data: dict[str, list[float]] = {}
    with csv_path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        for row in reader:
            for k, v in row.items():
                data.setdefault(k, []).append(float(v))
    return {k: np.array(v) for k, v in data.items()}


def plot_orbital_geometry(csv_data: dict[str, np.ndarray], output_path: Path) -> None:
    x_km = csv_data["x_m"] / 1000.0
    y_km = csv_data["y_m"] / 1000.0
    z_km = csv_data["z_m"] / 1000.0

    fig = plt.figure(figsize=(11, 9))
    ax = fig.add_subplot(111, projection="3d")

    # Plot Earth as a blue-green wireframe sphere
    r_earth_km = 6378.137
    u = np.linspace(0, 2 * np.pi, 40)
    v = np.linspace(0, np.pi, 20)
    xs = r_earth_km * np.outer(np.cos(u), np.sin(v))
    ys = r_earth_km * np.outer(np.sin(u), np.sin(v))
    zs = r_earth_km * np.outer(np.ones(np.size(u)), np.cos(v))
    ax.plot_surface(xs, ys, zs, color="#1f77b4", alpha=0.25, edgecolor="#0d47a1", linewidth=0.3)

    # Plot equatorial reference circle
    theta_eq = np.linspace(0, 2 * np.pi, 100)
    r_eq_grid = 12000.0
    ax.plot(
        r_eq_grid * np.cos(theta_eq),
        r_eq_grid * np.sin(theta_eq),
        np.zeros_like(theta_eq),
        color="gray",
        linestyle="--",
        linewidth=0.8,
        alpha=0.5,
        label="Equatorial Reference Plane (XY)",
    )

    # Plot full 3D orbit trajectory
    # One orbit slice
    ax.plot(
        x_km,
        y_km,
        z_km,
        color="#d32f2f",
        linewidth=2.0,
        label="Elliptical Inclined Orbit (a=10,000 km, e=0.2, i=45°)",
    )

    # Spacecraft position at t=0
    ax.scatter(
        [x_km[0]],
        [y_km[0]],
        [z_km[0]],
        color="#ffd600",
        edgecolors="black",
        s=90,
        zorder=10,
        label="Spacecraft Position (ν = 30°)",
    )

    # Coordinate axes
    axis_len = 14000.0
    ax.quiver(0, 0, 0, axis_len, 0, 0, color="red", arrow_length_ratio=0.08, linewidth=1.2)
    ax.quiver(0, 0, 0, 0, axis_len, 0, color="green", arrow_length_ratio=0.08, linewidth=1.2)
    ax.quiver(0, 0, 0, 0, 0, axis_len, color="blue", arrow_length_ratio=0.08, linewidth=1.2)
    ax.text(axis_len * 1.05, 0, 0, "+X (Vernal Equinox)", color="red", fontsize=9, weight="bold")
    ax.text(0, axis_len * 1.05, 0, "+Y_ECI", color="green", fontsize=9, weight="bold")
    ax.text(0, 0, axis_len * 1.05, "+Z (North Pole)", color="blue", fontsize=9, weight="bold")

    # Ascending node vector (n)
    raan_rad = float(csv_data["raan_rad"][0])
    node_len = 11000.0
    nx = node_len * math.cos(raan_rad)
    ny = node_len * math.sin(raan_rad)
    ax.quiver(
        0, 0, 0, nx, ny, 0, color="#ff9800", arrow_length_ratio=0.08, linewidth=2.0, linestyle="-"
    )
    ax.text(
        nx * 1.05,
        ny * 1.05,
        0,
        f"Line of Nodes (Ω = {math.degrees(raan_rad):.0f}°)",
        color="#ff9800",
        fontsize=9,
        weight="bold",
    )

    ax.set_xlabel("X_ECI [km]")
    ax.set_ylabel("Y_ECI [km]")
    ax.set_zlabel("Z_ECI [km]")
    ax.set_title(
        "AstraDock M07: Classical Keplerian Orbit Geometry in ECI\n"
        "(a = 10,000 km, e = 0.2, i = 45°, Ω = 120°, ω = 60°, ν₀ = 30°)",
        fontsize=12,
        weight="bold",
        pad=15,
    )

    limit = 14000.0
    ax.set_xlim([-limit, limit])
    ax.set_ylim([-limit, limit])
    ax.set_zlim([-limit, limit])
    ax.view_init(elev=28, azim=45)
    ax.legend(loc="upper left", fontsize=8)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.tight_layout()
    plt.savefig(output_path, dpi=180)
    plt.close()
    print(f"Saved Plot A: {output_path}")


def plot_elements_time_evolution(csv_data: dict[str, np.ndarray], output_path: Path) -> None:
    time_hr = csv_data["time_s"] / 3600.0
    a_km = csv_data["a_m"] / 1000.0
    e = csv_data["eccentricity"]
    inc_deg = np.degrees(csv_data["inclination_rad"])
    raan_deg = np.degrees(csv_data["raan_rad"])
    argp_deg = np.degrees(csv_data["argument_of_periapsis_rad"])
    nu_deg = np.degrees(csv_data["true_anomaly_rad"])

    fig, axes = plt.subplots(3, 2, figsize=(13, 10), sharex=True)

    # 1. Semi-major axis
    axes[0, 0].plot(time_hr, a_km, color="#1976d2", linewidth=1.5)
    axes[0, 0].set_ylabel("Semi-Major Axis $a$ [km]")
    axes[0, 0].set_title("Semi-Major Axis $a(t)$ (Invariant)")
    axes[0, 0].grid(True, linestyle=":", alpha=0.6)
    axes[0, 0].ticklabel_format(useOffset=False)

    # 2. Eccentricity
    axes[0, 1].plot(time_hr, e, color="#388e3c", linewidth=1.5)
    axes[0, 1].set_ylabel("Eccentricity $e$ [-]")
    axes[0, 1].set_title("Eccentricity $e(t)$ (Invariant)")
    axes[0, 1].grid(True, linestyle=":", alpha=0.6)
    axes[0, 1].ticklabel_format(useOffset=False)

    # 3. Inclination
    axes[1, 0].plot(time_hr, inc_deg, color="#7b1fa2", linewidth=1.5)
    axes[1, 0].set_ylabel("Inclination $i$ [deg]")
    axes[1, 0].set_title("Inclination $i(t)$ (Invariant)")
    axes[1, 0].grid(True, linestyle=":", alpha=0.6)
    axes[1, 0].ticklabel_format(useOffset=False)

    # 4. RAAN
    axes[1, 1].plot(time_hr, raan_deg, color="#f57c00", linewidth=1.5)
    axes[1, 1].set_ylabel("RAAN $\\Omega$ [deg]")
    axes[1, 1].set_title("RAAN $\\Omega(t)$ (Invariant)")
    axes[1, 1].grid(True, linestyle=":", alpha=0.6)
    axes[1, 1].ticklabel_format(useOffset=False)

    # 5. Argument of Periapsis
    axes[2, 0].plot(time_hr, argp_deg, color="#c2185b", linewidth=1.5)
    axes[2, 0].set_xlabel("Elapsed Time [hours]")
    axes[2, 0].set_ylabel("Arg of Periapsis $\\omega$ [deg]")
    axes[2, 0].set_title("Argument of Periapsis $\\omega(t)$ (Invariant)")
    axes[2, 0].grid(True, linestyle=":", alpha=0.6)
    axes[2, 0].ticklabel_format(useOffset=False)

    # 6. True Anomaly
    axes[2, 1].plot(time_hr, nu_deg, color="#d32f2f", linewidth=1.2, label="$\\nu(t)$")
    axes[2, 1].set_xlabel("Elapsed Time [hours]")
    axes[2, 1].set_ylabel("True Anomaly $\\nu$ [deg]")
    axes[2, 1].set_title("True Anomaly $\\nu(t)$ (Keplerian In-Orbit Motion)")
    axes[2, 1].grid(True, linestyle=":", alpha=0.6)

    legend_elements = [
        Line2D([0], [0], color="#1976d2", lw=2, label="Constant Elements (a, e, i, $\\Omega$, $\\omega$)"),
        Line2D([0], [0], color="#d32f2f", lw=2, label="Evolving Element ($\\nu$)"),
    ]
    fig.legend(handles=legend_elements, loc="upper center", ncol=2, fontsize=10, frameon=True)

    fig.suptitle(
        "AstraDock M07: Orbital Element Invariance & True Anomaly Evolution (3 Orbits)",
        fontsize=13,
        weight="bold",
        y=0.98,
    )
    plt.tight_layout(rect=(0, 0, 1, 0.95))
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=180)
    plt.close()
    print(f"Saved Plot B: {output_path}")


def main() -> None:
    csv_file = Path("artifacts/data/m07_elements_propagation.csv")
    if not csv_file.is_file():
        print(f"CSV file not found: {csv_file}. Please run astradock_elements_demo first.")
        return

    data = load_csv(csv_file)
    plot_orbital_geometry(data, Path("artifacts/figures/m07_orbital_geometry.png"))
    plot_elements_time_evolution(data, Path("artifacts/figures/m07_elements_time_evolution.png"))
    print("All M07 orbital element figures generated successfully.")


if __name__ == "__main__":
    main()
