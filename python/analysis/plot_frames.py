"""Plot authoritative AstraDock C++ coordinate frames CSV output.

This module performs display-unit conversions and plotting only. It contains no
orbital mechanics, basis construction, or DCM calculations. All plotted vectors
and metrics originate directly from C++.
"""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

plt.switch_backend("Agg")

EARTH_RADIUS_KM = 6378.137

# Color palette matching AstraDock theme
BLUE = "#2563EB"
RED = "#DC2626"
GREEN = "#059669"
PURPLE = "#7C3AED"
INK = "#1F2937"
MUTED = "#6B7280"
GRID = "#E5E7EB"
EARTH_FILL = "#DBEAFE"
EARTH_EDGE = "#3B82F6"

REQUIRED_COLUMNS = {
    "time_s",
    "position_eci_x_m",
    "position_eci_y_m",
    "position_eci_z_m",
    "velocity_eci_x_m_per_s",
    "velocity_eci_y_m_per_s",
    "velocity_eci_z_m_per_s",
    "lvlh_x_eci_x",
    "lvlh_x_eci_y",
    "lvlh_x_eci_z",
    "lvlh_y_eci_x",
    "lvlh_y_eci_y",
    "lvlh_y_eci_z",
    "lvlh_z_eci_x",
    "lvlh_z_eci_y",
    "lvlh_z_eci_z",
    "r_lvlh_x_m",
    "r_lvlh_y_m",
    "r_lvlh_z_m",
    "v_lvlh_x_m_per_s",
    "v_lvlh_y_m_per_s",
    "v_lvlh_z_m_per_s",
    "det_c",
    "orthonormality_max_err",
    "norm_diff_m",
    "round_trip_err_m",
}


def read_frames_csv(path: Path) -> dict[str, np.ndarray]:
    """Read numeric series from C++ exported frame data."""
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"CSV is missing required columns ({missing_list}): {path}")

        columns: dict[str, list[float]] = {name: [] for name in reader.fieldnames}
        for row in reader:
            for name in reader.fieldnames:
                columns[name].append(float(row[name]))

    return {name: np.array(vals) for name, vals in columns.items()}


def plot_single_point_frame(data: dict[str, np.ndarray], output_path: Path) -> None:
    """Plot A: ECI and LVLH coordinate triads at a single orbital position."""
    fig, ax = plt.subplots(figsize=(8, 8))

    # Plot Earth
    earth = plt.Circle((0, 0), EARTH_RADIUS_KM, color=EARTH_FILL, ec=EARTH_EDGE, lw=1.5, zorder=2)
    ax.add_patch(earth)
    ax.text(0, 0, "Earth\n(Inertial Center)", ha="center", va="center", color=INK, fontsize=10, weight="bold")

    # Plot Orbit track
    r_x_km = data["position_eci_x_m"] / 1000.0
    r_y_km = data["position_eci_y_m"] / 1000.0
    ax.plot(r_x_km, r_y_km, "--", color=MUTED, lw=1.2, label="500 km Orbit Track", zorder=3)

    # Spacecraft position at theta = 45 deg (sample ~ 70)
    idx = len(data["time_s"]) // 8
    sc_x = r_x_km[idx]
    sc_y = r_y_km[idx]

    # Spacecraft marker
    ax.plot(sc_x, sc_y, "o", color=INK, markersize=8, zorder=5, label="Spacecraft")

    # Scale factor for basis vector arrows (in km)
    arrow_len = 2500.0

    # LVLH basis vectors at this point
    er_x = data["lvlh_x_eci_x"][idx] * arrow_len
    er_y = data["lvlh_x_eci_y"][idx] * arrow_len
    et_x = data["lvlh_y_eci_x"][idx] * arrow_len
    et_y = data["lvlh_y_eci_y"][idx] * arrow_len

    # Plot LVLH triad
    ax.quiver(sc_x, sc_y, er_x, er_y, angles="xy", scale_units="xy", scale=1,
              color=RED, width=0.007, zorder=6, label=r"LVLH Radial $\mathbf{e}_r$ (X)")
    ax.quiver(sc_x, sc_y, et_x, et_y, angles="xy", scale_units="xy", scale=1,
              color=GREEN, width=0.007, zorder=6, label=r"LVLH Along-Track $\mathbf{e}_t$ (Y)")

    # Plot ECI reference axes at Earth center
    eci_arrow_len = 8000.0
    ax.quiver(0, 0, eci_arrow_len, 0, angles="xy", scale_units="xy", scale=1,
              color=BLUE, width=0.005, zorder=4, label=r"ECI $+X_{ECI}$")
    ax.quiver(0, 0, 0, eci_arrow_len, angles="xy", scale_units="xy", scale=1,
              color=PURPLE, width=0.005, zorder=4, label=r"ECI $+Y_{ECI}$")

    ax.set_aspect("equal")
    limit = 10500.0
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("ECI X Position (km)", fontsize=11)
    ax.set_ylabel("ECI Y Position (km)", fontsize=11)
    ax.set_title("AstraDock — ECI vs Local Orbital (LVLH) Reference Frames", fontsize=13, weight="bold", pad=12)
    ax.grid(True, linestyle=":", color=GRID, alpha=0.8)
    ax.legend(loc="upper right", framealpha=0.9, fontsize=9)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_orbit_evolution(data: dict[str, np.ndarray], output_path: Path) -> None:
    """Plot B: LVLH basis rotating as spacecraft traverses circular orbit."""
    fig, ax = plt.subplots(figsize=(9, 9))

    # Plot Earth
    earth = plt.Circle((0, 0), EARTH_RADIUS_KM, color=EARTH_FILL, ec=EARTH_EDGE, lw=1.5, zorder=2)
    ax.add_patch(earth)
    ax.text(0, 0, "Earth\n(ECI Origin)", ha="center", va="center", color=INK, fontsize=10, weight="bold")

    # Orbit path
    r_x_km = data["position_eci_x_m"] / 1000.0
    r_y_km = data["position_eci_y_m"] / 1000.0
    ax.plot(r_x_km, r_y_km, "--", color=MUTED, lw=1.2, label="500 km Orbit Path", zorder=3)

    # ECI fixed axes
    eci_len = 8500.0
    ax.quiver(0, 0, eci_len, 0, angles="xy", scale_units="xy", scale=1,
              color=BLUE, width=0.005, zorder=4, label=r"Fixed ECI $X$")
    ax.quiver(0, 0, 0, eci_len, angles="xy", scale_units="xy", scale=1,
              color=PURPLE, width=0.005, zorder=4, label=r"Fixed ECI $Y$")

    # 8 sample stations around the orbit
    num_samples = len(data["time_s"])
    indices = [int(i * (num_samples - 1) / 8) for i in range(8)]

    arrow_len = 1600.0
    for idx in indices:
        sc_x = r_x_km[idx]
        sc_y = r_y_km[idx]

        er_x = data["lvlh_x_eci_x"][idx] * arrow_len
        er_y = data["lvlh_x_eci_y"][idx] * arrow_len
        et_x = data["lvlh_y_eci_x"][idx] * arrow_len
        et_y = data["lvlh_y_eci_y"][idx] * arrow_len

        ax.plot(sc_x, sc_y, "o", color=INK, markersize=6, zorder=5)
        ax.quiver(sc_x, sc_y, er_x, er_y, angles="xy", scale_units="xy", scale=1,
                  color=RED, width=0.005, zorder=6)
        ax.quiver(sc_x, sc_y, et_x, et_y, angles="xy", scale_units="xy", scale=1,
                  color=GREEN, width=0.005, zorder=6)

    # Custom legend entries for LVLH vectors using Line2D proxy artists
    from matplotlib.lines import Line2D
    custom_lines = [
        Line2D([0], [0], color=MUTED, linestyle="--", lw=1.2, label="500 km Orbit Path"),
        Line2D([0], [0], color=BLUE, lw=2.0, label=r"Fixed ECI $X$"),
        Line2D([0], [0], color=PURPLE, lw=2.0, label=r"Fixed ECI $Y$"),
        Line2D([0], [0], color=RED, lw=2.5, label=r"LVLH $\mathbf{e}_r$ (Radial)"),
        Line2D([0], [0], color=GREEN, lw=2.5, label=r"LVLH $\mathbf{e}_t$ (Along-Track)"),
    ]

    ax.set_aspect("equal")
    limit = 10500.0
    ax.set_xlim(-limit, limit)
    ax.set_ylim(-limit, limit)
    ax.set_xlabel("ECI X (km)", fontsize=11)
    ax.set_ylabel("ECI Y (km)", fontsize=11)
    ax.set_title("AstraDock — LVLH Triad Rotation Over Full Orbit", fontsize=13, weight="bold", pad=12)
    ax.grid(True, linestyle=":", color=GRID, alpha=0.8)
    ax.legend(handles=custom_lines, loc="upper right", framealpha=0.9, fontsize=9)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close()


def plot_transform_sanity(data: dict[str, np.ndarray], output_path: Path) -> None:
    """Plot C: Vector coordinate transformation comparison and numerical invariants."""
    time_min = data["time_s"] / 60.0

    fig, axes = plt.subplots(3, 1, figsize=(10, 10), sharex=True)

    # Subplot 1: Position components (ECI vs LVLH)
    ax1 = axes[0]
    ax1.plot(time_min, data["position_eci_x_m"] / 1000.0, color=BLUE, label=r"$r_{ECI, x}$")
    ax1.plot(time_min, data["position_eci_y_m"] / 1000.0, color=PURPLE, label=r"$r_{ECI, y}$")
    ax1.plot(time_min, data["r_lvlh_x_m"] / 1000.0, "--", color=RED, lw=2.0, label=r"$r_{LVLH, x}$ (Radial)")
    ax1.plot(time_min, data["r_lvlh_y_m"] / 1000.0, ":", color=GREEN, lw=2.0, label=r"$r_{LVLH, y}$ (Along-track)")
    ax1.set_ylabel("Position (km)", fontsize=10)
    ax1.set_title("1. Position Coordinates: Oscillating in ECI vs Constant Radial in LVLH", fontsize=11, weight="bold")
    ax1.grid(True, linestyle=":", color=GRID)
    ax1.legend(loc="upper right", ncol=4, fontsize=9)

    # Subplot 2: Velocity components (ECI vs LVLH)
    ax2 = axes[1]
    ax2.plot(time_min, data["velocity_eci_x_m_per_s"] / 1000.0, color=BLUE, label=r"$v_{ECI, x}$")
    ax2.plot(time_min, data["velocity_eci_y_m_per_s"] / 1000.0, color=PURPLE, label=r"$v_{ECI, y}$")
    ax2.plot(time_min, data["v_lvlh_x_m_per_s"] / 1000.0, ":", color=RED, lw=2.0, label=r"$v_{LVLH, x}$")
    ax2.plot(time_min, data["v_lvlh_y_m_per_s"] / 1000.0, "--", color=GREEN, lw=2.0, label=r"$v_{LVLH, y}$ (Along-track)")
    ax2.set_ylabel("Velocity (km/s)", fontsize=10)
    ax2.set_title("2. Velocity Coordinates: Oscillating in ECI vs Constant Along-Track in LVLH", fontsize=11, weight="bold")
    ax2.grid(True, linestyle=":", color=GRID)
    ax2.legend(loc="upper right", ncol=4, fontsize=9)

    # Subplot 3: Norm difference & Round-trip transformation error
    ax3 = axes[2]
    r_eci_norm = np.sqrt(data["position_eci_x_m"]**2 + data["position_eci_y_m"]**2 + data["position_eci_z_m"]**2)
    r_lvlh_norm = np.sqrt(data["r_lvlh_x_m"]**2 + data["r_lvlh_y_m"]**2 + data["r_lvlh_z_m"]**2)
    norm_diff_m = np.abs(r_lvlh_norm - r_eci_norm)

    ax3.semilogy(time_min, np.maximum(norm_diff_m, 1e-15), color=BLUE, label=r"Norm Invariance $|\|\mathbf{r}_{LVLH}\| - \|\mathbf{r}_{ECI}\||$ (m)")
    ax3.semilogy(time_min, np.maximum(data["round_trip_err_m"], 1e-15), "--", color=RED, label=r"Round-Trip Error $\|C^T C \mathbf{r} - \mathbf{r}\|$ (m)")
    ax3.set_xlabel("Orbital Time (minutes)", fontsize=10)
    ax3.set_ylabel("Error Metric (m)", fontsize=10)
    ax3.set_title("3. Coordinate Transformation Numerical Invariants & Round-Trip Consistency", fontsize=11, weight="bold")
    ax3.grid(True, linestyle=":", color=GRID)
    ax3.legend(loc="upper right", fontsize=9)

    plt.tight_layout()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(output_path, dpi=200)
    plt.close()


def main() -> None:
    parser = argparse.ArgumentParser(description="Plot AstraDock coordinate frames analysis.")
    parser.add_argument("--csv", type=Path, default=Path("artifacts/data/m06_frames.csv"), help="Path to frames CSV")
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/figures"), help="Output directory for plots")
    args = parser.parse_args()

    if not args.csv.is_file():
        raise FileNotFoundError(f"Frames CSV file not found: {args.csv}")

    data = read_frames_csv(args.csv)

    single_point_png = args.output_dir / "m06_frame_single_point.png"
    evolution_png = args.output_dir / "m06_frame_orbit_evolution.png"
    transform_png = args.output_dir / "m06_frame_transform_sanity.png"

    print(f"Generating Plot A: {single_point_png}...")
    plot_single_point_frame(data, single_point_png)

    print(f"Generating Plot B: {evolution_png}...")
    plot_orbit_evolution(data, evolution_png)

    print(f"Generating Plot C: {transform_png}...")
    plot_transform_sanity(data, transform_png)

    print("All M06 coordinate frame plots generated successfully.")


if __name__ == "__main__":
    main()
