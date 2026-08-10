"""Plot authoritative AstraDock C++ orbit CSV output.

This module performs display-unit conversions and charting only. It contains no
gravity model, state derivative, or numerical integration implementation.
"""

from __future__ import annotations

import argparse
import csv
from collections.abc import Sequence
from pathlib import Path

import matplotlib.pyplot as plt

plt.switch_backend("Agg")

BLUE = "#2563EB"
ORANGE = "#D97706"
INK = "#1F2937"
MUTED = "#6B7280"
GRID = "#D1D5DB"
EARTH_FILL = "#E5E7EB"
EARTH_EDGE = "#4B5563"

REQUIRED_COLUMNS = {
    "integrator",
    "time_s",
    "x_m",
    "y_m",
    "radius_m",
    "altitude_m",
    "relative_energy_error",
    "relative_angular_momentum_error",
}

Series = dict[str, list[float]]


def read_cpp_trajectory(path: Path) -> dict[str, Series]:
    """Read numeric series grouped by the C++-reported integrator name."""
    with path.open(newline="", encoding="utf-8") as source:
        reader = csv.DictReader(source)
        if reader.fieldnames is None:
            raise ValueError(f"CSV has no header: {path}")
        missing = REQUIRED_COLUMNS.difference(reader.fieldnames)
        if missing:
            missing_list = ", ".join(sorted(missing))
            raise ValueError(f"CSV is missing required columns ({missing_list}): {path}")

        numeric_columns = [name for name in reader.fieldnames if name != "integrator"]
        grouped: dict[str, Series] = {}
        for row in reader:
            method = row["integrator"]
            series = grouped.setdefault(method, {name: [] for name in numeric_columns})
            for name in numeric_columns:
                series[name].append(float(row[name]))

    for required_method in ("euler", "rk4"):
        if required_method not in grouped or not grouped[required_method]["time_s"]:
            raise ValueError(f"CSV contains no {required_method} samples: {path}")
    return grouped


def style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.7)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_color(MUTED)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(colors=INK)
    ax.xaxis.label.set_color(INK)
    ax.yaxis.label.set_color(INK)


def make_figure() -> tuple[plt.Figure, plt.Axes]:
    fig, ax = plt.subplots(figsize=(10, 6.25))
    fig.subplots_adjust(left=0.11, right=0.98, bottom=0.12, top=0.82)
    return fig, ax


def add_title(fig: plt.Figure, title: str, subtitle: str) -> None:
    fig.text(0.08, 0.95, title, ha="left", fontsize=16, color=INK, weight="bold")
    fig.text(0.08, 0.90, subtitle, ha="left", fontsize=10, color=MUTED)


def save_figure(fig: plt.Figure, path: Path) -> None:
    fig.savefig(path, dpi=160, facecolor="white")
    plt.close(fig)


def plot_trajectory(one_orbit: dict[str, Series], output_path: Path) -> None:
    rk4 = one_orbit["rk4"]
    x_km = [value / 1_000.0 for value in rk4["x_m"]]
    y_km = [value / 1_000.0 for value in rk4["y_m"]]
    earth_radius_km = (rk4["radius_m"][0] - rk4["altitude_m"][0]) / 1_000.0
    dt_s = rk4["time_s"][1] - rk4["time_s"][0]

    fig, ax = make_figure()
    earth = plt.Circle(
        (0.0, 0.0),
        earth_radius_km,
        facecolor=EARTH_FILL,
        edgecolor=EARTH_EDGE,
        linewidth=1.5,
        label="Earth reference radius",
        zorder=1,
    )
    ax.add_patch(earth)
    ax.plot(x_km, y_km, color=BLUE, linewidth=2.0, label="RK4 trajectory", zorder=2)
    ax.scatter(
        [x_km[0]],
        [y_km[0]],
        color=ORANGE,
        edgecolor=INK,
        linewidth=0.8,
        s=55,
        label="Starting point",
        zorder=3,
    )
    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("Earth-centered inertial X (km)")
    ax.set_ylabel("Earth-centered inertial Y (km)")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    add_title(
        fig,
        "500 km circular-orbit trajectory",
        f"C++ two-body output · classical RK4 · nominal dt = {dt_s:g} s · equal axes",
    )
    save_figure(fig, output_path)


def plot_altitude_comparison(multi_orbit: dict[str, Series], output_path: Path) -> None:
    fig, ax = make_figure()
    for method, color, linestyle, label in (
        ("euler", ORANGE, "-", "Forward Euler"),
        ("rk4", BLUE, "--", "Classical RK4"),
    ):
        series = multi_orbit[method]
        time_min = [value / 60.0 for value in series["time_s"]]
        altitude_km = [value / 1_000.0 for value in series["altitude_m"]]
        ax.plot(
            time_min,
            altitude_km,
            color=color,
            linestyle=linestyle,
            linewidth=2.0,
            label=label,
        )

    ax.axhline(500.0, color=EARTH_EDGE, linewidth=1.0, alpha=0.8, label="Analytical 500 km")
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Altitude above reference radius (km)")
    ax.legend(frameon=False, loc="upper left")
    style_axes(ax)
    add_title(
        fig,
        "Orbital altitude: Euler vs RK4",
        "Five ideal two-body periods · same initial state and nominal 10 s step",
    )
    save_figure(fig, output_path)


def plot_relative_drift(
    multi_orbit: dict[str, Series],
    column: str,
    title: str,
    y_label: str,
    output_path: Path,
) -> None:
    fig, ax = make_figure()
    plotted_values: list[float] = []
    for method, color, linestyle, label in (
        ("euler", ORANGE, "-", "Forward Euler"),
        ("rk4", BLUE, "--", "Classical RK4"),
    ):
        series = multi_orbit[method]
        time_min = [value / 60.0 for value in series["time_s"]]
        ax.plot(
            time_min,
            series[column],
            color=color,
            linestyle=linestyle,
            linewidth=2.0,
            label=label,
        )
        plotted_values.extend(series[column])

    ax.axhline(0.0, color=EARTH_EDGE, linewidth=1.0, alpha=0.8)
    ax.set_yscale("symlog", linthresh=1.0e-12)
    negative_values = [value for value in plotted_values if value < 0.0]
    positive_values = [value for value in plotted_values if value > 0.0]
    lower_limit = min(negative_values) * 5.0 if negative_values else -1.0e-12
    upper_limit = max(positive_values) * 1.25 if positive_values else 1.0e-12
    ax.set_ylim(lower_limit, upper_limit)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel(y_label)
    ax.legend(frameon=False, loc="best")
    style_axes(ax)
    add_title(
        fig,
        title,
        "Five ideal two-body periods · signed relative error · symmetric-log vertical scale",
    )
    save_figure(fig, output_path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Plot C++-generated AstraDock two-body orbit CSV files."
    )
    parser.add_argument("one_orbit_csv", type=Path)
    parser.add_argument("five_orbit_csv", type=Path)
    parser.add_argument("output_directory", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    one_orbit = read_cpp_trajectory(args.one_orbit_csv)
    five_orbit = read_cpp_trajectory(args.five_orbit_csv)
    args.output_directory.mkdir(parents=True, exist_ok=True)

    outputs = {
        "trajectory": args.output_directory / "orbit_500km_rk4.png",
        "altitude": args.output_directory / "euler_vs_rk4_altitude.png",
        "energy": args.output_directory / "euler_vs_rk4_energy_drift.png",
        "angular_momentum": (
            args.output_directory / "euler_vs_rk4_angular_momentum_drift.png"
        ),
    }
    plot_trajectory(one_orbit, outputs["trajectory"])
    plot_altitude_comparison(five_orbit, outputs["altitude"])
    plot_relative_drift(
        five_orbit,
        "relative_energy_error",
        "Relative specific orbital energy error",
        r"$(\epsilon - \epsilon_0) / |\epsilon_0|$",
        outputs["energy"],
    )
    plot_relative_drift(
        five_orbit,
        "relative_angular_momentum_error",
        "Relative specific angular-momentum drift",
        r"$(|h| - |h_0|) / |h_0|$",
        outputs["angular_momentum"],
    )

    for name, path in outputs.items():
        print(f"{name}={path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
