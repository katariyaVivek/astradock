"""M05 propagator validation analysis and convergence plots.

Reads the authoritative C++-generated M05 CSV files and produces:
- log-log convergence plots (dt vs error)
- Euler vs RK4 comparison plots
- Convergence order summary

This module performs analysis and visualization only. It contains no gravity
model, state derivative, or numerical integration implementation.
"""

from __future__ import annotations

import argparse
import csv
import math
from collections import defaultdict
from collections.abc import Sequence
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Color palette matching existing plot_orbit.py
BLUE = "#2563EB"
ORANGE = "#D97706"
INK = "#1F2937"
MUTED = "#6B7280"
GRID = "#D1D5DB"


def read_csv(path: Path) -> list[dict[str, str]]:
    """Read a CSV into a list of row dicts."""
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        return list(reader)


def group_by_integrator(
    rows: list[dict[str, str]],
) -> dict[str, dict[str, list[float]]]:
    """Group numeric CSV columns by integrator name."""
    grouped: dict[str, dict[str, list[float]]] = defaultdict(lambda: defaultdict(list))
    for row in rows:
        name = row["integrator"]
        for key, val in row.items():
            if key == "integrator":
                continue
            try:
                grouped[name][key].append(float(val))
            except ValueError:
                grouped[name][key].append(math.nan)
    return dict(grouped)


def empirical_order(h1: float, e1: float, h2: float, e2: float) -> float:
    """Compute log(E1/E2) / log(h1/h2)."""
    if e1 <= 0 or e2 <= 0 or h1 <= h2 or h2 <= 0:
        return math.nan
    return math.log(e1 / e2) / math.log(h1 / h2)


def style_axes(ax: plt.Axes) -> None:
    ax.set_facecolor("white")
    ax.grid(True, color=GRID, linewidth=0.8, alpha=0.7, which="both")
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


# ---------------------------------------------------------------------------
# Plot A — RK4 timestep vs position error (log-log)
# ---------------------------------------------------------------------------
def plot_rk4_convergence(
    data: dict[str, dict[str, list[float]]], output_path: Path
) -> None:
    rk4 = data["rk4"]
    dts = rk4["dt_s"]
    pos_errs = rk4["final_position_error_m"]

    fig, ax = make_figure()
    ax.loglog(dts, pos_errs, "o-", color=BLUE, markersize=8, linewidth=2, label="RK4")

    # Overlay ideal 4th-order slope reference
    if len(dts) >= 2 and pos_errs[0] > 0:
        ref = [pos_errs[0] * ((dt / dts[0]) ** 4) for dt in dts]
        ax.loglog(
            dts,
            ref,
            "--",
            color=MUTED,
            linewidth=1.5,
            alpha=0.7,
            label=r"Ideal $O(h^4)$ reference",
        )

    ax.set_xlabel("Timestep $\\Delta t$ (s)")
    ax.set_ylabel("Final position closure error (m)")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    add_title(
        fig,
        "RK4 convergence: timestep vs position error",
        "500 km circular orbit · one analytical period · log-log scale",
    )
    save_figure(fig, output_path)


# ---------------------------------------------------------------------------
# Plot B — RK4 timestep vs energy drift (log-log)
# ---------------------------------------------------------------------------
def plot_rk4_energy_convergence(
    data: dict[str, dict[str, list[float]]], output_path: Path
) -> None:
    rk4 = data["rk4"]
    dts = rk4["dt_s"]
    energy_errs = rk4["max_relative_energy_error"]

    fig, ax = make_figure()
    ax.loglog(dts, energy_errs, "s-", color=BLUE, markersize=8, linewidth=2, label="RK4")

    ax.set_xlabel("Timestep $\\Delta t$ (s)")
    ax.set_ylabel("Maximum relative energy drift")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    add_title(
        fig,
        "RK4 convergence: timestep vs energy drift",
        "500 km circular orbit · one analytical period · log-log scale",
    )
    save_figure(fig, output_path)


# ---------------------------------------------------------------------------
# Plot C — RK4 timestep vs angular-momentum drift (log-log)
# ---------------------------------------------------------------------------
def plot_rk4_angular_momentum_convergence(
    data: dict[str, dict[str, list[float]]], output_path: Path
) -> None:
    rk4 = data["rk4"]
    dts = rk4["dt_s"]
    h_errs = rk4["max_relative_angular_momentum_error"]

    fig, ax = make_figure()
    ax.loglog(dts, h_errs, "D-", color=BLUE, markersize=8, linewidth=2, label="RK4")

    ax.set_xlabel("Timestep $\\Delta t$ (s)")
    ax.set_ylabel("Maximum relative angular-momentum drift")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    add_title(
        fig,
        "RK4 convergence: timestep vs angular-momentum drift",
        "500 km circular orbit · one analytical period · log-log scale",
    )
    save_figure(fig, output_path)


# ---------------------------------------------------------------------------
# Plot D — Euler vs RK4 error scaling (log-log)
# ---------------------------------------------------------------------------
def plot_euler_vs_rk4_convergence(
    data: dict[str, dict[str, list[float]]], output_path: Path
) -> None:
    fig, ax = make_figure()

    # RK4
    rk4 = data["rk4"]
    ax.loglog(
        rk4["dt_s"],
        rk4["final_position_error_m"],
        "o-",
        color=BLUE,
        markersize=8,
        linewidth=2,
        label="RK4",
    )

    # Euler
    euler = data["euler"]
    ax.loglog(
        euler["dt_s"],
        euler["final_position_error_m"],
        "s-",
        color=ORANGE,
        markersize=8,
        linewidth=2,
        label="Forward Euler",
    )

    # Ideal slope references
    if len(rk4["dt_s"]) >= 2 and rk4["final_position_error_m"][0] > 0:
        ref4 = [
            rk4["final_position_error_m"][0] * (dt / rk4["dt_s"][0]) ** 4
            for dt in rk4["dt_s"]
        ]
        ax.loglog(
            rk4["dt_s"],
            ref4,
            "--",
            color=BLUE,
            linewidth=1,
            alpha=0.4,
            label=r"Ideal $O(h^4)$",
        )

    if len(euler["dt_s"]) >= 2 and euler["final_position_error_m"][0] > 0:
        ref1 = [
            euler["final_position_error_m"][0] * (dt / euler["dt_s"][0]) ** 1
            for dt in euler["dt_s"]
        ]
        ax.loglog(
            euler["dt_s"],
            ref1,
            "--",
            color=ORANGE,
            linewidth=1,
            alpha=0.4,
            label=r"Ideal $O(h^1)$",
        )

    ax.set_xlabel("Timestep $\\Delta t$ (s)")
    ax.set_ylabel("Final position closure error (m)")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    add_title(
        fig,
        "Euler vs RK4: convergence comparison",
        "500 km circular orbit · one analytical period · log-log scale",
    )
    save_figure(fig, output_path)


# ---------------------------------------------------------------------------
# Plot E — Phase error comparison
# ---------------------------------------------------------------------------
def plot_phase_error(
    data: dict[str, dict[str, list[float]]], output_path: Path
) -> None:
    fig, ax = make_figure()

    rk4 = data["rk4"]
    ax.semilogx(
        rk4["dt_s"],
        rk4["phase_error_rad"],
        "o-",
        color=BLUE,
        markersize=8,
        linewidth=2,
        label="RK4",
    )

    euler = data["euler"]
    ax.semilogx(
        euler["dt_s"],
        euler["phase_error_rad"],
        "s-",
        color=ORANGE,
        markersize=8,
        linewidth=2,
        label="Forward Euler",
    )

    ax.set_xlabel("Timestep $\\Delta t$ (s)")
    ax.set_ylabel("Maximum phase error (rad)")
    ax.legend(frameon=False, loc="upper right")
    style_axes(ax)
    add_title(
        fig,
        "Phase error: Euler vs RK4",
        "500 km circular orbit · one analytical period",
    )
    save_figure(fig, output_path)


# ---------------------------------------------------------------------------
# Summary table
# ---------------------------------------------------------------------------
def print_summary(data: dict[str, dict[str, list[float]]]) -> None:
    """Print a convergence summary to stdout."""
    print("\n" + "=" * 80)
    print("M05 CONVERGENCE SUMMARY")
    print("=" * 80)

    for method in ("rk4", "euler"):
        d = data[method]
        label = "RK4" if method == "rk4" else "Euler"
        print(f"\n--- {label} ---")
        print(
            f"{'dt (s)':>8} | {'pos_err (m)':>16} | {'dE/E':>14} | {'dh/h':>14} | {'phase_err':>14}"
        )
        print("-" * 80)
        for i in range(len(d["dt_s"])):
            print(
                f"{d['dt_s'][i]:8.1f} | {d['final_position_error_m'][i]:16.6e} | "
                f"{d['max_relative_energy_error'][i]:14.6e} | "
                f"{d['max_relative_angular_momentum_error'][i]:14.6e} | "
                f"{d['phase_error_rad'][i]:14.6e}"
            )

        # Convergence orders
        if len(d["dt_s"]) >= 2:
            print("\nEmpirical convergence orders (position error):")
            for i in range(len(d["dt_s"]) - 1):
                h1 = d["dt_s"][i]
                h2 = d["dt_s"][i + 1]
                e1 = d["final_position_error_m"][i]
                e2 = d["final_position_error_m"][i + 1]
                p = empirical_order(h1, e1, h2, e2)
                print(f"  {h1}s -> {h2}s: p = {p:.4f}")

    print("\n" + "=" * 80)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="M05 propagator validation analysis and convergence plots."
    )
    parser.add_argument("convergence_csv", type=Path)
    parser.add_argument("output_directory", type=Path)
    return parser


def main(argv: Sequence[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    args.output_directory.mkdir(parents=True, exist_ok=True)

    rows = read_csv(args.convergence_csv)
    data = group_by_integrator(rows)

    # Print summary table
    print_summary(data)

    # Generate plots
    outputs: dict[str, Path] = {
        "rk4_convergence": args.output_directory / "m05_rk4_convergence.png",
        "rk4_energy": args.output_directory / "m05_rk4_energy_convergence.png",
        "rk4_angular_momentum": args.output_directory
        / "m05_rk4_angular_momentum_convergence.png",
        "euler_vs_rk4": args.output_directory / "m05_euler_vs_rk4_convergence.png",
        "phase_error": args.output_directory / "m05_phase_error.png",
    }

    plot_rk4_convergence(data, outputs["rk4_convergence"])
    plot_rk4_energy_convergence(data, outputs["rk4_energy"])
    plot_rk4_angular_momentum_convergence(data, outputs["rk4_angular_momentum"])
    plot_euler_vs_rk4_convergence(data, outputs["euler_vs_rk4"])
    plot_phase_error(data, outputs["phase_error"])

    for name, path in outputs.items():
        print(f"{name}={path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())