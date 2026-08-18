"""Plotting and visualization tools for M08 Attitude Representation & Quaternion Mathematics."""

from __future__ import annotations

import math
from pathlib import Path

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401


def quat_from_axis_angle(ax: tuple[float, float, float], angle_rad: float) -> tuple[float, float, float, float]:
    n = math.sqrt(ax[0] ** 2 + ax[1] ** 2 + ax[2] ** 2)
    ux, uy, uz = ax[0] / n, ax[1] / n, ax[2] / n
    half = 0.5 * angle_rad
    return (math.cos(half), math.sin(half) * ux, math.sin(half) * uy, math.sin(half) * uz)


def quat_rotate(q: tuple[float, float, float, float], v: tuple[float, float, float]) -> tuple[float, float, float]:
    w, x, y, z = q
    vx, vy, vz = v
    # u x v
    cx = y * vz - z * vy
    cy = z * vx - x * vz
    cz = x * vy - y * vx
    # u x (u x v)
    ccx = y * cz - z * cy
    ccy = z * cx - x * cz
    ccz = x * cy - y * cx
    return (
        vx + 2.0 * w * cx + 2.0 * ccx,
        vy + 2.0 * w * cy + 2.0 * ccy,
        vz + 2.0 * w * cz + 2.0 * ccz,
    )


def quat_mult(q1: tuple[float, float, float, float], q2: tuple[float, float, float, float]) -> tuple[float, float, float, float]:
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return (
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    )


def quat_to_dcm(q: tuple[float, float, float, float]) -> list[list[float]]:
    w, x, y, z = q
    return [
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
        [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
        [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)],
    ]


def plot_body_axes_in_eci(output_path: Path) -> None:
    """Plot A: Reference frame (ECI) vs rotated Spacecraft Body Frame."""
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Reference ECI axes (dashed)
    ax.quiver(0, 0, 0, 1, 0, 0, color="gray", linestyle="--", linewidth=1.5, label="ECI X (Reference)")
    ax.quiver(0, 0, 0, 0, 1, 0, color="gray", linestyle="--", linewidth=1.5, label="ECI Y (Reference)")
    ax.quiver(0, 0, 0, 0, 0, 1, color="gray", linestyle="--", linewidth=1.5, label="ECI Z (Reference)")

    # 45 deg rotation about [1, 1, 1]
    q = quat_from_axis_angle((1.0, 1.0, 1.0), math.radians(45))
    bx = quat_rotate(q, (1.0, 0.0, 0.0))
    by = quat_rotate(q, (0.0, 1.0, 0.0))
    bz = quat_rotate(q, (0.0, 0.0, 1.0))

    # Rotation Axis
    ax.quiver(0, 0, 0, 1 / math.sqrt(3), 1 / math.sqrt(3), 1 / math.sqrt(3), color="black", linewidth=2.5, label="Rotation Axis [1, 1, 1]")

    # Body axes (solid primary colors)
    ax.quiver(0, 0, 0, bx[0], bx[1], bx[2], color="red", linewidth=2.5, label="Body X (Rotated)")
    ax.quiver(0, 0, 0, by[0], by[1], by[2], color="green", linewidth=2.5, label="Body Y (Rotated)")
    ax.quiver(0, 0, 0, bz[0], bz[1], bz[2], color="blue", linewidth=2.5, label="Body Z (Rotated)")

    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_zlim([-1.2, 1.2])
    ax.set_xlabel("X (ECI)")
    ax.set_ylabel("Y (ECI)")
    ax.set_zlabel("Z (ECI)")
    ax.set_title("AstraDock M08 — Spacecraft Body Frame in ECI\n(Rotation: 45° about [1, 1, 1])", fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1.0), fontsize=9)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot A: {output_path}")


def plot_quaternion_dcm_consistency(output_path: Path) -> None:
    """Plot B: Visual & numerical consistency between Quaternion and DCM rotation."""
    fig = plt.figure(figsize=(10, 5))

    q = quat_from_axis_angle((1.0, 2.0, 3.0), math.radians(60))
    dcm = quat_to_dcm(q)

    # Subplot 1: Quaternion Rotated Triad
    ax1 = fig.add_subplot(121, projection="3d")
    bx_q = quat_rotate(q, (1, 0, 0))
    by_q = quat_rotate(q, (0, 1, 0))
    bz_q = quat_rotate(q, (0, 0, 1))

    ax1.quiver(0, 0, 0, bx_q[0], bx_q[1], bx_q[2], color="red", linewidth=2.5, label="X_q")
    ax1.quiver(0, 0, 0, by_q[0], by_q[1], by_q[2], color="green", linewidth=2.5, label="Y_q")
    ax1.quiver(0, 0, 0, bz_q[0], bz_q[1], bz_q[2], color="blue", linewidth=2.5, label="Z_q")
    ax1.set_title("Quaternion Rotation ($v' = q \\otimes v \\otimes q^*$)", fontsize=11, fontweight="bold")
    ax1.set_xlim([-1.1, 1.1])
    ax1.set_ylim([-1.1, 1.1])
    ax1.set_zlim([-1.1, 1.1])
    ax1.legend(loc="lower right")

    # Subplot 2: DCM Rotated Triad
    ax2 = fig.add_subplot(122, projection="3d")
    bx_dcm = (dcm[0][0], dcm[1][0], dcm[2][0])
    by_dcm = (dcm[0][1], dcm[1][1], dcm[2][1])
    bz_dcm = (dcm[0][2], dcm[1][2], dcm[2][2])

    ax2.quiver(0, 0, 0, bx_dcm[0], bx_dcm[1], bx_dcm[2], color="red", linestyle="--", linewidth=2.5, label="X_dcm")
    ax2.quiver(0, 0, 0, by_dcm[0], by_dcm[1], by_dcm[2], color="green", linestyle="--", linewidth=2.5, label="Y_dcm")
    ax2.quiver(0, 0, 0, bz_dcm[0], bz_dcm[1], bz_dcm[2], color="blue", linestyle="--", linewidth=2.5, label="Z_dcm")
    ax2.set_title("DCM Rotation ($v' = C(q) v$)", fontsize=11, fontweight="bold")
    ax2.set_xlim([-1.1, 1.1])
    ax2.set_ylim([-1.1, 1.1])
    ax2.set_zlim([-1.1, 1.1])
    ax2.legend(loc="lower right")

    fig.suptitle("AstraDock M08 — Quaternion vs Direction Cosine Matrix Exact Consistency", fontsize=12, fontweight="bold")
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot B: {output_path}")


def plot_composition(output_path: Path) -> None:
    """Plot C: Sequential rotations and composite attitude orientation."""
    fig = plt.figure(figsize=(9, 8))
    ax = fig.add_subplot(111, projection="3d")

    # Rotation 1: 45 deg about +Z
    q1 = quat_from_axis_angle((0, 0, 1), math.radians(45))
    # Rotation 2: 45 deg about +X'
    q2 = quat_from_axis_angle((1, 0, 0), math.radians(45))
    # Composite: q3 = q1 * q2
    q3 = quat_mult(q1, q2)

    # Initial frame (dashed)
    ax.quiver(0, 0, 0, 1, 0, 0, color="gray", linestyle=":", label="Initial Frame")
    ax.quiver(0, 0, 0, 0, 1, 0, color="gray", linestyle=":")
    ax.quiver(0, 0, 0, 0, 0, 1, color="gray", linestyle=":")

    # After R1
    r1_x = quat_rotate(q1, (1, 0, 0))
    r1_y = quat_rotate(q1, (0, 1, 0))
    r1_z = quat_rotate(q1, (0, 0, 1))
    ax.quiver(0, 0, 0, r1_x[0], r1_x[1], r1_x[2], color="orange", linestyle="--", label="After Rotation 1 (45° Z)")
    ax.quiver(0, 0, 0, r1_y[0], r1_y[1], r1_y[2], color="orange", linestyle="--")
    ax.quiver(0, 0, 0, r1_z[0], r1_z[1], r1_z[2], color="orange", linestyle="--")

    # After Composite R3 = R1 * R2
    r3_x = quat_rotate(q3, (1, 0, 0))
    r3_y = quat_rotate(q3, (0, 1, 0))
    r3_z = quat_rotate(q3, (0, 0, 1))
    ax.quiver(0, 0, 0, r3_x[0], r3_x[1], r3_x[2], color="red", linewidth=2.5, label="Composite Body X")
    ax.quiver(0, 0, 0, r3_y[0], r3_y[1], r3_y[2], color="green", linewidth=2.5, label="Composite Body Y")
    ax.quiver(0, 0, 0, r3_z[0], r3_z[1], r3_z[2], color="blue", linewidth=2.5, label="Composite Body Z")

    ax.set_xlim([-1.2, 1.2])
    ax.set_ylim([-1.2, 1.2])
    ax.set_zlim([-1.2, 1.2])
    ax.set_xlabel("X")
    ax.set_ylabel("Y")
    ax.set_zlabel("Z")
    ax.set_title("AstraDock M08 — Rotation Composition ($q_{composite} = q_1 \\otimes q_2$)", fontsize=12, fontweight="bold")
    ax.legend(loc="upper left", bbox_to_anchor=(1.05, 1.0), fontsize=9)
    plt.tight_layout()

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)
    print(f"Saved Plot C: {output_path}")


def main() -> None:
    figures_dir = Path("artifacts/figures")
    plot_body_axes_in_eci(figures_dir / "m08_body_frame_in_eci.png")
    plot_quaternion_dcm_consistency(figures_dir / "m08_quaternion_dcm_consistency.png")
    plot_composition(figures_dir / "m08_rotation_composition.png")
    print("All M08 attitude figures generated successfully.")


if __name__ == "__main__":
    main()
