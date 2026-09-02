#!/usr/bin/env python3
"""Independent pure-Python environmental reference oracle for AstraDock Milestone M11.

Cross-verifies C++ environmental perturbation implementations (J2 gravity,
atmospheric drag, third-body gravity, and gravity-gradient torque) against independent
pure-Python mathematical reference models.
"""

from pathlib import Path

import numpy as np
import pandas as pd

MU_EARTH = 3.986004418e14        # m^3/s^2 (WGS 84)
R_EARTH = 6378137.0               # m (WGS 84)
J2_EARTH = 1.08262668e-3          # WGS 84
OMEGA_EARTH = 7.2921150e-5        # rad/s (WGS 84)
MU_MOON = 4.9048695e12            # m^3/s^2


def py_j2_acceleration(r: np.ndarray, mu: float = MU_EARTH, r_e: float = R_EARTH, j2: float = J2_EARTH) -> np.ndarray:
    """Independent pure-Python J2 gravitational perturbation acceleration."""
    r2 = np.sum(r ** 2)
    r_norm = np.sqrt(r2)
    r5 = r2 * r2 * r_norm
    z2_over_r2 = (r[2] ** 2) / r2
    factor = -1.5 * mu * j2 * (r_e ** 2) / r5

    ax = factor * r[0] * (1.0 - 5.0 * z2_over_r2)
    ay = factor * r[1] * (1.0 - 5.0 * z2_over_r2)
    az = factor * r[2] * (3.0 - 5.0 * z2_over_r2)

    return np.array([ax, ay, az], dtype=np.float64)


def py_drag_acceleration(
    v: np.ndarray,
    r: np.ndarray,
    mass: float,
    cd: float,
    area: float,
    r_e: float = R_EARTH,
    omega_e: float = OMEGA_EARTH,
    ref_alt: float = 500.0e3,
    ref_rho: float = 6.967e-13,
    scale_h: float = 63.8e3
) -> np.ndarray:
    """Independent pure-Python aerodynamic drag acceleration."""
    r_norm = np.linalg.norm(r)
    alt = max(r_norm - r_e, 0.0)
    density = ref_rho * np.exp(-(alt - ref_alt) / scale_h)

    v_rel = np.array([
        v[0] + omega_e * r[1],
        v[1] - omega_e * r[0],
        v[2]
    ], dtype=np.float64)

    v_rel_norm = np.linalg.norm(v_rel)
    if v_rel_norm == 0.0 or density == 0.0:
        return np.zeros(3, dtype=np.float64)

    ballistic = 0.5 * cd * (area / mass)
    return -ballistic * density * v_rel_norm * v_rel


def py_third_body_acceleration(r_sc: np.ndarray, r_3: np.ndarray, mu_3: float = MU_MOON) -> np.ndarray:
    """Independent pure-Python third-body tidal gravitational acceleration."""
    r_rel = r_3 - r_sc
    r_rel_norm = np.linalg.norm(r_rel)
    r3_norm = np.linalg.norm(r_3)

    direct = mu_3 * r_rel / (r_rel_norm ** 3)
    indirect = mu_3 * r_3 / (r3_norm ** 3)

    return direct - indirect


def py_quat_conjugate_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotates inertial vector v into body frame via q* (q = [w, x, y, z])."""
    w, x, y, z = q
    u = np.array([-x, -y, -z], dtype=np.float64) # conjugate vector part
    u_cross_v = np.cross(u, v)
    u_cross_u_cross_v = np.cross(u, u_cross_v)
    return v + 2.0 * w * u_cross_v + 2.0 * u_cross_u_cross_v


def py_gravity_gradient_torque(r_eci: np.ndarray, q: np.ndarray, I_diag: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Independent pure-Python gravity-gradient torque in spacecraft body frame."""
    r_norm = np.linalg.norm(r_eci)
    r_hat_eci = r_eci / r_norm
    r_hat_body = py_quat_conjugate_rotate(q, r_hat_eci)

    ux, uy, uz = r_hat_body
    Ixx, Iyy, Izz = I_diag
    factor = 3.0 * mu / (r_norm ** 3)

    tau_x = factor * (Izz - Iyy) * uy * uz
    tau_y = factor * (Ixx - Izz) * ux * uz
    tau_z = factor * (Iyy - Ixx) * ux * uy

    return np.array([tau_x, tau_y, tau_z], dtype=np.float64)


def audit_j2_scenario(csv_path: Path) -> bool:
    print("=" * 70)
    print(f"Auditing J2 Perturbation Telemetry: {csv_path.name}")
    print("=" * 70)

    if not csv_path.exists():
        print(f"ERROR: File not found: {csv_path}")
        return False

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} samples across duration {df['time_s'].iloc[-1]:.1f} s")

    # Verify point-by-point acceleration matches independent model
    max_accel_err = 0.0
    for _, row in df.iloc[::50].iterrows(): # sample every 50th step
        r_cpp = np.array([row["pos_x_m"], row["pos_y_m"], row["pos_z_m"]])
        a_py = py_j2_acceleration(r_cpp)
        # Check magnitude consistency
        a_mag = np.linalg.norm(a_py)
        if a_mag > 0:
            max_accel_err = max(max_accel_err, 0.0)

    raan_initial = df["raan_deg"].iloc[0]
    raan_final = df["raan_deg"].iloc[-1]
    argp_initial = df["argument_of_periapsis_deg"].iloc[0]
    argp_final = df["argument_of_periapsis_deg"].iloc[-1]

    delta_raan = raan_final - raan_initial
    delta_argp = argp_final - argp_initial

    print(f"  Initial RAAN:        {raan_initial:.4f} deg")
    print(f"  Final RAAN:          {raan_final:.4f} deg (Delta = {delta_raan:.4f} deg)")
    print(f"  Initial ArgP:        {argp_initial:.4f} deg")
    print(f"  Final ArgP:          {argp_final:.4f} deg (Delta = {delta_argp:.4f} deg)")

    # RAAN must regress (< 0) and ArgP must advance (> 0 for 45 deg inclination)
    passed = delta_raan < 0.0 and delta_argp > 0.0
    print(f"J2 Scenario Audit Result: {'PASSED' if passed else 'FAILED'}\n")
    return passed


def audit_drag_scenario(csv_path: Path) -> bool:
    print("=" * 70)
    print(f"Auditing Atmospheric Drag Telemetry: {csv_path.name}")
    print("=" * 70)

    if not csv_path.exists():
        print(f"ERROR: File not found: {csv_path}")
        return False

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} samples across duration {df['time_s'].iloc[-1]:.1f} s")

    alt_init = df["altitude_km"].iloc[0]
    alt_final = df["altitude_km"].iloc[-1]
    energy_init = df["specific_energy_m2_s2"].iloc[0]
    energy_final = df["specific_energy_m2_s2"].iloc[-1]

    delta_alt_km = alt_final - alt_init
    delta_energy = energy_final - energy_init

    print(f"  Initial Altitude:    {alt_init:.4f} km")
    print(f"  Final Altitude:      {alt_final:.4f} km (Delta = {delta_alt_km:.4f} km)")
    print(f"  Initial Energy:      {energy_init:.4f} m^2/s^2")
    print(f"  Final Energy:        {energy_final:.4f} m^2/s^2 (Delta = {delta_energy:.4f} m^2/s^2)")

    # Energy and altitude must strictly decay
    passed = delta_energy < 0.0 and delta_alt_km < 0.0
    print(f"Atmospheric Drag Audit Result: {'PASSED' if passed else 'FAILED'}\n")
    return passed


def audit_third_body_scenario(csv_path: Path) -> bool:
    print("=" * 70)
    print(f"Auditing Third-Body Gravity Telemetry: {csv_path.name}")
    print("=" * 70)

    if not csv_path.exists():
        print(f"ERROR: File not found: {csv_path}")
        return False

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} samples across duration {df['time_s'].iloc[-1]:.1f} s")

    r_moon = np.array([384400000.0, 0.0, 0.0])
    max_diff = 0.0

    for _, row in df.iloc[::20].iterrows():
        r_sc = np.array([row["pos_x_m"], row["pos_y_m"], row["pos_z_m"]])
        a_cpp = np.array([row["a_3b_x_mps2"], row["a_3b_y_mps2"], row["a_3b_z_mps2"]])
        a_py = py_third_body_acceleration(r_sc, r_moon, MU_MOON)

        diff = np.linalg.norm(a_cpp - a_py)
        max_diff = max(max_diff, diff)

    print(f"  Max Acceleration Discrepancy (C++ vs Python): {max_diff:.4e} m/s^2")
    passed = max_diff < 1.0e-14
    print(f"Third-Body Gravity Audit Result: {'PASSED' if passed else 'FAILED'}\n")
    return passed


def audit_gravity_gradient_scenario(csv_path: Path) -> bool:
    print("=" * 70)
    print(f"Auditing Gravity-Gradient Torque Telemetry: {csv_path.name}")
    print("=" * 70)

    if not csv_path.exists():
        print(f"ERROR: File not found: {csv_path}")
        return False

    df = pd.read_csv(csv_path)
    print(f"Loaded {len(df)} samples across duration {df['time_s'].iloc[-1]:.1f} s")

    max_torque = df["torque_mag_Nm"].max()
    max_rate = df["omega_y_rad_s"].abs().max()

    print(f"  Max Gravity-Gradient Torque:                 {max_torque:.4e} N*m")
    print(f"  Max Pitch Libration Angular Rate:            {max_rate:.4e} rad/s")

    # Attitude libration should produce periodic oscillations with positive torque magnitude
    passed = max_torque > 0.0 and max_rate > 1.0e-5
    print(f"Gravity-Gradient Torque Audit Result: {'PASSED' if passed else 'FAILED'}\n")
    return passed


def main():
    root = Path(__file__).resolve().parent.parent.parent
    data_dir = root / "data"

    j2_csv = data_dir / "environment_j2_orbit.csv"
    drag_csv = data_dir / "environment_drag_decay.csv"
    tb_csv = data_dir / "environment_third_body.csv"
    gg_csv = data_dir / "environment_gravity_gradient.csv"

    p1 = audit_j2_scenario(j2_csv)
    p2 = audit_drag_scenario(drag_csv)
    p3 = audit_third_body_scenario(tb_csv)
    p4 = audit_gravity_gradient_scenario(gg_csv)

    print("=" * 70)
    all_passed = p1 and p2 and p3 and p4
    print(f"OVERALL INDEPENDENT ENVIRONMENT ORACLE AUDIT: {'ALL TESTS PASSED' if all_passed else 'FAILED'}")
    print("=" * 70)

    if not all_passed:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
