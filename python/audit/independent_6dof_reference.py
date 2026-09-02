#!/usr/bin/env python3
"""Independent pure-Python 6-DOF reference oracle for AstraDock Milestone M10.

Cross-verifies the C++ integrated 6-DOF simulation against an independent
implementation of two-body gravity orbital propagation and rigid-body Euler dynamics.
"""

from pathlib import Path

import numpy as np
import pandas as pd

MU_WGS84 = 3.986004418e14  # m^3/s^2
R_EARTH = 6378137.0         # m


def quat_mult(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Hamilton product for scalar-first quaternions q = [w, x, y, z]."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2,
    ], dtype=np.float64)


def quat_conjugate(q: np.ndarray) -> np.ndarray:
    """Quaternion conjugate q* = [w, -x, -y, -z]."""
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=np.float64)


def quat_angle_err_rad(q_ref: np.ndarray, q_test: np.ndarray) -> float:
    """Shortest physical orientation angle between two unit quaternions."""
    q_err = quat_mult(quat_conjugate(q_ref), q_test)
    clamped_w = np.clip(np.abs(q_err[0]), 0.0, 1.0)
    return 2.0 * float(np.arccos(clamped_w))


def six_dof_deriv(
    r: np.ndarray,
    v: np.ndarray,
    q: np.ndarray,
    omega: np.ndarray,
    I_diag: np.ndarray,
    mu: float,
    tau: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Evaluates the 13-component first-order ODE derivatives."""
    # Translation
    r_norm = np.linalg.norm(r)
    a_grav = -mu * r / (r_norm ** 3)

    # Rotation kinematics: dq/dt = 0.5 * q * [0, omega]
    omega_quat = np.array([0.0, omega[0], omega[1], omega[2]], dtype=np.float64)
    dq_dt = 0.5 * quat_mult(q, omega_quat)

    # Rotation dynamics: I * d(omega)/dt + omega x (I * omega) = tau
    I_omega = I_diag * omega
    omega_cross_I_omega = np.cross(omega, I_omega)
    domega_dt = (tau - omega_cross_I_omega) / I_diag

    return v, a_grav, dq_dt, domega_dt


def rk4_step_6dof(
    r: np.ndarray,
    v: np.ndarray,
    q: np.ndarray,
    omega: np.ndarray,
    dt: float,
    I_diag: np.ndarray,
    mu: float,
    tau: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Advances 6-DOF state by one RK4 step with post-step quaternion normalization."""
    # Stage 1
    dr1, dv1, dq1, dw1 = six_dof_deriv(r, v, q, omega, I_diag, mu, tau)

    # Stage 2
    r2 = r + 0.5 * dt * dr1
    v2 = v + 0.5 * dt * dv1
    q2 = q + 0.5 * dt * dq1
    w2 = omega + 0.5 * dt * dw1
    dr2, dv2, dq2, dw2 = six_dof_deriv(r2, v2, q2, w2, I_diag, mu, tau)

    # Stage 3
    r3 = r + 0.5 * dt * dr2
    v3 = v + 0.5 * dt * dv2
    q3 = q + 0.5 * dt * dq2
    w3 = omega + 0.5 * dt * dw2
    dr3, dv3, dq3, dw3 = six_dof_deriv(r3, v3, q3, w3, I_diag, mu, tau)

    # Stage 4
    r4 = r + dt * dr3
    v4 = v + dt * dv3
    q4 = q + dt * dq3
    w4 = omega + dt * dw3
    dr4, dv4, dq4, dw4 = six_dof_deriv(r4, v4, q4, w4, I_diag, mu, tau)

    # Combine weighted slopes
    r_next = r + (dt / 6.0) * (dr1 + 2.0 * dr2 + 2.0 * dr3 + dr4)
    v_next = v + (dt / 6.0) * (dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4)
    q_next = q + (dt / 6.0) * (dq1 + 2.0 * dq2 + 2.0 * dq3 + dq4)
    w_next = omega + (dt / 6.0) * (dw1 + 2.0 * dw2 + 2.0 * dw3 + dw4)

    # Reproject quaternion onto S^3
    q_next = q_next / np.linalg.norm(q_next)

    return r_next, v_next, q_next, w_next


def audit_canonical_scenario(csv_path: Path) -> bool:
    """Audits the canonical 500 km orbit + asymmetric tumbling simulation."""
    print("=" * 70)
    print(f"Auditing Canonical 6-DOF Telemetry: {csv_path.name}")
    print("=" * 70)

    if not csv_path.exists():
        print(f"ERROR: CSV file not found at {csv_path}")
        return False

    df = pd.read_csv(csv_path)
    times = df["time_s"].to_numpy()
    n_samples = len(times)
    print(f"Loaded {n_samples} telemetry samples across duration {times[-1]:.2f} s")

    # Initial state from CSV
    r_py = np.array([df["position_eci_x_m"].iloc[0], df["position_eci_y_m"].iloc[0], df["position_eci_z_m"].iloc[0]])
    v_py = np.array([df["velocity_eci_x_mps"].iloc[0], df["velocity_eci_y_mps"].iloc[0], df["velocity_eci_z_mps"].iloc[0]])
    q_py = np.array([df["quaternion_w"].iloc[0], df["quaternion_x"].iloc[0], df["quaternion_y"].iloc[0], df["quaternion_z"].iloc[0]])
    w_py = np.array([df["angular_velocity_body_x_rad_s"].iloc[0], df["angular_velocity_body_y_rad_s"].iloc[0], df["angular_velocity_body_z_rad_s"].iloc[0]])

    I_diag = np.array([10.0, 20.0, 30.0], dtype=np.float64)
    tau = np.zeros(3, dtype=np.float64)

    max_pos_err = 0.0
    max_vel_err = 0.0
    max_q_err = 0.0
    max_w_err = 0.0

    for i in range(1, n_samples):
        dt = times[i] - times[i - 1]
        r_py, v_py, q_py, w_py = rk4_step_6dof(r_py, v_py, q_py, w_py, dt, I_diag, MU_WGS84, tau)

        r_cpp = np.array([df["position_eci_x_m"].iloc[i], df["position_eci_y_m"].iloc[i], df["position_eci_z_m"].iloc[i]])
        v_cpp = np.array([df["velocity_eci_x_mps"].iloc[i], df["velocity_eci_y_mps"].iloc[i], df["velocity_eci_z_mps"].iloc[i]])
        q_cpp = np.array([df["quaternion_w"].iloc[i], df["quaternion_x"].iloc[i], df["quaternion_y"].iloc[i], df["quaternion_z"].iloc[i]])
        w_cpp = np.array([df["angular_velocity_body_x_rad_s"].iloc[i], df["angular_velocity_body_y_rad_s"].iloc[i], df["angular_velocity_body_z_rad_s"].iloc[i]])

        pos_err = np.linalg.norm(r_py - r_cpp)
        vel_err = np.linalg.norm(v_py - v_cpp)
        q_err = quat_angle_err_rad(q_py, q_cpp)
        w_err = np.linalg.norm(w_py - w_cpp)

        max_pos_err = max(max_pos_err, pos_err)
        max_vel_err = max(max_vel_err, vel_err)
        max_q_err = max(max_q_err, q_err)
        max_w_err = max(max_w_err, w_err)

    print(f"  Max Position Error (C++ vs Python):          {max_pos_err:.4e} m")
    print(f"  Max Velocity Error (C++ vs Python):          {max_vel_err:.4e} m/s")
    print(f"  Max Orientation Error (C++ vs Python):       {max_q_err:.4e} rad ({np.degrees(max_q_err):.4e} deg)")
    print(f"  Max Angular Velocity Error (C++ vs Python):  {max_w_err:.4e} rad/s")

    # Invariant audits from CSV
    init_orb_e = df["specific_orbital_energy_m2_s2"].iloc[0]
    final_orb_e = df["specific_orbital_energy_m2_s2"].iloc[-1]
    orb_e_drift = abs(final_orb_e - init_orb_e) / abs(init_orb_e)

    init_rot_e = df["rotational_kinetic_energy_J"].iloc[0]
    final_rot_e = df["rotational_kinetic_energy_J"].iloc[-1]
    rot_e_drift = abs(final_rot_e - init_rot_e) / init_rot_e

    max_q_norm_err = np.max(np.abs(df["quaternion_norm"].to_numpy() - 1.0))

    print(f"  Orbital Specific Energy Drift:               {orb_e_drift:.4e}")
    print(f"  Rotational Kinetic Energy Drift:             {rot_e_drift:.4e}")
    print(f"  Max Quaternion Norm Error |norm(q) - 1|:     {max_q_norm_err:.4e}")

    passed = (
        max_pos_err < 1.0e-3
        and max_vel_err < 1.0e-6
        and max_q_err < 1.0e-6
        and max_w_err < 1.0e-6
        and orb_e_drift < 1.0e-10
        and rot_e_drift < 1.0e-5
        and max_q_norm_err < 1.0e-14
    )

    print(f"Canonical Scenario Audit Result: {'PASSED' if passed else 'FAILED'}\n")
    return passed


def quat_to_dcm_matrix(q: np.ndarray) -> np.ndarray:
    """Converts a scalar-first unit quaternion into a 3x3 DCM (C_ECI_B)."""
    w, x, y, z = q
    xx, yy, zz = x * x, y * y, z * z
    wx, wy, wz = w * x, w * y, w * z
    xy, xz, yz = x * y, x * z, y * z

    return np.array([
        [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz),       2.0 * (xz + wy)],
        [2.0 * (xy + wz),       1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
        [2.0 * (xz - wy),       2.0 * (yz + wx),       1.0 - 2.0 * (xx + yy)],
    ], dtype=np.float64)


def audit_m10_refinement_study() -> tuple[bool, dict]:
    """Performs Stage A M10 timestep refinement study and convergence order measurement."""
    print("=" * 70)
    print("STAGE A: M10 Numerical Entry Gate — Timestep Refinement Study")
    print("=" * 70)

    # Initial state (500 km circular orbit + 45 deg Z attitude + tumbling rates)
    r_orbit = R_EARTH + 500.0e3
    v_circ = np.sqrt(MU_WGS84 / r_orbit)
    t_period = 2.0 * np.pi * np.sqrt((r_orbit ** 3) / MU_WGS84)

    r0 = np.array([r_orbit, 0.0, 0.0], dtype=np.float64)
    v0 = np.array([0.0, v_circ, 0.0], dtype=np.float64)
    q0 = np.array([np.cos(np.pi / 8.0), 0.0, 0.0, np.sin(np.pi / 8.0)], dtype=np.float64)
    w0 = np.array([0.05, 0.08, 0.02], dtype=np.float64)
    I_diag = np.array([10.0, 20.0, 30.0], dtype=np.float64)
    tau = np.zeros(3, dtype=np.float64)

    # Initial invariants
    e_rot0 = 0.5 * np.sum(I_diag * (w0 ** 2))
    h_body0 = I_diag * w0
    c0 = quat_to_dcm_matrix(q0)
    h_inertial0 = c0 @ h_body0
    h_inertial0_norm = np.linalg.norm(h_inertial0)

    def run_sim(dt: float, renormalize: bool = True):
        steps = int(np.round(t_period / dt))
        r, v, q, w = r0.copy(), v0.copy(), q0.copy(), w0.copy()
        for _ in range(steps):
            if renormalize:
                r, v, q, w = rk4_step_6dof(r, v, q, w, dt, I_diag, MU_WGS84, tau)
            else:
                # Unnormalized RK4
                dr1, dv1, dq1, dw1 = six_dof_deriv(r, v, q, w, I_diag, MU_WGS84, tau)
                r2, v2, q2, w2 = r + 0.5 * dt * dr1, v + 0.5 * dt * dv1, q + 0.5 * dt * dq1, w + 0.5 * dt * dw1
                dr2, dv2, dq2, dw2 = six_dof_deriv(r2, v2, q2, w2, I_diag, MU_WGS84, tau)
                r3, v3, q3, w3 = r + 0.5 * dt * dr2, v + 0.5 * dt * dv2, q + 0.5 * dt * dq2, w + 0.5 * dt * dw2
                dr3, dv3, dq3, dw3 = six_dof_deriv(r3, v3, q3, w3, I_diag, MU_WGS84, tau)
                r4, v4, q4, w4 = r + dt * dr3, v + dt * dv3, q + dt * dq3, w + dt * dw3
                dr4, dv4, dq4, dw4 = six_dof_deriv(r4, v4, q4, w4, I_diag, MU_WGS84, tau)
                r = r + (dt / 6.0) * (dr1 + 2.0 * dr2 + 2.0 * dr3 + dr4)
                v = v + (dt / 6.0) * (dv1 + 2.0 * dv2 + 2.0 * dv3 + dv4)
                q = q + (dt / 6.0) * (dq1 + 2.0 * dq2 + 2.0 * dq3 + dq4)
                w = w + (dt / 6.0) * (dw1 + 2.0 * dw2 + 2.0 * dw3 + dw4)
        return r, v, q, w

    # High fidelity reference run
    dt_ref = 0.05
    r_ref, v_ref, q_ref, w_ref = run_sim(dt_ref, renormalize=True)

    dt_list = [1.0, 0.5, 0.25]
    results = {}

    for dt in dt_list:
        r, v, q, w = run_sim(dt, renormalize=True)
        w_err = np.linalg.norm(w - w_ref)
        q_err = quat_angle_err_rad(q_ref, q)
        pos_err = np.linalg.norm(r - r_ref)

        e_rot = 0.5 * np.sum(I_diag * (w ** 2))
        e_rot_drift = abs(e_rot - e_rot0) / e_rot0

        c_mat = quat_to_dcm_matrix(q)
        h_inertial = c_mat @ (I_diag * w)
        h_drift = np.linalg.norm(h_inertial - h_inertial0) / h_inertial0_norm

        results[dt] = {
            "w_err": w_err,
            "q_err": q_err,
            "pos_err": pos_err,
            "e_rot_drift": e_rot_drift,
            "h_drift": h_drift,
            "q_norm_err": abs(np.linalg.norm(q) - 1.0),
        }

    # Calculate convergence orders
    p_w_12 = np.log(results[1.0]["w_err"] / results[0.5]["w_err"]) / np.log(2.0)
    p_w_23 = np.log(results[0.5]["w_err"] / results[0.25]["w_err"]) / np.log(2.0)
    p_q_12 = np.log(results[1.0]["q_err"] / results[0.5]["q_err"]) / np.log(2.0)
    p_q_23 = np.log(results[0.5]["q_err"] / results[0.25]["q_err"]) / np.log(2.0)

    # Invariant reduction factor
    ratio_e_12 = results[1.0]["e_rot_drift"] / results[0.5]["e_rot_drift"]
    ratio_e_23 = results[0.5]["e_rot_drift"] / results[0.25]["e_rot_drift"]
    ratio_h_12 = results[1.0]["h_drift"] / results[0.5]["h_drift"]
    ratio_h_23 = results[0.5]["h_drift"] / results[0.25]["h_drift"]

    print(f"Refinement Table (Duration = {t_period:.1f} s, 1 Orbit):")
    print(f"{'dt (s)':<8} | {'E_rot Drift':<14} | {'H_inertial Drift':<16} | {'w Error (rad/s)':<16} | {'Att Error (rad)':<16}")
    print("-" * 78)
    for dt in dt_list:
        res = results[dt]
        print(f"{dt:<8.2f} | {res['e_rot_drift']:<14.4e} | {res['h_drift']:<16.4e} | {res['w_err']:<16.4e} | {res['q_err']:<16.4e}")

    print("\nConvergence Order & Invariant Reduction Ratios:")
    print(f"  Angular Velocity Order p(1.0 -> 0.5):      {p_w_12:.2f}")
    print(f"  Angular Velocity Order p(0.5 -> 0.25):     {p_w_23:.2f}")
    print(f"  Attitude Orientation Order p(1.0 -> 0.5):  {p_q_12:.2f}")
    print(f"  Attitude Orientation Order p(0.5 -> 0.25): {p_q_23:.2f}")
    print(f"  Rotational Energy Drift Reduction Ratio:   {ratio_e_12:.2f} (dt/2), {ratio_e_23:.2f} (dt/4)")
    print(f"  Inertial Momentum Drift Reduction Ratio:   {ratio_h_12:.2f} (dt/2), {ratio_h_23:.2f} (dt/4)")

    # Gate decision criteria:
    # 1. Orders p ~ 4 (within 3.8 to 4.8)
    # 2. Invariant reduction ratio: Momentum ~ 16 (2^4), Energy >= 16 (2^5 ~ 32 due to quadratic E ~ w^2)
    # 3. Normalized quaternion error == 0
    passed = (
        3.8 <= p_w_12 <= 4.8
        and 3.8 <= p_w_23 <= 4.8
        and 3.8 <= p_q_12 <= 4.8
        and 3.8 <= p_q_23 <= 4.8
        and ratio_e_12 >= 16.0
        and ratio_e_23 >= 16.0
        and 14.0 <= ratio_h_12 <= 18.0
        and 14.0 <= ratio_h_23 <= 18.0
    )

    gate_decision = "PASS" if passed else "BLOCK"
    print(f"\nM10 NUMERICAL GATE DECISION: {gate_decision}\n")

    return passed, results


def audit_constant_torque_scenario(csv_path: Path) -> bool:
    """Audits the constant torque verification run."""
    print("=" * 70)
    print(f"Auditing Constant Torque Telemetry: {csv_path.name}")
    print("=" * 70)

    if not csv_path.exists():
        print(f"ERROR: CSV file not found at {csv_path}")
        return False

    df = pd.read_csv(csv_path)
    t_final = df["time_s"].iloc[-1]
    w_z_final = df["angular_velocity_body_z_rad_s"].iloc[-1]

    # In six_dof_demo.cpp: I_zz = 30.0 kg*m^2, tau_z = 0.2, t = 200s => omega_z = (0.2 / 30.0) * 200 = 1.333333 rad/s
    expected_w_z = (0.2 / 30.0) * t_final
    diff_w_z = abs(w_z_final - expected_w_z)

    print(f"  Duration:                                    {t_final:.1f} s")
    print(f"  Final omega_z (simulated):                   {w_z_final:.6f} rad/s")
    print(f"  Final omega_z (analytical):                  {expected_w_z:.6f} rad/s")
    print(f"  Difference:                                  {diff_w_z:.4e} rad/s")

    passed = diff_w_z < 1.0e-10
    print(f"Constant Torque Scenario Audit Result: {'PASSED' if passed else 'FAILED'}\n")
    return passed


def main():
    root = Path(__file__).resolve().parent.parent.parent
    data_dir = root / "data"

    canonical_csv = data_dir / "six_dof_canonical_orbit.csv"
    torque_csv = data_dir / "six_dof_constant_torque.csv"

    pass_gate, _ = audit_m10_refinement_study()
    pass1 = audit_canonical_scenario(canonical_csv)
    pass2 = audit_constant_torque_scenario(torque_csv)

    print("=" * 70)
    print(f"OVERALL INDEPENDENT 6-DOF ORACLE AUDIT: {'ALL TESTS PASSED' if pass_gate and pass1 and pass2 else 'FAILED'}")
    print("=" * 70)

    if not (pass_gate and pass1 and pass2):
        raise SystemExit(1)


if __name__ == "__main__":
    main()
