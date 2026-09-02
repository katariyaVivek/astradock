"""AstraDock -- Milestone M13B: Independent Python IMU-EKF Oracle.

Implements an IMU-aided translational EKF from scratch using only the
mathematical specification (not the C++ code structure) and compares its
estimates against the C++ production telemetry CSV.

Validation layers:
  1. Replays measurements from the C++ full-rate oracle-segment CSV
     (data/m13b_oracle_segment.csv; no C++ import).
  2. Runs its own inertial-navigation reconstruction
        a_I = C(q_hat) (f_m - b_a) + g(r_hat),
     its own RK4 propagator over that ODE, its own quaternion rotation,
     gravity Jacobian, Cholesky solve, Joseph update, and GNSS model.
  3. Verifies the BODY -> ECI rotation direction against a hand-computed DCM
     for a nontrivial attitude (catches conjugation errors).
  4. Verifies the gravity Jacobian against central finite differences.
  5. Asserts C++ vs Python estimate agreement within tolerance and
     recomputes NIS statistics.

Usage:
    python python/audit/independent_imu_ekf_reference.py [--segment PATH]
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

# ----------------------------------------------------------------------
# 1. Constants (must match tools/imu_ekf_demo.cpp exactly)
# ----------------------------------------------------------------------
MU = 3.986004418e14          # m^3/s^2
SIGMA_R = 10.0               # m   GNSS position noise
SIGMA_V = 0.05               # m/s GNSS velocity noise
SIGMA_F = 1.0e-3             # m/s^2 accelerometer noise (drives Q)
ACCEL_BIAS = (0.0, 0.0, 0.0) # configured bias removal (nominal run)
INIT_POS_ERR = (50_000.0, -30_000.0, 20_000.0)
INIT_VEL_ERR = (5.0, -3.0, 2.0)
P0_POS = 25_000.0 ** 2
P0_VEL = 5.0 ** 2


def main() -> None:
    data_dir = Path("data")
    segment_path = data_dir / "m13b_oracle_segment.csv"

    argv = sys.argv[1:]
    i = 0
    while i < len(argv):
        if argv[i] == "--segment" and i + 1 < len(argv):
            segment_path = Path(argv[i + 1])
            i += 2
        else:
            i += 1

    if not segment_path.exists():
        print(f"ERROR: {segment_path} not found")
        print("Run astradock_imu_ekf_demo first.")
        sys.exit(1)

    import numpy as np

    data = np.genfromtxt(segment_path, delimiter=",", names=True)
    n = len(data)
    print(f"Loaded {n} rows from {segment_path.name}")

    # ----------------------------------------------------------------------
    # 2. Core physics + EKF implementation (from specification, not C++)
    # ----------------------------------------------------------------------

    def quat_to_dcm(w, x, y, z):
        """C such that v_ECI = C @ v_BODY (standard textbook formula)."""
        return np.array([
            [1 - 2 * (y * y + z * z), 2 * (x * y - w * z), 2 * (x * z + w * y)],
            [2 * (x * y + w * z), 1 - 2 * (x * x + z * z), 2 * (y * z - w * x)],
            [2 * (x * z - w * y), 2 * (y * z + w * x), 1 - 2 * (x * x + y * y)],
        ])

    def two_body_accel(r):
        r_norm = np.linalg.norm(r)
        return -MU * r / r_norm**3

    def imu_accel(r, f_body, bias, q):
        c = quat_to_dcm(*q)
        return c @ (np.asarray(f_body) - np.asarray(bias)) + two_body_accel(r)

    def rk4_imu_step(r, v, dt, f_body, q):
        """RK4 across dr/dt = v, dv/dt = C(q)(f - b) + g(r), ZOH on f and q."""
        def deriv(ri, vi):
            return vi.copy(), imu_accel(ri, f_body, ACCEL_BIAS, q)

        k1r, k1v = deriv(r, v)
        k2r, k2v = deriv(r + 0.5 * dt * k1r, v + 0.5 * dt * k1v)
        k3r, k3v = deriv(r + 0.5 * dt * k2r, v + 0.5 * dt * k2v)
        k4r, k4v = deriv(r + dt * k3r, v + dt * k3v)
        r_new = r + dt / 6.0 * (k1r + 2 * k2r + 2 * k3r + k4r)
        v_new = v + dt / 6.0 * (k1v + 2 * k2v + 2 * k3v + k4v)
        return r_new, v_new

    def gravity_jacobian(r):
        r_norm = np.linalg.norm(r)
        rhat = r / r_norm
        return (-MU / r_norm**3) * (np.eye(3) - 3.0 * np.outer(rhat, rhat))

    def discrete_transition(r, dt):
        big_f = np.zeros((6, 6))
        big_f[0:3, 3:6] = np.eye(3)
        big_f[3:6, 0:3] = gravity_jacobian(r)
        return np.eye(6) + big_f * dt

    def process_noise(dt):
        variance = SIGMA_F**2
        q_mat = np.zeros((6, 6))
        for axis in range(3):
            q_mat[axis, axis] = variance * dt**3 / 3.0
            q_mat[axis + 3, axis + 3] = variance * dt
            q_mat[axis, axis + 3] = variance * dt**2 / 2.0
            q_mat[axis + 3, axis] = q_mat[axis, axis + 3]
        return q_mat

    def cholesky(a):
        size = a.shape[0]
        lower = np.zeros_like(a)
        for row in range(size):
            for col in range(row + 1):
                s = a[row, col] - np.dot(lower[row, :col], lower[col, :col])
                if row == col:
                    if s <= 0:
                        raise RuntimeError("Cholesky: not positive definite")
                    lower[row, col] = math.sqrt(s)
                else:
                    lower[row, col] = s / lower[col, col]
        return lower

    def solve_spd(a, b):
        lower = cholesky(a)
        size = len(b)
        y = np.zeros(size)
        for row in range(size):
            s = float(b[row])
            for col in range(row):
                s -= lower[row, col] * y[col]
            y[row] = s / lower[row, row]
        x = np.zeros(size)
        for row in range(size - 1, -1, -1):
            s = float(y[row])
            for col in range(row + 1, size):
                s -= lower[col, row] * x[col]
            x[row] = s / lower[row, row]
        return x

    def ekf_update(x, p, z, h_mat, r_mat):
        innovation = z - h_mat @ x
        s_mat = h_mat @ p @ h_mat.T + r_mat
        s_mat = 0.5 * (s_mat + s_mat.T)
        hp = h_mat @ p
        gain_t = np.zeros_like(hp)
        for col in range(hp.shape[1]):
            gain_t[:, col] = solve_spd(s_mat, hp[:, col])
        gain = gain_t.T
        x_new = x + gain @ innovation
        ikh = np.eye(len(x)) - gain @ h_mat
        p_new = ikh @ p @ ikh.T + gain @ r_mat @ gain.T
        p_new = 0.5 * (p_new + p_new.T)
        nis = float(innovation @ solve_spd(s_mat, innovation))
        return x_new, p_new, nis

    # ----------------------------------------------------------------------
    # 3. Independent rotation-direction audit (nontrivial attitude)
    # ----------------------------------------------------------------------
    print("\nRotation-direction audit (q: 90 deg about z, f = e_x):")
    half_sqrt = math.sin(math.pi / 4.0)
    q_audit = (half_sqrt, 0.0, 0.0, half_sqrt)
    rotated = quat_to_dcm(*q_audit) @ np.array([1.0, 0.0, 0.0])
    print(f"  C(q) @ [1,0,0] = [{rotated[0]:.6f}, {rotated[1]:.6f}, "
          f"{rotated[2]:.6f}] (expect ~[0, 1, 0])")
    rotation_ok = bool(np.allclose(rotated, [0.0, 1.0, 0.0], atol=1e-12))
    print("  PASS" if rotation_ok else "  FAIL")

    # ----------------------------------------------------------------------
    # 4. Gravity Jacobian finite-difference audit
    # ----------------------------------------------------------------------
    print("\nGravity Jacobian finite-difference audit:")
    positions = [
        np.array([7e6, 0.0, 0.0]),
        np.array([-3e6, 5e6, 8e6]),
        np.array([1.234e6, -6.789e6, 2.345e6]),
    ]
    steps = [1e6, 1e5, 1e4, 1e3, 1e2, 10.0, 1.0, 0.1]
    worst_best = 0.0
    for r in positions:
        g_analytic = gravity_jacobian(r)
        norm_analytic = np.linalg.norm(g_analytic)
        best_rel = np.inf
        for h in steps:
            g_num = np.zeros((3, 3))
            for j in range(3):
                rp = r.copy()
                rp[j] += h
                rm = r.copy()
                rm[j] -= h
                g_num[:, j] = (two_body_accel(rp) - two_body_accel(rm)) / (2 * h)
            rel = np.linalg.norm(g_analytic - g_num) / norm_analytic
            best_rel = min(best_rel, rel)
        worst_best = max(worst_best, best_rel)
    jacobian_ok = worst_best < 1e-10
    print(f"  Worst best relative error: {worst_best:.2e}")
    print("  PASS" if jacobian_ok else "  FAIL")

    # ----------------------------------------------------------------------
    # 5. Replay the segment through the independent filter
    # ----------------------------------------------------------------------
    truth_pos0 = None  # built from documented constants instead
    _ = truth_pos0
    x_state = np.zeros(6)
    # Initial TRUE state is the first row's implied orbit start; the oracle
    # reconstructs it from the documented circular-orbit constants rather
    # than reading truth columns (the segment carries no truth columns).
    r0 = 6378137.0 + 500.0e3
    speed = math.sqrt(MU / r0)
    true_pos0 = np.array([r0, 0.0, 0.0])
    true_vel0 = np.array([0.0, speed, 0.0])
    x_state[0:3] = true_pos0 + np.array(INIT_POS_ERR)
    x_state[3:6] = true_vel0 + np.array(INIT_VEL_ERR)
    p_cov = np.diag([P0_POS, P0_POS, P0_POS, P0_VEL, P0_VEL, P0_VEL])
    h_mat = np.eye(6)
    r_mat = np.diag([SIGMA_R**2] * 3 + [SIGMA_V**2] * 3)

    pos_diffs = []
    vel_diffs = []
    nis_values = []

    previous_time = None
    for k in range(n):
        row = data[k]
        t_now = float(row["time_s"])
        if previous_time is not None:
            dt = t_now - previous_time
            if dt > 0.0:
                f_body = np.array([
                    row["imu_specific_force_body_x_mps2"],
                    row["imu_specific_force_body_y_mps2"],
                    row["imu_specific_force_body_z_mps2"],
                ])
                q_attitude = (
                    row["attitude_w"], row["attitude_x"],
                    row["attitude_y"], row["attitude_z"])
                q_norm = math.sqrt(sum(component * component
                                       for component in q_attitude))
                q_attitude = tuple(c / q_norm for c in q_attitude)

                r_est = x_state[0:3].copy()
                v_est = x_state[3:6].copy()
                r_new, v_new = rk4_imu_step(r_est, v_est, dt, f_body, q_attitude)
                x_state[0:3] = r_new
                x_state[3:6] = v_new

                phi = discrete_transition(r_est, dt)
                p_cov = phi @ p_cov @ phi.T + process_noise(dt)
                p_cov = 0.5 * (p_cov + p_cov.T)
        previous_time = t_now

        if row["gnss_valid"] > 0.5:
            z_meas = np.array([
                row["gnss_pos_x_m"], row["gnss_pos_y_m"], row["gnss_pos_z_m"],
                row["gnss_vel_x_mps"], row["gnss_vel_y_mps"], row["gnss_vel_z_mps"],
            ])
            x_state, p_cov, nis = ekf_update(x_state, p_cov, z_meas, h_mat, r_mat)
            nis_values.append(nis)

        cpp_estimate = np.array([
            row["est_pos_x_m"], row["est_pos_y_m"], row["est_pos_z_m"],
            row["est_vel_x_mps"], row["est_vel_y_mps"], row["est_vel_z_mps"],
        ])
        diff = np.abs(x_state - cpp_estimate)
        pos_diffs.append(diff[0:3].max())
        vel_diffs.append(diff[3:6].max())

    # ----------------------------------------------------------------------
    # 6. Report
    # ----------------------------------------------------------------------
    pos_tolerance = 1e-6
    vel_tolerance = 1e-8
    max_pos_diff = max(pos_diffs[-100:])
    max_vel_diff = max(vel_diffs[-100:])
    mean_nis = (float(np.mean(nis_values))
                if nis_values else float("nan"))

    print(f"\n{'=' * 60}")
    print("Independent Python IMU-EKF Oracle Results")
    print(f"{'=' * 60}")
    print(f"Max position discrepancy (last 100 steps): {max_pos_diff:.3e} m")
    print(f"  Tolerance: {pos_tolerance:.0e} m")
    print(f"Max velocity discrepancy (last 100 steps): {max_vel_diff:.3e} m/s")
    print(f"  Tolerance: {vel_tolerance:.0e} m/s")
    print(f"Python mean NIS ({len(nis_values)} updates): {mean_nis:.4f}")

    overall_ok = rotation_ok and jacobian_ok
    if max_pos_diff > pos_tolerance:
        print(f"FAIL: position discrepancy {max_pos_diff:.3e} > {pos_tolerance}")
        overall_ok = False
    if max_vel_diff > vel_tolerance:
        print(f"FAIL: velocity discrepancy {max_vel_diff:.3e} > {vel_tolerance}")
        overall_ok = False

    print(f"\n{'=' * 60}")
    if overall_ok:
        print("OVERALL: PASS")
    else:
        print("OVERALL: FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
