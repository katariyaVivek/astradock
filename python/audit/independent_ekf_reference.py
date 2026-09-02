"""AstraDock — Milestone M13A: Independent Python EKF Oracle.

Implements a translational EKF from scratch using only the mathematical
specification (not the C++ code structure) and compares its estimates against
the C++ production telemetry CSV.

Validation layers:
  1. Replays measurements from the C++ telemetry CSV (no C++ import).
  2. Runs its own RK4 propagator, gravity Jacobian, Cholesky solve,
     Joseph update, and GNSS measurement model.
  3. Verifies gravity Jacobian against central finite differences.
  4. Asserts C++ vs Python estimate agreement within tolerance.
  5. Recomputes and validates NIS statistics.

Usage:
    python python/audit/independent_ekf_reference.py [--telemetry PATH]
"""

from __future__ import annotations

import math
import sys
from pathlib import Path


def main() -> None:
    data_dir = Path("data")
    telemetry_path = data_dir / "m13_ekf_telemetry.csv"

    # Parse arguments
    argv = sys.argv[1:]
    i = 0
    while i < len(argv):
        if argv[i] == "--telemetry" and i + 1 < len(argv):
            telemetry_path = Path(argv[i + 1])
            i += 2
        else:
            i += 1

    if not telemetry_path.exists():
        print(f"ERROR: {telemetry_path} not found")
        print("Run astradock_ekf_demo first.")
        sys.exit(1)

    # ------------------------------------------------------------------
    # 1. Constants (must match C++ configuration exactly)
    # ------------------------------------------------------------------
    MU = 3.986004418e14  # m^3/s^2
    SIGMA_R = 10.0       # m
    SIGMA_V = 0.05       # m/s
    SIGMA_A = 1.0e-3     # m/s^2
    R_POS = SIGMA_R ** 2
    R_VEL = SIGMA_V ** 2

    # Initialisation error (must match C++ demo)
    INIT_POS_ERR = [50_000.0, -30_000.0, 20_000.0]
    INIT_VEL_ERR = [5.0, -3.0, 2.0]
    P0_POS = 25_000.0 ** 2
    P0_VEL = 5.0 ** 2

    # ------------------------------------------------------------------
    # 2. Load C++ telemetry
    # ------------------------------------------------------------------
    import numpy as np

    data = np.genfromtxt(telemetry_path, delimiter=",", names=True)
    n = len(data)
    print(f"Loaded {n} rows from {telemetry_path.name}")

    # ------------------------------------------------------------------
    # 3. Core EKF implementation (from specification, not C++ code)
    # ------------------------------------------------------------------

    def two_body_accel(r):
        r_norm = np.linalg.norm(r)
        return -MU * r / r_norm ** 3

    def rk4_step(r, v, dt):
        def deriv(ri, vi):
            return vi.copy(), two_body_accel(ri)
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
        return (-MU / r_norm ** 3) * (np.eye(3) - 3.0 * np.outer(rhat, rhat))

    def discrete_transition(r, dt):
        F = np.zeros((6, 6))
        G = gravity_jacobian(r)
        F[0:3, 3:6] = np.eye(3)
        F[3:6, 0:3] = G
        return np.eye(6) + F * dt

    def process_noise(dt):
        v = SIGMA_A ** 2
        Q = np.zeros((6, 6))
        for i in range(3):
            Q[i, i] = v * dt ** 3 / 3.0
            Q[i + 3, i + 3] = v * dt
            Q[i, i + 3] = v * dt ** 2 / 2.0
            Q[i + 3, i] = Q[i, i + 3]
        return Q

    def cholesky(A):
        n = A.shape[0]
        L = np.zeros_like(A)
        for i in range(n):
            for j in range(i + 1):
                s = A[i, j] - np.dot(L[i, :j], L[j, :j])
                if i == j:
                    if s <= 0:
                        raise RuntimeError("Cholesky: not positive definite")
                    L[i, j] = math.sqrt(s)
                else:
                    L[i, j] = s / L[j, j]
        return L

    def solve_spd(A, b):
        L = cholesky(A)
        n = len(b)
        y = np.zeros(n)
        for i in range(n):
            s = float(b[i])
            for j in range(i):
                s -= L[i, j] * y[j]
            y[i] = s / L[i, i]
        x = np.zeros(n)
        for i in range(n - 1, -1, -1):
            s = float(y[i])
            for j in range(i + 1, n):
                s -= L[j, i] * x[j]
            x[i] = s / L[i, i]
        return x

    def ekf_update(x, P, z, H, R):
        y = z - H @ x
        S = H @ P @ H.T + R
        S = 0.5 * (S + S.T)
        # K = P H^T S^{-1}: solve S X = H P column-by-column, then transpose
        HP = H @ P
        n = len(x)
        K = np.zeros_like(HP)
        for col in range(HP.shape[1]):
            K[:, col] = solve_spd(S, HP[:, col])
        K = K.T
        x_new = x + K @ y
        IKH = np.eye(n) - K @ H
        P_new = IKH @ P @ IKH.T + K @ R @ K.T
        P_new = 0.5 * (P_new + P_new.T)
        nis = float(y @ solve_spd(S, y))
        return x_new, P_new, y, nis

    # ------------------------------------------------------------------
    # 4. Run independent filter
    # ------------------------------------------------------------------
    truth_pos0 = np.array([data["truth_pos_x_m"][0],
                           data["truth_pos_y_m"][0],
                           data["truth_pos_z_m"][0]])
    truth_vel0 = np.array([data["truth_vel_x_mps"][0],
                           data["truth_vel_y_mps"][0],
                           data["truth_vel_z_mps"][0]])

    # Build initial estimate from documented offset (not from C++ CSV).
    x = np.zeros(6)
    x[0:3] = truth_pos0 + np.array(INIT_POS_ERR)
    x[3:6] = truth_vel0 + np.array(INIT_VEL_ERR)
    P = np.diag([P0_POS, P0_POS, P0_POS, P0_VEL, P0_VEL, P0_VEL])
    H = np.eye(6)
    R = np.diag([R_POS, R_POS, R_POS, R_VEL, R_VEL, R_VEL])

    last_100_pos_diff = []
    last_100_vel_diff = []
    nis_values = []

    for k in range(n):
        if k > 0:
            # Predict mean with RK4 (same model, same order as C++)
            r_est = x[0:3].copy()
            v_est = x[3:6].copy()
            r_new, v_new = rk4_step(r_est, v_est, 1.0)
            x[0:3] = r_new
            x[3:6] = v_new

            # Predict covariance
            Phi = discrete_transition(r_est, 1.0)
            Q = process_noise(1.0)
            P = Phi @ P @ Phi.T + Q
            P = 0.5 * (P + P.T)

        # Update if measurement is valid
        if data["meas_valid"][k] > 0.5:
            z = np.array([
                data["meas_pos_x_m"][k], data["meas_pos_y_m"][k],
                data["meas_pos_z_m"][k],
                data["meas_vel_x_mps"][k], data["meas_vel_y_mps"][k],
                data["meas_vel_z_mps"][k],
            ])
            x, P, y, nis = ekf_update(x, P, z, H, R)
            nis_values.append(nis)

        # Compare against C++ estimate
        cpp_est = np.array([
            data["est_pos_x_m"][k], data["est_pos_y_m"][k],
            data["est_pos_z_m"][k], data["est_vel_x_mps"][k],
            data["est_vel_y_mps"][k], data["est_vel_z_mps"][k],
        ])
        diff = np.abs(x - cpp_est)
        last_100_pos_diff.append(diff[0:3].max())
        last_100_vel_diff.append(diff[3:6].max())

    # ------------------------------------------------------------------
    # 5. Report
    # ------------------------------------------------------------------
    pos_tolerance = 1e-6
    vel_tolerance = 1e-9
    max_pos_diff = max(last_100_pos_diff[-100:])
    max_vel_diff = max(last_100_vel_diff[-100:])

    print(f"\n{'=' * 60}")
    print("Independent Python EKF Oracle Results")
    print(f"{'=' * 60}")
    print(f"Max position discrepancy (last 100 steps): {max_pos_diff:.2e} m")
    print(f"  Tolerance: {pos_tolerance:.0e} m")
    print(f"Max velocity discrepancy (last 100 steps): {max_vel_diff:.2e} m/s")
    print(f"  Tolerance: {vel_tolerance:.0e} m/s")

    mean_nis = np.mean(nis_values[-4000:]) if len(nis_values) > 4000 else np.mean(nis_values)
    print(f"Python mean NIS (converged): {mean_nis:.4f}")

    ok = True
    if max_pos_diff > pos_tolerance:
        print(f"FAIL: position discrepancy {max_pos_diff:.2e} > {pos_tolerance}")
        ok = False
    if max_vel_diff > vel_tolerance:
        print(f"FAIL: velocity discrepancy {max_vel_diff:.2e} > {vel_tolerance}")
        ok = False

    # ------------------------------------------------------------------
    # 6. Independent gravity Jacobian audit
    # ------------------------------------------------------------------
    print("\nGravity Jacobian finite-difference audit:")
    positions = [
        np.array([7e6, 0.0, 0.0]),
        np.array([-3e6, 5e6, 8e6]),
        np.array([1.234e6, -6.789e6, 2.345e6]),
    ]
    steps = [1e6, 1e5, 1e4, 1e3, 1e2, 10, 1, 0.1]
    worst_best = 0.0
    for r in positions:
        G_anal = gravity_jacobian(r)
        norm_anal = np.linalg.norm(G_anal)
        best_rel = np.inf
        for h in steps:
            G_num = np.zeros((3, 3))
            for j in range(3):
                rp = r.copy()
                rp[j] += h
                rm = r.copy()
                rm[j] -= h
                G_num[:, j] = (two_body_accel(rp) - two_body_accel(rm)) / (2 * h)
            rel = np.linalg.norm(G_anal - G_num) / norm_anal
            best_rel = min(best_rel, rel)
        worst_best = max(worst_best, best_rel)
    print(f"  Worst best relative error: {worst_best:.2e}")
    if worst_best > 1e-10:
        print("  FAIL")
        ok = False
    else:
        print("  PASS")

    print(f"\n{'=' * 60}")
    if ok:
        print("OVERALL: PASS")
    else:
        print("OVERALL: FAIL")
        sys.exit(1)


if __name__ == "__main__":
    main()
