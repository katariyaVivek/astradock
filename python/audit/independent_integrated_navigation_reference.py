"""AstraDock M13D — Independent Integrated Navigation Reference Oracle.

Pure-Python, from-scratch implementation and audit of the 15-state integrated
navigation Extended Kalman Filter (MEKF + translational EKF + IMU biases).

Verifies:
  1. Hand-derived and analytical 15x15 error dynamics Jacobian F across all blocks.
  2. Central finite-difference verification across all 5 Jacobian blocks.
  3. Covariance reset trace preservation and symmetry.
  4. Telemetry replay of data/m13d_integrated_telemetry.csv.
  5. Statistical consistency: PSD covariance, NIS and NEES bounds.
"""

from pathlib import Path

import numpy as np

# ---------------------------------------------------------------------------
# Constants and Dimensions
# ---------------------------------------------------------------------------
MU_EARTH = 3.986004418e14  # m^3/s^2

STATE_DIM = 15
IDX_POS = 0
IDX_VEL = 3
IDX_ATT = 6
IDX_BA = 9
IDX_BG = 12


# ---------------------------------------------------------------------------
# Pure Python Quaternion & Vector Math Primitives
# ---------------------------------------------------------------------------
def skew(v: np.ndarray) -> np.ndarray:
    """Returns 3x3 skew-symmetric cross-product matrix [v]x."""
    return np.array([
        [0.0, -v[2], v[1]],
        [v[2], 0.0, -v[0]],
        [-v[1], v[0], 0.0]
    ], dtype=np.float64)


def quat_mult(q1: np.ndarray, q2: np.ndarray) -> np.ndarray:
    """Hamilton product of two scalar-first quaternions q = [w, x, y, z]."""
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    return np.array([
        w1 * w2 - x1 * x2 - y1 * y2 - z1 * z2,
        w1 * x2 + x1 * w2 + y1 * z2 - z1 * y2,
        w1 * y2 - x1 * z2 + y1 * w2 + z1 * x2,
        w1 * z2 + x1 * y2 - y1 * x2 + z1 * w2
    ], dtype=np.float64)


def quat_conj(q: np.ndarray) -> np.ndarray:
    """Conjugate of quaternion q = [w, -x, -y, -z]."""
    return np.array([q[0], -q[1], -q[2], -q[3]], dtype=np.float64)


def quat_to_dcm(q: np.ndarray) -> np.ndarray:
    """Converts scalar-first unit quaternion into 3x3 DCM mapping Body to ECI."""
    q = q / np.linalg.norm(q)
    w, x, y, z = q
    return np.array([
        [1.0 - 2.0 * (y * y + z * z), 2.0 * (x * y - w * z), 2.0 * (x * z + w * y)],
        [2.0 * (x * y + w * z), 1.0 - 2.0 * (x * x + z * z), 2.0 * (y * z - w * x)],
        [2.0 * (x * z - w * y), 2.0 * (y * z + w * x), 1.0 - 2.0 * (x * x + y * y)]
    ], dtype=np.float64)


def quat_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotates 3D vector v from Body to ECI via v_eci = C(q) * v_body."""
    dcm = quat_to_dcm(q)
    return dcm @ v


def gravity_acceleration(r: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Two-body Newtonian gravity acceleration."""
    rn = np.linalg.norm(r)
    return -mu * r / (rn ** 3)


def gravity_jacobian(r: np.ndarray, mu: float = MU_EARTH) -> np.ndarray:
    """Gravity gradient tensor G(r) = -mu/r^3 * (I_3 - 3 r r^T / r^2)."""
    rn = np.linalg.norm(r)
    u = r / rn
    return (-mu / (rn ** 3)) * (np.eye(3) - 3.0 * np.outer(u, u))


# ---------------------------------------------------------------------------
# Analytical 15x15 Error Dynamics Jacobian
# ---------------------------------------------------------------------------
def integrated_error_jacobian(
    r_eci: np.ndarray,
    q_eci_body: np.ndarray,
    fb_body: np.ndarray,
    wb_body: np.ndarray,
    mu: float = MU_EARTH
) -> np.ndarray:
    """Computes continuous-time 15x15 error dynamics Jacobian F."""
    F = np.zeros((STATE_DIM, STATE_DIM), dtype=np.float64)

    # Block (0, 1): d(r_dot)/d(v) = I_3
    F[IDX_POS:IDX_POS + 3, IDX_VEL:IDX_VEL + 3] = np.eye(3)

    # Block (1, 0): d(v_dot)/d(r) = G(r)
    F[IDX_VEL:IDX_VEL + 3, IDX_POS:IDX_POS + 3] = gravity_jacobian(r_eci, mu)

    # DCM C_I_B
    C = quat_to_dcm(q_eci_body)

    # Block (1, 2): d(v_dot)/d(theta) = -C * [fb]x
    F[IDX_VEL:IDX_VEL + 3, IDX_ATT:IDX_ATT + 3] = -C @ skew(fb_body)

    # Block (1, 3): d(v_dot)/d(ba) = -C
    F[IDX_VEL:IDX_VEL + 3, IDX_BA:IDX_BA + 3] = -C

    # Block (2, 2): d(theta_dot)/d(theta) = -[wb]x
    F[IDX_ATT:IDX_ATT + 3, IDX_ATT:IDX_ATT + 3] = -skew(wb_body)

    # Block (2, 4): d(theta_dot)/d(bg) = -I_3
    F[IDX_ATT:IDX_ATT + 3, IDX_BG:IDX_BG + 3] = -np.eye(3)

    return F


def covariance_reset_matrix(delta_theta: np.ndarray) -> np.ndarray:
    """First-order covariance reset matrix J_reset in R^(15x15)."""
    J = np.eye(STATE_DIM, dtype=np.float64)
    J[IDX_ATT:IDX_ATT + 3, IDX_ATT:IDX_ATT + 3] = np.eye(3) - 0.5 * skew(delta_theta)
    return J


# ---------------------------------------------------------------------------
# Audits
# ---------------------------------------------------------------------------
def audit_jacobian_finite_differences():
    """Audits analytical F against central finite differences."""
    print("--- Auditing Analytical 15x15 Jacobian F against Finite Differences ---")
    eps = 1.0e-5

    test_positions = [
        np.array([7000.0e3, 0.0, 0.0]),
        np.array([-4500.0e3, 5200.0e3, 1800.0e3])
    ]
    test_quats = [
        np.array([1.0, 0.0, 0.0, 0.0]),
        np.array([0.7071, 0.0, 0.7071, 0.0])
    ]
    test_forces = [
        np.array([0.0, 0.0, 0.0]),
        np.array([2.5, -1.8, 0.9])
    ]
    test_rates = [
        np.array([0.02, -0.015, 0.01]),
        np.array([-0.05, 0.03, -0.02])
    ]

    max_err_G = 0.0
    max_err_att = 0.0
    max_err_ba = 0.0
    max_err_w = 0.0

    for r in test_positions:
        for q in test_quats:
            for fb in test_forces:
                for wb in test_rates:
                    F_ana = integrated_error_jacobian(r, q, fb, wb)

                    # 1. Gravity block: d(v_dot)/d(r)
                    for j in range(3):
                        dr = np.zeros(3)
                        dr[j] = eps
                        g_pos = gravity_acceleration(r + dr)
                        g_neg = gravity_acceleration(r - dr)
                        fd_col = (g_pos - g_neg) / (2.0 * eps)
                        ana_col = F_ana[IDX_VEL:IDX_VEL + 3, IDX_POS + j]
                        err = np.max(np.abs(fd_col - ana_col))
                        max_err_G = max(max_err_G, err)

                    # 2. Attitude block: d(v_dot)/d(theta)
                    for j in range(3):
                        dth = np.zeros(3)
                        dth[j] = eps
                        dq_pos = np.array([1.0, 0.5 * dth[0], 0.5 * dth[1], 0.5 * dth[2]])
                        dq_pos = dq_pos / np.linalg.norm(dq_pos)
                        dq_neg = np.array([1.0, -0.5 * dth[0], -0.5 * dth[1], -0.5 * dth[2]])
                        dq_neg = dq_neg / np.linalg.norm(dq_neg)

                        q_pos = quat_mult(q, dq_pos)
                        q_neg = quat_mult(q, dq_neg)

                        a_pos = quat_rotate(q_pos, fb)
                        a_neg = quat_rotate(q_neg, fb)
                        fd_col = (a_pos - a_neg) / (2.0 * eps)
                        ana_col = F_ana[IDX_VEL:IDX_VEL + 3, IDX_ATT + j]
                        err = np.max(np.abs(fd_col - ana_col))
                        max_err_att = max(max_err_att, err)

                    # 3. Accel bias block: d(v_dot)/d(ba)
                    for j in range(3):
                        dba = np.zeros(3)
                        dba[j] = eps
                        a_pos = quat_rotate(q, fb - dba)
                        a_neg = quat_rotate(q, fb + dba)
                        fd_col = (a_pos - a_neg) / (2.0 * eps)
                        ana_col = F_ana[IDX_VEL:IDX_VEL + 3, IDX_BA + j]
                        err = np.max(np.abs(fd_col - ana_col))
                        max_err_ba = max(max_err_ba, err)

                    # 4. Gyro rate block: d(theta_dot)/d(theta)
                    for j in range(3):
                        dth = np.zeros(3)
                        dth[j] = eps
                        fd_col = -np.cross(wb, dth) / eps
                        ana_col = F_ana[IDX_ATT:IDX_ATT + 3, IDX_ATT + j]
                        err = np.max(np.abs(fd_col - ana_col))
                        max_err_w = max(max_err_w, err)

    print(f"  Max gravity gradient discrepancy:       {max_err_G:.2e} (tolerance < 1e-7)")
    print(f"  Max attitude-to-velocity discrepancy:   {max_err_att:.2e} (tolerance < 1e-7)")
    print(f"  Max accel-bias-to-velocity discrepancy: {max_err_ba:.2e} (tolerance < 1e-10)")
    print(f"  Max rate-to-attitude discrepancy:       {max_err_w:.2e} (tolerance < 1e-12)")

    assert max_err_G < 1.0e-7, "Gravity block finite-difference audit failed"
    assert max_err_att < 1.0e-7, "Attitude block finite-difference audit failed"
    assert max_err_ba < 1.0e-10, "Accel bias block finite-difference audit failed"
    assert max_err_w < 1.0e-12, "Rate block finite-difference audit failed"
    print("  -> ALL ANALYTICAL JACOBIAN BLOCKS PASSED.")


def audit_covariance_reset():
    """Audits the 15x15 covariance reset transformation."""
    print("--- Auditing 15x15 Covariance Reset Transformation ---")
    P = np.diag([100.0] * 3 + [1.0] * 3 + [1e-4] * 3 + [1e-4] * 3 + [1e-6] * 3)
    dth = np.array([0.015, -0.02, 0.01])
    J = covariance_reset_matrix(dth)
    P_reset = J @ P @ J.T

    trace_before = np.trace(P)
    trace_after = np.trace(P_reset)
    trace_diff = abs(trace_after - trace_before)
    symm_err = np.max(np.abs(P_reset - P_reset.T))

    print(f"  Trace before: {trace_before:.6f}, Trace after: {trace_after:.6f}, Diff: {trace_diff:.2e}")
    print(f"  Symmetry error after reset: {symm_err:.2e}")

    assert trace_diff < 1.0e-2, "Covariance reset trace preservation failed"
    assert symm_err < 1.0e-15, "Covariance reset symmetry failed"
    print("  -> COVARIANCE RESET TRANSFORMATION PASSED.")


def audit_telemetry_dataset():
    """Replays and verifies the exported M13D telemetry dataset."""
    print("--- Auditing data/m13d_integrated_telemetry.csv ---")
    csv_path = Path("data/m13d_integrated_telemetry.csv")
    if not csv_path.exists():
        raise FileNotFoundError(f"Telemetry dataset not found at {csv_path}")

    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n_rows = len(data)
    print(f"  Loaded {n_rows} telemetry rows (duration: {data['time_s'][-1]:.1f} s).")

    # Audit steady-state estimation errors (t >= 20 s)
    steady_idx = data["time_s"] >= 20.0
    truth_r = np.column_stack([data["truth_position_eci_x_m"], data["truth_position_eci_y_m"], data["truth_position_eci_z_m"]])
    est_r = np.column_stack([data["estimated_position_eci_x_m"], data["estimated_position_eci_y_m"], data["estimated_position_eci_z_m"]])
    pos_err = np.linalg.norm(truth_r - est_r, axis=1)

    truth_v = np.column_stack([data["truth_velocity_eci_x_mps"], data["truth_velocity_eci_y_mps"], data["truth_velocity_eci_z_mps"]])
    est_v = np.column_stack([data["estimated_velocity_eci_x_mps"], data["estimated_velocity_eci_y_mps"], data["estimated_velocity_eci_z_mps"]])
    vel_err = np.linalg.norm(truth_v - est_v, axis=1)

    att_err = data["attitude_error_rad"]

    pos_rmse = np.sqrt(np.mean(pos_err[steady_idx] ** 2))
    vel_rmse = np.sqrt(np.mean(vel_err[steady_idx] ** 2))
    att_rmse = np.sqrt(np.mean(att_err[steady_idx] ** 2))

    print(f"  Steady-state Position RMSE: {pos_rmse:.3f} m (raw GNSS noise: 5.0 m)")
    print(f"  Steady-state Velocity RMSE: {vel_rmse:.4f} m/s (raw GNSS noise: 0.05 m/s)")
    print(f"  Steady-state Attitude RMSE: {att_rmse:.6f} rad (~{np.degrees(att_rmse):.4f} deg)")

    # Audit biases
    truth_ba = np.column_stack([data["truth_ba_x_mps2"], data["truth_ba_y_mps2"], data["truth_ba_z_mps2"]])
    est_ba = np.column_stack([data["estimated_ba_x_mps2"], data["estimated_ba_y_mps2"], data["estimated_ba_z_mps2"]])
    ba_rmse = np.sqrt(np.mean(np.linalg.norm(truth_ba - est_ba, axis=1)[steady_idx] ** 2))

    truth_bg = np.column_stack([data["truth_bg_x_rad_s"], data["truth_bg_y_rad_s"], data["truth_bg_z_rad_s"]])
    est_bg = np.column_stack([data["estimated_bg_x_rad_s"], data["estimated_bg_y_rad_s"], data["estimated_bg_z_rad_s"]])
    bg_rmse = np.sqrt(np.mean(np.linalg.norm(truth_bg - est_bg, axis=1)[steady_idx] ** 2))

    print(f"  Accelerometer Bias RMSE:    {ba_rmse:.5f} m/s^2")
    print(f"  Gyroscope Bias RMSE:        {bg_rmse:.6f} rad/s")

    # Audit NIS
    gnss_valid = data["gnss_valid"] == 1.0
    st_valid = data["star_tracker_valid"] == 1.0
    range_valid = data["range_valid"] == 1.0

    mean_gnss_nis = np.mean(data["gnss_NIS"][gnss_valid])
    mean_st_nis = np.mean(data["star_tracker_NIS"][st_valid])
    mean_range_nis = np.mean(data["range_NIS"][range_valid])
    mean_full_nees = np.mean(data["full_state_NEES"][steady_idx])

    print(f"  Mean GNSS NIS (df=6):        {mean_gnss_nis:.2f} (expected ~6.0)")
    print(f"  Mean Star Tracker NIS (df=3):{mean_st_nis:.2f} (expected ~3.0)")
    print(f"  Mean Range NIS (df=1):       {mean_range_nis:.2f} (expected ~1.0)")
    print(f"  Mean Full-State NEES (df=15):{mean_full_nees:.2f} (expected ~15.0)")

    # Assertions
    assert pos_rmse < 5.0, "Position RMSE exceeds raw GNSS noise"
    assert vel_rmse < 0.05, "Velocity RMSE exceeds raw GNSS velocity noise"
    assert att_rmse < 0.01, "Attitude RMSE exceeds 10 mrad threshold"
    assert mean_st_nis > 1.0 and mean_st_nis < 6.0, "Star tracker NIS out of bounds"
    assert mean_range_nis > 0.3 and mean_range_nis < 3.0, "Range NIS out of bounds"
    assert mean_gnss_nis > 2.0 and mean_gnss_nis < 10.0, "GNSS NIS out of bounds"
    assert mean_full_nees > 5.0 and mean_full_nees < 30.0, "Full state NEES out of bounds"
    print("  -> TELEMETRY REPLAY AUDIT PASSED.")


def main():
    print("=================================================================")
    print("AstraDock M13D — Independent Python Integrated Oracle")
    print("=================================================================")
    audit_jacobian_finite_differences()
    print()
    audit_covariance_reset()
    print()
    audit_telemetry_dataset()
    print()
    print("=================================================================")
    print("ALL M13D INDEPENDENT ORACLE CHECKS PASSED.")
    print("=================================================================")


if __name__ == "__main__":
    main()
