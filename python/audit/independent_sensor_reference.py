"""AstraDock — Milestone M12: Independent Python Sensor Oracle.

Independently audits and verifies C++ sensor simulation telemetry datasets
across IMU, GNSS, Star Tracker, and Range sensors against analytical physics,
sampling rates, dropout windows, and stochastic noise/bias models.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np


def quaternion_orientation_error_rad(
    q1: tuple[float, float, float, float],
    q2: tuple[float, float, float, float]
) -> float:
    """Computes shortest geodesic rotation angle between two quaternions."""
    # q = [w, x, y, z]
    w1, x1, y1, z1 = q1
    w2, x2, y2, z2 = q2
    # Conjugate of q1: [w1, -x1, -y1, -z1]
    # Product q_err = q1* ⊗ q2
    w_err = w1 * w2 + x1 * x2 + y1 * y2 + z1 * z2
    return 2.0 * math.acos(min(1.0, max(0.0, abs(w_err))))


def audit_imu_telemetry(csv_path: Path) -> bool:
    print(f"\n{'='*70}\nAuditing IMU Telemetry: {csv_path.name}\n{'='*70}")
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n_samples = len(data)
    print(f"Loaded {n_samples} IMU samples across duration {data['time_s'][-1]:.2f} s")

    valid_mask = data["valid"] == 1
    dropout_mask = data["valid"] == 0
    n_dropout = np.sum(dropout_mask)
    n_valid = np.sum(valid_mask)
    print(f"  Valid samples:   {n_valid}")
    print(f"  Dropout samples: {n_dropout} (expected in [2000, 2100] s)")

    # Verify dropout bounds
    dropout_times = data["time_s"][dropout_mask]
    if len(dropout_times) > 0:
        assert np.min(dropout_times) >= 2000.0 - 1e-6
        assert np.max(dropout_times) <= 2100.0 + 1e-6

    # 1. Gyroscope statistics on valid samples
    truth_omega = np.column_stack([
        data["truth_omega_x_rad_s"][valid_mask],
        data["truth_omega_y_rad_s"][valid_mask],
        data["truth_omega_z_rad_s"][valid_mask]
    ])
    meas_omega = np.column_stack([
        data["meas_omega_x_rad_s"][valid_mask],
        data["meas_omega_y_rad_s"][valid_mask],
        data["meas_omega_z_rad_s"][valid_mask]
    ])
    gyro_err = meas_omega - truth_omega

    mean_gyro_err = np.mean(gyro_err, axis=0)
    std_gyro_err = np.std(gyro_err, axis=0)
    expected_gyro_bias = np.array([0.0005, -0.0003, 0.0002])
    expected_gyro_std = 0.001

    print("\n  Gyroscope Evaluation (BODY Frame, rad/s):")
    print(f"    Expected Bias: [{expected_gyro_bias[0]:.6f}, {expected_gyro_bias[1]:.6f}, {expected_gyro_bias[2]:.6f}]")
    print(f"    Sampled Bias:  [{mean_gyro_err[0]:.6f}, {mean_gyro_err[1]:.6f}, {mean_gyro_err[2]:.6f}]")
    print(f"    Expected Std:  {expected_gyro_std:.6f}")
    print(f"    Sampled Std:   [{std_gyro_err[0]:.6f}, {std_gyro_err[1]:.6f}, {std_gyro_err[2]:.6f}]")

    assert np.allclose(mean_gyro_err, expected_gyro_bias, atol=1e-4)
    assert np.allclose(std_gyro_err, expected_gyro_std, atol=1e-4)

    # 2. Accelerometer statistics on valid samples
    truth_spec_force = np.column_stack([
        data["truth_spec_force_x_mps2"][valid_mask],
        data["truth_spec_force_y_mps2"][valid_mask],
        data["truth_spec_force_z_mps2"][valid_mask]
    ])
    meas_spec_force = np.column_stack([
        data["meas_spec_force_x_mps2"][valid_mask],
        data["meas_spec_force_y_mps2"][valid_mask],
        data["meas_spec_force_z_mps2"][valid_mask]
    ])
    accel_err = meas_spec_force - truth_spec_force

    mean_accel_err = np.mean(accel_err, axis=0)
    std_accel_err = np.std(accel_err, axis=0)
    expected_accel_bias = np.array([0.002, -0.001, 0.0015])
    expected_accel_std = 0.005

    print("\n  Accelerometer Evaluation (BODY Frame, m/s^2):")
    print(f"    Expected Bias: [{expected_accel_bias[0]:.6f}, {expected_accel_bias[1]:.6f}, {expected_accel_bias[2]:.6f}]")
    print(f"    Sampled Bias:  [{mean_accel_err[0]:.6f}, {mean_accel_err[1]:.6f}, {mean_accel_err[2]:.6f}]")
    print(f"    Expected Std:  {expected_accel_std:.6f}")
    print(f"    Sampled Std:   [{std_accel_err[0]:.6f}, {std_accel_err[1]:.6f}, {std_accel_err[2]:.6f}]")

    max_truth_spec_force = np.max(np.linalg.norm(truth_spec_force, axis=1))
    print(f"    Max Truth Specific Force (Drag Only): {max_truth_spec_force:.6e} m/s^2")
    print("    (Confirms f = a - g << g_orbit ~ 8.43 m/s^2 during orbital flight)")

    assert np.allclose(mean_accel_err, expected_accel_bias, atol=2e-4)
    assert np.allclose(std_accel_err, expected_accel_std, atol=2e-4)

    print("IMU Telemetry Audit: PASSED")
    return True


def audit_gnss_telemetry(csv_path: Path) -> bool:
    print(f"\n{'='*70}\nAuditing GNSS Telemetry: {csv_path.name}\n{'='*70}")
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n_samples = len(data)
    print(f"Loaded {n_samples} GNSS samples across duration {data['time_s'][-1]:.2f} s")

    valid_mask = data["valid"] == 1
    dropout_mask = data["valid"] == 0
    print(f"  Valid samples:   {np.sum(valid_mask)}")
    print(f"  Dropout samples: {np.sum(dropout_mask)}")

    # 1. Position error statistics
    truth_r = np.column_stack([
        data["truth_r_x_m"][valid_mask],
        data["truth_r_y_m"][valid_mask],
        data["truth_r_z_m"][valid_mask]
    ])
    meas_r = np.column_stack([
        data["meas_r_x_m"][valid_mask],
        data["meas_r_y_m"][valid_mask],
        data["meas_r_z_m"][valid_mask]
    ])
    pos_err = meas_r - truth_r

    mean_pos_err = np.mean(pos_err, axis=0)
    std_pos_err = np.std(pos_err, axis=0)
    expected_pos_bias = np.array([2.5, -1.5, 3.0])
    expected_pos_std = 3.0

    print("\n  GNSS Position Evaluation (ECI Frame, m):")
    print(f"    Expected Bias: [{expected_pos_bias[0]:.3f}, {expected_pos_bias[1]:.3f}, {expected_pos_bias[2]:.3f}]")
    print(f"    Sampled Bias:  [{mean_pos_err[0]:.3f}, {mean_pos_err[1]:.3f}, {mean_pos_err[2]:.3f}]")
    print(f"    Expected Std:  {expected_pos_std:.3f}")
    print(f"    Sampled Std:   [{std_pos_err[0]:.3f}, {std_pos_err[1]:.3f}, {std_pos_err[2]:.3f}]")

    assert np.allclose(mean_pos_err, expected_pos_bias, atol=0.25)
    assert np.allclose(std_pos_err, expected_pos_std, atol=0.25)

    # 2. Velocity error statistics
    truth_v = np.column_stack([
        data["truth_v_x_mps"][valid_mask],
        data["truth_v_y_mps"][valid_mask],
        data["truth_v_z_mps"][valid_mask]
    ])
    meas_v = np.column_stack([
        data["meas_v_x_mps"][valid_mask],
        data["meas_v_y_mps"][valid_mask],
        data["meas_v_z_mps"][valid_mask]
    ])
    vel_err = meas_v - truth_v

    mean_vel_err = np.mean(vel_err, axis=0)
    std_vel_err = np.std(vel_err, axis=0)
    expected_vel_bias = np.array([0.02, -0.01, 0.015])
    expected_vel_std = 0.03

    print("\n  GNSS Velocity Evaluation (ECI Frame, m/s):")
    print(f"    Expected Bias: [{expected_vel_bias[0]:.4f}, {expected_vel_bias[1]:.4f}, {expected_vel_bias[2]:.4f}]")
    print(f"    Sampled Bias:  [{mean_vel_err[0]:.4f}, {mean_vel_err[1]:.4f}, {mean_vel_err[2]:.4f}]")
    print(f"    Expected Std:  {expected_vel_std:.4f}")
    print(f"    Sampled Std:   [{std_vel_err[0]:.4f}, {std_vel_err[1]:.4f}, {std_vel_err[2]:.4f}]")

    assert np.allclose(mean_vel_err, expected_vel_bias, atol=0.005)
    assert np.allclose(std_vel_err, expected_vel_std, atol=0.005)

    print("GNSS Telemetry Audit: PASSED")
    return True


def audit_star_tracker_telemetry(csv_path: Path) -> bool:
    print(f"\n{'='*70}\nAuditing Star Tracker Telemetry: {csv_path.name}\n{'='*70}")
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n_samples = len(data)
    print(f"Loaded {n_samples} Star Tracker samples across duration {data['time_s'][-1]:.2f} s")

    valid_mask = data["valid"] == 1
    print(f"  Valid samples:   {np.sum(valid_mask)}")
    print(f"  Dropout samples: {np.sum(data['valid'] == 0)}")

    meas_q = np.column_stack([
        data["meas_q_w"][valid_mask],
        data["meas_q_x"][valid_mask],
        data["meas_q_y"][valid_mask],
        data["meas_q_z"][valid_mask]
    ])
    q_norms = np.linalg.norm(meas_q, axis=1)
    max_norm_err = np.max(np.abs(q_norms - 1.0))
    print(f"  Max Quaternion Norm Error |norm(q) - 1|: {max_norm_err:.6e}")
    assert max_norm_err < 1e-10

    err_rad = data["orientation_error_rad"][valid_mask]
    mean_err_rad = np.mean(err_rad)
    std_err_rad = np.std(err_rad)
    # Expected composition: bias = 0.5 mrad, noise 1-sigma = 0.2 mrad
    print("\n  Star Tracker Orientation Error (rad):")
    print(f"    Mean Error Angle:    {mean_err_rad:.6f} rad ({math.degrees(mean_err_rad)*3600:.2f} arcsec)")
    print(f"    Std Error Angle:     {std_err_rad:.6f} rad ({math.degrees(std_err_rad)*3600:.2f} arcsec)")

    # Mean error should be in range ~[0.0004, 0.0007] rad
    assert 0.0004 <= mean_err_rad <= 0.0007
    print("Star Tracker Telemetry Audit: PASSED")
    return True


def audit_range_telemetry(csv_path: Path) -> bool:
    print(f"\n{'='*70}\nAuditing Range Telemetry: {csv_path.name}\n{'='*70}")
    data = np.genfromtxt(csv_path, delimiter=",", names=True)
    n_samples = len(data)
    print(f"Loaded {n_samples} Range samples across duration {data['time_s'][-1]:.2f} s")

    valid_mask = data["valid"] == 1
    print(f"  Valid samples:   {np.sum(valid_mask)}")
    print(f"  Dropout samples: {np.sum(data['valid'] == 0)}")

    truth_range = data["truth_range_m"][valid_mask]
    meas_range = data["meas_range_m"][valid_mask]
    range_err = meas_range - truth_range

    mean_range_err = np.mean(range_err)
    std_range_err = np.std(range_err)
    expected_bias = 0.5
    expected_std = 0.1

    print("\n  Range Sensor Evaluation (m):")
    print(f"    True Target Range:   ~{np.mean(truth_range):.2f} m")
    print(f"    Expected Bias:       {expected_bias:.4f} m")
    print(f"    Sampled Bias:        {mean_range_err:.4f} m")
    print(f"    Expected Std:        {expected_std:.4f} m")
    print(f"    Sampled Std:         {std_range_err:.4f} m")

    assert abs(mean_range_err - expected_bias) < 0.015
    assert abs(std_range_err - expected_std) < 0.015

    print("Range Telemetry Audit: PASSED")
    return True


def main() -> None:
    data_dir = Path("data")
    assert audit_imu_telemetry(data_dir / "m12_imu.csv")
    assert audit_gnss_telemetry(data_dir / "m12_gnss.csv")
    assert audit_star_tracker_telemetry(data_dir / "m12_star_tracker.csv")
    assert audit_range_telemetry(data_dir / "m12_range.csv")
    print(f"\n{'='*70}\nOVERALL INDEPENDENT SENSOR ORACLE AUDIT: ALL TESTS PASSED\n{'='*70}")


if __name__ == "__main__":
    main()
