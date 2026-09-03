#!/usr/bin/env python3
"""Independent pure-Python control-law reference oracle for AstraDock M16.

Cross-verifies quaternion-error PD, rate damping, CW-feedforward translation
control, the double-integrator LQR closed form, and per-axis saturation against
independent scalar models. Validation-only: never imports C++ code.
"""

import numpy as np


def quat_conj(q: np.ndarray) -> np.ndarray:
    q = np.asarray(q, dtype=np.float64)
    return np.array([q[0], -q[1], -q[2], -q[3]])


def quat_mul(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    aw, ax, ay, az = (float(v) for v in a)
    bw, bx, by, bz = (float(v) for v in b)
    return np.array([
        aw * bw - ax * bx - ay * by - az * bz,
        aw * bx + ax * bw + ay * bz - az * by,
        aw * by - ax * bz + ay * bw + az * bx,
        aw * bz + ax * by - ay * bx + az * bw,
    ])


def py_attitude_error_vector(q_cur: np.ndarray, q_des: np.ndarray) -> np.ndarray:
    """Shortest-arc error vector: vector part of q_des* ⊗ q_cur, w >= 0."""
    q_err = quat_mul(quat_conj(q_des), q_cur)
    if q_err[0] < 0.0:
        q_err = -q_err
    return q_err[1:4]


def py_attitude_pd(q_cur, q_des, w, w_des, kp: np.ndarray, kd: np.ndarray) -> np.ndarray:
    """tau = -Kp .* e_att - Kd .* (w - w_des)."""
    e = py_attitude_error_vector(np.asarray(q_cur), np.asarray(q_des))
    return -np.asarray(kp) * e - np.asarray(kd) * (np.asarray(w) - np.asarray(w_des))


def py_rate_damping(w, w_cmd, kw: np.ndarray) -> np.ndarray:
    """tau = -Kw .* (w - w_cmd)."""
    return -np.asarray(kw) * (np.asarray(w) - np.asarray(w_cmd))


def py_cw_feedforward(rho: np.ndarray, vel: np.ndarray, n: float) -> np.ndarray:
    """Linearized natural acceleration [3n²x + 2nvy, -2nvx, -n²z]."""
    rho = np.asarray(rho, dtype=np.float64)
    vel = np.asarray(vel, dtype=np.float64)
    return np.array([3 * n * n * rho[0] + 2 * n * vel[1], -2 * n * vel[0], -n * n * rho[2]])


def py_lqr_double_integrator(q_pos: float, q_vel: float, r: float) -> tuple:
    """Closed-form optimal gains (k_pos, k_vel) for the double integrator."""
    p12 = np.sqrt(q_pos * r)
    p22 = np.sqrt(r * (2 * p12 + q_vel))
    return p12 / r, p22 / r


def py_saturate(desired: np.ndarray, limit: np.ndarray) -> tuple:
    """Per-axis clamp returning (achieved, saturated, deficit)."""
    d = np.asarray(desired, dtype=np.float64)
    lim = np.asarray(limit, dtype=np.float64)
    achieved = np.clip(d, -lim, lim)
    deficit = d - achieved
    return achieved, bool(np.any(deficit != 0.0)), deficit
