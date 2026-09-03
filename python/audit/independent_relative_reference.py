#!/usr/bin/env python3
"""Independent pure-Python relative-dynamics reference oracle for AstraDock M15.

Cross-verifies C++ LVLH rate, transport theorem, relative-state, and CW
implementations against independent pure-Python models. Validation-only.
"""

import numpy as np


def py_lvlh_rate(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    """LVLH angular velocity: omega = (r x v) / |r|^2."""
    r = np.asarray(r, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    return np.cross(r, v) / float(np.dot(r, r))


def py_transport_first(v: np.ndarray, omega: np.ndarray, d_rot: np.ndarray) -> np.ndarray:
    """First-order transport: (dv/dt)_inertial = d_rot + omega x v."""
    return np.asarray(d_rot, dtype=np.float64) + np.cross(
        np.asarray(omega, dtype=np.float64), np.asarray(v, dtype=np.float64)
    )


def py_lvlh_basis(r: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rows of C_LVLH_ECI: [e_r; e_t; e_h] (AstraDock x-radial/y-along-track/z-normal)."""
    r = np.asarray(r, dtype=np.float64)
    v = np.asarray(v, dtype=np.float64)
    e_r = r / np.linalg.norm(r)
    h = np.cross(r, v)
    e_h = h / np.linalg.norm(h)
    e_t = np.cross(e_h, e_r)
    return np.array([e_r, e_t, e_h])


def py_relative_from_eci(r_tgt, v_tgt, r_ch, v_ch) -> tuple:
    """LVLH relative state from ECI pairs (rotating-frame velocity)."""
    c = py_lvlh_basis(r_tgt, v_tgt)
    rho = np.asarray(r_ch) - np.asarray(r_tgt)
    rho_lvlh = c @ rho
    omega_lvlh = c @ py_lvlh_rate(r_tgt, v_tgt)
    v_rel = c @ (np.asarray(v_ch) - np.asarray(v_tgt)) - np.cross(omega_lvlh, rho_lvlh)
    return rho_lvlh, v_rel


def py_cw_predict(rel0: np.ndarray, n: float, t: float) -> np.ndarray:
    """Exact unforced CW prediction for X = [x y z vx vy vz]."""
    x0, y0, z0, vx0, vy0, vz0 = (float(v) for v in rel0)
    c, s = np.cos(n * t), np.sin(n * t)
    x = (4 - 3 * c) * x0 + (s / n) * vx0 + (2 / n) * (1 - c) * vy0
    y = 6 * (s - n * t) * x0 + y0 - (2 / n) * (1 - c) * vx0 + ((4 * s - 3 * n * t) / n) * vy0
    z = c * z0 + (s / n) * vz0
    vx = 3 * n * s * x0 + c * vx0 + 2 * s * vy0
    vy = -6 * n * (1 - c) * x0 - 2 * s * vx0 + (4 * c - 3) * vy0
    vz = -n * s * z0 + c * vz0
    return np.array([x, y, z, vx, vy, vz])
