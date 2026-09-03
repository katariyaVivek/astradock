#!/usr/bin/env python3
"""Independent pure-Python docking reference oracle for AstraDock M19.

Cross-verifies port geometry, alignment, penalty contact, and acceptance logic
against independent models. Validation-only: never imports C++ code.
"""

import numpy as np


def quat_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Active rotation v' = q (*) v, scalar-first quaternion."""
    w, x, y, z = (float(c) for c in q)
    u = np.array([x, y, z])
    v = np.asarray(v, dtype=np.float64)
    return v + 2.0 * w * np.cross(u, v) + 2.0 * np.cross(u, np.cross(u, v))


def py_port_position(cm: np.ndarray, q: np.ndarray, offset: np.ndarray) -> np.ndarray:
    """Port world position P = R + q (*) r_body."""
    return np.asarray(cm, dtype=np.float64) + quat_rotate(np.asarray(q), np.asarray(offset))


def py_alignment_angle(q_t: np.ndarray, q_c: np.ndarray) -> float:
    """Relative attitude angle 2 acos(|w(q_t* x q_c)|)."""
    qt = np.asarray(q_t, dtype=np.float64)
    qc = np.asarray(q_c, dtype=np.float64)
    qc_conj = np.array([qt[0], -qt[1], -qt[2], -qt[3]])
    w = qc_conj[0] * qc[0] - qc_conj[1] * qc[1] - qc_conj[2] * qc[2] - qc_conj[3] * qc[3]
    return float(2.0 * np.arccos(min(abs(w), 1.0)))


def py_contact_force(pen: float, rate: float, k: float, c: float) -> float:
    """Penalty contact: k*pen + c*rate for pen > 0 (push only), else 0."""
    if pen <= 0.0:
        return 0.0
    return max(k * pen + c * rate, 0.0)


def py_acceptance(lateral: float, axial: float, closing: float, align: float, rate: float,
                  env: dict) -> bool:
    """All-criteria acceptance gate."""
    return (lateral <= env["lateral_limit"] and abs(axial) <= env["axial_limit"]
            and closing <= env["closing_limit"] and align <= env["align_limit"]
            and rate <= env["rate_limit"])
