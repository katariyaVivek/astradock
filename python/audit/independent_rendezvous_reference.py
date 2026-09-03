#!/usr/bin/env python3
"""Independent pure-Python rendezvous guidance reference oracle for AstraDock M18.

Cross-verifies safety geometry, leg sequencing, and the closing-speed profile
against independent scalar models. Validation-only: never imports C++ code.
"""

import numpy as np


def py_lateral_error(rho: np.ndarray) -> float:
    """Off-axis distance sqrt(x^2 + z^2)."""
    rho = np.asarray(rho, dtype=np.float64)
    return float(np.sqrt(rho[0] ** 2 + rho[2] ** 2))


def py_closing_speed(rho: np.ndarray, vel: np.ndarray) -> float:
    """Inward closing speed -(rho.v)/|rho|."""
    rho = np.asarray(rho, dtype=np.float64)
    vel = np.asarray(vel, dtype=np.float64)
    r = float(np.linalg.norm(rho))
    if r == 0.0:
        return 0.0
    return float(-np.dot(rho, vel) / r)


def py_profile_speed(distance: float, cruise: float, brake: float) -> float:
    """Closing-speed profile min(cruise, sqrt(2*brake*distance))."""
    return min(cruise, np.sqrt(2.0 * brake * max(distance, 0.0)))


def py_leg_advance(distance: float, capture: float) -> bool:
    """Capture predicate for leg advance."""
    return distance < capture
