#!/usr/bin/env python3
"""Independent pure-Python vision reference oracle for AstraDock M22.

Cross-verifies pinhole projection, distortion, bearing inversion, and the
reprojection cost against independent models. Validation-only.
"""

import numpy as np


def py_project(p_cam: np.ndarray, fx: float, fy: float, cx: float, cy: float,
               k1: float = 0.0) -> tuple:
    """Pinhole + radial distortion projection to pixels."""
    p = np.asarray(p_cam, dtype=np.float64)
    xn, yn = p[0] / p[2], p[1] / p[2]
    s = 1.0 + k1 * (xn ** 2 + yn ** 2)
    return fx * s * xn + cx, fy * s * yn + cy


def py_bearing(u: float, v: float, fx: float, fy: float, cx: float, cy: float) -> np.ndarray:
    """Back-project pixel to normalized bearing (undistorted path)."""
    xn, yn = (u - cx) / fx, (v - cy) / fy
    b = np.array([xn, yn, 1.0])
    return b / np.linalg.norm(b)


def py_reprojection_cost(pixels: list, points: list, R: np.ndarray, t: np.ndarray,
                         fx: float, fy: float, cx: float, cy: float) -> float:
    """Sum of squared reprojection errors for a candidate pose."""
    total = 0.0
    for (u, v), p in zip(pixels, points):
        q = R @ (np.asarray(p) - np.asarray(t))
        pu, pv = py_project(q, fx, fy, cx, cy)
        total += (u - pu) ** 2 + (v - pv) ** 2
    return float(total)
