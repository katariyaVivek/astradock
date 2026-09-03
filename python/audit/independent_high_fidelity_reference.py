#!/usr/bin/env python3
"""Independent pure-Python high-fidelity reference oracle for AstraDock M24.

Cross-verifies rotation, geodetic, zonal, and ephemeris models. Validation-only.
"""

import numpy as np

OMEGA = 7.2921150e-5
A = 6378137.0
F = 1.0 / 298.257223563


def py_dcm_ecef_eci(t: float) -> np.ndarray:
    """Principal Z-rotation by omega*t."""
    c, s = np.cos(OMEGA * t), np.sin(OMEGA * t)
    return np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]])


def py_geodetic_to_ecef(lat: float, lon: float, alt: float) -> np.ndarray:
    """Closed-form geodetic to ECEF."""
    e2 = 2 * F - F * F
    prime = A / np.sqrt(1 - e2 * np.sin(lat) ** 2)
    return np.array([
        (prime + alt) * np.cos(lat) * np.cos(lon),
        (prime + alt) * np.cos(lat) * np.sin(lon),
        ((1 - e2) * prime + alt) * np.sin(lat),
    ])


def py_zonal_potential(r: np.ndarray, mu: float, R: float, j2: float, j3: float, j4: float) -> float:
    """Unnormalized zonal potential through J4."""
    rn = float(np.linalg.norm(r))
    s = r[2] / rn
    rho = R / rn
    p2 = 0.5 * (3 * s * s - 1)
    p3 = 0.5 * (5 * s ** 3 - 3 * s)
    p4 = 0.125 * (35 * s ** 4 - 30 * s * s + 3)
    return mu / rn * (1 - j2 * rho ** 2 * p2 - j3 * rho ** 3 * p3 - j4 * rho ** 4 * p4)
