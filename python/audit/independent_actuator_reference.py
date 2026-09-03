#!/usr/bin/env python3
"""Independent pure-Python actuator reference oracle for AstraDock Milestone M14.

Cross-verifies C++ reaction-wheel, thruster, and propellant implementations
against independent scalar Python models. Validation-only: never imports C++ code.
"""

import numpy as np

G0 = 9.80665  # exact definitional standard gravity (m/s^2)


def py_wheel_speed(speed0: float, torque: float, dt: float, inertia: float) -> float:
    """Exact constant-torque wheel speed update: w += tau * dt / I."""
    return speed0 + torque * dt / inertia


def py_motor_lag(prev: float, target: float, dt: float, tau_lag: float) -> float:
    """Exact discrete first-order motor lag step."""
    if tau_lag == 0.0:
        return target
    return prev + (target - prev) * (1.0 - np.exp(-dt / tau_lag))


def py_saturate(cmd: float, limit: float) -> tuple:
    """Hard torque limit returning (limited, saturated)."""
    if cmd > limit:
        return limit, True
    if cmd < -limit:
        return -limit, True
    return cmd, False


def py_thruster_force(direction: np.ndarray, thrust: float) -> np.ndarray:
    """Body-frame thruster force: F = dir * T."""
    return np.asarray(direction, dtype=np.float64) * thrust


def py_thruster_torque(mount: np.ndarray, direction: np.ndarray, thrust: float) -> np.ndarray:
    """Body-frame thruster torque: tau = r x F."""
    return np.cross(np.asarray(mount, dtype=np.float64), py_thruster_force(direction, thrust))


def py_mass_flow(thrust: float, isp: float) -> float:
    """Propellant mass flow: m_dot = -T / (Isp * g0)."""
    return -thrust / (isp * G0)


def quat_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Active quaternion rotation v' = q (*) v (scalar-first [w, x, y, z])."""
    w, x, y, z = q
    u = np.array([x, y, z])
    return v + 2.0 * w * np.cross(u, v) + 2.0 * np.cross(u, np.cross(u, v))
