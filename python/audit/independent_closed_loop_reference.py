#!/usr/bin/env python3
"""Independent pure-Python closed-loop reference oracle for AstraDock M17.

Cross-verifies the estimate-based control tick, timeline selection, and metrics
summaries against independent scalar models. Validation-only.
"""

import numpy as np
from independent_control_reference import py_attitude_error_vector, py_attitude_pd, py_saturate


def py_closed_loop_tick(q_est, w_est, q_ref, w_ref, kp, kd, limit) -> dict:
    """One estimate-based PD tick with saturation accounting."""
    desired = py_attitude_pd(q_est, q_ref, w_est, w_ref, kp, kd)
    achieved, saturated, deficit = py_saturate(desired, limit)
    return {
        "desired": desired,
        "achieved": achieved,
        "saturated": saturated,
        "deficit": deficit,
        "error_vector": py_attitude_error_vector(np.asarray(q_est), np.asarray(q_ref)),
    }


def py_reference_at(segments: list, t: float):
    """Piecewise timeline selection: latest segment with start <= t."""
    active = segments[0]
    for seg in segments:
        if seg[0] <= t:
            active = seg
    return active


def py_settle_time(times: np.ndarray, errors: np.ndarray, threshold: float) -> float:
    """Latching settle: first index from which all errors hold the bound."""
    for i in range(len(times)):
        if bool(np.all(np.asarray(errors[i:]) < threshold)):
            return float(times[i])
    return -1.0
