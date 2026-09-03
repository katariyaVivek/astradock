#!/usr/bin/env python3
"""Independent pure-Python mission-framework reference oracle for AstraDock M21.

Cross-verifies seed derivation, percentile interpolation, Welford statistics,
and the recovery-free policy table (pass/fail directions). Validation-only.
"""

import numpy as np


def py_derive_seed(base: int, index: int, stride: int = 1000003) -> int:
    """Per-run seed = base * stride + index."""
    return base * stride + index


def py_percentile(sorted_vals: list, p: float) -> float:
    """Linear-interpolation percentile matching summarize_runs."""
    n = len(sorted_vals)
    rank = p * (n - 1)
    lo = int(rank)
    hi = min(lo + 1, n - 1)
    frac = rank - lo
    return sorted_vals[lo] * (1.0 - frac) + sorted_vals[hi] * frac


def py_welford(values: list) -> tuple:
    """Online mean + sample std (Welford)."""
    mean = 0.0
    m2 = 0.0
    for i, v in enumerate(values):
        delta = v - mean
        mean += delta / (i + 1)
        m2 += delta * (v - mean)
    std = np.sqrt(m2 / (len(values) - 1)) if len(values) > 1 else 0.0
    return mean, std


def py_passes(value: float, bound: float, higher_is_better: bool) -> bool:
    """Regression pass/fail in either direction."""
    return value >= bound if higher_is_better else value <= bound
