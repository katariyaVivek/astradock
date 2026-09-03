#!/usr/bin/env python3
"""Independent pure-Python FDIR reference oracle for AstraDock M20.

Cross-verifies M-of-N voting, chi-square gating, isolation margins, and the
recovery policy. Validation-only: never imports C++ code.
"""

from collections import deque


def py_persist(history: list, exceed: bool, n: int, m: int) -> bool:
    """M-of-N trigger state after pushing one sample."""
    history.append(1 if exceed else 0)
    window = history[-n:]
    return len(history) >= m and sum(window) >= m


def py_margin(nis: float, gate: float) -> float:
    """Normalized exceedance margin (nis - gate) / gate."""
    return (nis - gate) / gate


def py_isolate(margins: dict, triggered: dict) -> tuple:
    """Worst-margin suspect + ambiguity flag."""
    active = {k: margins[k] for k in margins if triggered.get(k, False)}
    if not active:
        return "none", False
    suspect = max(active, key=lambda k: active[k])
    return suspect, len(active) > 1


def py_recovery(suspect: str, ambiguous: bool, actuator: bool, critical: bool) -> str:
    """Recovery policy table."""
    if critical or actuator:
        return "safe_mode"
    if suspect == "none":
        return "none"
    if ambiguous:
        return "coast_predict"
    return "exclude_channel"


class PyVoter:
    """Stateful M-of-N voter mirroring PersistenceVoter."""

    def __init__(self, n: int, m: int):
        self.n = n
        self.m = m
        self.history: deque = deque(maxlen=n)

    def push(self, exceed: bool) -> bool:
        self.history.append(exceed)
        return len(self.history) >= self.m and sum(self.history) >= self.m
