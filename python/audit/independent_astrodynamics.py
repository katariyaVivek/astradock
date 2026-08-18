"""Independent pure-Python astrodynamics reference implementation for AstraDock audit.

This module is completely independent from AstraDock C++ source code.
It is an audit oracle used to verify C++ calculations, physical equations,
and numerical integration against independent reference implementations.
"""

from __future__ import annotations

import csv
import math
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Vec3:
    x: float
    y: float
    z: float

    def __add__(self, other: Vec3) -> Vec3:
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vec3) -> Vec3:
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, scalar: float) -> Vec3:
        return Vec3(self.x * scalar, self.y * scalar, self.z * scalar)

    def __rmul__(self, scalar: float) -> Vec3:
        return self.__mul__(scalar)

    def __truediv__(self, scalar: float) -> Vec3:
        if scalar == 0.0:
            raise ZeroDivisionError("Division by zero in Vec3")
        return Vec3(self.x / scalar, self.y / scalar, self.z / scalar)

    def dot(self, other: Vec3) -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def cross(self, other: Vec3) -> Vec3:
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def norm(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def normalized(self) -> Vec3:
        n = self.norm()
        if n == 0.0:
            raise ValueError("Cannot normalize zero vector")
        return self / n


@dataclass(frozen=True)
class State6:
    pos: Vec3
    vel: Vec3

    def __add__(self, other: State6) -> State6:
        return State6(self.pos + other.pos, self.vel + other.vel)

    def __mul__(self, scalar: float) -> State6:
        return State6(self.pos * scalar, self.vel * scalar)

    def __rmul__(self, scalar: float) -> State6:
        return self.__mul__(scalar)


# WGS-84 / Standard physical constants
MU_EARTH_M3_S2 = 3.986004418e14
R_EARTH_M = 6378137.0


def independent_two_body_accel(r: Vec3, mu: float = MU_EARTH_M3_S2) -> Vec3:
    """Compute gravitational acceleration: a = -mu * r / |r|^3."""
    dist = r.norm()
    if dist == 0.0:
        raise ValueError("Gravity undefined at origin")
    scale = -mu / (dist * dist * dist)
    return r * scale


def independent_state_deriv(t: float, state: State6, mu: float = MU_EARTH_M3_S2) -> State6:
    """dr/dt = v, dv/dt = a(r)."""
    del t  # Autonomous ODE
    return State6(pos=state.vel, vel=independent_two_body_accel(state.pos, mu))


def independent_rk4_step(
    t: float, state: State6, dt: float, mu: float = MU_EARTH_M3_S2
) -> State6:
    """One classical 4th-order Runge-Kutta step."""
    k1 = independent_state_deriv(t, state, mu)
    k2 = independent_state_deriv(t + 0.5 * dt, state + k1 * (0.5 * dt), mu)
    k3 = independent_state_deriv(t + 0.5 * dt, state + k2 * (0.5 * dt), mu)
    k4 = independent_state_deriv(t + dt, state + k3 * dt, mu)

    weighted = k1 + k2 * 2.0 + k3 * 2.0 + k4
    return state + weighted * (dt / 6.0)


def independent_euler_step(
    t: float, state: State6, dt: float, mu: float = MU_EARTH_M3_S2
) -> State6:
    """One explicit Forward Euler step."""
    slope = independent_state_deriv(t, state, mu)
    return state + slope * dt


def independent_propagate(
    state0: State6,
    t_end: float,
    dt: float,
    integrator: str = "rk4",
    mu: float = MU_EARTH_M3_S2,
) -> list[tuple[float, State6]]:
    """Propagates state from t=0 to t_end using fixed step dt."""
    t = 0.0
    state = state0
    trajectory = [(t, state)]

    step_fn = independent_rk4_step if integrator.lower() == "rk4" else independent_euler_step

    while t < t_end:
        step_dt = min(dt, t_end - t)
        state = step_fn(t, state, step_dt, mu)
        t += step_dt
        trajectory.append((t, state))

    return trajectory


def independent_energy(state: State6, mu: float = MU_EARTH_M3_S2) -> float:
    """Specific orbital energy: v^2/2 - mu/r."""
    r = state.pos.norm()
    v = state.vel.norm()
    return 0.5 * v * v - mu / r


def independent_angular_momentum(state: State6) -> Vec3:
    """Specific angular momentum vector: h = r x v."""
    return state.pos.cross(state.vel)


def independent_lvlh_basis(state: State6) -> tuple[Vec3, Vec3, Vec3]:
    """Compute LVLH basis (e_r, e_t, e_h) from orbital state.

    e_r = r / |r|
    e_h = (r x v) / |r x v|
    e_t = e_h x e_r
    """
    e_r = state.pos.normalized()
    h = state.pos.cross(state.vel)
    e_h = h.normalized()
    e_t = e_h.cross(e_r)
    return e_r, e_t, e_h


def independent_dcm_lvlh_from_eci(e_r: Vec3, e_t: Vec3, e_h: Vec3) -> list[list[float]]:
    """DCM whose rows are the basis vectors expressed in ECI."""
    return [
        [e_r.x, e_r.y, e_r.z],
        [e_t.x, e_t.y, e_t.z],
        [e_h.x, e_h.y, e_h.z],
    ]


def independent_dcm_mult_vec(dcm: list[list[float]], v: Vec3) -> Vec3:
    """Matrix-vector product."""
    return Vec3(
        dcm[0][0] * v.x + dcm[0][1] * v.y + dcm[0][2] * v.z,
        dcm[1][0] * v.x + dcm[1][1] * v.y + dcm[1][2] * v.z,
        dcm[2][0] * v.x + dcm[2][1] * v.y + dcm[2][2] * v.z,
    )


def independent_dcm_transpose_mult_vec(dcm: list[list[float]], v: Vec3) -> Vec3:
    """Transpose matrix-vector product (LVLH to ECI)."""
    return Vec3(
        dcm[0][0] * v.x + dcm[1][0] * v.y + dcm[2][0] * v.z,
        dcm[0][1] * v.x + dcm[1][1] * v.y + dcm[2][1] * v.z,
        dcm[0][2] * v.x + dcm[1][2] * v.y + dcm[2][2] * v.z,
    )


def run_cross_check(cpp_csv_path: Path | None = None) -> dict[str, float]:
    """Runs independent analytical calculations and compares against C++ if CSV is available."""
    alt_m = 500000.0
    r_orbit_m = R_EARTH_M + alt_m
    v_circ_m_s = math.sqrt(MU_EARTH_M3_S2 / r_orbit_m)
    period_s = 2.0 * math.pi * math.sqrt((r_orbit_m**3) / MU_EARTH_M3_S2)
    energy_analytical = -MU_EARTH_M3_S2 / (2.0 * r_orbit_m)
    h_analytical = r_orbit_m * v_circ_m_s

    print("=================================================================")
    print(" AstraDock — Independent Python Astrodynamics Reference & Audit ")
    print("=================================================================\n")
    print(f"Analytical Orbital Radius:      {r_orbit_m:.6f} m")
    print(f"Analytical Circular Speed:      {v_circ_m_s:.6f} m/s")
    print(f"Analytical Orbital Period:      {period_s:.6f} s ({period_s / 60.0:.4f} min)")
    print(f"Analytical Specific Energy:     {energy_analytical:.6f} m^2/s^2")
    print(f"Analytical Angular Momentum:    {h_analytical:.6f} m^2/s\n")

    state0 = State6(
        pos=Vec3(r_orbit_m, 0.0, 0.0),
        vel=Vec3(0.0, v_circ_m_s, 0.0),
    )

    # Propagate with independent RK4 (dt = 10.0 s)
    traj = independent_propagate(state0, period_s, dt=10.0, integrator="rk4")
    final_t, final_state = traj[-1]

    pos_err = (final_state.pos - state0.pos).norm()
    vel_err = (final_state.vel - state0.vel).norm()
    final_energy = independent_energy(final_state)
    energy_err = abs(final_energy - energy_analytical) / abs(energy_analytical)
    h_final = independent_angular_momentum(final_state).norm()
    h_err = abs(h_final - h_analytical) / h_analytical

    print("Independent Python RK4 (dt = 10 s, 1 Orbit):")
    print(f"  Samples:                      {len(traj)}")
    print(f"  Final Time:                   {final_t:.6f} s")
    print(f"  Position Closure Error:       {pos_err:.6f} m")
    print(f"  Velocity Closure Error:       {vel_err:.6e} m/s")
    print(f"  Relative Energy Error:        {energy_err:.6e}")
    print(f"  Relative Angular Momentum:    {h_err:.6e}\n")

    # If C++ CSV is available, compare sample by sample
    max_cpp_pos_diff = 0.0
    max_cpp_vel_diff = 0.0
    max_cpp_lvlh_diff = 0.0

    if cpp_csv_path and cpp_csv_path.is_file():
        print(f"Cross-checking against C++ CSV: {cpp_csv_path}")
        with cpp_csv_path.open(newline="", encoding="utf-8") as f:
            reader = csv.DictReader(f)
            for idx, row in enumerate(reader):
                cpp_pos = Vec3(
                    float(row["position_eci_x_m"]),
                    float(row["position_eci_y_m"]),
                    float(row["position_eci_z_m"]),
                )
                cpp_vel = Vec3(
                    float(row["velocity_eci_x_m_per_s"]),
                    float(row["velocity_eci_y_m_per_s"]),
                    float(row["velocity_eci_z_m_per_s"]),
                )
                cpp_er = Vec3(
                    float(row["lvlh_x_eci_x"]),
                    float(row["lvlh_x_eci_y"]),
                    float(row["lvlh_x_eci_z"]),
                )

                py_t, py_state = traj[idx]
                py_er, _, _ = independent_lvlh_basis(py_state)

                max_cpp_pos_diff = max(max_cpp_pos_diff, (py_state.pos - cpp_pos).norm())
                max_cpp_vel_diff = max(max_cpp_vel_diff, (py_state.vel - cpp_vel).norm())
                max_cpp_lvlh_diff = max(max_cpp_lvlh_diff, (py_er - cpp_er).norm())

        print(f"  Max C++ vs Python Position Diff: {max_cpp_pos_diff:.6e} m")
        print(f"  Max C++ vs Python Velocity Diff: {max_cpp_vel_diff:.6e} m/s")
        print(f"  Max C++ vs Python LVLH Diff:     {max_cpp_lvlh_diff:.6e}\n")

    return {
        "r_orbit_m": r_orbit_m,
        "v_circ_m_s": v_circ_m_s,
        "period_s": period_s,
        "pos_err_m": pos_err,
        "vel_err_m_s": vel_err,
        "energy_rel_err": energy_err,
        "h_rel_err": h_err,
        "max_cpp_pos_diff": max_cpp_pos_diff,
        "max_cpp_vel_diff": max_cpp_vel_diff,
    }


if __name__ == "__main__":
    csv_file = Path("artifacts/data/m06_frames.csv")
    run_cross_check(csv_file if csv_file.is_file() else None)
