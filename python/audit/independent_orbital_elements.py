"""Independent pure-Python orbital elements reference implementation for AstraDock audit.

This module is completely independent from AstraDock C++ source code.
It is an audit oracle used to verify C++ calculations, Keplerian orbital element
transformations, perifocal frame conversions, and round-trip state reconstruction.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

MU_EARTH_M3_S2 = 3.986004418e14
R_EARTH_M = 6378137.0


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
class KeplerianElements:
    semi_major_axis_m: float
    eccentricity: float
    inclination_rad: float
    raan_rad: float
    argument_of_periapsis_rad: float
    true_anomaly_rad: float


def normalize_2pi(angle_rad: float) -> float:
    two_pi = 2.0 * math.pi
    val = angle_rad % two_pi
    if val < 0.0:
        val += two_pi
    if val >= two_pi or abs(val) < 1e-15:
        val = 0.0
    return val


def independent_state_to_elements(
    r_vec: Vec3,
    v_vec: Vec3,
    mu: float = MU_EARTH_M3_S2,
    tol_e: float = 1e-11,
    tol_n: float = 1e-11,
) -> KeplerianElements:
    r = r_vec.norm()
    v = v_vec.norm()

    h_vec = r_vec.cross(v_vec)
    h = h_vec.norm()
    h_hat = h_vec.normalized()

    energy = 0.5 * v * v - mu / r
    if energy >= 0.0:
        raise ValueError("Only bound elliptic orbits supported (energy < 0)")
    a = -mu / (2.0 * energy)

    e_vec = (v_vec.cross(h_vec) / mu) - (r_vec / r)
    e = e_vec.norm()

    cos_i = max(-1.0, min(1.0, h_vec.z / h))
    i = math.acos(cos_i)

    n_vec = Vec3(-h_vec.y, h_vec.x, 0.0)
    n = n_vec.norm()

    is_circular = e < tol_e
    is_equatorial = n < tol_n

    if not is_equatorial:
        raan = normalize_2pi(math.atan2(n_vec.y, n_vec.x))
    else:
        raan = 0.0

    if not is_circular and not is_equatorial:
        sin_omega = n_vec.cross(e_vec).dot(h_hat)
        cos_omega = n_vec.dot(e_vec)
        omega = normalize_2pi(math.atan2(sin_omega, cos_omega))
    elif not is_circular and is_equatorial:
        omega = normalize_2pi(math.atan2(e_vec.y, e_vec.x))
    else:
        omega = 0.0

    if not is_circular:
        sin_nu = e_vec.cross(r_vec).dot(h_hat)
        cos_nu = e_vec.dot(r_vec)
        nu = normalize_2pi(math.atan2(sin_nu, cos_nu))
    elif not is_equatorial:
        sin_u = n_vec.cross(r_vec).dot(h_hat)
        cos_u = n_vec.dot(r_vec)
        nu = normalize_2pi(math.atan2(sin_u, cos_u))
    else:
        nu = normalize_2pi(math.atan2(r_vec.y, r_vec.x))

    return KeplerianElements(
        semi_major_axis_m=a,
        eccentricity=e,
        inclination_rad=i,
        raan_rad=raan,
        argument_of_periapsis_rad=omega,
        true_anomaly_rad=nu,
    )


def independent_elements_to_state(
    elem: KeplerianElements, mu: float = MU_EARTH_M3_S2
) -> tuple[Vec3, Vec3]:
    p = elem.semi_major_axis_m * (1.0 - elem.eccentricity * elem.eccentricity)
    cos_nu = math.cos(elem.true_anomaly_rad)
    sin_nu = math.sin(elem.true_anomaly_rad)

    r_mag = p / (1.0 + elem.eccentricity * cos_nu)
    v_coeff = math.sqrt(mu / p)

    r_p = r_mag * cos_nu
    r_q = r_mag * sin_nu

    v_p = v_coeff * (-sin_nu)
    v_q = v_coeff * (elem.eccentricity + cos_nu)

    cos_raan = math.cos(elem.raan_rad)
    sin_raan = math.sin(elem.raan_rad)
    cos_inc = math.cos(elem.inclination_rad)
    sin_inc = math.sin(elem.inclination_rad)
    cos_argp = math.cos(elem.argument_of_periapsis_rad)
    sin_argp = math.sin(elem.argument_of_periapsis_rad)

    # Unit vectors P and Q expressed in ECI:
    p_vec = Vec3(
        cos_raan * cos_argp - sin_raan * sin_argp * cos_inc,
        sin_raan * cos_argp + cos_raan * sin_argp * cos_inc,
        sin_argp * sin_inc,
    )

    q_vec = Vec3(
        -cos_raan * sin_argp - sin_raan * cos_argp * cos_inc,
        -sin_raan * sin_argp + cos_raan * cos_argp * cos_inc,
        cos_argp * sin_inc,
    )

    pos_eci = p_vec * r_p + q_vec * r_q
    vel_eci = p_vec * v_p + q_vec * v_q

    return pos_eci, vel_eci


def run_scenarios() -> None:
    scenarios = [
        (
            "Case A: Circular Equatorial",
            KeplerianElements(
                semi_major_axis_m=7000000.0,
                eccentricity=0.0,
                inclination_rad=0.0,
                raan_rad=0.0,
                argument_of_periapsis_rad=0.0,
                true_anomaly_rad=math.radians(45.0),
            ),
        ),
        (
            "Case B: Circular Inclined (i = 45 deg)",
            KeplerianElements(
                semi_major_axis_m=7000000.0,
                eccentricity=0.0,
                inclination_rad=math.radians(45.0),
                raan_rad=math.radians(60.0),
                argument_of_periapsis_rad=0.0,
                true_anomaly_rad=math.radians(90.0),
            ),
        ),
        (
            "Case C: Elliptical Equatorial (e = 0.2)",
            KeplerianElements(
                semi_major_axis_m=10000000.0,
                eccentricity=0.2,
                inclination_rad=0.0,
                raan_rad=0.0,
                argument_of_periapsis_rad=math.radians(30.0),
                true_anomaly_rad=math.radians(60.0),
            ),
        ),
        (
            "Case D: Elliptical Inclined (a = 10,000 km, e = 0.2, i = 45 deg, RAAN = 120 deg, argp = 60 deg, nu = 30 deg)",
            KeplerianElements(
                semi_major_axis_m=10000000.0,
                eccentricity=0.2,
                inclination_rad=math.radians(45.0),
                raan_rad=math.radians(120.0),
                argument_of_periapsis_rad=math.radians(60.0),
                true_anomaly_rad=math.radians(30.0),
            ),
        ),
        (
            "Case E: High-Inclination Elliptical (i = 98 deg Sun-Sync, e = 0.15)",
            KeplerianElements(
                semi_major_axis_m=8000000.0,
                eccentricity=0.15,
                inclination_rad=math.radians(98.0),
                raan_rad=math.radians(45.0),
                argument_of_periapsis_rad=math.radians(270.0),
                true_anomaly_rad=math.radians(150.0),
            ),
        ),
        (
            "Case F: Near-Circular (e = 1e-5)",
            KeplerianElements(
                semi_major_axis_m=7000000.0,
                eccentricity=1e-5,
                inclination_rad=math.radians(30.0),
                raan_rad=math.radians(15.0),
                argument_of_periapsis_rad=math.radians(45.0),
                true_anomaly_rad=math.radians(120.0),
            ),
        ),
    ]

    print("=================================================================")
    print(" AstraDock — Independent Python Orbital Elements Oracle ")
    print("=================================================================\n")

    for name, elem in scenarios:
        pos, vel = independent_elements_to_state(elem)
        elem_reconstructed = independent_state_to_elements(pos, vel)
        pos_roundtrip, vel_roundtrip = independent_elements_to_state(elem_reconstructed)

        pos_err = (pos_roundtrip - pos).norm()
        vel_err = (vel_roundtrip - vel).norm()

        print(f"Scenario: {name}")
        print(f"  Input Elements: a={elem.semi_major_axis_m:.1f} m, e={elem.eccentricity:.6f}, i={math.degrees(elem.inclination_rad):.2f} deg, RAAN={math.degrees(elem.raan_rad):.2f} deg, w={math.degrees(elem.argument_of_periapsis_rad):.2f} deg, nu={math.degrees(elem.true_anomaly_rad):.2f} deg")
        print(f"  ECI Position:   [{pos.x:.3f}, {pos.y:.3f}, {pos.z:.3f}] m (norm: {pos.norm():.3f} m)")
        print(f"  ECI Velocity:   [{vel.x:.3f}, {vel.y:.3f}, {vel.z:.3f}] m/s (norm: {vel.norm():.3f} m/s)")
        print(f"  Reconstructed:  a={elem_reconstructed.semi_major_axis_m:.1f} m, e={elem_reconstructed.eccentricity:.6f}, i={math.degrees(elem_reconstructed.inclination_rad):.2f} deg, RAAN={math.degrees(elem_reconstructed.raan_rad):.2f} deg, w={math.degrees(elem_reconstructed.argument_of_periapsis_rad):.2f} deg, nu={math.degrees(elem_reconstructed.true_anomaly_rad):.2f} deg")
        print(f"  Round-Trip Pos Error: {pos_err:.6e} m")
        print(f"  Round-Trip Vel Error: {vel_err:.6e} m/s\n")


if __name__ == "__main__":
    run_scenarios()
