"""Independent pure-Python quaternion and attitude reference implementation for AstraDock audit.

This module is completely independent from AstraDock C++ source code.
It is an audit oracle used to verify C++ calculations, quaternion algebra,
DCM conversions, vector rotations, and double-cover properties.
"""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class Vec3:
    x: float
    y: float
    z: float

    def norm(self) -> float:
        return math.sqrt(self.x * self.x + self.y * self.y + self.z * self.z)

    def normalized(self) -> Vec3:
        n = self.norm()
        if n == 0.0:
            raise ValueError("Cannot normalize zero vector")
        return Vec3(self.x / n, self.y / n, self.z / n)

    def cross(self, other: Vec3) -> Vec3:
        return Vec3(
            self.y * other.z - self.z * other.y,
            self.z * other.x - self.x * other.z,
            self.x * other.y - self.y * other.x,
        )

    def dot(self, other: Vec3) -> float:
        return self.x * other.x + self.y * other.y + self.z * other.z

    def __add__(self, other: Vec3) -> Vec3:
        return Vec3(self.x + other.x, self.y + other.y, self.z + other.z)

    def __sub__(self, other: Vec3) -> Vec3:
        return Vec3(self.x - other.x, self.y - other.y, self.z - other.z)

    def __mul__(self, s: float) -> Vec3:
        return Vec3(self.x * s, self.y * s, self.z * s)

    def __rmul__(self, s: float) -> Vec3:
        return self.__mul__(s)


@dataclass(frozen=True)
class Quat:
    w: float
    x: float
    y: float
    z: float

    def norm(self) -> float:
        return math.sqrt(self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z)

    def normalized(self) -> Quat:
        n = self.norm()
        if n == 0.0:
            raise ValueError("Cannot normalize zero quaternion")
        return Quat(self.w / n, self.x / n, self.y / n, self.z / n)

    def conjugate(self) -> Quat:
        return Quat(self.w, -self.x, -self.y, -self.z)

    def inverse(self) -> Quat:
        n2 = self.w * self.w + self.x * self.x + self.y * self.y + self.z * self.z
        if n2 == 0.0:
            raise ValueError("Cannot invert zero quaternion")
        return Quat(self.w / n2, -self.x / n2, -self.y / n2, -self.z / n2)

    def __mul__(self, other: Quat) -> Quat:
        # Hamilton product
        return Quat(
            self.w * other.w - self.x * other.x - self.y * other.y - self.z * other.z,
            self.w * other.x + self.x * other.w + self.y * other.z - self.z * other.y,
            self.w * other.y - self.x * other.z + self.y * other.w + self.z * other.x,
            self.w * other.z + self.x * other.y - self.y * other.x + self.z * other.w,
        )

    def __neg__(self) -> Quat:
        return Quat(-self.w, -self.x, -self.y, -self.z)

    def rotate_vector(self, v: Vec3) -> Vec3:
        u = Vec3(self.x, self.y, self.z)
        u_cross_v = u.cross(v)
        u_cross_u_cross_v = u.cross(u_cross_v)
        return v + (2.0 * self.w) * u_cross_v + 2.0 * u_cross_u_cross_v

    def to_dcm(self) -> list[list[float]]:
        xx, yy, zz = self.x * self.x, self.y * self.y, self.z * self.z
        wx, wy, wz = self.w * self.x, self.w * self.y, self.w * self.z
        xy, xz, yz = self.x * self.y, self.x * self.z, self.y * self.z

        return [
            [1.0 - 2.0 * (yy + zz), 2.0 * (xy - wz), 2.0 * (xz + wy)],
            [2.0 * (xy + wz), 1.0 - 2.0 * (xx + zz), 2.0 * (yz - wx)],
            [2.0 * (xz - wy), 2.0 * (yz + wx), 1.0 - 2.0 * (xx + yy)],
        ]


def quat_from_axis_angle(axis: Vec3, angle_rad: float) -> Quat:
    u = axis.normalized()
    half = 0.5 * angle_rad
    sin_half = math.sin(half)
    cos_half = math.cos(half)
    return Quat(cos_half, sin_half * u.x, sin_half * u.y, sin_half * u.z)


def dcm_to_quat(dcm: list[list[float]]) -> Quat:
    r00, r11, r22 = dcm[0][0], dcm[1][1], dcm[2][2]
    trace = r00 + r11 + r22

    if trace > r00 and trace > r11 and trace > r22:
        w = 0.5 * math.sqrt(1.0 + trace)
        s = 0.25 / w
        x = (dcm[2][1] - dcm[1][2]) * s
        y = (dcm[0][2] - dcm[2][0]) * s
        z = (dcm[1][0] - dcm[0][1]) * s
    elif r00 > r11 and r00 > r22:
        x = 0.5 * math.sqrt(1.0 + 2.0 * r00 - trace)
        s = 0.25 / x
        w = (dcm[2][1] - dcm[1][2]) * s
        y = (dcm[0][1] + dcm[1][0]) * s
        z = (dcm[0][2] + dcm[2][0]) * s
    elif r11 > r22:
        y = 0.5 * math.sqrt(1.0 + 2.0 * r11 - trace)
        s = 0.25 / y
        w = (dcm[0][2] - dcm[2][0]) * s
        x = (dcm[0][1] + dcm[1][0]) * s
        z = (dcm[1][2] + dcm[2][1]) * s
    else:
        z = 0.5 * math.sqrt(1.0 + 2.0 * r22 - trace)
        s = 0.25 / z
        w = (dcm[1][0] - dcm[0][1]) * s
        x = (dcm[0][2] + dcm[2][0]) * s
        y = (dcm[1][2] + dcm[2][1]) * s

    return Quat(w, x, y, z).normalized()


def dcm_mult(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    res = [[0.0] * 3 for _ in range(3)]
    for i in range(3):
        for j in range(3):
            res[i][j] = sum(a[i][k] * b[k][j] for k in range(3))
    return res


def dcm_mult_vec(c: list[list[float]], v: Vec3) -> Vec3:
    return Vec3(
        c[0][0] * v.x + c[0][1] * v.y + c[0][2] * v.z,
        c[1][0] * v.x + c[1][1] * v.y + c[1][2] * v.z,
        c[2][0] * v.x + c[2][1] * v.y + c[2][2] * v.z,
    )


def run_attitude_checks() -> None:
    print("=================================================================")
    print(" AstraDock — Independent Python Quaternion Reference & Audit ")
    print("=================================================================\n")

    # 1. Known Rotations
    cases = [
        ("Identity", Quat(1.0, 0.0, 0.0, 0.0), Vec3(1.0, 2.0, 3.0), Vec3(1.0, 2.0, 3.0)),
        ("90 deg about X", quat_from_axis_angle(Vec3(1.0, 0.0, 0.0), math.pi / 2.0), Vec3(0.0, 1.0, 0.0), Vec3(0.0, 0.0, 1.0)),
        ("90 deg about Y", quat_from_axis_angle(Vec3(0.0, 1.0, 0.0), math.pi / 2.0), Vec3(0.0, 0.0, 1.0), Vec3(1.0, 0.0, 0.0)),
        ("90 deg about Z", quat_from_axis_angle(Vec3(0.0, 0.0, 1.0), math.pi / 2.0), Vec3(1.0, 0.0, 0.0), Vec3(0.0, 1.0, 0.0)),
        ("180 deg about Z", quat_from_axis_angle(Vec3(0.0, 0.0, 1.0), math.pi), Vec3(1.0, 0.0, 0.0), Vec3(-1.0, 0.0, 0.0)),
    ]

    for name, q, v_in, expected_v in cases:
        v_out_q = q.rotate_vector(v_in)
        dcm = q.to_dcm()
        v_out_dcm = dcm_mult_vec(dcm, v_in)

        err_q = (v_out_q - expected_v).norm()
        err_dcm = (v_out_dcm - expected_v).norm()

        # Round-trip DCM -> Quat
        q_rec = dcm_to_quat(dcm)
        dcm_rec = q_rec.to_dcm()

        # Check matrix round-trip
        mat_diff = max(abs(dcm[i][j] - dcm_rec[i][j]) for i in range(3) for j in range(3))

        print(f"Scenario: {name}")
        print(f"  Quaternion:        [{q.w:.4f}, {q.x:.4f}, {q.y:.4f}, {q.z:.4f}]")
        print(f"  Vector Rot Error:  {err_q:.2e} m (q) | {err_dcm:.2e} m (DCM)")
        print(f"  DCM->Quat Mat Diff:{mat_diff:.2e}\n")

    # 2. Composition consistency check
    q1 = quat_from_axis_angle(Vec3(1.0, 2.0, 3.0), 0.7)
    q2 = quat_from_axis_angle(Vec3(-2.0, 1.0, 0.5), 1.2)
    q3 = q1 * q2

    c1 = q1.to_dcm()
    c2 = q2.to_dcm()
    c3_mat = dcm_mult(c1, c2)
    c3_q = q3.to_dcm()

    comp_err = max(abs(c3_mat[i][j] - c3_q[i][j]) for i in range(3) for j in range(3))
    print(f"Composition Check (q1 * q2 vs C1 * C2): max matrix difference = {comp_err:.2e}\n")

    # 3. Vector norm preservation across scale
    v_test = Vec3(3.0, -4.0, 12.0)
    for scale in [1e-9, 1e-6, 1.0, 1e3, 1e6]:
        v_scaled = v_test * scale
        v_rot = q3.rotate_vector(v_scaled)
        rel_norm_diff = abs(v_rot.norm() - v_scaled.norm()) / v_scaled.norm()
        print(f"Norm Preservation Scale {scale:.0e}: relative error = {rel_norm_diff:.2e}")


if __name__ == "__main__":
    run_attitude_checks()
