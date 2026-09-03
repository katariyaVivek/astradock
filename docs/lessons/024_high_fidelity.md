# Lesson 024 — High-Fidelity Environment, Time & Ephemeris

## 1. Physical Problem & Motivation

Through M23 the Earth is a spherical point mass in an idealized frame with a
fixed Moon. The autonomy stack (estimation, GNC, FDIR, vision) was built on
that clean room deliberately — upgrading gravity mid-way would have
contaminated every validation. M24 upgrades the WORLD the stack flies in,
each model verified independently before composition, and documents exactly
what is still NOT claimed (no UTC leap seconds, no EOP, no operational
atmosphere).

## 2. Time First (M24B)

One paragraph that everything else depends on: the clock `t_sim` (s) maps to
Earth angle `θ = ω_E · t` with `θ_0 = 0` (Greenwich ≡ inertial +X at t = 0)
and to Julian date `JD = 2451545.0 + t/86400`. A simulation time system — not
UTC/TT. Solar-vs-sidereal day difference is asserted in tests as the living
proof the distinction is understood.

## 3. Rotation, Geodetic, Zonals, Ephemeris

- **M24A**: ECI↔ECEF by principal Z-rotation; identity at epoch, sidereal-day
  periodicity, handedness regression (+X sweeps toward +Y ECEF).
- **M24C**: WGS-84 Bowring iteration with polar-axis closed form; landmarks
  (equator origin, pole, 500 km alt) + round-trips; ABSOLUTE altitude checks
  (relative tolerance is meaningless at 0 m).
- **M24D**: Legendre chain-rule gradient through J4. Returns TOTAL
  acceleration (documented contrast: legacy M11 J2 returns perturbation).
  J2-only reproduces legacy+central to 1e-9; J3 hand-derived pole push
  `4μJ3ρ³/r²`; finite-difference potential audit to 1e-8.
- **M24E**: circular coplanar Moon/Sun (inclination/eccentricity documented
  as NOT modeled); radius exact, period closure, quadrature phase.
- **M24F**: exponential atmosphere RETAINED with a written validity box —
  upgrading to NRLMSISE without a justifiable source would be fake fidelity.

## 4. Development Corrections (kept as findings)

1. "J3 antisymmetry" is wrong physics — odd zonals push both poles the same
   way; even zonals mirror vectors. Re-derived by hand, tests assert pushes.
2. Total-vs-perturbation API mismatch with legacy J2 — documented loudly,
   tests subtract central.
3. Relative-tolerance altitude check at 0 m — absolute checks.
4. M11 demo preserved untouched (`environment_demo.cpp` restored from git);
   M24 gets its own `high_fidelity_demo.cpp` + target. History is not
   overwritten.

## What you should now understand

1. Why upgrade the environment only AFTER the autonomy stack?
2. What is the simulator's time system, and what is it NOT?
3. Why does J3 push both poles the same direction?
4. Total vs perturbation: why does the distinction matter for propagators?
5. What makes an atmosphere model "operational" vs educational?
6. Why preserve the M11 demo instead of extending it?
