# M24 High-Fidelity Environment / Time / Ephemeris — Validation Report

## Scope

`cpp/environment/high_fidelity.hpp` (rotation, JD time, Bowring geodetic,
J2–J4 gradient, Moon/Sun ephemeris), `tests/cpp/test_high_fidelity.cpp`
(7 cases / 45 assertions), demo (`tools/high_fidelity_demo.cpp`: 100-min LEO
telemetry), oracle (`python/audit/independent_high_fidelity_reference.py`),
telemetry (`data/m24_environment.csv`), figures
(`artifacts/figures/m24_*`), lesson (`docs/lessons/024_high_fidelity.md`).

## Analytical checks (all pass)

- Rotation identity/periodicity/orthonormality/round-trip/sign.
- JD at J2000 + one day; solar-vs-sidereal documented by assertion.
- Geodetic landmarks + 3-regime round-trips; polar axis; invalid rejected.
- Zonal J2-only = legacy + central to 1e-9 (4 positions); J3 pole push
  hand value; J4 mirrored vectors; FD potential audit 1e-8.
- Ephemeris radius/period/quadrature; equatorial-plane simplification stated.

## Independent verification

Pure-Python oracle (rotation, geodetic, potential); FD audit cross-checks.

## Scenario validation

100-min LEO: ground track spans longitudes, altitude holds 500 km, J2/J234
magnitudes tracked, fixed-Moon error grows (ephemeris justification).

## Development corrections

J3-antisymmetry fallacy; total-vs-perturbation API; altitude tolerance;
M11 demo preservation (restored from git, M24 separate target).

## Regression

M11 demo byte-identical to HEAD; new files only otherwise. Full suite
308/308. Ruff clean.

## Limitations

No leap seconds/EOP (not UTC); no lunar inclination/eccentricity; J3/J4
reference values (documented EGM-class); educational atmosphere retained.
