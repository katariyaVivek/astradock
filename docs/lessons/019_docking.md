# Lesson 019 — Autonomous Docking

## 1. Physical Problem & Motivation

Rendezvous parks the chaser at 50 m. Docking flies the last meters into
physical contact: ports must meet within centimeters, alignment within degrees,
closing speed within cm/s — then a compliant mechanism absorbs residual energy
and latches. The frame problem dominates: ports live on the vehicle bodies,
the approach flies in LVLH, and a sign flip in the approach axis is a mission
failure. Three development bugs carried the teaching weight:

1. **Axial sign**: separation-positive vs penetration-positive. With the axis
   pointing target→chaser, `axial = e·u > 0` means APART. Contact at `s ≤ 0`,
   penetration `−s`, closing `−ds/dt`. One axis choice, three sign sites —
   all pinned by hand-derived regressions.
2. **Port-offset lateral**: body port offsets (±1 m) read as 2 m of port-gap
   lateral at long range on a perfectly nominal approach. The lateral abort is
   therefore range-gated (≤ 10 m); far-field lateral safety belongs to the M18
   corridor. Envelopes apply where their geometry lives.
3. **Contact-unaware guidance + coarse nav**: 5 m GNSS noise vs a 0.25 m
   envelope limit-cycles past the plane into crush (500 N). Docking declares a
   sensor requirement — 0.05 m relative-grade sensing — and freezes the
   approach command inside contact range so the spring settles.

## 2. Docking Frames and State

Target/chaser ports: body offset + body approach axis. World positions via
`P = R + q(*)r`. Port gap in LVLH; axial/lateral split about the LVLH approach
axis; relative attitude `q_t* ⊗ q_c` with angle `φ`; relative rate as a
conservative norm bound (frames differ — documented limitation).

## 3. Contact, Acceptance, Abort

Axial spring-damper penalty (`k = 1000 N/m`, `c = 100 N·s/m`), push-only.
Acceptance is a five-way AND (lateral, axial, closing, alignment, rate) plus
crush check; latch needs 5 s of continuous acceptance. Abort on lateral
(gated), 3×-envelope closing, 15° attitude, crush force, timeout — priority
ordered with reason codes. Retreat trajectories are M21 scope.

## 4. Verification

7 cases / 44 assertions: port frame hand cases (identity, 90° yaw), alignment
+ double cover, penalty push/no-pull, per-criterion acceptance rejection,
per-reason aborts, latch banking, axial sign. Oracle (5 functions). Demo:
nominal latches; 3× cruise latches (envelope-scaled, correctly); 2 m offset
recovered; 20° tilt aborts `attitude_exceeded` on tick 1; dropout latches.
Two figures.

## What you should now understand

1. Why does separation-positive force three coordinated sign choices?
2. Why is the lateral abort range-gated?
3. What sensor grade does a 0.25 m envelope demand, and why?
4. Why freeze guidance inside contact range?
5. Why is relative rate a norm bound, and when is that conservative?
6. Why can a hot-but-tracked approach latch while a runaway aborts?
