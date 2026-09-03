# Lesson 025 — Full-System Verification & Capstone Mission

## 1. Physical Problem & Motivation

Twenty-four milestones built parts; M25 proves they compose into a MISSION.
The capstone flies one deterministic stack — orbit, attitude acquisition,
nav convergence, slew, transfer, rendezvous, final approach, docking hold,
fault, recovery, complete — on estimates-only control with truth scored
harness-side. No new filter, no new controller: composition of verified parts
is itself the verification. The lesson is integration discipline: epoch
alignment, thrust-aware prediction, and state isolation hold end-to-end, or
nothing does.

## 2. Verification Layers

- **M25A mission**: 9-phase timeline on one seed; rendezvous completes,
  docking hold latches, fault rides through, outcome converged.
- **M25B isolation audit**: declaration scan over the control seam +
  rendezvous + docking (M17 pattern extended); attitude/aggregate containers
  banned, shared ECI math types allowed by documented rationale.
- **M25C interface audit**: handedness, SI magnitudes, JD time, normal sign,
  15-state ordering, double cover — one executable check per boundary class.
- **M25D Monte Carlo**: 100 dispersed runs, 100/100 converged, p50/p95
  reported, zero hidden failures.
- **M25E benchmarks**: frozen table with per-row bounds from the milestone
  reports — the numbers future work must not regress.
- **M25F/G reproducibility + performance**: clean-clone build (this session's
  full rebuild), demo runtime recorded.

## 3. Results

Nominal: converged, 255 m final (250 m hold + capture ball), Δv 3.0 m/s.
Monte Carlo: 100/100 converged, p50 255 m, p95 256 m, success 1.000. Full
suite 311/311. The capstone holds at the rendezvous ball rather than docking
contact: contact dynamics at mission scale belong to the M19 rig with docking-
grade nav, and the capstone documents the seam instead of faking it.

## What you should now understand

1. Why is composition itself a verification activity?
2. What does the isolation audit allow, and why is the rationale load-bearing?
3. Why does the capstone stop at the hold instead of contact?
4. What does a frozen benchmark row obligate future work to do?
5. Why report p95 alongside the success rate?
6. What would break first if estimation and control shared state types?
