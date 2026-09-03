# Lesson 020 — Fault Detection, Isolation & Recovery

## 1. Physical Problem & Motivation

Sensors lie and actuators degrade mid-mission: GNSS jumps 50 m, star trackers
blind, gyros jump bias, thrusters deliver half thrust, wheels stick. Classical
FDIR catches these with signals the loop already computes — no ML (M23), no
new sensors. The chain is inject → detect → isolate → recover, with roles
never collapsed: the injector lives harness-side (between sensor and filter,
where hardware faults live); detection sees estimates + measurements only.

## 2. Detection: NIS Gating with M-of-N Teeth

A 95% chi-square gate alone false-alarms 5% (~1/20 samples — 15 per 300 s).
The M-of-N persistence voter (3-of-5) is what makes detection reliable: a
single spike never triggers; sustained exceedance triggers in ≤ 3 samples.
Thresholds are chi-square 95% bounds per channel DOF (GNSS-6: 14.449, ST-3:
9.348, range-1: 5.024) — physics numbers, not tuned constants. Actuator
tracking uses desired-vs-achieved residuals with the same discipline.

## 3. Isolation and Recovery

Per-channel monitors vote; worst normalized margin wins suspect; multiple
triggers report ambiguity honestly instead of forcing certainty. Recovery is a
pure policy table: single measurement fault → exclude the channel; ambiguous →
coast on predict; actuator fault or critical geometry → safe mode
(hold + null rates); clean → none.

## 4. Verification

7 cases / 62 assertions: exact voter window arithmetic, gate values,
single-spike-vs-sustained, schedule windows, worst-margin + ambiguity,
policy table, in-loop 50 m bias detection ≤ 10 s with zero false alarms.
Oracle (voter, margin, isolation, policy). Demo: bias detected in 2 s,
excluded, final error 3.6 m; dropout coasts cleanly; 300 s clean run has
ZERO false alarms (raw gating would give ~15); thruster degradation in 2 s.

## 5. Development Bugs (kept as findings)

1. Voter test expectations miscounted window occupancy — the voter was right,
   the test was wrong (window still holds old exceeds). Rewrote with exact
   window tableaux.
2. `CHECK_THROWS_AS(a && (b, c))` trips Catch2's `&&` static assert — split
   into two statements.
3. Static-truth-vs-orbital-filter scenario diverged to 1e6 m — scenarios must
   propagate truth on a real orbit (same class as the M18 epoch bug).

## What you should now understand

1. Why does a 95% gate alone false-alarm, and what fixes it?
2. Why M-of-N instead of a higher threshold?
3. Why does isolation report ambiguity instead of picking?
4. Why does exclusion beat retuning after a bias jump?
5. What does a single NIS spike mean — and what should you do about it?
6. Why must fault scenarios propagate truth dynamically?
