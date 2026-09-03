# M20 Fault Detection, Isolation & Recovery — Validation Report

## Scope

`cpp/fdir/fdir.hpp` (injector, voter, NIS/residual monitors, isolation,
recovery policy, metrics), `tests/cpp/test_fdir.cpp` (7 cases / 62
assertions), demo (`tools/fdir_demo.cpp`: bias jump, dropout, clean census,
thruster degradation), oracle (`python/audit/independent_fdir_reference.py`),
telemetry (`data/m20_*.csv`), figures (`artifacts/figures/m20_*`), lesson
(`docs/lessons/020_fdir.md`).

## Analytical checks (all pass)

- Voter: exact window tableaux incl. m == n boundary; invalid config throws.
- Gates equal M13 chi-square 95% bounds to 1e-4; single spike never triggers.
- Injector windows, bias/authority queries, stuck flags, name table.
- Isolation worst-margin + ambiguity + none-case.
- Recovery policy table (6 entries) incl. safe-mode dominance.
- In-loop 50 m bias: detected ≤ 10 s, zero pre-onset triggers.

## Independent verification

Pure-Python oracle (vote, margin, isolate, policy, stateful voter); 7 spot
checks pass.

## Scenario validation

- Bias jump: latency 2 s, suspect gnss, exclude_channel, 0 false alarms,
  final 3.6 m (coast on predict after exclusion).
- Dropout 60 s: no false trigger, coasts, reacquires, final 0 m-class.
- Clean 300 s: 0 false alarms (M-of-N vs ~15 raw).
- Thruster 50%: residual trigger 2 s after onset.

## Regression

New `cpp/fdir/` + tests/demo only. Full suite 284/284. Ruff clean.

## Limitations

Single-channel exclusion only (no multi-fault combinatorics); coast is
dynamics-only (no IMU bridging in this demo); safe mode declares but does not
execute a retreat (M21); no sensor-bias re-estimation post-exclusion.
