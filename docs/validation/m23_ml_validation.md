# M23 Machine Learning for Space Systems — Validation Report

## Scope

Dataset generator (`tools/ml_dataset_demo.cpp`: seeded EKF NIS streams,
2880/960/960 windows, seed-family splits), protocol
(`python/analysis/train_anomaly_classifier.py`: audit → baselines → selection
→ test → advisory table → latency), oracle (leakage + metric checks), results
(`data/m23_*.csv`), figure (`artifacts/figures/m23_metrics.png`), lesson
(`docs/lessons/023_machine_learning.md`).

## Results (held-out seeds, 10 m bias)

| Model | Precision | Recall | F1 | FAR | Miss |
|---|---|---|---|---|---|
| m20_rule_3of5 | 1.000 | 0.013 | 0.025 | 0.000 | 0.988 |
| logreg | 0.621 | 0.300 | 0.404 | 0.061 | 0.700 |
| ml_knn5 | 0.550 | 0.346 | 0.425 | 0.094 | 0.654 |

Latency: 20/20 fault scenarios detected, mean 3.5 s, median 0 s (first window).
Advisory agreement ML vs rule: 0.846; 148 ML-only alarms held advisory-only.

## Independent verification

Leakage audit passes (120/40/40 disjoint families); metric spot checks pass.
No production C++ imports ML; no ML output reaches actuators (no such path).

## Key methodological events (kept as findings)

1. The 50 m bias gave all-F1 = 1.000 (10σ/axis, trivially separable) —
   discarded as uninformative; redesigned at 10 m for overlapping classes.
2. sklearn install via pip stalled; conda-forge channel provided 1.9.0.
   Declared in `pyproject.toml [project.optional-dependencies] ml`.
3. Ruff flags C++ files (pre-existing scope: Python-only lint) — verified
   Python-only scope clean.

## Regression

New demo + Python only; no C++ library changes. Full suite 301/301. Ruff
(Python scope) clean.

## Limitations

Uncalibrated votes; single fault kind (bias jump); single sensor (GNSS NIS);
kNN stores training data (memory vs logreg); advisory-only (no closed-loop
ML); no physics+ML residual experiment (documented future work, needs a
residual with learnable structure — drag or vision, not NIS windows).
