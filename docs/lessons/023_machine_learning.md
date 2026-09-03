# Lesson 023 — Machine Learning for Space Systems

## 1. Physical Problem & Motivation

ML enters at M23 — after EKFs, PD/LQR, geometric PnP, and classical FDIR all
stand as validated baselines. The application is deliberately narrow and
aerospace-real: **GNSS bias-jump classification from EKF innovation windows**.
The question is not "can ML classify?" but "does ML add anything over the
chi-square monitor, and can it be integrated without endangering the loop?"
The answer, measured: at 10 m bias (2σ/axis, overlapping classes), kNN reaches
F1 = 0.425 vs the M20 rule's 0.025 and logreg's 0.404 — a real but modest gain
bought with 9.4% false alarms. That tradeoff, not accuracy, is the result.

## 2. Leakage Prevention (The Load-Bearing Wall)

Time-series rows from one trajectory share dynamics, noise seeds, and geometry:
random row splits leak trajectory identity into every split and report fantasy
scores. The dataset splits by SCENARIO (seed-disjoint families: train
1000/2000s, val 3000/4000s, test 5000/6000s), enforced by an audit that raises
on overlap. Windows touching the fault onset are excluded (ambiguous labels).
Features are window statistics (max/mean/over-gate/last) — no raw time index,
no seed, no scenario ID.

## 3. Classical Baselines First

Baseline 0 (no training): the M20 3-of-5 gate — precision 1.000, recall 0.013
(it never false-alarms and almost never fires at 10 m: honest and useless).
Baseline 1: logistic regression (F1 = 0.404). Only then the challenger:
RBF-SVM vs kNN selected on VAL (kNN wins, F1 = 0.372 val), evaluated ONCE on
test (F1 = 0.425). No deep learning: a 5-element innovation window does not
justify it, and the lesson says so explicitly.

## 4. Metrics and Safe Integration

Precision/recall/F1/false-alarm/miss/latency — never accuracy alone (75%
nominal base rate makes accuracy a liar). Latency: all 20 fault scenarios
detected, mean 3.5 s. Calibration is reported as a limitation (uncalibrated
votes, no probability claims). Integration is ADVISORY: the ML vote is
displayed beside the M20 monitor; agreement 0.846; 148 ML-only alarms would
NOT trigger exclusion — the validated monitor acts, ML advises. An
unconstrained model never touches actuators (architectural rule, audited by
the absence of any C++ ML-in-the-loop path: there is none).

## 5. Verification

Dataset generator: deterministic C++ (seeded EKF streams, 2880/960/960
windows). Oracle: leakage audit + metric spot checks. Demo: full protocol run
with results CSV. Figure: precision/recall + error-rate comparison.

## What you should now understand

1. Why do random row splits lie on trajectory data?
2. What does the M20 rule's P=1.000/R=0.013 tell you about thresholds?
3. Why was the 50 m bias experiment discarded (all-F1 = 1.000)?
4. What does kNN buy over logreg here, and what does it cost?
5. Why advisory-only, and what would change that?
6. Why is accuracy the wrong metric for this dataset?
