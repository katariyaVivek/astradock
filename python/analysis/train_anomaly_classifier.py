#!/usr/bin/env python3
"""M23 sensor-anomaly classification: classical baselines first, ML challenger second.

Protocol (mirrors the lesson):
  1. Leakage audit: scenario_seed families disjoint across train/val/test.
  2. Classical baselines: chi-square gate count (M20 rule) + logistic regression
     on window features. NO model selection on test.
  3. ML challenger: RBF-SVM and kNN selected on VAL, evaluated once on TEST.
  4. Metrics: precision/recall/F1/false-alarm/miss/latency + calibration note.
  5. Safe integration: advisory gate — ML votes, M20 monitor acts (agreement
     table); ML never commands exclusion alone.
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC

FEATURES = ["nis_0", "nis_1", "nis_2", "nis_3", "nis_4"]
GATE = 14.4494  # M20 6-DOF 95% bound


def window_features(df: pd.DataFrame) -> pd.DataFrame:
    """Hand features: max, mean, count-over-gate, last value."""
    out = pd.DataFrame()
    out["max"] = df[FEATURES].max(axis=1)
    out["mean"] = df[FEATURES].mean(axis=1)
    out["over"] = (df[FEATURES] > GATE).sum(axis=1)
    out["last"] = df["nis_4"]
    return out


def leakage_audit(train: pd.DataFrame, val: pd.DataFrame, test: pd.DataFrame) -> None:
    """Seed families must be disjoint across splits (raises on leak)."""
    tr, va, te = set(train["scenario_seed"]), set(val["scenario_seed"]), set(test["scenario_seed"])
    assert not (tr & va), f"train/val seed overlap: {tr & va}"
    assert not (tr & te), f"train/test seed overlap: {tr & te}"
    assert not (va & te), f"val/test seed overlap: {va & te}"
    print(f"Leakage audit PASS: {len(tr)}/{len(va)}/{len(te)} disjoint seed families")


def metrics(name: str, y_true: np.ndarray, y_pred: np.ndarray) -> dict:
    """Precision/recall/F1/false-alarm/miss with zero-division guards."""
    tn, fp, fn, tp = confusion_matrix(y_true, y_pred, labels=[0, 1]).ravel()
    out = {
        "model": name,
        "precision": float(precision_score(y_true, y_pred, zero_division=0)),
        "recall": float(recall_score(y_true, y_pred, zero_division=0)),
        "f1": float(f1_score(y_true, y_pred, zero_division=0)),
        "false_alarm_rate": float(fp / max(tn + fp, 1)),
        "miss_rate": float(fn / max(fn + tp, 1)),
        "fp": int(fp),
        "fn": int(fn),
    }
    print(f"{name:>22}: P={out['precision']:.3f} R={out['recall']:.3f} "
          f"F1={out['f1']:.3f} FAR={out['false_alarm_rate']:.4f} MISS={out['miss_rate']:.4f}")
    return out


def main():
    root = Path(__file__).resolve().parents[2]
    data = root / "data"
    train = pd.read_csv(data / "m23_train.csv")
    val = pd.read_csv(data / "m23_val.csv")
    test = pd.read_csv(data / "m23_test.csv")
    leakage_audit(train, val, test)

    X_train, y_train = window_features(train), train["label"].to_numpy()
    X_val, y_val = window_features(val), val["label"].to_numpy()
    X_test, y_test = window_features(test), test["label"].to_numpy()

    # Baseline 0: M20 rule — 3-of-5 over gate (no training at all).
    rule = lambda df: ((df[FEATURES] > GATE).sum(axis=1) >= 3).astype(int).to_numpy()  # noqa: E731
    results = [metrics("m20_rule_3of5", y_test, rule(test))]

    # Baseline 1: logistic regression (trained on train, untouched test).
    logreg = LogisticRegression(max_iter=2000)
    logreg.fit(X_train, y_train)
    results.append(metrics("logreg", y_test, logreg.predict(X_test)))

    # Challenger selection on VAL only.
    scaler = StandardScaler().fit(X_train)
    candidates = {
        "svm_rbf": SVC(kernel="rbf", random_state=0),
        "knn5": KNeighborsClassifier(n_neighbors=5),
    }
    best_name, best_f1, best_model = "", -1.0, None
    for name, model in candidates.items():
        model.fit(scaler.transform(X_train), y_train)
        f1 = f1_score(y_val, model.predict(scaler.transform(X_val)))
        print(f"val {name}: F1={f1:.3f}")
        if f1 > best_f1:
            best_name, best_f1, best_model = name, f1, model
    print(f"Selected challenger on VAL: {best_name}")
    y_ml = best_model.predict(scaler.transform(X_test))
    results.append(metrics(f"ml_{best_name}", y_test, y_ml))

    # Safe integration: advisory agreement table (ML votes, M20 acts).
    y_rule = rule(test)
    agree = (y_ml == y_rule).mean()
    ml_only = ((y_ml == 1) & (y_rule == 0)).sum()
    print(f"Advisory agreement ML vs M20 rule: {agree:.3f}; ML-only alarms: {ml_only} "
          f"(would NOT trigger exclusion — monitor acts, ML advises)")

    # Latency: first fault-window index flagged per faulty test scenario.
    faulty = test[test["label"] == 1].copy()
    faulty["pred"] = y_ml[test["label"].to_numpy() == 1]
    latencies = faulty.groupby("scenario_seed")["pred"].apply(
        lambda s: int(np.argmax(s.to_numpy())) if s.to_numpy().any() else -1)
    detected = latencies[latencies >= 0]
    print(f"Fault scenarios detected: {len(detected)}/{latencies.size}; "
          f"mean latency {detected.mean() * 5:.1f} s in 5 s windows "
          f"(median {detected.median() * 5:.1f} s)")

    summary = pd.DataFrame(results)
    summary.to_csv(data / "m23_results.csv", index=False)
    print("Wrote data/m23_results.csv")


if __name__ == "__main__":
    main()
