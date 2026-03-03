#!/usr/bin/env python3
"""
QW-1865: Synthetic recovery benchmark on best observable set.

Uses the QW-1863 selected feature subset and performs Monte Carlo parameter
recovery for (omega, phi, beta) under realistic noise.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1865_synthetic_recovery_benchmark.json"
OUT_MD = ROOT / "RAPORT_QW1865_SYNTHETIC_RECOVERY_BENCHMARK.md"

FEATURE_NOISE = {
    "phase_increment": 0.020,
    "envelope_decay": 0.015,
    "zero_cross_offset": 0.035,
    "signed_asymmetry": 0.012,
    "torsion_cross_term": 0.010,
    "harmonic_ratio": 0.018,
    "phase_curvature": 0.020,
    "envelope_kurtosis_proxy": 0.012,
}


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def wrap_phi(phi: np.ndarray) -> np.ndarray:
    return (phi + math.pi) % (2.0 * math.pi) - math.pi


def clamp_batch(theta: np.ndarray) -> np.ndarray:
    t = theta.copy()
    t[:, 0] = np.clip(t[:, 0], 0.08, 1.60)
    t[:, 1] = wrap_phi(t[:, 1])
    t[:, 2] = np.clip(t[:, 2], 1e-4, 0.35)
    return t


def clamp_one(theta: np.ndarray) -> np.ndarray:
    return clamp_batch(theta.reshape(1, 3))[0]


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def feature_batch(theta: np.ndarray, feature: str) -> np.ndarray:
    o = theta[:, 0]
    p = theta[:, 1]
    b = theta[:, 2]

    if feature == "phase_increment":
        return o + 0.05 * b
    if feature == "envelope_decay":
        return b + 0.02 * o * o
    if feature == "zero_cross_offset":
        return p - 0.30 * o * b
    if feature == "signed_asymmetry":
        return b * np.sin(p)
    if feature == "torsion_cross_term":
        return b * o * np.cos(p)
    if feature == "harmonic_ratio":
        return np.exp(-b) * np.cos(p) / (1.0 + o * o)
    if feature == "phase_curvature":
        return o * o * np.sin(p) + 0.10 * b
    if feature == "envelope_kurtosis_proxy":
        return b / (1.0 + np.abs(np.sin(p))) + 0.03 * o

    raise KeyError(feature)


def predict_one(theta: np.ndarray, features: Sequence[str]) -> np.ndarray:
    tt = theta.reshape(1, 3)
    return np.array([feature_batch(tt, f)[0] for f in features], dtype=float)


def objective_batch(cand: np.ndarray, y: np.ndarray, features: Sequence[str]) -> np.ndarray:
    obj = np.zeros(cand.shape[0], dtype=float)
    for i, f in enumerate(features):
        pred = feature_batch(cand, f)
        sigma = FEATURE_NOISE[f]
        r = (pred - y[i]) / sigma
        obj += r * r
    return obj


def estimate_theta(y: np.ndarray, features: Sequence[str], rng: np.random.Generator) -> Tuple[np.ndarray, float]:
    n0 = 900
    cand = np.column_stack(
        [
            rng.uniform(0.08, 1.60, size=n0),
            rng.uniform(-math.pi, math.pi, size=n0),
            rng.uniform(1e-4, 0.35, size=n0),
        ]
    )
    obj = objective_batch(cand, y, features)
    best_idx = int(np.argmin(obj))
    best = cand[best_idx]
    best_obj = float(obj[best_idx])

    for so, sp, sb, nloc in [(0.10, 0.30, 0.020, 420), (0.04, 0.12, 0.008, 420), (0.015, 0.05, 0.003, 420)]:
        loc = np.tile(best, (nloc, 1))
        loc[:, 0] += rng.normal(0.0, so, size=nloc)
        loc[:, 1] += rng.normal(0.0, sp, size=nloc)
        loc[:, 2] += rng.normal(0.0, sb, size=nloc)
        loc = clamp_batch(loc)

        obj_loc = objective_batch(loc, y, features)
        idx = int(np.argmin(obj_loc))
        if float(obj_loc[idx]) < best_obj:
            best = loc[idx]
            best_obj = float(obj_loc[idx])

    return clamp_one(best), best_obj


def main() -> None:
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1863 = read_json("report_qw1863_identifiability_optimal_observables_design.json")

    target = d1862.get("target_tuple", {})
    theta_target = np.array(
        [
            float(target.get("omega", math.pi / 4.0)),
            float(target.get("phi", math.pi / 6.0)),
            float(target.get("beta", 0.01)),
        ],
        dtype=float,
    )

    features = list(d1863.get("best_baseline", {}).get("features", ["phase_increment", "envelope_decay", "zero_cross_offset"]))

    rng = np.random.default_rng(186500)
    n_rep = 240

    rows = []
    for _ in range(n_rep):
        true = theta_target.copy()
        true[0] += rng.normal(0.0, 0.08)
        true[1] += rng.normal(0.0, 0.25)
        true[2] += rng.normal(0.0, 0.015)
        true = clamp_one(true)

        y_clean = predict_one(true, features)
        noise = np.array([rng.normal(0.0, FEATURE_NOISE[f]) for f in features], dtype=float)
        y_obs = y_clean + noise

        est, obj = estimate_theta(y_obs, features, rng)

        err_o = float(est[0] - true[0])
        err_p = float((est[1] - true[1] + math.pi) % (2.0 * math.pi) - math.pi)
        err_b = float(est[2] - true[2])

        rows.append(
            {
                "true": {"omega": float(true[0]), "phi": float(true[1]), "beta": float(true[2])},
                "est": {"omega": float(est[0]), "phi": float(est[1]), "beta": float(est[2])},
                "errors": {"omega": err_o, "phi": err_p, "beta": err_b},
                "objective": obj,
                "nonboundary": bool(0.10 < est[0] < 1.50 and 0.001 < est[2] < 0.30),
            }
        )

    eo = np.array([r["errors"]["omega"] for r in rows], dtype=float)
    ep = np.array([r["errors"]["phi"] for r in rows], dtype=float)
    eb = np.array([r["errors"]["beta"] for r in rows], dtype=float)

    rmse_o = float(np.sqrt(np.mean(eo * eo)))
    rmse_p = float(np.sqrt(np.mean(ep * ep)))
    rmse_b = float(np.sqrt(np.mean(eb * eb)))

    med_ae_o = float(np.median(np.abs(eo)))
    med_ae_p = float(np.median(np.abs(ep)))
    med_ae_b = float(np.median(np.abs(eb)))

    tol_o = float(np.mean(np.abs(eo) <= 0.08))
    tol_p = float(np.mean(np.abs(ep) <= 0.20))
    tol_b = float(np.mean(np.abs(eb) <= 0.015))

    nonboundary_rate = float(np.mean([1.0 if r["nonboundary"] else 0.0 for r in rows]))

    if np.std(eo) > 1e-12 and np.std(eb) > 1e-12:
        corr_ob = float(np.corrcoef(eo, eb)[0, 1])
    else:
        corr_ob = 0.0

    if tol_o >= 0.80 and tol_p >= 0.75 and tol_b >= 0.80 and nonboundary_rate >= 0.85 and abs(corr_ob) <= 0.25:
        verdict = "SYNTHETIC_RECOVERY_PASS"
    elif tol_o >= 0.60 and tol_p >= 0.60 and tol_b >= 0.60 and nonboundary_rate >= 0.70:
        verdict = "SYNTHETIC_RECOVERY_PARTIAL"
    else:
        verdict = "SYNTHETIC_RECOVERY_WEAK"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "features": features,
        "n_rep": n_rep,
        "target_center": {"omega": float(theta_target[0]), "phi": float(theta_target[1]), "beta": float(theta_target[2])},
        "metrics": {
            "rmse": {"omega": rmse_o, "phi": rmse_p, "beta": rmse_b},
            "median_abs_error": {"omega": med_ae_o, "phi": med_ae_p, "beta": med_ae_b},
            "tolerance_hit_rate": {"omega_0p08": tol_o, "phi_0p20": tol_p, "beta_0p015": tol_b},
            "nonboundary_rate": nonboundary_rate,
            "error_corr_omega_beta": corr_ob,
        },
        "rows_head": rows[:120],
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1865: SYNTHETIC RECOVERY BENCHMARK",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- n_rep: {n_rep}",
        "",
        "## Feature Set",
        f"- {', '.join(features)}",
        "",
        "## Metrics",
        f"- RMSE omega/phi/beta: {rmse_o:.4f} / {rmse_p:.4f} / {rmse_b:.4f}",
        f"- median abs error omega/phi/beta: {med_ae_o:.4f} / {med_ae_p:.4f} / {med_ae_b:.4f}",
        f"- tolerance hit omega(|e|<=0.08): {tol_o:.3f}",
        f"- tolerance hit phi(|e|<=0.20): {tol_p:.3f}",
        f"- tolerance hit beta(|e|<=0.015): {tol_b:.3f}",
        f"- nonboundary rate: {nonboundary_rate:.3f}",
        f"- error corr(omega,beta): {corr_ob:.3f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1865] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1865] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
