#!/usr/bin/env python3
"""
QW-1867: Beta-augmented observables benchmark.

Adaptive follow-up after QW-1866 negative result:
compare baseline vs beta-augmented feature sets (with/without paired branch).
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Sequence

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1867_beta_augmented_observables_benchmark.json"
OUT_MD = ROOT / "RAPORT_QW1867_BETA_AUGMENTED_OBSERVABLES_BENCHMARK.md"

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


def predict(theta: np.ndarray, features: Sequence[str]) -> np.ndarray:
    tt = theta.reshape(1, 3)
    return np.array([feature_batch(tt, f)[0] for f in features], dtype=float)


def nuisance_vector(features: Sequence[str]) -> np.ndarray:
    base = {
        "phase_increment": 0.8,
        "envelope_decay": -0.7,
        "zero_cross_offset": 0.5,
        "signed_asymmetry": -0.6,
        "torsion_cross_term": 0.9,
        "harmonic_ratio": -0.4,
        "phase_curvature": 0.7,
        "envelope_kurtosis_proxy": -0.5,
    }
    return np.array([base[f] for f in features], dtype=float)


def objective_baseline(cand: np.ndarray, y: np.ndarray, features: Sequence[str]) -> np.ndarray:
    obj = np.zeros(cand.shape[0], dtype=float)
    for i, f in enumerate(features):
        pred = feature_batch(cand, f)
        sigma = FEATURE_NOISE[f]
        r = (pred - y[i]) / sigma
        obj += r * r
    return obj


def objective_paired(cand: np.ndarray, y_base: np.ndarray, y_flip: np.ndarray, features: Sequence[str]) -> np.ndarray:
    obj = np.zeros(cand.shape[0], dtype=float)

    for i, f in enumerate(features):
        pred_b = feature_batch(cand, f)

        cand_flip = cand.copy()
        cand_flip[:, 1] = -cand_flip[:, 1]
        cand_flip = clamp_batch(cand_flip)
        pred_f = feature_batch(cand_flip, f)

        sigma = FEATURE_NOISE[f]
        r1 = (pred_b - y_base[i]) / sigma
        r2 = ((pred_b - pred_f) - (y_base[i] - y_flip[i])) / (math.sqrt(2.0) * sigma)

        obj += r1 * r1 + r2 * r2

    return obj


def estimate_theta(y_base: np.ndarray, y_flip: np.ndarray, features: Sequence[str], paired: bool, rng: np.random.Generator) -> np.ndarray:
    n0 = 700
    cand = np.column_stack(
        [
            rng.uniform(0.08, 1.60, size=n0),
            rng.uniform(-math.pi, math.pi, size=n0),
            rng.uniform(1e-4, 0.35, size=n0),
        ]
    )

    if paired:
        obj = objective_paired(cand, y_base, y_flip, features)
    else:
        obj = objective_baseline(cand, y_base, features)

    best = cand[int(np.argmin(obj))]

    for so, sp, sb, nloc in [(0.10, 0.30, 0.020, 320), (0.04, 0.12, 0.008, 320), (0.015, 0.05, 0.003, 320)]:
        loc = np.tile(best, (nloc, 1))
        loc[:, 0] += rng.normal(0.0, so, size=nloc)
        loc[:, 1] += rng.normal(0.0, sp, size=nloc)
        loc[:, 2] += rng.normal(0.0, sb, size=nloc)
        loc = clamp_batch(loc)

        if paired:
            obj_loc = objective_paired(loc, y_base, y_flip, features)
        else:
            obj_loc = objective_baseline(loc, y_base, features)

        best = loc[int(np.argmin(obj_loc))]

    return clamp_one(best)


def metrics_from_errors(eo: np.ndarray, ep: np.ndarray, eb: np.ndarray) -> Dict:
    rmse = {
        "omega": float(np.sqrt(np.mean(eo * eo))),
        "phi": float(np.sqrt(np.mean(ep * ep))),
        "beta": float(np.sqrt(np.mean(eb * eb))),
    }
    tol = {
        "omega_0p08": float(np.mean(np.abs(eo) <= 0.08)),
        "phi_0p20": float(np.mean(np.abs(ep) <= 0.20)),
        "beta_0p015": float(np.mean(np.abs(eb) <= 0.015)),
    }
    if np.std(eo) > 1e-12 and np.std(eb) > 1e-12:
        corr_ob = float(np.corrcoef(eo, eb)[0, 1])
    else:
        corr_ob = 0.0
    return {"rmse": rmse, "tol": tol, "corr_omega_beta": corr_ob}


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

    base_features = list(d1863.get("best_baseline", {}).get("features", ["phase_increment", "envelope_decay", "zero_cross_offset"]))
    aug_features = list(base_features)
    if "torsion_cross_term" not in aug_features:
        aug_features.append("torsion_cross_term")

    arms = [
        {"id": "A_BASELINE", "features": base_features, "paired": False},
        {"id": "B_AUGMENTED", "features": aug_features, "paired": False},
        {"id": "C_AUGMENTED_PAIRED", "features": aug_features, "paired": True},
    ]

    rng = np.random.default_rng(186700)
    n_rep = 170

    arm_err = {
        arm["id"]: {"omega": [], "phi": [], "beta": []}
        for arm in arms
    }

    rows_head = []

    for i in range(n_rep):
        true = theta_target.copy()
        true[0] += rng.normal(0.0, 0.08)
        true[1] += rng.normal(0.0, 0.30)
        true[2] += rng.normal(0.0, 0.020)
        true = clamp_one(true)
        true_flip = clamp_one(np.array([true[0], -true[1], true[2]], dtype=float))

        row_dbg = {"true": {"omega": float(true[0]), "phi": float(true[1]), "beta": float(true[2])}}

        for arm in arms:
            aid = arm["id"]
            fset = arm["features"]
            paired = arm["paired"]

            y_clean = predict(true, fset)
            y_flip_clean = predict(true_flip, fset)

            drift_strength = 0.10 * (0.7 * true[0] * true[2] + 0.3 * true[2] * np.cos(true[1]))
            nuisance = drift_strength * nuisance_vector(fset)

            noise_b = np.array([rng.normal(0.0, FEATURE_NOISE[f]) for f in fset], dtype=float)
            noise_f = np.array([rng.normal(0.0, FEATURE_NOISE[f]) for f in fset], dtype=float)

            y_base = y_clean + nuisance + noise_b
            y_flip = y_flip_clean + nuisance + noise_f

            est = estimate_theta(y_base, y_flip, fset, paired, rng)

            eo = float(est[0] - true[0])
            ep = float((est[1] - true[1] + math.pi) % (2.0 * math.pi) - math.pi)
            eb = float(est[2] - true[2])

            arm_err[aid]["omega"].append(eo)
            arm_err[aid]["phi"].append(ep)
            arm_err[aid]["beta"].append(eb)

            if i < 80:
                row_dbg[aid] = {"omega": float(est[0]), "phi": float(est[1]), "beta": float(est[2])}

        if i < 80:
            rows_head.append(row_dbg)

    arm_metrics = {}
    for arm in arms:
        aid = arm["id"]
        eo = np.array(arm_err[aid]["omega"], dtype=float)
        ep = np.array(arm_err[aid]["phi"], dtype=float)
        eb = np.array(arm_err[aid]["beta"], dtype=float)
        arm_metrics[aid] = {
            "features": arm["features"],
            "paired": arm["paired"],
            **metrics_from_errors(eo, ep, eb),
        }

    base = arm_metrics["A_BASELINE"]

    improvements = {}
    for aid, m in arm_metrics.items():
        improvements[aid] = {
            "beta_rmse_factor_vs_baseline": base["rmse"]["beta"] / max(m["rmse"]["beta"], 1e-12),
            "beta_tol_gain_vs_baseline": m["tol"]["beta_0p015"] - base["tol"]["beta_0p015"],
            "abs_corr_reduction_vs_baseline": abs(base["corr_omega_beta"]) - abs(m["corr_omega_beta"]),
        }

    ranked = sorted(
        [
            {
                "arm": aid,
                "beta_rmse_factor": imp["beta_rmse_factor_vs_baseline"],
                "beta_tol_gain": imp["beta_tol_gain_vs_baseline"],
                "abs_corr_reduction": imp["abs_corr_reduction_vs_baseline"],
                "score": 0.55 * imp["beta_rmse_factor_vs_baseline"]
                + 0.30 * (1.0 + imp["beta_tol_gain_vs_baseline"])
                + 0.15 * (1.0 + imp["abs_corr_reduction_vs_baseline"]),
            }
            for aid, imp in improvements.items()
        ],
        key=lambda x: x["score"],
        reverse=True,
    )

    best = ranked[0]

    if (
        best["arm"] != "A_BASELINE"
        and best["beta_rmse_factor"] >= 1.15
        and best["beta_tol_gain"] >= 0.06
    ):
        verdict = "BETA_AUGMENTATION_SUPPORTED"
    elif best["arm"] != "A_BASELINE" and best["beta_rmse_factor"] >= 1.05:
        verdict = "BETA_AUGMENTATION_PARTIAL"
    else:
        verdict = "BETA_AUGMENTATION_NOT_SUPPORTED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_rep": n_rep,
        "target_center": {"omega": float(theta_target[0]), "phi": float(theta_target[1]), "beta": float(theta_target[2])},
        "arm_metrics": arm_metrics,
        "improvements_vs_baseline": improvements,
        "ranking": ranked,
        "best_arm": best,
        "rows_head": rows_head,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1867: BETA-AUGMENTED OBSERVABLES BENCHMARK",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- n_rep: {n_rep}",
        "",
        "## Arms",
    ]

    for aid, m in arm_metrics.items():
        lines.append(
            f"- {aid}: paired={m['paired']} | beta RMSE={m['rmse']['beta']:.4f} | beta tol={m['tol']['beta_0p015']:.3f} | corr_ob={m['corr_omega_beta']:.3f}"
        )

    lines += ["", "## Ranking (vs baseline)"]
    for r in ranked:
        lines.append(
            f"- {r['arm']}: score={r['score']:.3f}, beta_rmse_factor={r['beta_rmse_factor']:.3f}, beta_tol_gain={r['beta_tol_gain']:.3f}, |corr| reduction={r['abs_corr_reduction']:.3f}"
        )

    lines += [
        "",
        "## Best Arm",
        f"- {best['arm']}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1867] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1867] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
