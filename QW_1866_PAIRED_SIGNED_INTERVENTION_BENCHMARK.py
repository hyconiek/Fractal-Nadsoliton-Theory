#!/usr/bin/env python3
"""
QW-1866: Paired signed-intervention benchmark.

Compares baseline recovery vs paired (baseline + sign-flip difference) recovery
under hidden nuisance drift that couples omega and beta.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Sequence, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1866_paired_signed_intervention_benchmark.json"
OUT_MD = ROOT / "RAPORT_QW1866_PAIRED_SIGNED_INTERVENTION_BENCHMARK.md"

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


def estimate_theta_baseline(y: np.ndarray, features: Sequence[str], rng: np.random.Generator) -> np.ndarray:
    n0 = 850
    cand = np.column_stack(
        [
            rng.uniform(0.08, 1.60, size=n0),
            rng.uniform(-math.pi, math.pi, size=n0),
            rng.uniform(1e-4, 0.35, size=n0),
        ]
    )
    obj = objective_baseline(cand, y, features)
    best = cand[int(np.argmin(obj))]

    for so, sp, sb, nloc in [(0.10, 0.30, 0.020, 380), (0.04, 0.12, 0.008, 380), (0.015, 0.05, 0.003, 380)]:
        loc = np.tile(best, (nloc, 1))
        loc[:, 0] += rng.normal(0.0, so, size=nloc)
        loc[:, 1] += rng.normal(0.0, sp, size=nloc)
        loc[:, 2] += rng.normal(0.0, sb, size=nloc)
        loc = clamp_batch(loc)
        obj_loc = objective_baseline(loc, y, features)
        best = loc[int(np.argmin(obj_loc))]

    return clamp_one(best)


def estimate_theta_paired(y_base: np.ndarray, y_flip: np.ndarray, features: Sequence[str], rng: np.random.Generator) -> np.ndarray:
    n0 = 850
    cand = np.column_stack(
        [
            rng.uniform(0.08, 1.60, size=n0),
            rng.uniform(-math.pi, math.pi, size=n0),
            rng.uniform(1e-4, 0.35, size=n0),
        ]
    )
    obj = objective_paired(cand, y_base, y_flip, features)
    best = cand[int(np.argmin(obj))]

    for so, sp, sb, nloc in [(0.10, 0.30, 0.020, 380), (0.04, 0.12, 0.008, 380), (0.015, 0.05, 0.003, 380)]:
        loc = np.tile(best, (nloc, 1))
        loc[:, 0] += rng.normal(0.0, so, size=nloc)
        loc[:, 1] += rng.normal(0.0, sp, size=nloc)
        loc[:, 2] += rng.normal(0.0, sb, size=nloc)
        loc = clamp_batch(loc)
        obj_loc = objective_paired(loc, y_base, y_flip, features)
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

    features = list(d1863.get("best_baseline", {}).get("features", ["phase_increment", "envelope_decay", "zero_cross_offset"]))
    nu_vec = nuisance_vector(features)

    rng = np.random.default_rng(186600)
    n_rep = 180

    err_base_o = []
    err_base_p = []
    err_base_b = []
    err_pair_o = []
    err_pair_p = []
    err_pair_b = []

    rows_head = []

    for i in range(n_rep):
        true = theta_target.copy()
        true[0] += rng.normal(0.0, 0.08)
        true[1] += rng.normal(0.0, 0.30)
        true[2] += rng.normal(0.0, 0.020)
        true = clamp_one(true)

        true_flip = clamp_one(np.array([true[0], -true[1], true[2]], dtype=float))

        y_clean = predict(true, features)
        y_flip_clean = predict(true_flip, features)

        # Hidden nuisance term confounding baseline branch.
        drift_strength = 0.10 * (0.7 * true[0] * true[2] + 0.3 * true[2] * np.cos(true[1]))
        nuisance = drift_strength * nu_vec

        noise_b = np.array([rng.normal(0.0, FEATURE_NOISE[f]) for f in features], dtype=float)
        noise_f = np.array([rng.normal(0.0, FEATURE_NOISE[f]) for f in features], dtype=float)

        y_base = y_clean + nuisance + noise_b
        y_flip = y_flip_clean + nuisance + noise_f

        est_base = estimate_theta_baseline(y_base, features, rng)
        est_pair = estimate_theta_paired(y_base, y_flip, features, rng)

        ebo = float(est_base[0] - true[0])
        ebp = float((est_base[1] - true[1] + math.pi) % (2.0 * math.pi) - math.pi)
        ebb = float(est_base[2] - true[2])

        epo = float(est_pair[0] - true[0])
        epp = float((est_pair[1] - true[1] + math.pi) % (2.0 * math.pi) - math.pi)
        epb = float(est_pair[2] - true[2])

        err_base_o.append(ebo)
        err_base_p.append(ebp)
        err_base_b.append(ebb)

        err_pair_o.append(epo)
        err_pair_p.append(epp)
        err_pair_b.append(epb)

        if i < 120:
            rows_head.append(
                {
                    "true": {"omega": float(true[0]), "phi": float(true[1]), "beta": float(true[2])},
                    "est_base": {"omega": float(est_base[0]), "phi": float(est_base[1]), "beta": float(est_base[2])},
                    "est_pair": {"omega": float(est_pair[0]), "phi": float(est_pair[1]), "beta": float(est_pair[2])},
                }
            )

    eo_b = np.array(err_base_o, dtype=float)
    ep_b = np.array(err_base_p, dtype=float)
    eb_b = np.array(err_base_b, dtype=float)

    eo_p = np.array(err_pair_o, dtype=float)
    ep_p = np.array(err_pair_p, dtype=float)
    eb_p = np.array(err_pair_b, dtype=float)

    m_base = metrics_from_errors(eo_b, ep_b, eb_b)
    m_pair = metrics_from_errors(eo_p, ep_p, eb_p)

    improve = {
        "rmse_omega_factor": m_base["rmse"]["omega"] / max(m_pair["rmse"]["omega"], 1e-12),
        "rmse_phi_factor": m_base["rmse"]["phi"] / max(m_pair["rmse"]["phi"], 1e-12),
        "rmse_beta_factor": m_base["rmse"]["beta"] / max(m_pair["rmse"]["beta"], 1e-12),
        "tol_beta_gain": m_pair["tol"]["beta_0p015"] - m_base["tol"]["beta_0p015"],
        "abs_corr_omega_beta_reduction": abs(m_base["corr_omega_beta"]) - abs(m_pair["corr_omega_beta"]),
    }

    if (
        improve["rmse_beta_factor"] >= 1.20
        and improve["rmse_omega_factor"] >= 1.10
        and improve["tol_beta_gain"] >= 0.08
        and improve["abs_corr_omega_beta_reduction"] >= 0.10
    ):
        verdict = "PAIRED_SIGNED_INTERVENTION_STRONGLY_SUPPORTED"
    elif improve["rmse_beta_factor"] >= 1.05 and improve["tol_beta_gain"] >= 0.03:
        verdict = "PAIRED_SIGNED_INTERVENTION_PARTIAL_SUPPORT"
    else:
        verdict = "PAIRED_SIGNED_INTERVENTION_NOT_SUPPORTED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_rep": n_rep,
        "features": features,
        "target_center": {"omega": float(theta_target[0]), "phi": float(theta_target[1]), "beta": float(theta_target[2])},
        "baseline": m_base,
        "paired": m_pair,
        "improvement": improve,
        "rows_head": rows_head,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1866: PAIRED SIGNED-INTERVENTION BENCHMARK",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- n_rep: {n_rep}",
        "",
        "## Feature Set",
        f"- {', '.join(features)}",
        "",
        "## Baseline vs Paired RMSE",
        f"- omega: {m_base['rmse']['omega']:.4f} -> {m_pair['rmse']['omega']:.4f}",
        f"- phi: {m_base['rmse']['phi']:.4f} -> {m_pair['rmse']['phi']:.4f}",
        f"- beta: {m_base['rmse']['beta']:.4f} -> {m_pair['rmse']['beta']:.4f}",
        "",
        "## Baseline vs Paired Tolerance Hit",
        f"- omega(|e|<=0.08): {m_base['tol']['omega_0p08']:.3f} -> {m_pair['tol']['omega_0p08']:.3f}",
        f"- phi(|e|<=0.20): {m_base['tol']['phi_0p20']:.3f} -> {m_pair['tol']['phi_0p20']:.3f}",
        f"- beta(|e|<=0.015): {m_base['tol']['beta_0p015']:.3f} -> {m_pair['tol']['beta_0p015']:.3f}",
        "",
        "## Coupling Diagnostic",
        f"- corr(omega_error, beta_error): {m_base['corr_omega_beta']:.3f} -> {m_pair['corr_omega_beta']:.3f}",
        "",
        "## Improvement Factors",
        f"- rmse omega factor: {improve['rmse_omega_factor']:.3f}",
        f"- rmse phi factor: {improve['rmse_phi_factor']:.3f}",
        f"- rmse beta factor: {improve['rmse_beta_factor']:.3f}",
        f"- beta tolerance gain: {improve['tol_beta_gain']:.3f}",
        f"- |corr| reduction: {improve['abs_corr_omega_beta_reduction']:.3f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1866] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1866] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
