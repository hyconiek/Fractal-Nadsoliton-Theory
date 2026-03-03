#!/usr/bin/env python3
"""
QW-1863: Optimal observable design for ansatz-free identifiability.

Builds a synthetic information-geometry benchmark and selects feature sets that
best separate (beta, omega, phi) in the signed dynamic micromodel regime.
"""

from __future__ import annotations

import itertools
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Sequence, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1863_identifiability_optimal_observables_design.json"
OUT_MD = ROOT / "RAPORT_QW1863_IDENTIFIABILITY_OPTIMAL_OBSERVABLES_DESIGN.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def target_tuple() -> Dict[str, float]:
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    t = d1862.get("target_tuple", {})
    return {
        "omega": float(t.get("omega", math.pi / 4.0)),
        "phi": float(t.get("phi", math.pi / 6.0)),
        "beta": float(t.get("beta", 0.01)),
    }


def clamp_theta(theta: np.ndarray) -> np.ndarray:
    o = float(np.clip(theta[0], 0.08, 1.60))
    p = float(((theta[1] + math.pi) % (2.0 * math.pi)) - math.pi)
    b = float(np.clip(theta[2], 1e-4, 0.35))
    return np.array([o, p, b], dtype=float)


def feature_map(theta: np.ndarray) -> Dict[str, float]:
    omega, phi, beta = float(theta[0]), float(theta[1]), float(theta[2])

    return {
        "phase_increment": omega + 0.05 * beta,
        "envelope_decay": beta + 0.02 * omega * omega,
        "zero_cross_offset": phi - 0.30 * omega * beta,
        "signed_asymmetry": beta * math.sin(phi),
        "torsion_cross_term": beta * omega * math.cos(phi),
        "harmonic_ratio": math.exp(-beta) * math.cos(phi) / (1.0 + omega * omega),
        "phase_curvature": omega * omega * math.sin(phi) + 0.10 * beta,
        "envelope_kurtosis_proxy": beta / (1.0 + abs(math.sin(phi))) + 0.03 * omega,
    }


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


def jacobian(theta: np.ndarray, features: Sequence[str]) -> np.ndarray:
    eps = np.array([1e-3, 1e-3, 1e-4], dtype=float)
    j = np.zeros((len(features), 3), dtype=float)

    for k in range(3):
        t_plus = theta.copy()
        t_minus = theta.copy()
        t_plus[k] += eps[k]
        t_minus[k] -= eps[k]
        t_plus = clamp_theta(t_plus)
        t_minus = clamp_theta(t_minus)

        f_plus = feature_map(t_plus)
        f_minus = feature_map(t_minus)

        for i, f_name in enumerate(features):
            deriv = (f_plus[f_name] - f_minus[f_name]) / (2.0 * eps[k])
            j[i, k] = deriv / FEATURE_NOISE[f_name]

    return j


def row_beta_omega_coupling(j: np.ndarray) -> float:
    v_o = j[:, 0]
    v_b = j[:, 2]
    n_o = float(np.linalg.norm(v_o))
    n_b = float(np.linalg.norm(v_b))
    if n_o <= 1e-12 or n_b <= 1e-12:
        return 1.0
    c = float(np.dot(v_o, v_b) / (n_o * n_b))
    return abs(max(-1.0, min(1.0, c)))


def evaluate_combo(features: Sequence[str], theta0: np.ndarray, n_samples: int, intervention: bool, seed: int) -> Dict:
    rng = np.random.default_rng(seed)
    conds = []
    logdets = []
    mineigs = []
    couplings = []

    for _ in range(n_samples):
        d_omega = rng.normal(0.0, 0.12)
        d_phi = rng.normal(0.0, 0.35)
        d_beta = rng.normal(0.0, 0.02)
        t = clamp_theta(theta0 + np.array([d_omega, d_phi, d_beta], dtype=float))

        j = jacobian(t, features)
        if intervention:
            # Paired signed intervention: compare baseline to sign-flipped phase branch.
            t_flip = clamp_theta(np.array([t[0], -t[1], t[2]], dtype=float))
            j_flip = jacobian(t_flip, features)
            j = np.vstack([j, j - j_flip])

        fim = j.T @ j + 1e-10 * np.eye(3)
        eig = np.linalg.eigvalsh(fim)

        e_min = float(max(eig[0], 1e-14))
        e_max = float(max(eig[-1], 1e-14))
        cond = e_max / e_min
        logdet = float(np.sum(np.log(np.clip(eig, 1e-14, None))))

        conds.append(cond)
        logdets.append(logdet)
        mineigs.append(e_min)
        couplings.append(row_beta_omega_coupling(j))

    cond_med = float(np.median(conds))
    cond_q90 = float(np.quantile(conds, 0.90))
    logdet_med = float(np.median(logdets))
    mineig_q10 = float(np.quantile(mineigs, 0.10))
    coupling_med = float(np.median(couplings))

    # Higher is better.
    cond_term = 1.0 / (1.0 + math.log10(max(cond_q90, 1.0)))
    det_term = max(0.0, min(1.0, (logdet_med + 12.0) / 18.0))
    coupling_term = 1.0 - coupling_med
    score = cond_term * det_term * max(0.0, coupling_term)

    return {
        "features": list(features),
        "n_features": len(features),
        "intervention": intervention,
        "metrics": {
            "cond_median": cond_med,
            "cond_q90": cond_q90,
            "logdet_median": logdet_med,
            "mineig_q10": mineig_q10,
            "beta_omega_coupling_median": coupling_med,
        },
        "score": float(score),
    }


def main() -> None:
    t = target_tuple()
    theta0 = np.array([t["omega"], t["phi"], t["beta"]], dtype=float)

    feature_names = list(FEATURE_NOISE.keys())

    baseline_rows = []
    seed = 186300

    for k in [3, 4, 5]:
        for comb in itertools.combinations(feature_names, k):
            row = evaluate_combo(comb, theta0, n_samples=320, intervention=False, seed=seed)
            baseline_rows.append(row)
            seed += 1

    baseline_rows_sorted = sorted(baseline_rows, key=lambda x: x["score"], reverse=True)
    best_base = baseline_rows_sorted[0]

    best_features = tuple(best_base["features"])
    inter_row = evaluate_combo(best_features, theta0, n_samples=400, intervention=True, seed=999001)

    gain_score = inter_row["score"] / max(best_base["score"], 1e-12)
    gain_cond_q90 = best_base["metrics"]["cond_q90"] / max(inter_row["metrics"]["cond_q90"], 1e-12)

    if best_base["score"] >= 0.12 and inter_row["score"] >= 0.20 and gain_cond_q90 >= 1.5:
        verdict = "IDENTIFIABILITY_DESIGN_READY_FOR_SYNTHETIC_VALIDATION"
    elif best_base["score"] >= 0.08:
        verdict = "IDENTIFIABILITY_DESIGN_PARTIAL"
    else:
        verdict = "IDENTIFIABILITY_DESIGN_WEAK"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_tuple": t,
        "baseline_top_rows": baseline_rows_sorted[:40],
        "best_baseline": best_base,
        "best_with_signed_intervention": inter_row,
        "intervention_gain": {
            "score_ratio": gain_score,
            "cond_q90_reduction_factor": gain_cond_q90,
        },
        "recommended_next_studies": [
            {
                "id": "QW-1865",
                "title": "Synthetic recovery benchmark on best observable set",
                "goal": "Verify unbiased recovery of beta, omega, phi under realistic noise.",
            },
            {
                "id": "QW-1866",
                "title": "Paired signed-intervention execution protocol",
                "goal": "Implement baseline vs sign-flipped branch measurements to break omega-beta coupling.",
            },
        ],
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1863: IDENTIFIABILITY OPTIMAL OBSERVABLES DESIGN",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Target Tuple",
        f"- omega={t['omega']:.6f}",
        f"- phi={t['phi']:.6f}",
        f"- beta={t['beta']:.6f}",
        "",
        "## Best Baseline Feature Set",
        f"- features: {', '.join(best_base['features'])}",
        f"- score: {best_base['score']:.4f}",
        f"- cond_q90: {best_base['metrics']['cond_q90']:.3e}",
        f"- beta-omega coupling (median): {best_base['metrics']['beta_omega_coupling_median']:.3f}",
        "",
        "## Signed Intervention Result",
        f"- score: {inter_row['score']:.4f}",
        f"- cond_q90: {inter_row['metrics']['cond_q90']:.3e}",
        f"- beta-omega coupling (median): {inter_row['metrics']['beta_omega_coupling_median']:.3f}",
        "",
        "## Gains",
        f"- score ratio (intervention/base): {gain_score:.3f}",
        f"- cond_q90 reduction factor: {gain_cond_q90:.3f}",
        "",
        "## Top Baseline Combos",
    ]

    for row in baseline_rows_sorted[:12]:
        lines.append(
            f"- k={row['n_features']}, score={row['score']:.4f}, cond_q90={row['metrics']['cond_q90']:.3e}, "
            f"coupling={row['metrics']['beta_omega_coupling_median']:.3f} | {', '.join(row['features'])}"
        )

    lines += ["", "## Recommended Next Studies"]
    for r in out["recommended_next_studies"]:
        lines.append(f"- {r['id']}: {r['title']} -> {r['goal']}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1863] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1863] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
