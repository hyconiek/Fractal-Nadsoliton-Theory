#!/usr/bin/env python3
"""
QW-1875: Canon-anchored constrained fit.

Fixes phi to canonical value and fits (omega, beta, node dynamics) to test whether
canonical compatibility can escape near-zero regime.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1875_canon_anchored_constrained_fit.json"
OUT_MD = ROOT / "RAPORT_QW1875_CANON_ANCHORED_CONSTRAINED_FIT.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def node_field(d: np.ndarray, nodes: List[int], sigma: float) -> np.ndarray:
    if not nodes:
        return np.zeros_like(d)
    terms = [np.exp(-0.5 * ((d - float(n)) / sigma) ** 2) for n in nodes]
    return np.clip(np.max(np.vstack(terms), axis=0), 0.0, 1.0)


def model_signal(d: np.ndarray, omega: float, phi: float, beta: float, amp: float, kappa: float, nf: np.ndarray) -> np.ndarray:
    env = 1.0 / (1.0 + beta * d)
    mod = np.clip(1.0 - kappa * nf, 0.12, 1.0)
    return amp * np.cos(omega * d + phi) * env * mod


def estimate_amp(y: np.ndarray, d: np.ndarray, omega: float, phi: float, beta: float, kappa: float, nf: np.ndarray) -> float:
    basis = np.cos(omega * d + phi) / (1.0 + beta * d) * np.clip(1.0 - kappa * nf, 0.12, 1.0)
    den = float(np.dot(basis, basis))
    if den <= 1e-12:
        return 0.0
    return float(np.dot(y, basis) / den)


def canon_score_omega_beta(omega: float, beta: float, t_omega: float, t_beta: float) -> float:
    z_o = abs(omega - t_omega) / 0.12
    z_b = abs(beta - t_beta) / 0.020
    z = min(20.0, math.sqrt((z_o * z_o + z_b * z_b) / 2.0))
    return float(math.exp(-0.5 * z * z))


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1871 = read_json("report_qw1871_primary_node_evidence_corpus.json")
    d1874 = read_json("report_qw1874_beta_omega_orthogonal_forcing.json")

    target = d1862.get("target_tuple", {})
    t_omega = float(target.get("omega", math.pi / 4.0))
    t_phi = float(target.get("phi", math.pi / 6.0))
    t_beta = float(target.get("beta", 0.01))

    post = d1871.get("model_posterior", {})
    models = {
        "M_A": {"nodes": [2, 5, 8, 11], "prior": float(post.get("M_A_2_5_8_11_or_2plus3n", 1/3))},
        "M_B": {"nodes": [2, 8, 14], "prior": float(post.get("M_B_2_8_14", 1/3))},
        "M_MIX": {"nodes": [2, 5, 8, 11, 14], "prior": float(post.get("M_MIX_A_and_B", 1/3))},
    }

    base_rows = {int(x["seed"]): x for x in d1874.get("rows", [])}

    rng = np.random.default_rng(187500)
    sigma_grid = [0.40, 0.70, 1.00, 1.40]

    rows_out = []

    for r in d1739.get("rows", []):
        seed = int(r["seed"])
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        base = base_rows.get(seed, {}).get("best", None)
        if base is None:
            base = {
                "omega": float(r["omega_hat"]),
                "beta": float(r["beta_hat"]),
                "rmse": float(r["rmse"]),
                "canon_score": canon_score_omega_beta(float(r["omega_hat"]), float(r["beta_hat"]), t_omega, t_beta),
            }
        base_score = float(base["canon_score"])
        base_rmse = float(base["rmse"])

        best = {
            "model": "BASELINE",
            "omega": float(base["omega"]),
            "phi_fixed": t_phi,
            "beta": float(base["beta"]),
            "amp": 0.0,
            "kappa": 0.0,
            "sigma": None,
            "rmse": base_rmse,
            "canon_score": base_score,
            "objective": base_rmse + 0.20 * (1.0 - base_score),
        }

        for mid, md in models.items():
            prior = max(1e-6, md["prior"])
            nodes = [n for n in md["nodes"] if n <= len(y) + 2]

            for sigma in sigma_grid:
                nf = node_field(d, nodes, sigma)

                # anchors near canonical + shifted regime
                anchor_params = [
                    (0.78, 0.010, 0.35),
                    (0.75, 0.020, 0.45),
                    (0.20, 0.250, 0.40),
                ]
                for om, be, ka in anchor_params:
                    amp = estimate_amp(y, d, om, t_phi, be, ka, nf)
                    yhat = model_signal(d, om, t_phi, be, amp, ka, nf)
                    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
                    cscore = canon_score_omega_beta(om, be, t_omega, t_beta)
                    obj = rmse + 0.20 * (1.0 - cscore) + 0.03 * (-math.log(prior)) + 0.01 * ka
                    if obj < best["objective"]:
                        best = {
                            "model": mid,
                            "omega": om,
                            "phi_fixed": t_phi,
                            "beta": be,
                            "amp": amp,
                            "kappa": ka,
                            "sigma": sigma,
                            "rmse": rmse,
                            "canon_score": cscore,
                            "objective": obj,
                        }

                for _ in range(280):
                    coin = rng.random()
                    if coin < 0.55:
                        om = rng.uniform(0.62, 0.92)
                        be = rng.uniform(0.004, 0.06)
                    elif coin < 0.85:
                        om = rng.uniform(0.10, 0.35)
                        be = rng.uniform(0.10, 0.30)
                    else:
                        om = rng.uniform(0.08, 0.95)
                        be = rng.uniform(0.004, 0.30)
                    ka = rng.uniform(0.05, 0.95)

                    amp = estimate_amp(y, d, om, t_phi, be, ka, nf)
                    yhat = model_signal(d, om, t_phi, be, amp, ka, nf)
                    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
                    cscore = canon_score_omega_beta(om, be, t_omega, t_beta)
                    obj = rmse + 0.20 * (1.0 - cscore) + 0.03 * (-math.log(prior)) + 0.01 * ka

                    if obj < best["objective"]:
                        best = {
                            "model": mid,
                            "omega": om,
                            "phi_fixed": t_phi,
                            "beta": be,
                            "amp": amp,
                            "kappa": ka,
                            "sigma": sigma,
                            "rmse": rmse,
                            "canon_score": cscore,
                            "objective": obj,
                        }

        rows_out.append(
            {
                "seed": seed,
                "baseline": {
                    "omega": float(base["omega"]),
                    "beta": float(base["beta"]),
                    "rmse": base_rmse,
                    "canon_score": base_score,
                },
                "anchored_best": best,
                "gains": {
                    "rmse_gain": base_rmse - best["rmse"],
                    "canon_gain": best["canon_score"] - base_score,
                },
            }
        )

    rmse_b = np.array([x["baseline"]["rmse"] for x in rows_out], dtype=float)
    rmse_a = np.array([x["anchored_best"]["rmse"] for x in rows_out], dtype=float)
    cs_b = np.array([x["baseline"]["canon_score"] for x in rows_out], dtype=float)
    cs_a = np.array([x["anchored_best"]["canon_score"] for x in rows_out], dtype=float)

    counts = {}
    for x in rows_out:
        m = x["anchored_best"]["model"]
        counts[m] = counts.get(m, 0) + 1

    summary = {
        "n_profiles": len(rows_out),
        "baseline_rmse_median": float(np.median(rmse_b)),
        "anchored_rmse_median": float(np.median(rmse_a)),
        "baseline_canon_median": float(np.median(cs_b)),
        "anchored_canon_median": float(np.median(cs_a)),
        "rmse_improved_fraction": float(np.mean(rmse_a < rmse_b)),
        "canon_improved_fraction": float(np.mean(cs_a > cs_b)),
        "canon_gain_median": float(np.median(cs_a - cs_b)),
        "best_model_counts": counts,
    }

    rmse_ratio = summary["anchored_rmse_median"] / max(summary["baseline_rmse_median"], 1e-12)

    if (
        summary["anchored_canon_median"] >= 1e-3
        and summary["canon_gain_median"] > 1e-4
        and summary["rmse_improved_fraction"] >= 0.60
        and rmse_ratio <= 1.10
    ):
        verdict = "CANON_ANCHORED_CONSTRAINED_FIT_STRONG_PROGRESS"
    elif (
        summary["canon_gain_median"] > 0
        and summary["canon_improved_fraction"] >= 0.60
        and rmse_ratio <= 1.25
    ):
        verdict = "CANON_ANCHORED_CONSTRAINED_FIT_PARTIAL_PROGRESS"
    elif summary["canon_gain_median"] > 0:
        verdict = "CANON_ANCHORED_CONSTRAINED_FIT_TRADEOFF_NOT_ACCEPTABLE"
    else:
        verdict = "CANON_ANCHORED_CONSTRAINED_FIT_WEAK"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target": {
            "omega": t_omega,
            "phi_fixed": t_phi,
            "beta": t_beta,
        },
        "summary": summary,
        "rows": rows_out,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1875: CANON-ANCHORED CONSTRAINED FIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Summary",
        f"- n_profiles: {summary['n_profiles']}",
        f"- rmse median: {summary['baseline_rmse_median']:.4f} -> {summary['anchored_rmse_median']:.4f}",
        f"- canon score median: {summary['baseline_canon_median']:.4e} -> {summary['anchored_canon_median']:.4e}",
        f"- rmse improved fraction: {summary['rmse_improved_fraction']:.3f}",
        f"- canon improved fraction: {summary['canon_improved_fraction']:.3f}",
        f"- canon gain median: {summary['canon_gain_median']:.4e}",
        f"- rmse ratio anchored/base: {rmse_ratio:.3f}",
        "",
        "## Best Model Counts",
    ]

    for k, v in sorted(counts.items(), key=lambda x: x[1], reverse=True):
        lines.append(f"- {k}: {v}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1875] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1875] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
