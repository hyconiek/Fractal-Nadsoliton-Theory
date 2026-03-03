#!/usr/bin/env python3
"""
QW-1874: Targeted beta-omega orthogonal forcing.

Fits node-augmented model with explicit gradient-orthogonality penalty between
omega and beta channels to reduce identifiability degeneracy.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1874_beta_omega_orthogonal_forcing.json"
OUT_MD = ROOT / "RAPORT_QW1874_BETA_OMEGA_ORTHOGONAL_FORCING.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return abs(d)


def node_field(d: np.ndarray, nodes: List[int], sigma: float) -> np.ndarray:
    if not nodes:
        return np.zeros_like(d)
    terms = [np.exp(-0.5 * ((d - float(n)) / sigma) ** 2) for n in nodes]
    return np.clip(np.max(np.vstack(terms), axis=0), 0.0, 1.0)


def model_signal(d: np.ndarray, omega: float, phi: float, beta: float, amp: float, kappa: float, nf: np.ndarray) -> np.ndarray:
    env = 1.0 / (1.0 + beta * d)
    mod = np.clip(1.0 - kappa * nf, 0.12, 1.0)
    return amp * np.cos(omega * d + phi) * env * mod


def grad_vectors(d: np.ndarray, omega: float, phi: float, beta: float, amp: float, kappa: float, nf: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    env = 1.0 / (1.0 + beta * d)
    mod = np.clip(1.0 - kappa * nf, 0.12, 1.0)

    g_omega = -amp * d * np.sin(omega * d + phi) * env * mod
    g_beta = -amp * np.cos(omega * d + phi) * d * (env ** 2) * mod
    return g_omega, g_beta


def corr_abs(x: np.ndarray, y: np.ndarray) -> float:
    sx = float(np.std(x))
    sy = float(np.std(y))
    if sx <= 1e-12 or sy <= 1e-12:
        return 1.0
    c = float(np.corrcoef(x, y)[0, 1])
    return abs(max(-1.0, min(1.0, c)))


def canon_score(omega: float, phi: float, beta: float, tgt: Dict[str, float]) -> float:
    z_o = abs(omega - tgt["omega"]) / 0.20
    z_p = circular_diff(phi, tgt["phi"]) / 0.30
    z_b = abs(beta - tgt["beta"]) / 0.030
    z = min(20.0, math.sqrt((z_o * z_o + z_p * z_p + z_b * z_b) / 3.0))
    return float(math.exp(-0.5 * z * z))


def estimate_amp(y: np.ndarray, d: np.ndarray, omega: float, phi: float, beta: float, kappa: float, nf: np.ndarray) -> float:
    basis = np.cos(omega * d + phi) / (1.0 + beta * d) * np.clip(1.0 - kappa * nf, 0.12, 1.0)
    num = float(np.dot(y, basis))
    den = float(np.dot(basis, basis))
    if den <= 1e-12:
        return 0.0
    return num / den


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1871 = read_json("report_qw1871_primary_node_evidence_corpus.json")

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    post = d1871.get("model_posterior", {})
    pA = float(post.get("M_A_2_5_8_11_or_2plus3n", 1.0 / 3.0))
    pB = float(post.get("M_B_2_8_14", 1.0 / 3.0))
    pM = float(post.get("M_MIX_A_and_B", 1.0 / 3.0))

    models = {
        "M_A": {"nodes": [2, 5, 8, 11], "prior": pA},
        "M_B": {"nodes": [2, 8, 14], "prior": pB},
        "M_MIX": {"nodes": [2, 5, 8, 11, 14], "prior": pM},
    }

    # Randomized targeted search (computationally tractable).
    sigma_grid = [0.40, 0.70, 1.00, 1.40]
    rng = np.random.default_rng(187400)

    rows_out = []

    for r in d1739.get("rows", []):
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        # Baseline point from 1739 row.
        base = {
            "omega": float(r["omega_hat"]),
            "phi": float(r["phi_hat"]),
            "beta": float(r["beta_hat"]),
            "amp": float(r["amp_hat"]),
        }

        yb = model_signal(d, base["omega"], base["phi"], base["beta"], base["amp"], 0.0, np.zeros_like(d))
        rmse_b = float(np.sqrt(np.mean((y - yb) ** 2)))
        gbo, gbb = grad_vectors(d, base["omega"], base["phi"], base["beta"], base["amp"], 0.0, np.zeros_like(d))
        corr_b = corr_abs(gbo, gbb)
        cscore_b = canon_score(base["omega"], base["phi"], base["beta"], t)

        best = {
            "model": "BASELINE",
            "omega": base["omega"],
            "phi": base["phi"],
            "beta": base["beta"],
            "amp": base["amp"],
            "kappa": 0.0,
            "sigma": None,
            "rmse": rmse_b,
            "corr_omega_beta_grad": corr_b,
            "canon_score": cscore_b,
            "objective": rmse_b + 0.10 * corr_b + 0.12 * (1.0 - cscore_b),
        }

        def eval_candidate(mid: str, prior: float, nodes: List[int], sigma: float, kappa: float, omega: float, phi: float, beta: float):
            nf = node_field(d, nodes, sigma)
            amp = estimate_amp(y, d, float(omega), float(phi), float(beta), float(kappa), nf)
            yhat = model_signal(d, float(omega), float(phi), float(beta), amp, float(kappa), nf)
            rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))

            go, gb = grad_vectors(d, float(omega), float(phi), float(beta), amp, float(kappa), nf)
            corr = corr_abs(go, gb)
            cscore = canon_score(float(omega), float(phi), float(beta), t)

            obj = rmse + 0.26 * corr + 0.12 * (1.0 - cscore) + 0.02 * (-math.log(prior))
            return {
                "model": mid,
                "omega": float(omega),
                "phi": float(phi),
                "beta": float(beta),
                "amp": float(amp),
                "kappa": float(kappa),
                "sigma": float(sigma),
                "rmse": rmse,
                "corr_omega_beta_grad": corr,
                "canon_score": cscore,
                "objective": obj,
            }

        for mid, md in models.items():
            nodes = [n for n in md["nodes"] if n <= len(y) + 2]
            prior = max(1e-6, md["prior"])

            # Directed anchor candidates.
            anchors = [
                (0.75, t["phi"], 0.01, 0.30),
                (0.78, t["phi"], 0.02, 0.45),
                (0.20, base["phi"], 0.25, 0.40),
                (base["omega"], base["phi"], base["beta"], 0.20),
            ]

            for sigma in sigma_grid:
                for omega, phi, beta, kappa in anchors:
                    cand = eval_candidate(mid, prior, nodes, sigma, kappa, omega, phi, beta)
                    if cand["objective"] < best["objective"]:
                        best = cand

                # Random sweep.
                for _ in range(160):
                    coin = rng.random()
                    if coin < 0.45:
                        omega = rng.uniform(0.60, 0.95)
                        beta = rng.uniform(0.005, 0.06)
                    elif coin < 0.80:
                        omega = rng.uniform(0.08, 0.30)
                        beta = rng.uniform(0.10, 0.30)
                    else:
                        omega = rng.uniform(0.08, 0.95)
                        beta = rng.uniform(0.005, 0.30)
                    phi = rng.uniform(-0.9, 0.9)
                    kappa = rng.uniform(0.08, 0.90)

                    cand = eval_candidate(mid, prior, nodes, sigma, kappa, omega, phi, beta)
                    if cand["objective"] < best["objective"]:
                        best = cand

            # Local refinement around best.
            for _ in range(220):
                omega = float(np.clip(rng.normal(best["omega"], 0.03), 0.08, 1.00))
                phi = float(np.clip(rng.normal(best["phi"], 0.10), -1.2, 1.2))
                beta = float(np.clip(rng.normal(best["beta"], 0.015), 0.003, 0.30))
                kappa = float(np.clip(rng.normal(best["kappa"], 0.08), 0.05, 0.95))
                sigma = float(np.clip(rng.normal(best["sigma"] if best["sigma"] is not None else 0.8, 0.20), 0.30, 1.60))

                cand = eval_candidate(mid, prior, nodes, sigma, kappa, omega, phi, beta)
                if cand["objective"] < best["objective"]:
                    best = cand

        rows_out.append(
            {
                "seed": int(r["seed"]),
                "baseline": {
                    "omega": base["omega"],
                    "phi": base["phi"],
                    "beta": base["beta"],
                    "rmse": rmse_b,
                    "corr_omega_beta_grad": corr_b,
                    "canon_score": cscore_b,
                },
                "best": best,
                "gains": {
                    "rmse_gain": rmse_b - best["rmse"],
                    "corr_reduction": corr_b - best["corr_omega_beta_grad"],
                    "canon_gain": best["canon_score"] - cscore_b,
                },
            }
        )

    rmse_b = np.array([x["baseline"]["rmse"] for x in rows_out], dtype=float)
    rmse_n = np.array([x["best"]["rmse"] for x in rows_out], dtype=float)
    corr_b = np.array([x["baseline"]["corr_omega_beta_grad"] for x in rows_out], dtype=float)
    corr_n = np.array([x["best"]["corr_omega_beta_grad"] for x in rows_out], dtype=float)
    cs_b = np.array([x["baseline"]["canon_score"] for x in rows_out], dtype=float)
    cs_n = np.array([x["best"]["canon_score"] for x in rows_out], dtype=float)

    cnt_models: Dict[str, int] = {}
    for x in rows_out:
        m = x["best"]["model"]
        cnt_models[m] = cnt_models.get(m, 0) + 1

    summary = {
        "n_profiles": len(rows_out),
        "baseline_rmse_median": float(np.median(rmse_b)),
        "orthogonal_rmse_median": float(np.median(rmse_n)),
        "baseline_corr_median": float(np.median(corr_b)),
        "orthogonal_corr_median": float(np.median(corr_n)),
        "baseline_canon_median": float(np.median(cs_b)),
        "orthogonal_canon_median": float(np.median(cs_n)),
        "rmse_improved_fraction": float(np.mean(rmse_n < rmse_b)),
        "corr_reduced_fraction": float(np.mean(corr_n < corr_b)),
        "canon_improved_fraction": float(np.mean(cs_n > cs_b)),
        "best_model_counts": cnt_models,
    }

    if (
        summary["corr_reduced_fraction"] >= 0.70
        and summary["rmse_improved_fraction"] >= 0.70
        and summary["orthogonal_canon_median"] > 5.0 * summary["baseline_canon_median"]
    ):
        verdict = "ORTHOGONAL_FORCING_STRONG_PROGRESS"
    elif summary["corr_reduced_fraction"] >= 0.55 and summary["rmse_improved_fraction"] >= 0.55:
        verdict = "ORTHOGONAL_FORCING_PARTIAL_PROGRESS"
    else:
        verdict = "ORTHOGONAL_FORCING_WEAK_PROGRESS"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "target_tuple": t,
        "node_model_prior": {"A": pA, "B": pB, "MIX": pM},
        "summary": summary,
        "rows": rows_out,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1874: BETA-OMEGA ORTHOGONAL FORCING",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Summary",
        f"- n_profiles: {summary['n_profiles']}",
        f"- rmse median: {summary['baseline_rmse_median']:.4f} -> {summary['orthogonal_rmse_median']:.4f}",
        f"- grad corr median: {summary['baseline_corr_median']:.4f} -> {summary['orthogonal_corr_median']:.4f}",
        f"- canonical score median: {summary['baseline_canon_median']:.4e} -> {summary['orthogonal_canon_median']:.4e}",
        f"- rmse improved fraction: {summary['rmse_improved_fraction']:.3f}",
        f"- corr reduced fraction: {summary['corr_reduced_fraction']:.3f}",
        f"- canon improved fraction: {summary['canon_improved_fraction']:.3f}",
        "",
        "## Best Model Counts",
    ]

    for k, v in sorted(cnt_models.items(), key=lambda x: x[1], reverse=True):
        lines.append(f"- {k}: {v}")

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1874] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1874] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
