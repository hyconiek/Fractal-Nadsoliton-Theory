#!/usr/bin/env python3
"""
QW-1881: Node-state identifiability stress under perturbations.

Uses locked configuration from QW-1880 and evaluates parameter stability under:
- additive noise perturbations,
- edge/boundary drift perturbations,
- combined perturbations.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1881_node_state_identifiability_stress.json"
OUT_MD = ROOT / "RAPORT_QW1881_NODE_STATE_IDENTIFIABILITY_STRESS.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def circular_diff(a: float, b: float) -> float:
    d = (a - b + math.pi) % (2.0 * math.pi) - math.pi
    return float(d)


def abs_circular_diff(a: float, b: float) -> float:
    return abs(circular_diff(a, b))


def canon_score(omega: float, phi: float, beta: float, t: Dict[str, float]) -> float:
    z_o = abs(omega - t["omega"]) / 0.20
    z_p = abs_circular_diff(phi, t["phi"]) / 0.30
    z_b = abs(beta - t["beta"]) / 0.030
    z = min(20.0, math.sqrt((z_o * z_o + z_p * z_p + z_b * z_b) / 3.0))
    return float(math.exp(-0.5 * z * z))


def node_indicator(d: np.ndarray, nodes: List[int]) -> np.ndarray:
    s = np.zeros_like(d, dtype=float)
    node_set = {int(x) for x in nodes}
    for i, x in enumerate(d):
        s[i] = 1.0 if int(x) in node_set else 0.0
    return s


def simulate_node_state(d: np.ndarray, omega: float, phi: float, rho: float, xi: float, zeta: float, nodes: List[int]) -> np.ndarray:
    ind = node_indicator(d, nodes)
    st = np.zeros_like(d, dtype=float)
    prev = 0.0
    for i, di in enumerate(d):
        drive = xi * math.sin(omega * float(di) + phi)
        node_term = -zeta * ind[i]
        cur = rho * prev + drive + node_term
        st[i] = cur
        prev = cur
    return st


def estimate_amp_gamma(y: np.ndarray, d: np.ndarray, omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, nodes: List[int]) -> Dict[str, float]:
    env_basis = np.cos(omega * d + phi) / (1.0 + beta * d)
    st_basis = simulate_node_state(d, omega, phi, rho, xi, zeta, nodes)
    X = np.column_stack([env_basis, st_basis])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    amp = float(coef[0])
    gamma = float(coef[1])
    yhat = X @ coef

    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    return {
        "amp": amp,
        "gamma": gamma,
        "rmse": rmse,
        "yhat": yhat,
    }


def fit_profile(
    y: np.ndarray,
    d: np.ndarray,
    nodes: List[int],
    prior_mean: Dict[str, float],
    prior_scale: Dict[str, float],
    target: Dict[str, float],
    lambda_c: float,
    lambda_p: float,
    rng: np.random.Generator,
    n_samples: int = 180,
) -> Dict:
    def prior_penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float) -> float:
        z1 = (omega - prior_mean["omega"]) / prior_scale["omega"]
        z2 = circular_diff(phi, prior_mean["phi"]) / prior_scale["phi"]
        z3 = (beta - prior_mean["beta"]) / prior_scale["beta"]
        z4 = (rho - prior_mean["rho"]) / prior_scale["rho"]
        z5 = (xi - prior_mean["xi"]) / prior_scale["xi"]
        z6 = (zeta - prior_mean["zeta"]) / prior_scale["zeta"]
        return float((z1*z1 + z2*z2 + z3*z3 + z4*z4 + z5*z5 + z6*z6) / 6.0)

    best = None

    def eval_one(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float):
        est = estimate_amp_gamma(y, d, omega, phi, beta, rho, xi, zeta, nodes)
        cscore = canon_score(omega, phi, beta, target)
        pp = prior_penalty(omega, phi, beta, rho, xi, zeta)
        obj = est["rmse"] + lambda_c * (1.0 - cscore) + lambda_p * pp
        return {
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "rho": rho,
            "xi": xi,
            "zeta": zeta,
            "amp": est["amp"],
            "gamma": est["gamma"],
            "rmse": est["rmse"],
            "canon_score": cscore,
            "objective": obj,
        }

    anchors = [
        (prior_mean["omega"], prior_mean["phi"], prior_mean["beta"], prior_mean["rho"], prior_mean["xi"], prior_mean["zeta"]),
        (target["omega"], target["phi"], target["beta"], 0.60, 0.20, 0.20),
    ]

    for a in anchors:
        cand = eval_one(*a)
        if best is None or cand["objective"] < best["objective"]:
            best = cand

    for _ in range(n_samples):
        omega = float(np.clip(rng.normal(best["omega"], 0.05), 0.08, 0.95))
        phi = float(np.clip(rng.normal(best["phi"], 0.18), -1.4, 1.4))
        beta = float(np.clip(rng.normal(best["beta"], 0.02), 0.003, 0.30))
        rho = float(np.clip(rng.normal(best["rho"], 0.10), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.06), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.06), 0.01, 0.70))
        cand = eval_one(omega, phi, beta, rho, xi, zeta)
        if cand["objective"] < best["objective"]:
            best = cand

    return best


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1880 = read_json("report_qw1880_node_state_strict_oos.json")

    target = d1880.get("test_rows", [{}])[0].get("fit", None)
    if target is None:
        # fallback to canonical from 1862 if no test rows, but should not happen.
        d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
        target_tuple = d1862.get("target_tuple", {})
    else:
        d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
        target_tuple = d1862.get("target_tuple", {})

    t = {
        "omega": float(target_tuple.get("omega", math.pi / 4.0)),
        "phi": float(target_tuple.get("phi", math.pi / 6.0)),
        "beta": float(target_tuple.get("beta", 0.01)),
    }

    prior_mean = d1880.get("prior_from_train", {}).get("mean", {})
    prior_scale = d1880.get("prior_from_train", {}).get("scale", {})

    lc = float(d1880.get("best_hyperparams", {}).get("lambda_c", 0.1))
    lp = float(d1880.get("best_hyperparams", {}).get("lambda_p", 0.1))

    nodes = d1880.get("protocol", {}).get("nodes_used", [2, 5, 8, 11])

    # Stress on combined val+test rows from 1880 protocol scope.
    split_seeds = set(int(x["seed"]) for x in d1880.get("test_rows", []))
    split_seeds.update(int(x["seed"]) for x in d1880.get("test_rows", []))

    # Fallback: if only few rows, include all rows.
    source_rows = [r for r in d1739.get("rows", []) if int(r["seed"]) in split_seeds]
    if len(source_rows) < 4:
        source_rows = d1739.get("rows", [])

    scenarios = [
        {"name": "noise_low", "noise_scale": 0.5, "edge": 0.0},
        {"name": "noise_mid", "noise_scale": 1.0, "edge": 0.0},
        {"name": "noise_high", "noise_scale": 1.5, "edge": 0.0},
        {"name": "edge_mid", "noise_scale": 0.7, "edge": 0.06},
        {"name": "combined_high", "noise_scale": 1.3, "edge": 0.08},
    ]

    rng = np.random.default_rng(188100)

    results = []

    for r in source_rows:
        seed = int(r["seed"])
        y0 = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y0) + 1, dtype=float)

        # baseline estimate under locked config
        fit0 = fit_profile(y0, d, nodes, prior_mean, prior_scale, t, lc, lp, rng=rng, n_samples=220)
        resid0 = y0 - estimate_amp_gamma(y0, d, fit0["omega"], fit0["phi"], fit0["beta"], fit0["rho"], fit0["xi"], fit0["zeta"], nodes)["yhat"]
        sigma0 = float(np.std(resid0))
        sigma0 = max(0.005, sigma0)

        for sc in scenarios:
            for rep in range(36):
                noise = rng.normal(0.0, sc["noise_scale"] * sigma0, size=len(y0))

                # boundary/edge perturbation: affine tilt + endpoint offsets
                edge = sc["edge"]
                x = (d - np.mean(d)) / max(np.std(d), 1e-9)
                tilt = 1.0 + edge * rng.normal(0.0, 1.0) * x
                endpoint = np.zeros_like(y0)
                endpoint[0] += edge * rng.normal(0.0, sigma0)
                endpoint[-1] += edge * rng.normal(0.0, sigma0)

                y = y0 * tilt + noise + endpoint

                fit = fit_profile(y, d, nodes, prior_mean, prior_scale, t, lc, lp, rng=rng, n_samples=140)

                nonboundary = bool(0.09 < fit["omega"] < 0.92 and 0.004 < fit["beta"] < 0.28 and 0.28 < fit["rho"] < 0.97)

                results.append(
                    {
                        "seed": seed,
                        "scenario": sc["name"],
                        "rep": rep,
                        "fit": fit,
                        "nonboundary": nonboundary,
                    }
                )

    omega = np.array([x["fit"]["omega"] for x in results], dtype=float)
    phi = np.array([x["fit"]["phi"] for x in results], dtype=float)
    beta = np.array([x["fit"]["beta"] for x in results], dtype=float)

    omega_iqr = float(np.quantile(omega, 0.75) - np.quantile(omega, 0.25))
    beta_iqr = float(np.quantile(beta, 0.75) - np.quantile(beta, 0.25))

    # circular spread proxy
    phi_cstd = float(np.sqrt(-2.0 * np.log(max(1e-12, np.sqrt(np.mean(np.sin(phi)) ** 2 + np.mean(np.cos(phi)) ** 2)))))

    nonboundary_rate = float(np.mean([1.0 if x["nonboundary"] else 0.0 for x in results]))

    if np.std(omega) > 1e-12 and np.std(beta) > 1e-12:
        corr_ob = float(np.corrcoef(omega, beta)[0, 1])
    else:
        corr_ob = 0.0

    # scenario-wise stability
    scenario_summary = {}
    for sc in scenarios:
        name = sc["name"]
        rr = [x for x in results if x["scenario"] == name]
        om = np.array([x["fit"]["omega"] for x in rr], dtype=float)
        be = np.array([x["fit"]["beta"] for x in rr], dtype=float)
        cs = np.array([x["fit"]["canon_score"] for x in rr], dtype=float)
        scenario_summary[name] = {
            "n": len(rr),
            "omega_iqr": float(np.quantile(om, 0.75) - np.quantile(om, 0.25)),
            "beta_iqr": float(np.quantile(be, 0.75) - np.quantile(be, 0.25)),
            "canon_median": float(np.median(cs)),
            "nonboundary_rate": float(np.mean([1.0 if x["nonboundary"] else 0.0 for x in rr])),
        }

    stability_index = (
        0.30 * max(0.0, min(1.0, 1.0 - omega_iqr / 0.12))
        + 0.30 * max(0.0, min(1.0, 1.0 - beta_iqr / 0.05))
        + 0.15 * max(0.0, min(1.0, 1.0 - phi_cstd / 0.8))
        + 0.15 * max(0.0, min(1.0, 1.0 - abs(corr_ob) / 0.6))
        + 0.10 * nonboundary_rate
    )

    if stability_index >= 0.70 and nonboundary_rate >= 0.70 and abs(corr_ob) <= 0.45:
        verdict = "NODE_STATE_IDENTIFIABILITY_STRESS_PASS"
    elif stability_index >= 0.50 and nonboundary_rate >= 0.50:
        verdict = "NODE_STATE_IDENTIFIABILITY_STRESS_PARTIAL"
    else:
        verdict = "NODE_STATE_IDENTIFIABILITY_STRESS_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_source_profiles": len(source_rows),
            "n_results": len(results),
            "locked_lambda_c": lc,
            "locked_lambda_p": lp,
            "locked_nodes": nodes,
            "scenarios": scenarios,
        },
        "global_metrics": {
            "omega_iqr": omega_iqr,
            "phi_circular_std": phi_cstd,
            "beta_iqr": beta_iqr,
            "corr_omega_beta": corr_ob,
            "nonboundary_rate": nonboundary_rate,
            "stability_index": float(stability_index),
        },
        "scenario_summary": scenario_summary,
        "rows_head": results[:500],
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    gm = out["global_metrics"]
    lines = [
        "# RAPORT QW-1881: NODE-STATE IDENTIFIABILITY STRESS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- source profiles: {out['protocol']['n_source_profiles']}",
        f"- stress samples: {out['protocol']['n_results']}",
        "",
        "## Global Metrics",
        f"- omega_iqr: {gm['omega_iqr']:.4f}",
        f"- phi_circular_std: {gm['phi_circular_std']:.4f}",
        f"- beta_iqr: {gm['beta_iqr']:.4f}",
        f"- corr(omega,beta): {gm['corr_omega_beta']:.4f}",
        f"- nonboundary_rate: {gm['nonboundary_rate']:.3f}",
        f"- stability_index: {gm['stability_index']:.3f}",
        "",
        "## Scenario Summary",
    ]

    for k, v in scenario_summary.items():
        lines.append(
            f"- {k}: n={v['n']}, omega_iqr={v['omega_iqr']:.4f}, beta_iqr={v['beta_iqr']:.4f}, canon_median={v['canon_median']:.4f}, nonboundary={v['nonboundary_rate']:.3f}"
        )

    lines += ["", "## Artifacts", f"- JSON: `{OUT_JSON.name}`"]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1881] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1881] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
