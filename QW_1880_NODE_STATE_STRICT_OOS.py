#!/usr/bin/env python3
"""
QW-1880: Strict OOS validation for node-state micromodel.

Protocol:
- Seed-level split (stratified by n_nodes): train/val/test.
- Priors estimated on train only.
- Hyperparameter tuning (lambda_c, lambda_p) on validation only.
- Single locked evaluation on test.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1880_node_state_strict_oos.json"
OUT_MD = ROOT / "RAPORT_QW1880_NODE_STATE_STRICT_OOS.md"


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
    ss_res = float(np.sum((y - yhat) ** 2))
    ss_tot = float(np.sum((y - np.mean(y)) ** 2))
    r2 = 1.0 if ss_tot <= 1e-15 else float(1.0 - ss_res / ss_tot)

    return {
        "amp": amp,
        "gamma": gamma,
        "rmse": rmse,
        "r2": r2,
        "yhat": yhat,
    }


def robust_phi_mean(phi_vals: List[float]) -> float:
    x = np.array(phi_vals, dtype=float)
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def split_dataset(rows: List[Dict], seed: int = 188000) -> Dict[str, List[Dict]]:
    rng = np.random.default_rng(seed)
    by_n = {}
    for r in rows:
        by_n.setdefault(int(r["n_nodes"]), []).append(r)

    train, val, test = [], [], []

    for n in sorted(by_n.keys()):
        group = list(by_n[n])
        idx = np.arange(len(group))
        rng.shuffle(idx)
        group = [group[i] for i in idx]

        # per group: 8 train, 3 val, rest test (for n=14 -> 3 test)
        n_train = min(8, len(group) - 4)
        n_val = min(3, len(group) - n_train - 1)

        train.extend(group[:n_train])
        val.extend(group[n_train : n_train + n_val])
        test.extend(group[n_train + n_val :])

    return {"train": train, "val": val, "test": test}


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
    n_samples: int = 280,
) -> Dict:
    def prior_penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float) -> float:
        z1 = (omega - prior_mean["omega"]) / prior_scale["omega"]
        z2 = circular_diff(phi, prior_mean["phi"]) / prior_scale["phi"]
        z3 = (beta - prior_mean["beta"]) / prior_scale["beta"]
        z4 = (rho - prior_mean["rho"]) / prior_scale["rho"]
        z5 = (xi - prior_mean["xi"]) / prior_scale["xi"]
        z6 = (zeta - prior_mean["zeta"]) / prior_scale["zeta"]
        return float((z1*z1 + z2*z2 + z3*z3 + z4*z4 + z5*z5 + z6*z6) / 6.0)

    anchors = [
        (prior_mean["omega"], prior_mean["phi"], prior_mean["beta"], prior_mean["rho"], prior_mean["xi"], prior_mean["zeta"]),
        (target["omega"], target["phi"], target["beta"], 0.60, 0.20, 0.20),
        (0.20, prior_mean["phi"], 0.25, 0.70, 0.12, 0.15),
    ]

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
            "r2": est["r2"],
            "canon_score": cscore,
            "prior_penalty": pp,
            "objective": obj,
        }

    for a in anchors:
        cand = eval_one(*a)
        if best is None or cand["objective"] < best["objective"]:
            best = cand

    # randomized guided search
    for _ in range(n_samples):
        coin = rng.random()
        if coin < 0.45:
            omega = float(np.clip(rng.normal(prior_mean["omega"], 0.08), 0.08, 0.95))
            beta = float(np.clip(rng.normal(prior_mean["beta"], 0.04), 0.005, 0.30))
        elif coin < 0.80:
            omega = float(rng.uniform(0.08, 0.30))
            beta = float(rng.uniform(0.10, 0.30))
        else:
            omega = float(rng.uniform(0.55, 0.95))
            beta = float(rng.uniform(0.005, 0.08))

        phi = float(np.clip(rng.normal(prior_mean["phi"], 0.35), -1.3, 1.3))
        rho = float(rng.uniform(0.30, 0.95))
        xi = float(rng.uniform(0.02, 0.55))
        zeta = float(rng.uniform(0.02, 0.60))

        cand = eval_one(omega, phi, beta, rho, xi, zeta)
        if cand["objective"] < best["objective"]:
            best = cand

    # local refinement
    for _ in range(180):
        omega = float(np.clip(rng.normal(best["omega"], 0.03), 0.08, 0.95))
        phi = float(np.clip(rng.normal(best["phi"], 0.12), -1.4, 1.4))
        beta = float(np.clip(rng.normal(best["beta"], 0.015), 0.003, 0.30))
        rho = float(np.clip(rng.normal(best["rho"], 0.08), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.05), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.05), 0.01, 0.70))

        cand = eval_one(omega, phi, beta, rho, xi, zeta)
        if cand["objective"] < best["objective"]:
            best = cand

    return best


def summarize_rows(rows: List[Dict]) -> Dict:
    rmse = np.array([r["fit"]["rmse"] for r in rows], dtype=float)
    cs = np.array([r["fit"]["canon_score"] for r in rows], dtype=float)
    nb = np.array([1.0 if r["nonboundary"] else 0.0 for r in rows], dtype=float)

    return {
        "n": len(rows),
        "rmse_median": float(np.median(rmse)) if len(rmse) else None,
        "rmse_mean": float(np.mean(rmse)) if len(rmse) else None,
        "canon_median": float(np.median(cs)) if len(cs) else None,
        "canon_mean": float(np.mean(cs)) if len(cs) else None,
        "nonboundary_rate": float(np.mean(nb)) if len(nb) else None,
        "canon_above_0p2_rate": float(np.mean(cs >= 0.2)) if len(cs) else None,
    }


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1871 = read_json("report_qw1871_primary_node_evidence_corpus.json")

    rows = d1739.get("rows", [])
    split = split_dataset(rows, seed=188000)

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    winner = d1871.get("winner", "M_A_2_5_8_11_or_2plus3n")
    if winner.startswith("M_A"):
        nodes = [2, 5, 8, 11]
    elif winner.startswith("M_B"):
        nodes = [2, 8, 14]
    else:
        nodes = [2, 5, 8, 11, 14]

    rng = np.random.default_rng(188001)

    # Stage A: train-only base fits for prior estimation.
    train_base = []
    for r in split["train"]:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        # weak broad priors in stage A
        prior_mean = {"omega": 0.30, "phi": 0.0, "beta": 0.12, "rho": 0.60, "xi": 0.20, "zeta": 0.20}
        prior_scale = {"omega": 0.35, "phi": 0.8, "beta": 0.15, "rho": 0.25, "xi": 0.20, "zeta": 0.20}

        fit = fit_profile(y, d, nodes, prior_mean, prior_scale, t, lambda_c=0.0, lambda_p=0.0, rng=rng, n_samples=220)
        train_base.append({"seed": int(r["seed"]), "fit": fit, "n_nodes": int(r["n_nodes"])})

    # Priors from train-only estimates.
    omega_vals = np.array([x["fit"]["omega"] for x in train_base], dtype=float)
    phi_vals = [x["fit"]["phi"] for x in train_base]
    beta_vals = np.array([x["fit"]["beta"] for x in train_base], dtype=float)
    rho_vals = np.array([x["fit"]["rho"] for x in train_base], dtype=float)
    xi_vals = np.array([x["fit"]["xi"] for x in train_base], dtype=float)
    zeta_vals = np.array([x["fit"]["zeta"] for x in train_base], dtype=float)

    prior_mean = {
        "omega": float(np.median(omega_vals)),
        "phi": robust_phi_mean(phi_vals),
        "beta": float(np.median(beta_vals)),
        "rho": float(np.median(rho_vals)),
        "xi": float(np.median(xi_vals)),
        "zeta": float(np.median(zeta_vals)),
    }

    def robust_scale(x: np.ndarray, floor: float) -> float:
        iqr = float(np.quantile(x, 0.75) - np.quantile(x, 0.25))
        s = iqr / 1.349 if iqr > 0 else float(np.std(x))
        return float(max(floor, s if np.isfinite(s) else floor))

    prior_scale = {
        "omega": robust_scale(omega_vals, 0.06),
        "phi": max(0.15, float(np.std(np.unwrap(np.array(phi_vals, dtype=float))))),
        "beta": robust_scale(beta_vals, 0.02),
        "rho": robust_scale(rho_vals, 0.08),
        "xi": robust_scale(xi_vals, 0.06),
        "zeta": robust_scale(zeta_vals, 0.06),
    }

    # Stage B: tune lambdas on val only.
    lambda_c_grid = [0.0, 0.05, 0.10, 0.20, 0.35]
    lambda_p_grid = [0.0, 0.05, 0.10, 0.20]

    tuning_rows = []
    best_combo = None

    for lc in lambda_c_grid:
        for lp in lambda_p_grid:
            val_fits = []
            for r in split["val"]:
                y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
                d = np.arange(1, len(y) + 1, dtype=float)
                fit = fit_profile(y, d, nodes, prior_mean, prior_scale, t, lc, lp, rng=rng, n_samples=240)

                nonboundary = bool(0.09 < fit["omega"] < 0.92 and 0.004 < fit["beta"] < 0.28 and 0.28 < fit["rho"] < 0.97)
                val_fits.append({"seed": int(r["seed"]), "fit": fit, "nonboundary": nonboundary})

            s = summarize_rows(val_fits)
            val_obj = float(s["rmse_median"] + 0.20 * (1.0 - s["canon_median"]) + 0.10 * (1.0 - s["nonboundary_rate"]))

            row = {
                "lambda_c": lc,
                "lambda_p": lp,
                "val_summary": s,
                "val_objective": val_obj,
            }
            tuning_rows.append(row)

            if best_combo is None or val_obj < best_combo["val_objective"]:
                best_combo = row

    lc_star = float(best_combo["lambda_c"])
    lp_star = float(best_combo["lambda_p"])

    # Stage C: locked evaluation on test.
    test_rows = []
    baseline_test_rows = []
    for r in split["test"]:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        fit = fit_profile(y, d, nodes, prior_mean, prior_scale, t, lc_star, lp_star, rng=rng, n_samples=280)

        nonboundary = bool(0.09 < fit["omega"] < 0.92 and 0.004 < fit["beta"] < 0.28 and 0.28 < fit["rho"] < 0.97)
        test_rows.append(
            {
                "seed": int(r["seed"]),
                "n_nodes": int(r["n_nodes"]),
                "fit": fit,
                "nonboundary": nonboundary,
            }
        )

        # baseline from QW-1739 row-level estimates
        b_omega = float(r["omega_hat"])
        b_phi = float(r["phi_hat"])
        b_beta = float(r["beta_hat"])
        b_rmse = float(r["rmse"])
        b_cs = canon_score(b_omega, b_phi, b_beta, t)
        b_nb = bool(0.09 < b_omega < 0.92 and 0.004 < b_beta < 0.28)
        baseline_test_rows.append(
            {
                "seed": int(r["seed"]),
                "fit": {
                    "omega": b_omega,
                    "phi": b_phi,
                    "beta": b_beta,
                    "rmse": b_rmse,
                    "canon_score": b_cs,
                },
                "nonboundary": b_nb,
            }
        )

    train_summary = summarize_rows([{"fit": x["fit"], "nonboundary": True} for x in train_base])
    val_summary_best = best_combo["val_summary"]
    test_summary = summarize_rows(test_rows)
    baseline_test_summary = summarize_rows(baseline_test_rows)

    delta_test = {
        "rmse_median_gain_vs_baseline": float(baseline_test_summary["rmse_median"] - test_summary["rmse_median"]),
        "canon_median_gain_vs_baseline": float(test_summary["canon_median"] - baseline_test_summary["canon_median"]),
        "nonboundary_gain_vs_baseline": float(test_summary["nonboundary_rate"] - baseline_test_summary["nonboundary_rate"]),
    }

    if (
        test_summary["canon_median"] >= 0.15
        and test_summary["nonboundary_rate"] >= 0.70
        and test_summary["rmse_median"] <= 1.12 * val_summary_best["rmse_median"]
        and delta_test["canon_median_gain_vs_baseline"] > 0
    ):
        verdict = "NODE_STATE_STRICT_OOS_PASS"
    elif (
        test_summary["canon_median"] >= 0.05
        and test_summary["nonboundary_rate"] >= 0.50
    ):
        verdict = "NODE_STATE_STRICT_OOS_PARTIAL"
    else:
        verdict = "NODE_STATE_STRICT_OOS_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "split_strategy": "stratified_seed_level_by_n_nodes",
            "split_counts": {k: len(v) for k, v in split.items()},
            "nodes_model_locked": winner,
            "nodes_used": nodes,
            "lambda_grid": {
                "lambda_c": lambda_c_grid,
                "lambda_p": lambda_p_grid,
            },
        },
        "prior_from_train": {
            "mean": prior_mean,
            "scale": prior_scale,
        },
        "tuning_rows": tuning_rows,
        "best_hyperparams": {
            "lambda_c": lc_star,
            "lambda_p": lp_star,
            "val_objective": best_combo["val_objective"],
        },
        "summaries": {
            "train": train_summary,
            "val_best": val_summary_best,
            "test": test_summary,
            "test_baseline": baseline_test_summary,
        },
        "test_delta_vs_baseline": delta_test,
        "test_rows": test_rows,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1880: NODE-STATE STRICT OOS",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- split train/val/test: {len(split['train'])}/{len(split['val'])}/{len(split['test'])}",
        f"- locked nodes model: {winner} -> {nodes}",
        f"- best lambdas: lambda_c={lc_star}, lambda_p={lp_star}",
        "",
        "## Summaries",
        f"- train rmse/canon/nonboundary: {train_summary['rmse_median']:.4f} / {train_summary['canon_median']:.4f} / {train_summary['nonboundary_rate']:.3f}",
        f"- val rmse/canon/nonboundary: {val_summary_best['rmse_median']:.4f} / {val_summary_best['canon_median']:.4f} / {val_summary_best['nonboundary_rate']:.3f}",
        f"- test rmse/canon/nonboundary: {test_summary['rmse_median']:.4f} / {test_summary['canon_median']:.4f} / {test_summary['nonboundary_rate']:.3f}",
        f"- test baseline rmse/canon/nonboundary: {baseline_test_summary['rmse_median']:.4f} / {baseline_test_summary['canon_median']:.4f} / {baseline_test_summary['nonboundary_rate']:.3f}",
        "",
        "## Test Delta vs Baseline",
        f"- rmse median gain: {delta_test['rmse_median_gain_vs_baseline']:.4e}",
        f"- canon median gain: {delta_test['canon_median_gain_vs_baseline']:.4e}",
        f"- nonboundary gain: {delta_test['nonboundary_gain_vs_baseline']:.4e}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1880] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1880] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
