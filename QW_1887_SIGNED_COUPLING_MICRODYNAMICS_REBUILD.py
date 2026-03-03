#!/usr/bin/env python3
"""
QW-1887: Signed-coupling micromodel rebuild (QW1739/QW1740 direction).

Model extension over node-state baseline:
- Adds signed coupling drive term to dynamic state equation:
    s_t = rho * s_{t-1} + xi * sin(omega*d + phi) + eta * C_signed(d) - zeta * I_node(d)
- C_signed(d) encodes direction-sensitive influence from kernel nodes.

Protocol:
- Locked split strategy as QW-1880 (seed-level stratified).
- Train-only prior estimation for extended parameter set.
- Validation tuning of lambda_c and lambda_b.
- One locked test evaluation and comparison with QW-1880 + QW-1884.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1887_signed_coupling_micromodel_rebuild.json"
OUT_MD = ROOT / "RAPORT_QW1887_SIGNED_COUPLING_MICRODYNAMICS_REBUILD.md"


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
    out = np.zeros_like(d, dtype=float)
    node_set = {int(x) for x in nodes}
    for i, x in enumerate(d):
        out[i] = 1.0 if int(x) in node_set else 0.0
    return out


def signed_coupling_drive(d: np.ndarray, nodes: List[int], tau: float = 2.5) -> np.ndarray:
    # Signed influence: points to the right of a node contribute +, left contribute -.
    out = np.zeros_like(d, dtype=float)
    for i, di in enumerate(d):
        s = 0.0
        for n in nodes:
            delta = float(di - n)
            sign = 1.0 if delta > 0 else (-1.0 if delta < 0 else 0.0)
            s += sign * math.exp(-abs(delta) / tau)
        out[i] = s

    m = float(np.max(np.abs(out)))
    if m > 1e-12:
        out = out / m
    return out


def simulate_state_signed(
    d: np.ndarray,
    omega: float,
    phi: float,
    rho: float,
    xi: float,
    zeta: float,
    eta: float,
    nodes: List[int],
    tau: float,
) -> np.ndarray:
    ind = node_indicator(d, nodes)
    coupl = signed_coupling_drive(d, nodes, tau=tau)

    st = np.zeros_like(d, dtype=float)
    prev = 0.0
    for i, di in enumerate(d):
        drive = xi * math.sin(omega * float(di) + phi)
        node_term = -zeta * ind[i]
        coup_term = eta * coupl[i]
        cur = rho * prev + drive + node_term + coup_term
        st[i] = cur
        prev = cur
    return st


def estimate_amp_gamma(
    y: np.ndarray,
    d: np.ndarray,
    omega: float,
    phi: float,
    beta: float,
    rho: float,
    xi: float,
    zeta: float,
    eta: float,
    nodes: List[int],
    tau: float,
) -> Dict[str, float]:
    env_basis = np.cos(omega * d + phi) / (1.0 + beta * d)
    st_basis = simulate_state_signed(d, omega, phi, rho, xi, zeta, eta, nodes, tau)
    X = np.column_stack([env_basis, st_basis])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ coef
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    return {
        "amp": float(coef[0]),
        "gamma": float(coef[1]),
        "rmse": rmse,
    }


def boundary_penalty(omega: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> float:
    prefs = [
        (omega, 0.08, 0.95, 0.12, 0.90),
        (beta, 0.003, 0.30, 0.008, 0.25),
        (rho, 0.25, 0.98, 0.35, 0.93),
        (xi, 0.01, 0.65, 0.03, 0.58),
        (zeta, 0.01, 0.70, 0.03, 0.55),
        (eta, -0.70, 0.70, -0.45, 0.45),
    ]
    pen = 0.0
    for v, L, U, l, u in prefs:
        if v < l:
            pen += ((l - v) / max(l - L, 1e-9)) ** 2
        elif v > u:
            pen += ((v - u) / max(U - u, 1e-9)) ** 2
    return float(pen / len(prefs))


def fit_profile(
    y: np.ndarray,
    d: np.ndarray,
    nodes: List[int],
    prior_mean: Dict[str, float],
    prior_scale: Dict[str, float],
    target: Dict[str, float],
    lambda_c: float,
    lambda_p: float,
    lambda_b: float,
    rng: np.random.Generator,
    tau: float,
    n_samples: int,
) -> Dict:
    def prior_penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> float:
        z1 = (omega - prior_mean["omega"]) / prior_scale["omega"]
        z2 = circular_diff(phi, prior_mean["phi"]) / prior_scale["phi"]
        z3 = (beta - prior_mean["beta"]) / prior_scale["beta"]
        z4 = (rho - prior_mean["rho"]) / prior_scale["rho"]
        z5 = (xi - prior_mean["xi"]) / prior_scale["xi"]
        z6 = (zeta - prior_mean["zeta"]) / prior_scale["zeta"]
        z7 = (eta - prior_mean["eta"]) / prior_scale["eta"]
        return float((z1 * z1 + z2 * z2 + z3 * z3 + z4 * z4 + z5 * z5 + z6 * z6 + z7 * z7) / 7.0)

    def eval_one(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> Dict:
        est = estimate_amp_gamma(y, d, omega, phi, beta, rho, xi, zeta, eta, nodes, tau=tau)
        cscore = canon_score(omega, phi, beta, target)
        pp = prior_penalty(omega, phi, beta, rho, xi, zeta, eta)
        bp = boundary_penalty(omega, beta, rho, xi, zeta, eta)
        obj = est["rmse"] + lambda_c * (1.0 - cscore) + lambda_p * pp + lambda_b * bp
        return {
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "rho": rho,
            "xi": xi,
            "zeta": zeta,
            "eta": eta,
            "amp": est["amp"],
            "gamma": est["gamma"],
            "rmse": est["rmse"],
            "canon_score": cscore,
            "objective": obj,
        }

    anchors = [
        (
            prior_mean["omega"],
            prior_mean["phi"],
            prior_mean["beta"],
            prior_mean["rho"],
            prior_mean["xi"],
            prior_mean["zeta"],
            prior_mean["eta"],
        ),
        (target["omega"], target["phi"], target["beta"], 0.60, 0.20, 0.20, 0.0),
        (0.70, prior_mean["phi"], 0.03, 0.75, 0.20, 0.20, 0.2),
        (0.20, prior_mean["phi"], 0.20, 0.75, 0.20, 0.20, -0.2),
    ]

    best = None
    for a in anchors:
        cand = eval_one(*a)
        if best is None or cand["objective"] < best["objective"]:
            best = cand

    for _ in range(n_samples):
        omega = float(np.clip(rng.normal(best["omega"], 0.045), 0.08, 0.95))
        phi = float(np.clip(rng.normal(best["phi"], 0.16), -1.4, 1.4))
        beta = float(np.clip(rng.normal(best["beta"], 0.018), 0.003, 0.30))
        rho = float(np.clip(rng.normal(best["rho"], 0.10), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.06), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.06), 0.01, 0.70))
        eta = float(np.clip(rng.normal(best["eta"], 0.10), -0.70, 0.70))

        cand = eval_one(omega, phi, beta, rho, xi, zeta, eta)
        if cand["objective"] < best["objective"]:
            best = cand

    return best


def split_dataset(rows: List[Dict], seed: int = 188000) -> Dict[str, List[Dict]]:
    rng = np.random.default_rng(seed)
    by_n = {}
    for r in rows:
        by_n.setdefault(int(r["n_nodes"]), []).append(r)

    split = {"train": [], "val": [], "test": []}
    for n in sorted(by_n.keys()):
        group = list(by_n[n])
        idx = np.arange(len(group))
        rng.shuffle(idx)
        group = [group[i] for i in idx]

        n_train = min(8, len(group) - 4)
        n_val = min(3, len(group) - n_train - 1)

        split["train"].extend(group[:n_train])
        split["val"].extend(group[n_train : n_train + n_val])
        split["test"].extend(group[n_train + n_val :])

    return split


def robust_phi_mean(phi_vals: List[float]) -> float:
    x = np.array(phi_vals, dtype=float)
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def robust_scale(x: np.ndarray, floor: float) -> float:
    iqr = float(np.quantile(x, 0.75) - np.quantile(x, 0.25))
    s = iqr / 1.349 if iqr > 0 else float(np.std(x))
    return float(max(floor, s if np.isfinite(s) else floor))


def nonboundary_flag(f: Dict) -> bool:
    return bool(
        0.12 < f["omega"] < 0.90
        and 0.008 < f["beta"] < 0.25
        and 0.35 < f["rho"] < 0.93
        and 0.03 < f["xi"] < 0.58
        and abs(f["eta"]) < 0.45
    )


def summarize(rows: List[Dict]) -> Dict:
    rmse = np.array([r["fit"]["rmse"] for r in rows], dtype=float)
    cs = np.array([r["fit"]["canon_score"] for r in rows], dtype=float)
    nb = np.array([1.0 if r["nonboundary"] else 0.0 for r in rows], dtype=float)
    eta_abs = np.array([abs(r["fit"]["eta"]) for r in rows], dtype=float)
    return {
        "n": len(rows),
        "rmse_median": float(np.median(rmse)) if len(rmse) else None,
        "canon_median": float(np.median(cs)) if len(cs) else None,
        "nonboundary_rate": float(np.mean(nb)) if len(nb) else None,
        "eta_abs_median": float(np.median(eta_abs)) if len(eta_abs) else None,
    }


def fit_partition(
    part: List[Dict],
    nodes: List[int],
    prior_mean: Dict[str, float],
    prior_scale: Dict[str, float],
    t: Dict[str, float],
    lambda_c: float,
    lambda_p: float,
    lambda_b: float,
    rng: np.random.Generator,
    tau: float,
    n_samples: int,
) -> Tuple[List[Dict], Dict]:
    rows = []
    for r in part:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        fit = fit_profile(
            y,
            d,
            nodes,
            prior_mean,
            prior_scale,
            t,
            lambda_c,
            lambda_p,
            lambda_b,
            rng,
            tau=tau,
            n_samples=n_samples,
        )
        rows.append({"seed": int(r["seed"]), "fit": fit, "nonboundary": nonboundary_flag(fit)})
    return rows, summarize(rows)


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1880 = read_json("report_qw1880_node_state_strict_oos.json")
    d1884 = read_json("report_qw1884_node_state_pareto_oos_rebalancing.json")

    rows = sorted(d1739.get("rows", []), key=lambda x: int(x["seed"]))
    split = split_dataset(rows, seed=188000)

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    nodes = list(d1880.get("protocol", {}).get("nodes_used", [2, 5, 8, 11]))
    tau = 2.5

    rng = np.random.default_rng(188700)

    # Stage A: train-only weak prior estimation for extended parameter set.
    weak_mean = {
        "omega": 0.30,
        "phi": 0.0,
        "beta": 0.12,
        "rho": 0.60,
        "xi": 0.20,
        "zeta": 0.20,
        "eta": 0.0,
    }
    weak_scale = {
        "omega": 0.35,
        "phi": 0.8,
        "beta": 0.15,
        "rho": 0.25,
        "xi": 0.20,
        "zeta": 0.20,
        "eta": 0.35,
    }

    train_base = []
    for r in split["train"]:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        fit = fit_profile(
            y,
            d,
            nodes,
            weak_mean,
            weak_scale,
            t,
            lambda_c=0.0,
            lambda_p=0.0,
            lambda_b=0.0,
            rng=rng,
            tau=tau,
            n_samples=180,
        )
        train_base.append(fit)

    omega_vals = np.array([x["omega"] for x in train_base], dtype=float)
    phi_vals = [x["phi"] for x in train_base]
    beta_vals = np.array([x["beta"] for x in train_base], dtype=float)
    rho_vals = np.array([x["rho"] for x in train_base], dtype=float)
    xi_vals = np.array([x["xi"] for x in train_base], dtype=float)
    zeta_vals = np.array([x["zeta"] for x in train_base], dtype=float)
    eta_vals = np.array([x["eta"] for x in train_base], dtype=float)

    prior_mean = {
        "omega": float(np.median(omega_vals)),
        "phi": robust_phi_mean(phi_vals),
        "beta": float(np.median(beta_vals)),
        "rho": float(np.median(rho_vals)),
        "xi": float(np.median(xi_vals)),
        "zeta": float(np.median(zeta_vals)),
        "eta": float(np.median(eta_vals)),
    }

    prior_scale = {
        "omega": robust_scale(omega_vals, 0.06),
        "phi": max(0.15, float(np.std(np.unwrap(np.array(phi_vals, dtype=float))))),
        "beta": robust_scale(beta_vals, 0.02),
        "rho": robust_scale(rho_vals, 0.08),
        "xi": robust_scale(xi_vals, 0.06),
        "zeta": robust_scale(zeta_vals, 0.06),
        "eta": robust_scale(eta_vals, 0.08),
    }

    # Baseline validation rmse from 1880-style setting.
    lc_ref = float(d1880.get("best_hyperparams", {}).get("lambda_c", 0.2))
    lp_ref = float(d1880.get("best_hyperparams", {}).get("lambda_p", 0.0))

    _, val_base = fit_partition(
        split["val"],
        nodes,
        prior_mean,
        prior_scale,
        t,
        lambda_c=lc_ref,
        lambda_p=lp_ref,
        lambda_b=0.0,
        rng=rng,
        tau=tau,
        n_samples=220,
    )

    val_rmse_base = float(val_base["rmse_median"])

    # Stage B: tune lambda_c/lambda_b on val only.
    lambda_c_grid = [0.10, 0.20, 0.30]
    lambda_b_grid = [0.20, 0.35, 0.50]

    tuning = []
    feasible = []

    for lc in lambda_c_grid:
        for lb in lambda_b_grid:
            _, s = fit_partition(
                split["val"],
                nodes,
                prior_mean,
                prior_scale,
                t,
                lambda_c=lc,
                lambda_p=lp_ref,
                lambda_b=lb,
                rng=rng,
                tau=tau,
                n_samples=220,
            )
            rmse_ratio = float(s["rmse_median"] / max(val_rmse_base, 1e-12))
            obj = float(s["rmse_median"] + 0.20 * (1.0 - s["canon_median"]) + 0.15 * (1.0 - s["nonboundary_rate"]))
            ok = bool(s["nonboundary_rate"] >= 0.50 and s["canon_median"] >= 0.75 and rmse_ratio <= 1.45)
            row = {
                "lambda_c": lc,
                "lambda_b": lb,
                "val_summary": s,
                "val_rmse_ratio_vs_base": rmse_ratio,
                "val_objective": obj,
                "feasible": ok,
            }
            tuning.append(row)
            if ok:
                feasible.append(row)

    if feasible:
        selected = min(feasible, key=lambda r: r["val_objective"])
        selection_mode = "feasible_min_objective"
    else:
        selected = min(tuning, key=lambda r: r["val_objective"])
        selection_mode = "no_feasible_fallback"

    lc_star = float(selected["lambda_c"])
    lb_star = float(selected["lambda_b"])

    test_rows, test_summary = fit_partition(
        split["test"],
        nodes,
        prior_mean,
        prior_scale,
        t,
        lambda_c=lc_star,
        lambda_p=lp_ref,
        lambda_b=lb_star,
        rng=rng,
        tau=tau,
        n_samples=260,
    )

    s1880 = d1880.get("summaries", {}).get("test", {})
    s1884 = d1884.get("summaries", {}).get("test", {})

    delta_1880 = {
        "rmse_median_gain": float(float(s1880.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1880.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1880.get("nonboundary_rate", 0.0))),
        "rmse_ratio": float(test_summary["rmse_median"] / max(float(s1880.get("rmse_median", 1.0)), 1e-12)),
    }
    delta_1884 = {
        "rmse_median_gain": float(float(s1884.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1884.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1884.get("nonboundary_rate", 0.0))),
    }

    if (
        test_summary["nonboundary_rate"] >= 0.50
        and test_summary["canon_median"] >= 0.80
        and delta_1880["rmse_ratio"] <= 1.45
        and delta_1884["rmse_median_gain"] >= 0.010
    ):
        verdict = "SIGNED_COUPLING_OOS_IMPROVEMENT"
    elif (
        test_summary["nonboundary_rate"] >= 0.50
        and test_summary["canon_median"] >= 0.80
        and delta_1884["rmse_median_gain"] >= -0.010
    ):
        verdict = "SIGNED_COUPLING_PARTIAL"
    elif test_summary["nonboundary_rate"] > 0.0:
        verdict = "SIGNED_COUPLING_WEAK"
    else:
        verdict = "SIGNED_COUPLING_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "split_seed": 188000,
            "split_counts": {k: len(v) for k, v in split.items()},
            "nodes": nodes,
            "tau": tau,
            "selection_mode": selection_mode,
        },
        "prior_from_train": {
            "mean": prior_mean,
            "scale": prior_scale,
        },
        "val_baseline": val_base,
        "tuning_rows": tuning,
        "selected_hyperparams": {
            "lambda_c": lc_star,
            "lambda_p": lp_ref,
            "lambda_b": lb_star,
        },
        "summaries": {
            "test": test_summary,
            "test_1880_reference": s1880,
            "test_1884_reference": s1884,
        },
        "delta_vs_1880": delta_1880,
        "delta_vs_1884": delta_1884,
        "test_rows": test_rows,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1887: SIGNED-COUPLING MICRODYNAMICS REBUILD",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- selection_mode: `{selection_mode}`",
        f"- selected lambda_c/lambda_b: `{lc_star}` / `{lb_star}`",
        "",
        "## Test Summary",
        f"- rmse median: {test_summary['rmse_median']:.4f}",
        f"- canon median: {test_summary['canon_median']:.4f}",
        f"- nonboundary rate: {test_summary['nonboundary_rate']:.3f}",
        f"- |eta| median: {test_summary['eta_abs_median']:.4f}",
        "",
        "## Delta vs QW-1880",
        f"- rmse gain: {delta_1880['rmse_median_gain']:.4e}",
        f"- canon gain: {delta_1880['canon_median_gain']:.4e}",
        f"- nonboundary gain: {delta_1880['nonboundary_gain']:.4e}",
        f"- rmse ratio: {delta_1880['rmse_ratio']:.4f}",
        "",
        "## Delta vs QW-1884",
        f"- rmse gain: {delta_1884['rmse_median_gain']:.4e}",
        f"- canon gain: {delta_1884['canon_median_gain']:.4e}",
        f"- nonboundary gain: {delta_1884['nonboundary_gain']:.4e}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1887] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1887] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
