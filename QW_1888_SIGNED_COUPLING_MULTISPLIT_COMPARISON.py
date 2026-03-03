#!/usr/bin/env python3
"""
QW-1888: Multi-split comparison of unsigned vs signed-coupling node-state models.

Control: unsigned (eta fixed to 0), hyperparams from QW-1884.
Treatment: signed-coupling (eta free), hyperparams from QW-1887.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1888_signed_coupling_multisplit_comparison.json"
OUT_MD = ROOT / "RAPORT_QW1888_SIGNED_COUPLING_MULTISPLIT_COMPARISON.md"


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


def simulate_state(
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
        cur = rho * prev + drive + node_term + eta * coupl[i]
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
    st_basis = simulate_state(d, omega, phi, rho, xi, zeta, eta, nodes, tau)
    X = np.column_stack([env_basis, st_basis])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ coef
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    return {"amp": float(coef[0]), "gamma": float(coef[1]), "rmse": rmse}


def robust_phi_mean(phi_vals: List[float]) -> float:
    x = np.array(phi_vals, dtype=float)
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def robust_scale(x: np.ndarray, floor: float) -> float:
    iqr = float(np.quantile(x, 0.75) - np.quantile(x, 0.25))
    s = iqr / 1.349 if iqr > 0 else float(np.std(x))
    return float(max(floor, s if np.isfinite(s) else floor))


def split_dataset(rows: List[Dict], seed: int) -> Dict[str, List[Dict]]:
    rng = np.random.default_rng(seed)
    by_n = {}
    for r in rows:
        by_n.setdefault(int(r["n_nodes"]), []).append(r)

    split = {"train": [], "val": [], "test": []}
    for n in sorted(by_n.keys()):
        g = list(by_n[n])
        idx = np.arange(len(g))
        rng.shuffle(idx)
        g = [g[i] for i in idx]
        n_train = min(8, len(g) - 4)
        n_val = min(3, len(g) - n_train - 1)
        split["train"].extend(g[:n_train])
        split["val"].extend(g[n_train : n_train + n_val])
        split["test"].extend(g[n_train + n_val :])

    return split


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
    eta_mode: str,
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
        est = estimate_amp_gamma(y, d, omega, phi, beta, rho, xi, zeta, eta, nodes, tau)
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

    eta0 = float(prior_mean["eta"]) if eta_mode == "free" else 0.0
    anchors = [
        (prior_mean["omega"], prior_mean["phi"], prior_mean["beta"], prior_mean["rho"], prior_mean["xi"], prior_mean["zeta"], eta0),
        (target["omega"], target["phi"], target["beta"], 0.60, 0.20, 0.20, 0.0),
        (0.70, prior_mean["phi"], 0.03, 0.75, 0.20, 0.20, 0.2 if eta_mode == "free" else 0.0),
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

        if eta_mode == "free":
            eta = float(np.clip(rng.normal(best["eta"], 0.10), -0.70, 0.70))
        else:
            eta = 0.0

        cand = eval_one(omega, phi, beta, rho, xi, zeta, eta)
        if cand["objective"] < best["objective"]:
            best = cand

    return best


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


def estimate_train_prior(
    train_rows: List[Dict],
    nodes: List[int],
    t: Dict[str, float],
    eta_mode: str,
    rng: np.random.Generator,
    tau: float,
) -> Tuple[Dict[str, float], Dict[str, float]]:
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

    fits = []
    for r in train_rows:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        f = fit_profile(
            y,
            d,
            nodes,
            weak_mean,
            weak_scale,
            t,
            lambda_c=0.0,
            lambda_p=0.0,
            lambda_b=0.0,
            eta_mode=eta_mode,
            rng=rng,
            tau=tau,
            n_samples=120,
        )
        fits.append(f)

    omega_vals = np.array([x["omega"] for x in fits], dtype=float)
    phi_vals = [x["phi"] for x in fits]
    beta_vals = np.array([x["beta"] for x in fits], dtype=float)
    rho_vals = np.array([x["rho"] for x in fits], dtype=float)
    xi_vals = np.array([x["xi"] for x in fits], dtype=float)
    zeta_vals = np.array([x["zeta"] for x in fits], dtype=float)
    eta_vals = np.array([x["eta"] for x in fits], dtype=float)

    mean = {
        "omega": float(np.median(omega_vals)),
        "phi": robust_phi_mean(phi_vals),
        "beta": float(np.median(beta_vals)),
        "rho": float(np.median(rho_vals)),
        "xi": float(np.median(xi_vals)),
        "zeta": float(np.median(zeta_vals)),
        "eta": float(np.median(eta_vals)),
    }
    scale = {
        "omega": robust_scale(omega_vals, 0.06),
        "phi": max(0.15, float(np.std(np.unwrap(np.array(phi_vals, dtype=float))))),
        "beta": robust_scale(beta_vals, 0.02),
        "rho": robust_scale(rho_vals, 0.08),
        "xi": robust_scale(xi_vals, 0.06),
        "zeta": robust_scale(zeta_vals, 0.06),
        "eta": robust_scale(eta_vals, 0.08),
    }
    return mean, scale


def eval_test(
    test_rows: List[Dict],
    nodes: List[int],
    prior_mean: Dict[str, float],
    prior_scale: Dict[str, float],
    t: Dict[str, float],
    eta_mode: str,
    lambda_c: float,
    lambda_p: float,
    lambda_b: float,
    rng: np.random.Generator,
    tau: float,
) -> Dict:
    rows = []
    for r in test_rows:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        fit = fit_profile(
            y,
            d,
            nodes,
            prior_mean,
            prior_scale,
            t,
            lambda_c=lambda_c,
            lambda_p=lambda_p,
            lambda_b=lambda_b,
            eta_mode=eta_mode,
            rng=rng,
            tau=tau,
            n_samples=220,
        )
        rows.append({"seed": int(r["seed"]), "fit": fit, "nonboundary": nonboundary_flag(fit)})
    return summarize(rows)


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1884 = read_json("report_qw1884_node_state_pareto_oos_rebalancing.json")
    d1887 = read_json("report_qw1887_signed_coupling_micromodel_rebuild.json")

    rows_all = sorted(d1739.get("rows", []), key=lambda x: int(x["seed"]))

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    nodes = [2, 5, 8, 11]
    tau = 2.5

    lc_ctrl = float(d1884.get("selected_hyperparams", {}).get("lambda_c", 0.2))
    lb_ctrl = float(d1884.get("selected_hyperparams", {}).get("lambda_b", 0.35))
    lc_trt = float(d1887.get("selected_hyperparams", {}).get("lambda_c", 0.2))
    lb_trt = float(d1887.get("selected_hyperparams", {}).get("lambda_b", 0.2))
    lp = 0.0

    split_seeds = list(range(188000, 188025))
    split_reports = []

    for i, s in enumerate(split_seeds):
        split = split_dataset(rows_all, seed=s)
        rng_ctrl = np.random.default_rng(188800 + i)
        rng_trt = np.random.default_rng(199800 + i)

        pm_ctrl, ps_ctrl = estimate_train_prior(split["train"], nodes, t, "fixed0", rng_ctrl, tau)
        pm_trt, ps_trt = estimate_train_prior(split["train"], nodes, t, "free", rng_trt, tau)

        ctrl = eval_test(
            split["test"],
            nodes,
            pm_ctrl,
            ps_ctrl,
            t,
            eta_mode="fixed0",
            lambda_c=lc_ctrl,
            lambda_p=lp,
            lambda_b=lb_ctrl,
            rng=rng_ctrl,
            tau=tau,
        )
        trt = eval_test(
            split["test"],
            nodes,
            pm_trt,
            ps_trt,
            t,
            eta_mode="free",
            lambda_c=lc_trt,
            lambda_p=lp,
            lambda_b=lb_trt,
            rng=rng_trt,
            tau=tau,
        )

        rmse_gain = float(ctrl["rmse_median"] - trt["rmse_median"])
        canon_gain = float(trt["canon_median"] - ctrl["canon_median"])
        nb_gain = float(trt["nonboundary_rate"] - ctrl["nonboundary_rate"])

        split_success = bool(rmse_gain >= 0.005 and canon_gain >= -0.05 and nb_gain >= -0.05)

        split_reports.append(
            {
                "split_seed": s,
                "counts": {k: len(v) for k, v in split.items()},
                "control_unsigned": ctrl,
                "treatment_signed": trt,
                "delta_signed_minus_unsigned": {
                    "rmse_gain": rmse_gain,
                    "canon_gain": canon_gain,
                    "nonboundary_gain": nb_gain,
                },
                "split_success": split_success,
            }
        )

    arr_rmse = np.array([r["delta_signed_minus_unsigned"]["rmse_gain"] for r in split_reports], dtype=float)
    arr_canon = np.array([r["delta_signed_minus_unsigned"]["canon_gain"] for r in split_reports], dtype=float)
    arr_nb = np.array([r["delta_signed_minus_unsigned"]["nonboundary_gain"] for r in split_reports], dtype=float)
    arr_succ = np.array([1.0 if r["split_success"] else 0.0 for r in split_reports], dtype=float)

    summary = {
        "n_splits": len(split_reports),
        "success_rate": float(np.mean(arr_succ)),
        "rmse_gain_median": float(np.median(arr_rmse)),
        "rmse_gain_q25": float(np.quantile(arr_rmse, 0.25)),
        "canon_gain_median": float(np.median(arr_canon)),
        "nonboundary_gain_median": float(np.median(arr_nb)),
        "nonboundary_gain_q25": float(np.quantile(arr_nb, 0.25)),
    }

    if (
        summary["success_rate"] >= 0.60
        and summary["rmse_gain_median"] >= 0.005
        and summary["canon_gain_median"] >= -0.05
    ):
        verdict = "SIGNED_COUPLING_MULTISPLIT_SUPERIOR"
    elif summary["success_rate"] >= 0.40 and summary["rmse_gain_median"] > 0.0:
        verdict = "SIGNED_COUPLING_MULTISPLIT_PARTIAL"
    else:
        verdict = "SIGNED_COUPLING_MULTISPLIT_NOT_SUPERIOR"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "split_seed_range": [split_seeds[0], split_seeds[-1]],
            "n_splits": len(split_seeds),
            "nodes": nodes,
            "tau": tau,
            "control": {"lambda_c": lc_ctrl, "lambda_p": lp, "lambda_b": lb_ctrl, "eta_mode": "fixed0"},
            "treatment": {"lambda_c": lc_trt, "lambda_p": lp, "lambda_b": lb_trt, "eta_mode": "free"},
            "split_success_rule": {
                "rmse_gain_min": 0.005,
                "canon_gain_min": -0.05,
                "nonboundary_gain_min": -0.05,
            },
        },
        "summary": summary,
        "split_reports": split_reports,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1888: SIGNED-COUPLING MULTISPLIT COMPARISON",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Splits: {summary['n_splits']} (seed {split_seeds[0]}..{split_seeds[-1]})",
        "",
        "## Global Summary (signed - unsigned)",
        f"- success_rate: {summary['success_rate']:.3f}",
        f"- rmse_gain_median: {summary['rmse_gain_median']:.4f}",
        f"- rmse_gain_q25: {summary['rmse_gain_q25']:.4f}",
        f"- canon_gain_median: {summary['canon_gain_median']:.4f}",
        f"- nonboundary_gain_median: {summary['nonboundary_gain_median']:.4f}",
        f"- nonboundary_gain_q25: {summary['nonboundary_gain_q25']:.4f}",
        "",
        "## Split Success Rule",
        "- rmse_gain >= 0.005",
        "- canon_gain >= -0.05",
        "- nonboundary_gain >= -0.05",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1888] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1888] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
