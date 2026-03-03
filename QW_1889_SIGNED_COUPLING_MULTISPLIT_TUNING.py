#!/usr/bin/env python3
"""
QW-1889: Multi-split tuning of signed-coupling model with holdout evaluation.

Goal:
- Find (lambda_c, lambda_b) for signed model that preserves canonical behavior
  while retaining RMSE/nonboundary gains vs unsigned control.
- Tune on split seeds 188000..188011, evaluate on holdout 188012..188024.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1889_signed_coupling_multisplit_tuning.json"
OUT_MD = ROOT / "RAPORT_QW1889_SIGNED_COUPLING_MULTISPLIT_TUNING.md"


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
    return {
        "n": len(rows),
        "rmse_median": float(np.median(rmse)) if len(rmse) else None,
        "canon_median": float(np.median(cs)) if len(cs) else None,
        "nonboundary_rate": float(np.mean(nb)) if len(nb) else None,
    }


def estimate_train_prior(
    train_rows: List[Dict],
    nodes: List[int],
    t: Dict[str, float],
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
            rng=rng,
            tau=tau,
            n_samples=100,
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


def eval_split_signed(
    test_rows: List[Dict],
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
) -> Dict:
    rows = []
    for r in test_rows:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        f = fit_profile(
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
            tau,
            n_samples=n_samples,
        )
        rows.append({"seed": int(r["seed"]), "fit": f, "nonboundary": nonboundary_flag(f)})
    return summarize(rows)


def aggregate_vs_control(reports: List[Dict], control_map: Dict[int, Dict]) -> Dict:
    deltas = []
    for r in reports:
        s = int(r["split_seed"])
        ctrl = control_map[s]
        trt = r["treatment"]
        deltas.append(
            {
                "split_seed": s,
                "rmse_gain": float(ctrl["rmse_median"] - trt["rmse_median"]),
                "canon_gain": float(trt["canon_median"] - ctrl["canon_median"]),
                "nonboundary_gain": float(trt["nonboundary_rate"] - ctrl["nonboundary_rate"]),
            }
        )

    arr_r = np.array([x["rmse_gain"] for x in deltas], dtype=float)
    arr_c = np.array([x["canon_gain"] for x in deltas], dtype=float)
    arr_n = np.array([x["nonboundary_gain"] for x in deltas], dtype=float)

    success = np.array(
        [
            1.0 if (x["rmse_gain"] >= 0.005 and x["canon_gain"] >= -0.10 and x["nonboundary_gain"] >= 0.0) else 0.0
            for x in deltas
        ],
        dtype=float,
    )

    return {
        "deltas": deltas,
        "summary": {
            "n_splits": int(len(deltas)),
            "success_rate": float(np.mean(success)),
            "rmse_gain_median": float(np.median(arr_r)),
            "canon_gain_median": float(np.median(arr_c)),
            "nonboundary_gain_median": float(np.median(arr_n)),
            "rmse_gain_q25": float(np.quantile(arr_r, 0.25)),
            "canon_gain_q25": float(np.quantile(arr_c, 0.25)),
            "nonboundary_gain_q25": float(np.quantile(arr_n, 0.25)),
        },
    }


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1888 = read_json("report_qw1888_signed_coupling_multisplit_comparison.json")

    rows_all = sorted(d1739.get("rows", []), key=lambda x: int(x["seed"]))

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    nodes = [2, 5, 8, 11]
    tau = 2.5
    lp = 0.0

    split_seeds_tune = list(range(188000, 188012))
    split_seeds_holdout = list(range(188012, 188025))

    control_map = {
        int(r["split_seed"]): r["control_unsigned"]
        for r in d1888.get("split_reports", [])
    }

    # Prepare split-specific priors once.
    split_cache = {}
    for i, s in enumerate(split_seeds_tune + split_seeds_holdout):
        split = split_dataset(rows_all, seed=s)
        rng = np.random.default_rng(188900 + i)
        pm, ps = estimate_train_prior(split["train"], nodes, t, rng, tau=tau)
        split_cache[s] = {
            "test_rows": split["test"],
            "prior_mean": pm,
            "prior_scale": ps,
        }

    lambda_c_grid = [0.20, 0.30, 0.40, 0.50]
    lambda_b_grid = [0.10, 0.20, 0.35]

    tuning_rows = []
    feasible = []

    for lc in lambda_c_grid:
        for lb in lambda_b_grid:
            reports = []
            for j, s in enumerate(split_seeds_tune):
                cache = split_cache[s]
                rng = np.random.default_rng(199900 + 100 * j + int(100 * lc) + int(100 * lb))
                trt = eval_split_signed(
                    cache["test_rows"],
                    nodes,
                    cache["prior_mean"],
                    cache["prior_scale"],
                    t,
                    lambda_c=lc,
                    lambda_p=lp,
                    lambda_b=lb,
                    rng=rng,
                    tau=tau,
                    n_samples=170,
                )
                reports.append({"split_seed": s, "treatment": trt})

            agg = aggregate_vs_control(reports, control_map)
            ssum = agg["summary"]

            # Preregistered feasibility for tuning phase.
            ok = bool(
                ssum["rmse_gain_median"] >= 0.005
                and ssum["canon_gain_median"] >= -0.10
                and ssum["nonboundary_gain_median"] >= 0.0
            )

            # maximize objective (bigger is better).
            obj = float(
                ssum["rmse_gain_median"]
                + 0.20 * ssum["nonboundary_gain_median"]
                + 0.15 * ssum["canon_gain_median"]
            )

            row = {
                "lambda_c": lc,
                "lambda_b": lb,
                "tune_summary": ssum,
                "tune_objective": obj,
                "feasible": ok,
            }
            tuning_rows.append(row)
            if ok:
                feasible.append(row)

    if feasible:
        selected = max(feasible, key=lambda r: r["tune_objective"])
        selection_mode = "feasible_max_objective"
    else:
        selected = max(tuning_rows, key=lambda r: r["tune_objective"])
        selection_mode = "no_feasible_fallback_max_objective"

    lc_star = float(selected["lambda_c"])
    lb_star = float(selected["lambda_b"])

    # Holdout evaluation on unseen split seeds.
    hold_reports = []
    for j, s in enumerate(split_seeds_holdout):
        cache = split_cache[s]
        rng = np.random.default_rng(288900 + 100 * j + int(100 * lc_star) + int(100 * lb_star))
        trt = eval_split_signed(
            cache["test_rows"],
            nodes,
            cache["prior_mean"],
            cache["prior_scale"],
            t,
            lambda_c=lc_star,
            lambda_p=lp,
            lambda_b=lb_star,
            rng=rng,
            tau=tau,
            n_samples=220,
        )
        hold_reports.append({"split_seed": s, "treatment": trt})

    hold_agg = aggregate_vs_control(hold_reports, control_map)
    hold_sum = hold_agg["summary"]

    if (
        hold_sum["success_rate"] >= 0.50
        and hold_sum["rmse_gain_median"] >= 0.005
        and hold_sum["canon_gain_median"] >= -0.10
    ):
        verdict = "SIGNED_COUPLING_TUNED_HOLDOUT_PARTIAL_SUCCESS"
    elif hold_sum["success_rate"] >= 0.30 and hold_sum["rmse_gain_median"] > 0.0:
        verdict = "SIGNED_COUPLING_TUNED_HOLDOUT_WEAK"
    else:
        verdict = "SIGNED_COUPLING_TUNED_HOLDOUT_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "split_seeds_tune": split_seeds_tune,
            "split_seeds_holdout": split_seeds_holdout,
            "nodes": nodes,
            "tau": tau,
            "lambda_p": lp,
            "selection_mode": selection_mode,
            "success_rule": {
                "rmse_gain_min": 0.005,
                "canon_gain_min": -0.10,
                "nonboundary_gain_min": 0.0,
            },
        },
        "tuning_rows": tuning_rows,
        "selected_hyperparams": {"lambda_c": lc_star, "lambda_b": lb_star},
        "selected_tune_row": selected,
        "holdout_summary": hold_sum,
        "holdout_deltas": hold_agg["deltas"],
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1889: SIGNED-COUPLING MULTISPLIT TUNING",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- selection_mode: `{selection_mode}`",
        f"- selected lambda_c/lambda_b: `{lc_star}` / `{lb_star}`",
        "",
        "## Holdout Summary (signed - unsigned)",
        f"- success_rate: {hold_sum['success_rate']:.3f}",
        f"- rmse_gain_median: {hold_sum['rmse_gain_median']:.4f}",
        f"- canon_gain_median: {hold_sum['canon_gain_median']:.4f}",
        f"- nonboundary_gain_median: {hold_sum['nonboundary_gain_median']:.4f}",
        f"- rmse_gain_q25: {hold_sum['rmse_gain_q25']:.4f}",
        f"- canon_gain_q25: {hold_sum['canon_gain_q25']:.4f}",
        f"- nonboundary_gain_q25: {hold_sum['nonboundary_gain_q25']:.4f}",
        "",
        "## Success Rule",
        "- rmse_gain >= 0.005",
        "- canon_gain >= -0.10",
        "- nonboundary_gain >= 0.0",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1889] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1889] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
