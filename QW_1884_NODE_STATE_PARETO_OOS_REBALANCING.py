#!/usr/bin/env python3
"""
QW-1884: Pareto OOS rebalancing for node-state model.

Goal:
- Keep strict split/protocol from QW-1880.
- Tune lambda_c and lambda_b on validation only.
- Enforce preregistered feasibility constraints to avoid arbitrary tradeoff selection.
- Run one locked test evaluation and compare with QW-1880 and QW-1883.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1884_node_state_pareto_oos_rebalancing.json"
OUT_MD = ROOT / "RAPORT_QW1884_NODE_STATE_PARETO_OOS_REBALANCING.md"


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


def simulate_node_state(
    d: np.ndarray,
    omega: float,
    phi: float,
    rho: float,
    xi: float,
    zeta: float,
    nodes: List[int],
) -> np.ndarray:
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


def estimate_amp_gamma(
    y: np.ndarray,
    d: np.ndarray,
    omega: float,
    phi: float,
    beta: float,
    rho: float,
    xi: float,
    zeta: float,
    nodes: List[int],
) -> Dict[str, float]:
    env_basis = np.cos(omega * d + phi) / (1.0 + beta * d)
    st_basis = simulate_node_state(d, omega, phi, rho, xi, zeta, nodes)
    X = np.column_stack([env_basis, st_basis])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    yhat = X @ coef
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))
    return {"amp": float(coef[0]), "gamma": float(coef[1]), "rmse": rmse}


def boundary_penalty(omega: float, beta: float, rho: float, xi: float, zeta: float) -> float:
    prefs = [
        (omega, 0.08, 0.95, 0.12, 0.90),
        (beta, 0.003, 0.30, 0.008, 0.25),
        (rho, 0.25, 0.98, 0.35, 0.93),
        (xi, 0.01, 0.65, 0.03, 0.58),
        (zeta, 0.01, 0.70, 0.03, 0.55),
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
    n_samples: int,
) -> Dict:
    def prior_penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float) -> float:
        z1 = (omega - prior_mean["omega"]) / prior_scale["omega"]
        z2 = circular_diff(phi, prior_mean["phi"]) / prior_scale["phi"]
        z3 = (beta - prior_mean["beta"]) / prior_scale["beta"]
        z4 = (rho - prior_mean["rho"]) / prior_scale["rho"]
        z5 = (xi - prior_mean["xi"]) / prior_scale["xi"]
        z6 = (zeta - prior_mean["zeta"]) / prior_scale["zeta"]
        return float((z1 * z1 + z2 * z2 + z3 * z3 + z4 * z4 + z5 * z5 + z6 * z6) / 6.0)

    def eval_one(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float) -> Dict:
        est = estimate_amp_gamma(y, d, omega, phi, beta, rho, xi, zeta, nodes)
        cscore = canon_score(omega, phi, beta, target)
        pp = prior_penalty(omega, phi, beta, rho, xi, zeta)
        bp = boundary_penalty(omega, beta, rho, xi, zeta)
        obj = est["rmse"] + lambda_c * (1.0 - cscore) + lambda_p * pp + lambda_b * bp
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
            "prior_penalty": pp,
            "boundary_penalty": bp,
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
        ),
        (target["omega"], target["phi"], target["beta"], 0.60, 0.20, 0.20),
        (0.70, prior_mean["phi"], 0.03, 0.75, 0.20, 0.20),
    ]

    best = None
    for a in anchors:
        cand = eval_one(*a)
        if best is None or cand["objective"] < best["objective"]:
            best = cand

    for _ in range(n_samples):
        # Keep search dynamics aligned with QW-1883 to avoid optimizer confound.
        omega = float(np.clip(rng.normal(best["omega"], 0.045), 0.08, 0.95))
        phi = float(np.clip(rng.normal(best["phi"], 0.16), -1.4, 1.4))
        beta = float(np.clip(rng.normal(best["beta"], 0.018), 0.003, 0.30))
        rho = float(np.clip(rng.normal(best["rho"], 0.10), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.06), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.06), 0.01, 0.70))

        cand = eval_one(omega, phi, beta, rho, xi, zeta)
        if cand["objective"] < best["objective"]:
            best = cand

    return best


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


def nonboundary_flag(f: Dict) -> bool:
    return bool(
        0.12 < f["omega"] < 0.90
        and 0.008 < f["beta"] < 0.25
        and 0.35 < f["rho"] < 0.93
        and 0.03 < f["xi"] < 0.58
    )


def make_locked_split(rows_all: List[Dict]) -> Dict[str, List[Dict]]:
    by_n = {}
    for r in rows_all:
        by_n.setdefault(int(r["n_nodes"]), []).append(r)

    rng_split = np.random.default_rng(188000)
    split = {"train": [], "val": [], "test": []}
    for n in sorted(by_n.keys()):
        g = list(by_n[n])
        idx = np.arange(len(g))
        rng_split.shuffle(idx)
        g = [g[i] for i in idx]
        n_train = min(8, len(g) - 4)
        n_val = min(3, len(g) - n_train - 1)
        split["train"].extend(g[:n_train])
        split["val"].extend(g[n_train : n_train + n_val])
        split["test"].extend(g[n_train + n_val :])

    return split


def eval_split(
    split_part: List[Dict],
    nodes: List[int],
    prior_mean: Dict[str, float],
    prior_scale: Dict[str, float],
    t: Dict[str, float],
    lambda_c: float,
    lambda_p: float,
    lambda_b: float,
    rng: np.random.Generator,
    n_samples: int,
) -> Tuple[List[Dict], Dict]:
    rows = []
    for r in split_part:
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
            rng=rng,
            n_samples=n_samples,
        )
        rows.append({"seed": int(r["seed"]), "fit": fit, "nonboundary": nonboundary_flag(fit)})

    return rows, summarize(rows)


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1880 = read_json("report_qw1880_node_state_strict_oos.json")
    d1883 = read_json("report_qw1883_node_state_boundary_aware_oos.json")

    rows_all = sorted(d1739.get("rows", []), key=lambda x: int(x["seed"]))
    split = make_locked_split(rows_all)

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    prior_mean = d1880.get("prior_from_train", {}).get("mean", {})
    prior_scale = d1880.get("prior_from_train", {}).get("scale", {})
    protocol = d1880.get("protocol", {})
    nodes = list(protocol.get("nodes_used", [2, 5, 8, 11]))
    lp = float(d1880.get("best_hyperparams", {}).get("lambda_p", 0.0))
    lc_ref = float(d1880.get("best_hyperparams", {}).get("lambda_c", 0.2))

    rng = np.random.default_rng(188400)

    lambda_c_grid = [0.02, 0.05, 0.10, 0.20, 0.30]
    lambda_b_grid = [0.00, 0.05, 0.10, 0.20, 0.35, 0.50]

    # Baseline: 1880 setting without boundary penalty.
    _, val_baseline = eval_split(
        split["val"],
        nodes,
        prior_mean,
        prior_scale,
        t,
        lambda_c=lc_ref,
        lambda_p=lp,
        lambda_b=0.0,
        rng=rng,
        n_samples=220,
    )
    val_rmse_baseline = float(val_baseline["rmse_median"])

    tuning = []
    feasible = []

    for lc in lambda_c_grid:
        for lb in lambda_b_grid:
            _, val_s = eval_split(
                split["val"],
                nodes,
                prior_mean,
                prior_scale,
                t,
                lambda_c=lc,
                lambda_p=lp,
                lambda_b=lb,
                rng=rng,
                n_samples=220,
            )

            rmse_ratio = float(val_s["rmse_median"] / max(val_rmse_baseline, 1e-12))
            obj = float(
                val_s["rmse_median"]
                + 0.18 * (1.0 - val_s["canon_median"])
                + 0.18 * (1.0 - val_s["nonboundary_rate"])
            )
            is_feasible = bool(
                val_s["nonboundary_rate"] >= 0.50
                and val_s["canon_median"] >= 0.80
                and rmse_ratio <= 1.35
            )

            row = {
                "lambda_c": lc,
                "lambda_b": lb,
                "val_summary": val_s,
                "val_rmse_ratio_vs_1880baseline": rmse_ratio,
                "val_objective": obj,
                "feasible": is_feasible,
            }
            tuning.append(row)
            if is_feasible:
                feasible.append(row)

    if feasible:
        selected = min(feasible, key=lambda r: r["val_objective"])
        selection_mode = "feasible_min_objective"
    else:
        # fallback only if no feasible point exists.
        selected = min(tuning, key=lambda r: r["val_objective"])
        selection_mode = "no_feasible_fallback_min_objective"

    lc_star = float(selected["lambda_c"])
    lb_star = float(selected["lambda_b"])

    test_rows, test_summary = eval_split(
        split["test"],
        nodes,
        prior_mean,
        prior_scale,
        t,
        lambda_c=lc_star,
        lambda_p=lp,
        lambda_b=lb_star,
        rng=rng,
        n_samples=260,
    )

    s1880 = d1880.get("summaries", {}).get("test", {})
    s1883 = d1883.get("summaries", {}).get("test", {})

    rmse_ratio_1880 = float(test_summary["rmse_median"] / max(float(s1880.get("rmse_median", 1.0)), 1e-12))

    delta_1880 = {
        "rmse_median_gain": float(float(s1880.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1880.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1880.get("nonboundary_rate", 0.0))),
        "rmse_ratio": rmse_ratio_1880,
    }

    delta_1883 = {
        "rmse_median_gain": float(float(s1883.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1883.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1883.get("nonboundary_rate", 0.0))),
    }

    if (
        test_summary["nonboundary_rate"] >= 0.50
        and test_summary["canon_median"] >= 0.80
        and rmse_ratio_1880 <= 1.35
    ):
        verdict = "PARETO_OOS_REBALANCING_SUCCESS"
    elif (
        test_summary["nonboundary_rate"] >= 0.50
        and test_summary["canon_median"] >= 0.80
        and rmse_ratio_1880 <= 1.60
    ):
        verdict = "PARETO_OOS_REBALANCING_PARTIAL_TRADEOFF"
    elif test_summary["nonboundary_rate"] > 0.0:
        verdict = "PARETO_OOS_REBALANCING_WEAK_PROGRESS"
    else:
        verdict = "PARETO_OOS_REBALANCING_FAIL"

    tuning_sorted = sorted(tuning, key=lambda r: (not r["feasible"], r["val_objective"]))

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol_locked_from_1880": {
            "split_counts": {k: len(v) for k, v in split.items()},
            "nodes": nodes,
            "lambda_p_locked": lp,
            "baseline_lambda_c": lc_ref,
        },
        "preregistered_selection": {
            "feasibility": {
                "val_nonboundary_min": 0.50,
                "val_canon_min": 0.80,
                "val_rmse_ratio_max_vs_1880baseline": 1.35,
            },
            "objective": "rmse_median + 0.18*(1-canon_median) + 0.18*(1-nonboundary_rate)",
            "selection_mode": selection_mode,
        },
        "val_baseline_1880style": val_baseline,
        "tuning_rows": tuning,
        "selected_hyperparams": {
            "lambda_c": lc_star,
            "lambda_b": lb_star,
        },
        "selected_val_row": selected,
        "summaries": {
            "test": test_summary,
            "test_1880_reference": s1880,
            "test_1883_reference": s1883,
        },
        "delta_vs_1880": delta_1880,
        "delta_vs_1883": delta_1883,
        "test_rows": test_rows,
        "verdict": verdict,
        "top5_val_rows": tuning_sorted[:5],
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1884: NODE-STATE PARETO OOS REBALANCING",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- selection_mode: `{selection_mode}`",
        f"- selected lambda_c/lambda_b: `{lc_star}` / `{lb_star}`",
        "",
        "## Val Feasibility Rule",
        "- nonboundary_rate >= 0.50",
        "- canon_median >= 0.80",
        "- rmse_ratio_vs_1880baseline <= 1.35",
        "",
        "## Test Summary",
        f"- rmse median: {test_summary['rmse_median']:.4f}",
        f"- canon median: {test_summary['canon_median']:.4f}",
        f"- nonboundary rate: {test_summary['nonboundary_rate']:.3f}",
        "",
        "## Delta vs QW-1880",
        f"- rmse gain: {delta_1880['rmse_median_gain']:.4e}",
        f"- canon gain: {delta_1880['canon_median_gain']:.4e}",
        f"- nonboundary gain: {delta_1880['nonboundary_gain']:.4e}",
        f"- rmse ratio: {delta_1880['rmse_ratio']:.4f}",
        "",
        "## Delta vs QW-1883",
        f"- rmse gain: {delta_1883['rmse_median_gain']:.4e}",
        f"- canon gain: {delta_1883['canon_median_gain']:.4e}",
        f"- nonboundary gain: {delta_1883['nonboundary_gain']:.4e}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1884] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1884] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
