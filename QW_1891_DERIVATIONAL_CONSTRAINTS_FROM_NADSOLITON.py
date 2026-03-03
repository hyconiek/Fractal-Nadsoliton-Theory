#!/usr/bin/env python3
"""
QW-1891: Derivational constraints from nadsoliton characteristics.

Imposes explicit physics-motivated constraints:
- omega near octave resonance: omega ~= pi/4
- phi near hexagonal family: phi ~= +/- pi/6
- beta_tors small positive damping: beta in [0.005, 0.05]

Then tests constrained signed-coupling model on strict OOS split.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1891_derivational_constraints_from_nadsoliton.json"
OUT_MD = ROOT / "RAPORT_QW1891_DERIVATIONAL_CONSTRAINTS_FROM_NADSOLITON.md"


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


def nonboundary_flag(f: Dict) -> bool:
    return bool(
        0.70 < f["omega"] < 0.87
        and (abs(f["phi"] - math.pi / 6.0) < 0.22 or abs(f["phi"] + math.pi / 6.0) < 0.22)
        and 0.005 <= f["beta"] <= 0.05
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


def fit_profile_constrained(
    y: np.ndarray,
    d: np.ndarray,
    nodes: List[int],
    target: Dict[str, float],
    lambda_c: float,
    lambda_b: float,
    rng: np.random.Generator,
    tau: float,
    n_samples: int,
) -> Dict:
    def boundary_penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> float:
        pen = 0.0
        # omega around pi/4 +/- 0.12
        if not (math.pi / 4 - 0.12 <= omega <= math.pi / 4 + 0.12):
            pen += 1.0
        # phi around +/- pi/6 within 0.22
        d1 = abs(phi - math.pi / 6.0)
        d2 = abs(phi + math.pi / 6.0)
        if min(d1, d2) > 0.22:
            pen += 1.0
        # beta small positive
        if not (0.005 <= beta <= 0.05):
            pen += 1.0
        # mild interior constraints
        if not (0.35 < rho < 0.93):
            pen += 0.5
        if not (0.03 < xi < 0.58):
            pen += 0.5
        if not (0.03 < zeta < 0.55):
            pen += 0.5
        if not (-0.45 < eta < 0.45):
            pen += 0.5
        return pen / 5.0

    def eval_one(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> Dict:
        est = estimate_amp_gamma(y, d, omega, phi, beta, rho, xi, zeta, eta, nodes, tau)
        cs = canon_score(omega, phi, beta, target)
        bp = boundary_penalty(omega, phi, beta, rho, xi, zeta, eta)
        obj = est["rmse"] + lambda_c * (1.0 - cs) + lambda_b * bp
        return {
            "omega": omega,
            "phi": phi,
            "beta": beta,
            "rho": rho,
            "xi": xi,
            "zeta": zeta,
            "eta": eta,
            "rmse": est["rmse"],
            "canon_score": cs,
            "objective": obj,
        }

    phi_family = [math.pi / 6.0, -math.pi / 6.0]
    anchors = [
        (math.pi / 4.0, phi_family[0], 0.01, 0.60, 0.20, 0.20, 0.0),
        (math.pi / 4.0, phi_family[1], 0.01, 0.60, 0.20, 0.20, 0.0),
        (math.pi / 4.0 + 0.05, phi_family[0], 0.02, 0.75, 0.20, 0.20, 0.2),
    ]

    best = None
    for a in anchors:
        cand = eval_one(*a)
        if best is None or cand["objective"] < best["objective"]:
            best = cand

    for _ in range(n_samples):
        omega = float(np.clip(rng.normal(best["omega"], 0.03), math.pi / 4 - 0.14, math.pi / 4 + 0.14))
        if rng.random() < 0.5:
            phi_center = math.pi / 6.0
        else:
            phi_center = -math.pi / 6.0
        phi = float(np.clip(rng.normal(phi_center, 0.11), -1.2, 1.2))
        beta = float(np.clip(rng.normal(best["beta"], 0.010), 0.003, 0.08))
        rho = float(np.clip(rng.normal(best["rho"], 0.10), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.06), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.06), 0.01, 0.70))
        eta = float(np.clip(rng.normal(best["eta"], 0.10), -0.70, 0.70))
        cand = eval_one(omega, phi, beta, rho, xi, zeta, eta)
        if cand["objective"] < best["objective"]:
            best = cand

    return best


def fit_part(
    part: List[Dict],
    nodes: List[int],
    t: Dict[str, float],
    lambda_c: float,
    lambda_b: float,
    rng: np.random.Generator,
    tau: float,
    n_samples: int,
) -> Tuple[List[Dict], Dict]:
    rows = []
    for r in part:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        fit = fit_profile_constrained(y, d, nodes, t, lambda_c, lambda_b, rng, tau, n_samples)
        rows.append({"seed": int(r["seed"]), "fit": fit, "nonboundary": nonboundary_flag(fit)})
    return rows, summarize(rows)


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1880 = read_json("report_qw1880_node_state_strict_oos.json")
    d1884 = read_json("report_qw1884_node_state_pareto_oos_rebalancing.json")
    d1887 = read_json("report_qw1887_signed_coupling_micromodel_rebuild.json")

    rows = sorted(d1739.get("rows", []), key=lambda x: int(x["seed"]))
    split = split_dataset(rows, seed=188000)

    target = d1862.get("target_tuple", {})
    t = {
        "omega": float(target.get("omega", math.pi / 4.0)),
        "phi": float(target.get("phi", math.pi / 6.0)),
        "beta": float(target.get("beta", 0.01)),
    }

    nodes = [2, 5, 8, 11]
    tau = 2.5
    rng = np.random.default_rng(189100)

    # tune only lambda_b under locked lambda_c to isolate derivational constraint impact.
    lc = 0.2
    lb_grid = [0.10, 0.20, 0.35, 0.50]
    tuning = []
    best = None

    for lb in lb_grid:
        _, s = fit_part(split["val"], nodes, t, lambda_c=lc, lambda_b=lb, rng=rng, tau=tau, n_samples=220)
        obj = float(s["rmse_median"] + 0.20 * (1.0 - s["canon_median"]) + 0.15 * (1.0 - s["nonboundary_rate"]))
        row = {"lambda_b": lb, "val_summary": s, "val_objective": obj}
        tuning.append(row)
        if best is None or obj < best["val_objective"]:
            best = row

    lb_star = float(best["lambda_b"])

    test_rows, test_summary = fit_part(split["test"], nodes, t, lambda_c=lc, lambda_b=lb_star, rng=rng, tau=tau, n_samples=260)

    s1880 = d1880.get("summaries", {}).get("test", {})
    s1884 = d1884.get("summaries", {}).get("test", {})
    s1887 = d1887.get("summaries", {}).get("test", {})

    delta_1880 = {
        "rmse_median_gain": float(float(s1880.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1880.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1880.get("nonboundary_rate", 0.0))),
    }
    delta_1884 = {
        "rmse_median_gain": float(float(s1884.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1884.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1884.get("nonboundary_rate", 0.0))),
    }
    delta_1887 = {
        "rmse_median_gain": float(float(s1887.get("rmse_median", 0.0)) - test_summary["rmse_median"]),
        "canon_median_gain": float(test_summary["canon_median"] - float(s1887.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_summary["nonboundary_rate"] - float(s1887.get("nonboundary_rate", 0.0))),
    }

    if (
        test_summary["canon_median"] >= 0.85
        and test_summary["nonboundary_rate"] >= 0.50
        and delta_1887["rmse_median_gain"] >= -0.01
    ):
        verdict = "DERIVATIONAL_CONSTRAINTS_PARTIAL_COMPATIBLE"
    elif test_summary["canon_median"] >= 0.70 and test_summary["nonboundary_rate"] >= 0.30:
        verdict = "DERIVATIONAL_CONSTRAINTS_WEAK_COMPATIBLE"
    else:
        verdict = "DERIVATIONAL_CONSTRAINTS_INCOMPATIBLE"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "constraints": {
            "omega": "pi/4 +/- 0.12",
            "phi": "{+/- pi/6} +/- 0.22",
            "beta_tors": "[0.005, 0.05]",
            "interpretation": "nadsoliton resonance + hex symmetry + torsion damping",
        },
        "protocol": {
            "split_seed": 188000,
            "split_counts": {k: len(v) for k, v in split.items()},
            "nodes": nodes,
            "tau": tau,
            "lambda_c_locked": lc,
            "selected_lambda_b": lb_star,
        },
        "tuning": tuning,
        "summaries": {
            "test": test_summary,
            "test_1880_reference": s1880,
            "test_1884_reference": s1884,
            "test_1887_reference": s1887,
        },
        "delta_vs_1880": delta_1880,
        "delta_vs_1884": delta_1884,
        "delta_vs_1887": delta_1887,
        "test_rows": test_rows,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1891: DERIVATIONAL CONSTRAINTS FROM NADSOLITON",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- selected lambda_b: {lb_star}",
        "",
        "## Imposed Constraints",
        "- omega: pi/4 +/- 0.12",
        "- phi: {+/- pi/6} +/- 0.22",
        "- beta_tors: [0.005, 0.05]",
        "",
        "## Test Summary",
        f"- rmse median: {test_summary['rmse_median']:.4f}",
        f"- canon median: {test_summary['canon_median']:.4f}",
        f"- nonboundary rate: {test_summary['nonboundary_rate']:.3f}",
        "",
        "## Delta vs QW-1887",
        f"- rmse gain: {delta_1887['rmse_median_gain']:.4e}",
        f"- canon gain: {delta_1887['canon_median_gain']:.4e}",
        f"- nonboundary gain: {delta_1887['nonboundary_gain']:.4e}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1891] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1891] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
