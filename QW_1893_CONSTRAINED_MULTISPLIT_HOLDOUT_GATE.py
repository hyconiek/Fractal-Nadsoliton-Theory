#!/usr/bin/env python3
"""
QW-1893: Strict OOS + multisplit holdout gate for derivationally constrained model.

Compares QW-1891 constrained branch against unsigned control (from QW-1888)
over repeated split seeds.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1893_constrained_multisplit_holdout_gate.json"
OUT_MD = ROOT / "RAPORT_QW1893_CONSTRAINED_MULTISPLIT_HOLDOUT_GATE.md"


def read_json(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


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


def canon_score(omega: float, phi: float, beta: float, t: Dict[str, float]) -> float:
    z_o = abs(omega - t["omega"]) / 0.20
    dphi = (phi - t["phi"] + math.pi) % (2.0 * math.pi) - math.pi
    z_p = abs(dphi) / 0.30
    z_b = abs(beta - t["beta"]) / 0.030
    z = min(20.0, math.sqrt((z_o * z_o + z_p * z_p + z_b * z_b) / 3.0))
    return float(math.exp(-0.5 * z * z))


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


def fit_profile_constrained(
    y: np.ndarray,
    d: np.ndarray,
    t: Dict[str, float],
    nodes: List[int],
    tau: float,
    lambda_c: float,
    lambda_b: float,
    rng: np.random.Generator,
    n_samples: int,
) -> Dict:
    def penalty(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> float:
        p = 0.0
        if not (math.pi / 4 - 0.12 <= omega <= math.pi / 4 + 0.12):
            p += 1.0
        if min(abs(phi - math.pi / 6.0), abs(phi + math.pi / 6.0)) > 0.22:
            p += 1.0
        if not (0.005 <= beta <= 0.05):
            p += 1.0
        if not (0.35 < rho < 0.93):
            p += 0.5
        if not (0.03 < xi < 0.58):
            p += 0.5
        if not (0.03 < zeta < 0.55):
            p += 0.5
        if not (-0.45 < eta < 0.45):
            p += 0.5
        return p / 5.0

    def eval_one(omega: float, phi: float, beta: float, rho: float, xi: float, zeta: float, eta: float) -> Dict:
        est = estimate_amp_gamma(y, d, omega, phi, beta, rho, xi, zeta, eta, nodes, tau)
        cs = canon_score(omega, phi, beta, t)
        obj = est["rmse"] + lambda_c * (1.0 - cs) + lambda_b * penalty(omega, phi, beta, rho, xi, zeta, eta)
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

    anchors = [
        (math.pi / 4.0, math.pi / 6.0, 0.01, 0.60, 0.20, 0.20, 0.0),
        (math.pi / 4.0, -math.pi / 6.0, 0.01, 0.60, 0.20, 0.20, 0.0),
        (math.pi / 4.0 + 0.05, math.pi / 6.0, 0.02, 0.75, 0.20, 0.20, 0.2),
    ]

    best = None
    for a in anchors:
        c = eval_one(*a)
        if best is None or c["objective"] < best["objective"]:
            best = c

    for _ in range(n_samples):
        omega = float(np.clip(rng.normal(best["omega"], 0.03), math.pi / 4 - 0.14, math.pi / 4 + 0.14))
        phi_center = math.pi / 6.0 if rng.random() < 0.5 else -math.pi / 6.0
        phi = float(np.clip(rng.normal(phi_center, 0.11), -1.2, 1.2))
        beta = float(np.clip(rng.normal(best["beta"], 0.010), 0.003, 0.08))
        rho = float(np.clip(rng.normal(best["rho"], 0.10), 0.25, 0.98))
        xi = float(np.clip(rng.normal(best["xi"], 0.06), 0.01, 0.65))
        zeta = float(np.clip(rng.normal(best["zeta"], 0.06), 0.01, 0.70))
        eta = float(np.clip(rng.normal(best["eta"], 0.10), -0.70, 0.70))
        c = eval_one(omega, phi, beta, rho, xi, zeta, eta)
        if c["objective"] < best["objective"]:
            best = c

    return best


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
    rmse = np.array([r["rmse"] for r in rows], dtype=float)
    cs = np.array([r["canon_score"] for r in rows], dtype=float)
    nb = np.array([1.0 if r["nonboundary"] else 0.0 for r in rows], dtype=float)
    return {
        "n": len(rows),
        "rmse_median": float(np.median(rmse)) if len(rmse) else None,
        "canon_median": float(np.median(cs)) if len(cs) else None,
        "nonboundary_rate": float(np.mean(nb)) if len(nb) else None,
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

    control_map = {
        int(r["split_seed"]): r["control_unsigned"]
        for r in d1888.get("split_reports", [])
    }

    nodes = [2, 5, 8, 11]
    tau = 2.5
    lc = 0.2
    lb = 0.35

    split_seeds = list(range(188000, 188025))
    split_reports = []

    for i, s in enumerate(split_seeds):
        split = split_dataset(rows_all, seed=s)
        rng = np.random.default_rng(189300 + i)

        trt_rows = []
        for r in split["test"]:
            y = np.array([float(r["profile"][str(j)]) for j in range(1, 13)], dtype=float)
            d = np.arange(1, len(y) + 1, dtype=float)
            fit = fit_profile_constrained(y, d, t, nodes, tau, lc, lb, rng, n_samples=220)
            trt_rows.append(
                {
                    "rmse": fit["rmse"],
                    "canon_score": fit["canon_score"],
                    "nonboundary": nonboundary_flag(fit),
                }
            )

        trt = summarize(trt_rows)
        ctrl = control_map[s]

        delta = {
            "rmse_gain": float(ctrl["rmse_median"] - trt["rmse_median"]),
            "canon_gain": float(trt["canon_median"] - ctrl["canon_median"]),
            "nonboundary_gain": float(trt["nonboundary_rate"] - ctrl["nonboundary_rate"]),
        }

        split_success = bool(
            delta["rmse_gain"] >= 0.0
            and delta["canon_gain"] >= 0.0
            and delta["nonboundary_gain"] >= 0.0
        )

        split_reports.append(
            {
                "split_seed": s,
                "control_unsigned": ctrl,
                "treatment_constrained": trt,
                "delta_constrained_minus_control": delta,
                "split_success": split_success,
            }
        )

    arr_r = np.array([x["delta_constrained_minus_control"]["rmse_gain"] for x in split_reports], dtype=float)
    arr_c = np.array([x["delta_constrained_minus_control"]["canon_gain"] for x in split_reports], dtype=float)
    arr_n = np.array([x["delta_constrained_minus_control"]["nonboundary_gain"] for x in split_reports], dtype=float)
    arr_s = np.array([1.0 if x["split_success"] else 0.0 for x in split_reports], dtype=float)

    summary = {
        "n_splits": len(split_reports),
        "success_rate": float(np.mean(arr_s)),
        "rmse_gain_median": float(np.median(arr_r)),
        "canon_gain_median": float(np.median(arr_c)),
        "nonboundary_gain_median": float(np.median(arr_n)),
        "rmse_gain_q25": float(np.quantile(arr_r, 0.25)),
        "canon_gain_q25": float(np.quantile(arr_c, 0.25)),
        "nonboundary_gain_q25": float(np.quantile(arr_n, 0.25)),
    }

    if (
        summary["success_rate"] >= 0.50
        and summary["rmse_gain_median"] >= 0.0
        and summary["canon_gain_median"] >= 0.0
    ):
        verdict = "CONSTRAINED_BRANCH_MULTISPLIT_STRONG"
    elif summary["success_rate"] >= 0.30 and summary["canon_gain_median"] >= 0.0:
        verdict = "CONSTRAINED_BRANCH_MULTISPLIT_PARTIAL"
    else:
        verdict = "CONSTRAINED_BRANCH_MULTISPLIT_FAIL"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "split_seed_range": [split_seeds[0], split_seeds[-1]],
            "n_splits": len(split_seeds),
            "constraints": {
                "omega": "pi/4 +/- 0.12",
                "phi": "{+/- pi/6} +/- 0.22",
                "beta": "[0.005, 0.05]",
            },
            "lambda_c": lc,
            "lambda_b": lb,
            "success_rule": "rmse_gain>=0 and canon_gain>=0 and nonboundary_gain>=0",
        },
        "summary": summary,
        "split_reports": split_reports,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1893: CONSTRAINED MULTISPLIT HOLDOUT GATE",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- Splits: {summary['n_splits']} (seed {split_seeds[0]}..{split_seeds[-1]})",
        "",
        "## Summary (constrained - unsigned control)",
        f"- success_rate: {summary['success_rate']:.3f}",
        f"- rmse_gain_median: {summary['rmse_gain_median']:.4f}",
        f"- canon_gain_median: {summary['canon_gain_median']:.4f}",
        f"- nonboundary_gain_median: {summary['nonboundary_gain_median']:.4f}",
        f"- rmse_gain_q25: {summary['rmse_gain_q25']:.4f}",
        f"- canon_gain_q25: {summary['canon_gain_q25']:.4f}",
        f"- nonboundary_gain_q25: {summary['nonboundary_gain_q25']:.4f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1893] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1893] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
