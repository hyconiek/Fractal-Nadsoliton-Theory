#!/usr/bin/env python3
"""
QW-1894: Amplitude-channel reformulation under derivational constraints.

Adds a separate linear amplitude baseline channel to constrained signed model:
  y_hat = b0 + b1*d + b2*(1/d) + A*env + G*state

Goal: recover RMSE while keeping constrained (omega, phi, beta_tors) regime.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1894_amplitude_channel_reformulation.json"
OUT_MD = ROOT / "RAPORT_QW1894_AMPLITUDE_CHANNEL_REFORMULATION.md"


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


def split_dataset(rows: List[Dict], seed: int = 188000) -> Dict[str, List[Dict]]:
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


def estimate_linear_components(
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
    lambda_amp: float,
) -> Dict:
    env_basis = np.cos(omega * d + phi) / (1.0 + beta * d)
    st_basis = simulate_state(d, omega, phi, rho, xi, zeta, eta, nodes, tau)

    # Amplitude baseline channel (independent of phase/torsion tuple).
    X = np.column_stack([
        np.ones_like(d),
        d,
        1.0 / d,
        env_basis,
        st_basis,
    ])

    # Ridge on baseline coefficients only (first 3 cols).
    R = np.diag([lambda_amp, lambda_amp, lambda_amp, 1e-9, 1e-9])
    XtX = X.T @ X + R
    Xty = X.T @ y
    coef = np.linalg.solve(XtX, Xty)
    yhat = X @ coef
    rmse = float(np.sqrt(np.mean((y - yhat) ** 2)))

    return {
        "coef": coef.tolist(),
        "rmse": rmse,
    }


def fit_profile(
    y: np.ndarray,
    d: np.ndarray,
    t: Dict[str, float],
    nodes: List[int],
    tau: float,
    lambda_c: float,
    lambda_b: float,
    lambda_amp: float,
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
        est = estimate_linear_components(
            y, d, omega, phi, beta, rho, xi, zeta, eta, nodes, tau, lambda_amp
        )
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


def fit_part(
    part: List[Dict],
    t: Dict[str, float],
    nodes: List[int],
    tau: float,
    lambda_c: float,
    lambda_b: float,
    lambda_amp: float,
    rng: np.random.Generator,
    n_samples: int,
) -> Dict:
    rows = []
    for r in part:
        y = np.array([float(r["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)
        f = fit_profile(y, d, t, nodes, tau, lambda_c, lambda_b, lambda_amp, rng, n_samples)
        rows.append(f)

    rmse = np.array([x["rmse"] for x in rows], dtype=float)
    cs = np.array([x["canon_score"] for x in rows], dtype=float)
    nb = np.array([1.0 if nonboundary_flag(x) else 0.0 for x in rows], dtype=float)

    return {
        "rows": rows,
        "summary": {
            "n": len(rows),
            "rmse_median": float(np.median(rmse)),
            "canon_median": float(np.median(cs)),
            "nonboundary_rate": float(np.mean(nb)),
        },
    }


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1862 = read_json("report_qw1862_micromodel_canonical_compatibility_audit.json")
    d1891 = read_json("report_qw1891_derivational_constraints_from_nadsoliton.json")
    d1887 = read_json("report_qw1887_signed_coupling_micromodel_rebuild.json")
    d1884 = read_json("report_qw1884_node_state_pareto_oos_rebalancing.json")

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
    lc = 0.2

    rng = np.random.default_rng(189400)

    lb_grid = [0.20, 0.35, 0.50]
    la_grid = [0.01, 0.05, 0.10, 0.20]

    tuning = []
    best = None

    for lb in lb_grid:
        for la in la_grid:
            res = fit_part(split["val"], t, nodes, tau, lc, lb, la, rng, n_samples=220)
            s = res["summary"]
            obj = float(s["rmse_median"] + 0.20 * (1.0 - s["canon_median"]) + 0.15 * (1.0 - s["nonboundary_rate"]))
            row = {"lambda_b": lb, "lambda_amp": la, "val_summary": s, "val_objective": obj}
            tuning.append(row)
            if best is None or obj < best["val_objective"]:
                best = row

    lb_star = float(best["lambda_b"])
    la_star = float(best["lambda_amp"])

    test_res = fit_part(split["test"], t, nodes, tau, lc, lb_star, la_star, rng, n_samples=260)
    test_s = test_res["summary"]

    s1891 = d1891.get("summaries", {}).get("test", {})
    s1887 = d1887.get("summaries", {}).get("test", {})
    s1884 = d1884.get("summaries", {}).get("test", {})

    delta_1891 = {
        "rmse_median_gain": float(float(s1891.get("rmse_median", 0.0)) - test_s["rmse_median"]),
        "canon_median_gain": float(test_s["canon_median"] - float(s1891.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_s["nonboundary_rate"] - float(s1891.get("nonboundary_rate", 0.0))),
    }
    delta_1887 = {
        "rmse_median_gain": float(float(s1887.get("rmse_median", 0.0)) - test_s["rmse_median"]),
        "canon_median_gain": float(test_s["canon_median"] - float(s1887.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_s["nonboundary_rate"] - float(s1887.get("nonboundary_rate", 0.0))),
    }
    delta_1884 = {
        "rmse_median_gain": float(float(s1884.get("rmse_median", 0.0)) - test_s["rmse_median"]),
        "canon_median_gain": float(test_s["canon_median"] - float(s1884.get("canon_median", 0.0))),
        "nonboundary_gain": float(test_s["nonboundary_rate"] - float(s1884.get("nonboundary_rate", 0.0))),
    }

    if (
        delta_1891["rmse_median_gain"] >= 0.01
        and test_s["canon_median"] >= 0.85
        and test_s["nonboundary_rate"] >= 0.50
    ):
        verdict = "AMPLITUDE_REFORMULATION_PROMISING"
    elif delta_1891["rmse_median_gain"] > 0.0:
        verdict = "AMPLITUDE_REFORMULATION_WEAK"
    else:
        verdict = "AMPLITUDE_REFORMULATION_NO_RECOVERY"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "split_seed": 188000,
            "split_counts": {k: len(v) for k, v in split.items()},
            "constraints": {
                "omega": "pi/4 +/- 0.12",
                "phi": "{+/- pi/6} +/- 0.22",
                "beta": "[0.005, 0.05]",
            },
            "lambda_c": lc,
            "selected_lambda_b": lb_star,
            "selected_lambda_amp": la_star,
        },
        "tuning_rows": tuning,
        "summaries": {
            "test": test_s,
            "test_1891_reference": s1891,
            "test_1887_reference": s1887,
            "test_1884_reference": s1884,
        },
        "delta_vs_1891": delta_1891,
        "delta_vs_1887": delta_1887,
        "delta_vs_1884": delta_1884,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1894: AMPLITUDE CHANNEL REFORMULATION",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        f"- selected lambda_b/lambda_amp: {lb_star}/{la_star}",
        "",
        "## Test Summary",
        f"- rmse median: {test_s['rmse_median']:.4f}",
        f"- canon median: {test_s['canon_median']:.4f}",
        f"- nonboundary rate: {test_s['nonboundary_rate']:.3f}",
        "",
        "## Delta vs QW-1891",
        f"- rmse gain: {delta_1891['rmse_median_gain']:.4e}",
        f"- canon gain: {delta_1891['canon_median_gain']:.4e}",
        f"- nonboundary gain: {delta_1891['nonboundary_gain']:.4e}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1894] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1894] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
