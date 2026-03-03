#!/usr/bin/env python3
"""
QW-1892: Identifiability-first audit for constrained nadsoliton branch.

Compares local Jacobian/Fisher proxies between:
- QW-1887 signed-coupling (less constrained)
- QW-1891 constrained-by-derivation branch
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1892_identifiability_first_audit.json"
OUT_MD = ROOT / "RAPORT_QW1892_IDENTIFIABILITY_FIRST_AUDIT.md"


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


def predict(
    y: np.ndarray,
    d: np.ndarray,
    p: Dict[str, float],
    nodes: List[int],
    tau: float,
) -> np.ndarray:
    omega = float(p["omega"])
    phi = float(p["phi"])
    beta = float(p["beta"])
    rho = float(p["rho"])
    xi = float(p["xi"])
    zeta = float(p["zeta"])
    eta = float(p["eta"])

    env_basis = np.cos(omega * d + phi) / (1.0 + beta * d)
    st_basis = simulate_state(d, omega, phi, rho, xi, zeta, eta, nodes, tau)
    X = np.column_stack([env_basis, st_basis])
    coef, *_ = np.linalg.lstsq(X, y, rcond=None)
    return X @ coef


def local_jacobian(
    y: np.ndarray,
    d: np.ndarray,
    p0: Dict[str, float],
    nodes: List[int],
    tau: float,
) -> Dict:
    params = ["omega", "phi", "beta", "rho", "xi", "zeta", "eta"]
    steps = {
        "omega": 0.005,
        "phi": 0.010,
        "beta": 0.003,
        "rho": 0.008,
        "xi": 0.006,
        "zeta": 0.006,
        "eta": 0.008,
    }

    y0 = predict(y, d, p0, nodes, tau)
    J = np.zeros((len(y0), len(params)), dtype=float)

    for j, k in enumerate(params):
        h = steps[k]
        p_plus = dict(p0)
        p_minus = dict(p0)
        p_plus[k] = float(p_plus[k] + h)
        p_minus[k] = float(p_minus[k] - h)
        y_plus = predict(y, d, p_plus, nodes, tau)
        y_minus = predict(y, d, p_minus, nodes, tau)
        J[:, j] = (y_plus - y_minus) / (2.0 * h)

    F = J.T @ J + 1e-8 * np.eye(J.shape[1])
    eigvals = np.linalg.eigvalsh(F)

    cond = float(np.linalg.cond(F))
    min_eig = float(np.min(eigvals))
    max_eig = float(np.max(eigvals))

    col_w = J[:, 0]
    col_b = J[:, 2]
    if np.std(col_w) < 1e-12 or np.std(col_b) < 1e-12:
        corr_wb = 1.0
    else:
        corr_wb = float(abs(np.corrcoef(col_w, col_b)[0, 1]))

    return {
        "cond": cond,
        "log10_cond": float(np.log10(max(cond, 1.0))),
        "min_eig": min_eig,
        "max_eig": max_eig,
        "corr_abs_omega_beta": corr_wb,
    }


def summarize(rows: List[Dict]) -> Dict:
    arr_cond = np.array([r["diag"]["log10_cond"] for r in rows], dtype=float)
    arr_min = np.array([r["diag"]["min_eig"] for r in rows], dtype=float)
    arr_corr = np.array([r["diag"]["corr_abs_omega_beta"] for r in rows], dtype=float)
    return {
        "n": len(rows),
        "log10_cond_median": float(np.median(arr_cond)),
        "log10_cond_q75": float(np.quantile(arr_cond, 0.75)),
        "min_eig_median": float(np.median(arr_min)),
        "corr_abs_omega_beta_median": float(np.median(arr_corr)),
    }


def evaluate_branch(test_rows: List[Dict], profile_map: Dict[int, Dict], nodes: List[int], tau: float) -> Dict:
    out_rows = []
    for row in test_rows:
        seed = int(row["seed"])
        prof = profile_map[seed]
        y = np.array([float(prof["profile"][str(i)]) for i in range(1, 13)], dtype=float)
        d = np.arange(1, len(y) + 1, dtype=float)

        fit = row["fit"]
        p = {
            "omega": float(fit["omega"]),
            "phi": float(fit["phi"]),
            "beta": float(fit["beta"]),
            "rho": float(fit["rho"]),
            "xi": float(fit["xi"]),
            "zeta": float(fit["zeta"]),
            "eta": float(fit.get("eta", 0.0)),
        }

        diag = local_jacobian(y, d, p, nodes, tau)
        out_rows.append({"seed": seed, "params": p, "diag": diag})

    return {"rows": out_rows, "summary": summarize(out_rows)}


def main() -> None:
    d1739 = read_json("report_qw1739_signed_dynamic_micromodel_derivation.json")
    d1887 = read_json("report_qw1887_signed_coupling_micromodel_rebuild.json")
    d1891 = read_json("report_qw1891_derivational_constraints_from_nadsoliton.json")

    profile_map = {int(r["seed"]): r for r in d1739.get("rows", [])}
    nodes = [2, 5, 8, 11]
    tau = 2.5

    b1887 = evaluate_branch(d1887.get("test_rows", []), profile_map, nodes, tau)
    b1891 = evaluate_branch(d1891.get("test_rows", []), profile_map, nodes, tau)

    s87 = b1887["summary"]
    s91 = b1891["summary"]

    delta = {
        "log10_cond_median_gain": float(s87["log10_cond_median"] - s91["log10_cond_median"]),
        "min_eig_median_gain": float(s91["min_eig_median"] - s87["min_eig_median"]),
        "corr_abs_omega_beta_gain": float(s87["corr_abs_omega_beta_median"] - s91["corr_abs_omega_beta_median"]),
    }

    if delta["log10_cond_median_gain"] >= 0.40 and delta["corr_abs_omega_beta_gain"] >= 0.10:
        verdict = "IDENTIFIABILITY_IMPROVED_BY_DERIVATIONAL_CONSTRAINTS"
    elif delta["log10_cond_median_gain"] >= 0.15:
        verdict = "IDENTIFIABILITY_WEAKLY_IMPROVED"
    else:
        verdict = "IDENTIFIABILITY_NOT_IMPROVED"

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "branch_a": "QW-1887 signed-coupling baseline",
            "branch_b": "QW-1891 constrained by nadsoliton derivation",
            "metric": "local Jacobian/Fisher proxy on test seeds",
            "nodes": nodes,
            "tau": tau,
        },
        "branch_1887": b1887,
        "branch_1891": b1891,
        "delta_1891_minus_1887": delta,
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1892: IDENTIFIABILITY-FIRST AUDIT",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{verdict}**",
        "",
        "## Branch 1887",
        f"- log10(cond) median: {s87['log10_cond_median']:.3f}",
        f"- min eig median: {s87['min_eig_median']:.3e}",
        f"- |corr(omega,beta)| median: {s87['corr_abs_omega_beta_median']:.3f}",
        "",
        "## Branch 1891",
        f"- log10(cond) median: {s91['log10_cond_median']:.3f}",
        f"- min eig median: {s91['min_eig_median']:.3e}",
        f"- |corr(omega,beta)| median: {s91['corr_abs_omega_beta_median']:.3f}",
        "",
        "## Delta (1891 - 1887)",
        f"- log10(cond) gain: {delta['log10_cond_median_gain']:.3f}",
        f"- min eig gain: {delta['min_eig_median_gain']:.3e}",
        f"- |corr(omega,beta)| gain: {delta['corr_abs_omega_beta_gain']:.3f}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1892] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1892] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
