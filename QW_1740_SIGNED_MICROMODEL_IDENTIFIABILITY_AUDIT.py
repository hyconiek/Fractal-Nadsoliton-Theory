#!/usr/bin/env python3
"""
QW-1740: Identifiability audit for (omega, phi, beta) from QW-1739 profiles.

Methods:
1) Joint weighted SSE objective with per-run nuisance amplitude eliminated analytically.
2) Multistart coordinate search.
3) Local Hessian-based uncertainty and parameter correlation diagnostics.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1740_signed_micromodel_identifiability_audit.json"
OUT_MD = ROOT / "RAPORT_QW1740_SIGNED_MICROMODEL_IDENTIFIABILITY_AUDIT.md"


def load_1739() -> Dict:
    p = ROOT / "report_qw1739_signed_dynamic_micromodel_derivation.json"
    return json.loads(p.read_text(encoding="utf-8"))


def wrap_phi(phi: float) -> float:
    return float((phi + math.pi) % (2.0 * math.pi) - math.pi)


def build_dataset(rows: List[Dict[str, object]]) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    d = np.arange(1, 13, dtype=float)
    y = []
    w = []
    for r in rows:
        prof = r["profile"]
        vec = np.array([float(prof[str(i)]) for i in range(1, 13)], dtype=float)
        y.append(vec)
        var = float(np.var(vec))
        w.append(1.0 / max(var, 1e-5))
    return d, np.array(y, dtype=float), np.array(w, dtype=float)


def objective(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    omega, phi, beta = theta
    if omega < 0.05 or omega > 1.9 or beta <= 0.0005 or beta > 0.4:
        return float("inf")
    b = np.cos(omega * d + phi) / (1.0 + beta * d)
    bb = float(np.dot(b, b))
    if bb < 1e-12:
        return float("inf")
    # Per-run best amplitude nuisance:
    a = (y @ b) / bb  # shape [R]
    pred = a[:, None] * b[None, :]
    res = y - pred
    sse = np.sum(weights[:, None] * (res ** 2))
    return float(sse)


def safe_objective(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    f = objective(theta, d, y, weights)
    if not np.isfinite(f):
        return 1e12
    return float(f)


def coordinate_refine(start: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, w: np.ndarray) -> Tuple[Tuple[float, float, float], float]:
    omega, phi, beta = start
    omega = float(np.clip(omega, 0.05, 1.9))
    phi = wrap_phi(phi)
    beta = float(np.clip(beta, 0.001, 0.4))
    cur = (omega, phi, beta)
    fcur = safe_objective(cur, d, y, w)

    step_plan = [
        (0.20, 0.80, 0.060),
        (0.08, 0.35, 0.025),
        (0.03, 0.14, 0.010),
        (0.012, 0.06, 0.004),
    ]
    for so, sp, sb in step_plan:
        improved = True
        while improved:
            improved = False

            # omega scan
            cand_omega = np.linspace(max(0.05, cur[0] - so), min(1.9, cur[0] + so), 9)
            best_local = (cur, fcur)
            for om in cand_omega:
                th = (float(om), cur[1], cur[2])
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

            # phi scan
            cand_phi = np.linspace(cur[1] - sp, cur[1] + sp, 9)
            best_local = (cur, fcur)
            for ph in cand_phi:
                th = (cur[0], wrap_phi(float(ph)), cur[2])
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

            # beta scan
            cand_beta = np.linspace(max(0.001, cur[2] - sb), min(0.4, cur[2] + sb), 9)
            best_local = (cur, fcur)
            for be in cand_beta:
                th = (cur[0], cur[1], float(be))
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

    return cur, fcur


def hessian_numeric(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, w: np.ndarray) -> np.ndarray:
    x = np.array(theta, dtype=float)
    h = np.array([0.004, 0.020, 0.0025], dtype=float)
    f0 = safe_objective(tuple(x), d, y, w)
    hess = np.zeros((3, 3), dtype=float)

    for i in range(3):
        xp = x.copy()
        xm = x.copy()
        xp[i] += h[i]
        xm[i] -= h[i]
        if i == 1:
            xp[i] = wrap_phi(float(xp[i]))
            xm[i] = wrap_phi(float(xm[i]))
        fp = safe_objective(tuple(xp), d, y, w)
        fm = safe_objective(tuple(xm), d, y, w)
        hess[i, i] = (fp - 2.0 * f0 + fm) / (h[i] ** 2)

    for i in range(3):
        for j in range(i + 1, 3):
            xpp = x.copy()
            xpm = x.copy()
            xmp = x.copy()
            xmm = x.copy()
            xpp[i] += h[i]
            xpp[j] += h[j]
            xpm[i] += h[i]
            xpm[j] -= h[j]
            xmp[i] -= h[i]
            xmp[j] += h[j]
            xmm[i] -= h[i]
            xmm[j] -= h[j]
            if i == 1:
                xpp[i] = wrap_phi(float(xpp[i]))
                xpm[i] = wrap_phi(float(xpm[i]))
                xmp[i] = wrap_phi(float(xmp[i]))
                xmm[i] = wrap_phi(float(xmm[i]))
            if j == 1:
                xpp[j] = wrap_phi(float(xpp[j]))
                xpm[j] = wrap_phi(float(xpm[j]))
                xmp[j] = wrap_phi(float(xmp[j]))
                xmm[j] = wrap_phi(float(xmm[j]))
            fpp = safe_objective(tuple(xpp), d, y, w)
            fpm = safe_objective(tuple(xpm), d, y, w)
            fmp = safe_objective(tuple(xmp), d, y, w)
            fmm = safe_objective(tuple(xmm), d, y, w)
            val = (fpp - fpm - fmp + fmm) / (4.0 * h[i] * h[j])
            hess[i, j] = val
            hess[j, i] = val

    return hess


def unique_modes(params: List[Tuple[float, float, float]], tol: Tuple[float, float, float] = (0.03, 0.10, 0.01)) -> List[Tuple[float, float, float]]:
    modes: List[Tuple[float, float, float]] = []
    for p in params:
        add = True
        for m in modes:
            d0 = abs(p[0] - m[0]) < tol[0]
            d1 = abs(((p[1] - m[1] + math.pi) % (2 * math.pi)) - math.pi) < tol[1]
            d2 = abs(p[2] - m[2]) < tol[2]
            if d0 and d1 and d2:
                add = False
                break
        if add:
            modes.append(p)
    return modes


def main() -> None:
    d1739 = load_1739()
    rows = d1739.get("rows", [])
    if not rows:
        raise RuntimeError("No rows in QW-1739 report.")

    d, y, weights = build_dataset(rows)

    omega_med = float(d1739["summary"]["omega_median"])
    phi_med = float(d1739["summary"]["phi_circular_mean"])
    beta_med = float(d1739["summary"]["beta_median"])

    starts = [
        (omega_med, phi_med, beta_med),
        (max(0.05, omega_med * 0.8), wrap_phi(phi_med - 0.5), min(0.4, beta_med * 1.2)),
        (min(1.9, omega_med * 1.2), wrap_phi(phi_med + 0.6), max(0.001, beta_med * 0.8)),
        (0.55, 0.0, 0.02),
        (0.95, 0.7, 0.05),
        (0.35, -0.9, 0.12),
    ]

    sols = []
    for st in starts:
        th, f = coordinate_refine(st, d, y, weights)
        sols.append((th, f))

    sols = sorted(sols, key=lambda x: x[1])
    theta_best = sols[0][0]
    f_best = float(sols[0][1])

    param_solutions = [s[0] for s in sols]
    modes = unique_modes(param_solutions)

    hess = hessian_numeric(theta_best, d, y, weights)
    # Stabilize Hessian for inversion.
    hess_sym = 0.5 * (hess + hess.T)
    hess_sym = np.nan_to_num(hess_sym, nan=1e6, posinf=1e6, neginf=-1e6)
    eigvals = np.linalg.eigvalsh(hess_sym)
    min_eig = float(np.min(eigvals))
    reg = 0.0
    if min_eig <= 1e-8:
        reg = (1e-8 - min_eig) + 1e-6
    hess_reg = hess_sym + reg * np.eye(3)
    cov = np.linalg.inv(hess_reg)

    std = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    ci95 = 1.96 * std
    cond = float(np.linalg.cond(hess_reg))

    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    corr = np.clip(corr, -1.0, 1.0)

    omega_hat, phi_hat, beta_hat = theta_best
    ci_omega, ci_phi, ci_beta = float(ci95[0]), float(ci95[1]), float(ci95[2])

    pass_cond = cond <= 1e5
    pass_ci = ci_omega <= 0.18 and ci_phi <= 0.55 and ci_beta <= 0.06
    pass_modes = len(modes) <= 2
    pass_corr = float(np.max(np.abs(corr - np.eye(3)))) <= 0.95

    if pass_cond and pass_ci and pass_modes and pass_corr:
        verdict = "IDENTIFIABILITY_STRONG"
    elif pass_modes and (pass_cond or pass_ci):
        verdict = "IDENTIFIABILITY_MODERATE"
    else:
        verdict = "IDENTIFIABILITY_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "source_report": "report_qw1739_signed_dynamic_micromodel_derivation.json",
        "dataset": {
            "n_runs": int(y.shape[0]),
            "n_points_per_profile": int(y.shape[1]),
        },
        "optimum": {
            "omega": float(omega_hat),
            "phi": float(phi_hat),
            "beta": float(beta_hat),
            "objective_sse_weighted": f_best,
        },
        "hessian": {
            "matrix": hess_sym.tolist(),
            "regularization_added": reg,
            "condition_number": cond,
            "eigenvalues": [float(v) for v in eigvals],
        },
        "uncertainty": {
            "std": {"omega": float(std[0]), "phi": float(std[1]), "beta": float(std[2])},
            "ci95_halfwidth": {"omega": ci_omega, "phi": ci_phi, "beta": ci_beta},
            "correlation_matrix": corr.tolist(),
        },
        "multistart": {
            "n_starts": len(starts),
            "solutions": [{"theta": [float(x) for x in th], "objective": float(f)} for th, f in sols],
            "n_unique_modes": len(modes),
            "unique_modes": [[float(x) for x in m] for m in modes],
        },
        "pass_flags": {
            "conditioning": pass_cond,
            "ci_width": pass_ci,
            "mode_count": pass_modes,
            "correlation": pass_corr,
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1740: SIGNED MICROMODEL IDENTIFIABILITY AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Dataset runs: {y.shape[0]}",
        f"- Optimum: omega={omega_hat:.6f}, phi={phi_hat:.6f}, beta={beta_hat:.6f}",
        f"- Weighted SSE: {f_best:.6f}",
        f"- Hessian cond: {cond:.3e}",
        f"- CI95 half-widths: domega={ci_omega:.4f}, dphi={ci_phi:.4f}, dbeta={ci_beta:.4f}",
        f"- Unique modes: {len(modes)}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- conditioning: {pass_cond}",
        f"- ci_width: {pass_ci}",
        f"- mode_count: {pass_modes}",
        f"- correlation: {pass_corr}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1740] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1740] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
