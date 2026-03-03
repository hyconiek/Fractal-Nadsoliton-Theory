#!/usr/bin/env python3
"""
QW-1917: Triad derivation (beta_tors, omega, phi) from micromodel profiles, without ansatz.

Design:
- uses raw profile rows from QW-1739 (signed dynamic micromodel),
- does not inject canonical targets (pi/4, pi/6, 0.01),
- performs extended-bound global fitting with multistart and bootstrap stability.

Model:
    y_r(d) = a_r * cos(omega*d + phi) / (1 + beta*d) + eps
with per-run nuisance amplitude a_r eliminated analytically.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1917_triad_derivation_no_ansatz_profile.json"
OUT_MD = ROOT / "RAPORT_QW1917_TRIAD_DERIVATION_NO_ANSATZ_PROFILE.md"


def load_1739() -> Dict:
    return json.loads((ROOT / "report_qw1739_signed_dynamic_micromodel_derivation.json").read_text(encoding="utf-8"))


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
    if omega < 0.02 or omega > 1.90 or beta <= 1e-4 or beta > 2.0:
        return float("inf")

    b = np.cos(omega * d + phi) / (1.0 + beta * d)
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")

    # Per-run optimal nuisance amplitude a_r.
    a = (y @ b) / bb
    pred = a[:, None] * b[None, :]
    res = y - pred
    sse = np.sum(weights[:, None] * (res ** 2))
    return float(sse)


def safe_objective(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, weights: np.ndarray) -> float:
    f = objective(theta, d, y, weights)
    if not np.isfinite(f):
        return 1e30
    return float(f)


def coordinate_refine(
    start: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    quick: bool = False,
) -> Tuple[Tuple[float, float, float], float]:
    omega, phi, beta = start
    omega = float(np.clip(omega, 0.02, 1.90))
    phi = wrap_phi(phi)
    beta = float(np.clip(beta, 1e-4, 2.0))

    cur = (omega, phi, beta)
    fcur = safe_objective(cur, d, y, w)

    if quick:
        step_plan = [
            (0.18, 0.65, 0.35),
            (0.07, 0.25, 0.14),
            (0.025, 0.09, 0.05),
        ]
    else:
        step_plan = [
            (0.30, 1.00, 0.50),
            (0.12, 0.45, 0.22),
            (0.045, 0.16, 0.08),
            (0.015, 0.07, 0.03),
        ]

    for so, sp, sb in step_plan:
        improved = True
        while improved:
            improved = False

            cand_omega = np.linspace(max(0.02, cur[0] - so), min(1.90, cur[0] + so), 9)
            best_local = (cur, fcur)
            for om in cand_omega:
                th = (float(om), cur[1], cur[2])
                f = safe_objective(th, d, y, w)
                if f < best_local[1]:
                    best_local = (th, f)
            if best_local[1] < fcur:
                cur, fcur = best_local
                improved = True

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

            cand_beta = np.linspace(max(1e-4, cur[2] - sb), min(2.0, cur[2] + sb), 9)
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


def unique_modes(params: List[Tuple[float, float, float]], tol: Tuple[float, float, float]) -> List[Tuple[float, float, float]]:
    modes: List[Tuple[float, float, float]] = []
    for p in params:
        add = True
        for m in modes:
            d0 = abs(p[0] - m[0]) < tol[0]
            d1 = abs(((p[1] - m[1] + math.pi) % (2.0 * math.pi)) - math.pi) < tol[1]
            d2 = abs(p[2] - m[2]) < tol[2]
            if d0 and d1 and d2:
                add = False
                break
        if add:
            modes.append(p)
    return modes


def hessian_numeric(theta: Tuple[float, float, float], d: np.ndarray, y: np.ndarray, w: np.ndarray) -> np.ndarray:
    x = np.array(theta, dtype=float)
    h = np.array([0.006, 0.030, 0.018], dtype=float)
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


def bootstrap_stability(
    theta0: Tuple[float, float, float],
    d: np.ndarray,
    y: np.ndarray,
    w: np.ndarray,
    n_boot: int,
    seed: int,
) -> Dict[str, object]:
    rng = np.random.default_rng(seed)
    n = y.shape[0]

    vals = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        yb = y[idx]
        wb = w[idx]
        th, _ = coordinate_refine(theta0, d, yb, wb, quick=True)
        vals.append(th)

    arr = np.array(vals, dtype=float)
    om = arr[:, 0]
    ph = arr[:, 1]
    be = arr[:, 2]

    ph_unwrap = np.unwrap(ph)
    ph_q05 = float(np.quantile(ph_unwrap, 0.05))
    ph_q95 = float(np.quantile(ph_unwrap, 0.95))

    out = {
        "n_boot": int(n_boot),
        "omega": {
            "median": float(np.median(om)),
            "q05": float(np.quantile(om, 0.05)),
            "q95": float(np.quantile(om, 0.95)),
        },
        "phi": {
            "median": float(np.median(ph_unwrap)),
            "q05": ph_q05,
            "q95": ph_q95,
        },
        "beta": {
            "median": float(np.median(be)),
            "q05": float(np.quantile(be, 0.05)),
            "q95": float(np.quantile(be, 0.95)),
        },
        "boundary_fraction": {
            "omega_low": float(np.mean(om <= 0.03)),
            "omega_high": float(np.mean(om >= 1.88)),
            "beta_low": float(np.mean(be <= 0.005)),
            "beta_high": float(np.mean(be >= 1.80)),
        },
    }
    return out


def main() -> None:
    d1739 = load_1739()
    rows = d1739.get("rows", [])
    if not rows:
        raise RuntimeError("No rows in QW-1739 report.")

    d, y, w = build_dataset(rows)

    # Multistart: include previous shifted solutions and broad random starts.
    starts: List[Tuple[float, float, float]] = [
        (0.118, 0.475, 0.300),
        (0.088, 0.890, 0.400),
        (0.35, 0.00, 0.08),
        (0.70, 0.50, 0.02),
        (1.10, -0.30, 0.20),
        (0.20, -0.90, 0.80),
    ]

    rng = np.random.default_rng(191700)
    for _ in range(22):
        starts.append(
            (
                float(rng.uniform(0.02, 1.90)),
                float(rng.uniform(-math.pi, math.pi)),
                float(10 ** rng.uniform(-3.0, math.log10(2.0))),
            )
        )

    sols: List[Tuple[Tuple[float, float, float], float]] = []
    for st in starts:
        th, f = coordinate_refine(st, d, y, w, quick=False)
        sols.append((th, f))

    sols = sorted(sols, key=lambda x: x[1])
    theta_best = sols[0][0]
    f_best = float(sols[0][1])

    params = [s[0] for s in sols]
    modes = unique_modes(params, tol=(0.05, 0.18, 0.08))

    hess = hessian_numeric(theta_best, d, y, w)
    hsym = 0.5 * (hess + hess.T)
    hsym = np.nan_to_num(hsym, nan=1e6, posinf=1e6, neginf=-1e6)
    eigvals = np.linalg.eigvalsh(hsym)
    min_eig = float(np.min(eigvals))
    reg = 0.0
    if min_eig <= 1e-8:
        reg = (1e-8 - min_eig) + 1e-6
    hreg = hsym + reg * np.eye(3)
    cov = np.linalg.inv(hreg)

    std = np.sqrt(np.clip(np.diag(cov), 0.0, None))
    ci95 = 1.96 * std
    cond = float(np.linalg.cond(hreg))

    corr = cov / np.sqrt(np.outer(np.diag(cov), np.diag(cov)))
    corr = np.clip(corr, -1.0, 1.0)

    omega_hat, phi_hat, beta_hat = theta_best

    bstrap = bootstrap_stability(theta_best, d, y, w, n_boot=120, seed=191701)
    om_w = float(bstrap["omega"]["q95"] - bstrap["omega"]["q05"])
    ph_w = float(bstrap["phi"]["q95"] - bstrap["phi"]["q05"])
    be_w = float(bstrap["beta"]["q95"] - bstrap["beta"]["q05"])
    bf = bstrap["boundary_fraction"]
    bf_max = max(float(v) for v in bf.values())

    pass_cond = cond <= 1e7
    pass_ci = (om_w <= 0.40) and (ph_w <= 1.60) and (be_w <= 0.90)
    pass_modes = len(modes) <= 4
    pass_boundary = bf_max <= 0.55
    pass_corr = float(np.max(np.abs(corr - np.eye(3)))) <= 0.98

    if pass_cond and pass_ci and pass_modes and pass_boundary and pass_corr:
        verdict = "TRIAD_NO_ANSATZ_DERIVATION_STRONG"
    elif pass_modes and (pass_ci or pass_cond):
        verdict = "TRIAD_NO_ANSATZ_DERIVATION_PARTIAL"
    else:
        verdict = "TRIAD_NO_ANSATZ_DERIVATION_WEAK"

    out = {
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
            "matrix": hsym.tolist(),
            "regularization_added": reg,
            "condition_number": cond,
            "eigenvalues": [float(v) for v in eigvals],
        },
        "uncertainty": {
            "std": {"omega": float(std[0]), "phi": float(std[1]), "beta": float(std[2])},
            "ci95_halfwidth": {"omega": float(ci95[0]), "phi": float(ci95[1]), "beta": float(ci95[2])},
            "correlation_matrix": corr.tolist(),
        },
        "multistart": {
            "n_starts": len(starts),
            "solutions": [{"theta": [float(x) for x in th], "objective": float(f)} for th, f in sols],
            "n_unique_modes": len(modes),
            "unique_modes": [[float(x) for x in m] for m in modes],
        },
        "bootstrap": bstrap,
        "pass_flags": {
            "conditioning": bool(pass_cond),
            "bootstrap_ci_width": bool(pass_ci),
            "mode_count": bool(pass_modes),
            "boundary_pressure": bool(pass_boundary),
            "correlation": bool(pass_corr),
        },
        "verdict": verdict,
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1917: TRIAD DERIVATION NO ANSATZ (PROFILE)",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Dataset runs: {y.shape[0]}",
        f"- Optimum: omega={omega_hat:.6f}, phi={phi_hat:.6f}, beta={beta_hat:.6f}",
        f"- Weighted SSE: {f_best:.6f}",
        f"- Hessian cond: {cond:.3e}",
        f"- Unique modes: {len(modes)}",
        f"- Bootstrap widths (q95-q05): domega={om_w:.4f}, dphi={ph_w:.4f}, dbeta={be_w:.4f}",
        f"- Boundary pressure max fraction: {bf_max:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- conditioning: {pass_cond}",
        f"- bootstrap_ci_width: {pass_ci}",
        f"- mode_count: {pass_modes}",
        f"- boundary_pressure: {pass_boundary}",
        f"- correlation: {pass_corr}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1917] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1917] Saved MD:   {OUT_MD.name}")
    print(f"[QW-1917] verdict={verdict}")


if __name__ == "__main__":
    main()
