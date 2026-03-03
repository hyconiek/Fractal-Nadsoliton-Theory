#!/usr/bin/env python3
"""
QW-1743: Derivation on oscillatory cohort only.

Rationale:
- Parameter identifiability of spatial oscillation needs actual oscillatory data.
- Build large signed-dynamic ensemble, retain runs with clear sign changes in C(d).
- Fit shared (omega, phi, beta) on selected cohort.
"""

from __future__ import annotations

import json
import math
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1743_oscillatory_cohort_derivation.json"
OUT_MD = ROOT / "RAPORT_QW1743_OSCILLATORY_COHORT_DERIVATION.md"


@dataclass
class Cfg:
    n: int
    seed: int
    xi: float
    p_short: float
    eta: float
    tau: float
    rho: float
    c2: float
    damp: float
    n_steps: int = 300
    burn: int = 120


def ring_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def smooth_periodic(x: np.ndarray, it: int = 7) -> np.ndarray:
    y = x.copy()
    for _ in range(it):
        y = 0.25 * np.roll(y, -1) + 0.5 * y + 0.25 * np.roll(y, 1)
    return y


def build_w(cfg: Cfg, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n
    th = math.pi * np.tanh(smooth_periodic(rng.normal(size=n)) / 1.2)
    q = rng.integers(-2, 3, size=n).astype(float)

    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = math.exp(-d / cfg.xi)
            if rng.random() < cfg.p_short / (d ** cfg.eta):
                amp += 0.9 * (d ** -0.95) * (1.0 + 0.18 * rng.normal())

            dt = abs(th[i] - th[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.4 * dt

            sgn = 1.0 if (math.cos(th[i] - th[j]) + 0.25 * (q[i] - q[j])) >= 0 else -1.0
            sym = sgn * amp * math.exp(-cfg.tau * tors)
            anti = cfg.rho * math.sin(th[i] - th[j]) / (d ** 1.1)
            w[i, j] = sym + anti
            w[j, i] = sym - anti

    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w /= smax
    return w


def simulate_profile(cfg: Cfg) -> np.ndarray:
    rng = np.random.default_rng(cfg.seed)
    w = build_w(cfg, rng)
    n = cfg.n

    x_prev = 0.05 * rng.normal(size=n)
    x_curr = 0.04 * rng.normal(size=n)
    x_curr[0] += 1.0

    lap = w - np.eye(n)
    snaps = []
    for _ in range(cfg.n_steps):
        x_next = 2 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev)
        x_next = np.tanh(1.4 * x_next)
        x_prev, x_curr = x_curr, x_next
        snaps.append(x_curr.copy())

    arr = np.array(snaps[cfg.burn :], dtype=float)
    dmax = min(12, n // 2)
    y = []
    for d in range(1, dmax + 1):
        c = arr * np.roll(arr, -d, axis=1)
        y.append(float(np.mean(c)))
    y = np.array(y, dtype=float)
    if abs(y[0]) > 1e-12:
        y = y / abs(y[0])
    return y


def sign_changes(y: np.ndarray) -> int:
    s = np.sign(y)
    s[s == 0] = 1
    return int(np.sum(s[1:] != s[:-1]))


def oscillation_strength(y: np.ndarray) -> float:
    return float(np.std(y) / max(np.mean(np.abs(y)), 1e-8))


def robust_loss(r: np.ndarray, delta: float = 0.08) -> np.ndarray:
    a = np.abs(r)
    return np.where(a <= delta, 0.5 * a * a, delta * (a - 0.5 * delta))


def pred_basis(d: np.ndarray, omega: float, phi: float, beta: float) -> np.ndarray:
    return np.cos(omega * d + phi) / (1.0 + beta * d)


def boundary_penalty(x: float, lo: float, hi: float, band: float) -> float:
    p = 0.0
    if x < lo + band:
        p += ((lo + band - x) / band) ** 2
    if x > hi - band:
        p += ((x - (hi - band)) / band) ** 2
    return float(p)


def obj(theta: Tuple[float, float, float], y: np.ndarray, w: np.ndarray) -> float:
    om, ph, be = theta
    if not (0.20 <= om <= 1.5 and -math.pi <= ph <= math.pi and 0.001 <= be <= 0.20):
        return float("inf")
    d = np.arange(1, y.shape[1] + 1, dtype=float)
    b = pred_basis(d, om, ph, be)
    bb = float(np.dot(b, b))
    if bb <= 1e-12:
        return float("inf")
    a = (y @ b) / bb
    p = a[:, None] * b[None, :]
    r = y - p
    core = np.sum(w[:, None] * robust_loss(r))
    # Soft non-boundary pressure for identifiability
    pen = 0.7 * boundary_penalty(om, 0.20, 1.5, 0.10) + 0.7 * boundary_penalty(be, 0.001, 0.20, 0.015)
    return float(core + pen)


def refine(start: Tuple[float, float, float], y: np.ndarray, w: np.ndarray) -> Tuple[Tuple[float, float, float], float]:
    cur = (start[0], start[1], start[2])
    fcur = obj(cur, y, w)
    steps = [(0.20, 0.70, 0.03), (0.08, 0.25, 0.012), (0.03, 0.10, 0.005), (0.012, 0.04, 0.002)]
    for so, sp, sb in steps:
        improved = True
        while improved:
            improved = False
            # omega
            best = (cur, fcur)
            for om in np.linspace(max(0.20, cur[0] - so), min(1.5, cur[0] + so), 9):
                th = (float(om), cur[1], cur[2])
                f = obj(th, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # phi
            best = (cur, fcur)
            for ph in np.linspace(cur[1] - sp, cur[1] + sp, 9):
                ph = float((ph + math.pi) % (2 * math.pi) - math.pi)
                th = (cur[0], ph, cur[2])
                f = obj(th, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True

            # beta
            best = (cur, fcur)
            for be in np.linspace(max(0.001, cur[2] - sb), min(0.20, cur[2] + sb), 9):
                th = (cur[0], cur[1], float(be))
                f = obj(th, y, w)
                if f < best[1]:
                    best = (th, f)
            if best[1] < fcur:
                cur, fcur = best
                improved = True
    return cur, fcur


def bootstrap(y: np.ndarray, w: np.ndarray, th0: Tuple[float, float, float], n_boot: int = 120, seed: int = 1743) -> np.ndarray:
    rng = np.random.default_rng(seed)
    n = y.shape[0]
    out = []
    for _ in range(n_boot):
        idx = rng.integers(0, n, size=n)
        th, _ = refine(th0, y[idx], w[idx])
        out.append(th)
    return np.array(out, dtype=float)


def circular_mean(x: np.ndarray) -> float:
    return float(np.arctan2(np.mean(np.sin(x)), np.mean(np.cos(x))))


def circular_std(x: np.ndarray) -> float:
    r = np.sqrt(np.mean(np.sin(x)) ** 2 + np.mean(np.cos(x)) ** 2)
    if r <= 1e-12:
        return float("inf")
    return float(np.sqrt(-2 * np.log(r)))


def main() -> None:
    # Large candidate ensemble.
    candidates: List[Tuple[Cfg, np.ndarray, int, float]] = []
    for n in [72, 96, 120]:
        for k in range(20):
            seed = 174300 + 100 * n + k
            cfg = Cfg(
                n=n,
                seed=seed,
                xi=1.5 + 0.25 * ((k % 3) - 1),
                p_short=0.10 + 0.05 * (k % 2),
                eta=1.45 + 0.20 * ((k + 1) % 3),
                tau=0.45 + 0.18 * (k % 3),
                rho=0.20 + 0.12 * (k % 2),
                c2=0.30 + 0.08 * ((k + 2) % 3),
                damp=0.03 + 0.02 * (k % 2),
            )
            y = simulate_profile(cfg)
            sc = sign_changes(y)
            os = oscillation_strength(y)
            candidates.append((cfg, y, sc, os))

    # Cohort selection: force informative oscillatory regime.
    selected = [x for x in candidates if (x[2] >= 2 and x[3] >= 0.35)]
    if len(selected) < 12:
        # fallback with relaxed threshold if signal is sparse.
        selected = [x for x in candidates if (x[2] >= 1 and x[3] >= 0.28)]

    y = np.array([x[1] for x in selected], dtype=float)
    w = 1.0 / np.clip(np.var(y, axis=1), 1e-5, None)

    starts = [
        (0.35, 0.2, 0.02),
        (0.55, -0.7, 0.04),
        (0.80, 0.5, 0.06),
        (1.05, -1.0, 0.03),
        (1.25, 0.9, 0.09),
    ]
    sols = sorted([refine(st, y, w) for st in starts], key=lambda z: z[1])
    th, fbest = sols[0]

    d = np.arange(1, y.shape[1] + 1, dtype=float)
    b = pred_basis(d, th[0], th[1], th[2])
    bb = float(np.dot(b, b))
    a = (y @ b) / bb
    pred = a[:, None] * b[None, :]
    rmse = np.sqrt(np.mean((y - pred) ** 2, axis=1))
    ss_tot = np.sum((y - np.mean(y, axis=1, keepdims=True)) ** 2, axis=1)
    ss_res = np.sum((y - pred) ** 2, axis=1)
    r2 = np.where(ss_tot > 1e-12, 1.0 - ss_res / ss_tot, 1.0)

    boot = bootstrap(y, w, th, n_boot=120, seed=1743)
    om = boot[:, 0]
    ph = boot[:, 1]
    be = boot[:, 2]

    om_ci = [float(np.quantile(om, 0.025)), float(np.quantile(om, 0.975))]
    ph_std = circular_std(ph)
    ph_mean = circular_mean(ph)
    be_ci = [float(np.quantile(be, 0.025)), float(np.quantile(be, 0.975))]

    om_med = float(np.median(om))
    be_med = float(np.median(be))

    pass_fit = float(np.median(r2)) >= 0.72
    pass_info = len(selected) >= 14
    pass_nonboundary = 0.24 < om_med < 1.45 and 0.003 < be_med < 0.18
    pass_ci = (om_ci[1] - om_ci[0]) <= 0.30 and (be_ci[1] - be_ci[0]) <= 0.09 and ph_std <= 1.0

    if pass_fit and pass_info and pass_nonboundary and pass_ci:
        verdict = "OSCILLATORY_COHORT_DERIVATION_IDENTIFIABLE"
    elif pass_fit and pass_info and pass_nonboundary:
        verdict = "OSCILLATORY_COHORT_DERIVATION_PARTIAL"
    else:
        verdict = "OSCILLATORY_COHORT_DERIVATION_NOT_CLOSED"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_candidates": len(candidates),
            "selection_rule_primary": "sign_changes>=2 and osc_strength>=0.35",
            "selection_rule_fallback": "sign_changes>=1 and osc_strength>=0.28",
            "n_selected": len(selected),
            "n_bootstrap": 120,
        },
        "optimum": {"omega": th[0], "phi": th[1], "beta": th[2], "objective": fbest},
        "cohort_metrics": {
            "median_r2": float(np.median(r2)),
            "q10_r2": float(np.quantile(r2, 0.1)),
            "median_rmse": float(np.median(rmse)),
            "median_sign_changes_selected": float(np.median([x[2] for x in selected])) if selected else None,
            "median_osc_strength_selected": float(np.median([x[3] for x in selected])) if selected else None,
        },
        "bootstrap_summary": {
            "omega_median": om_med,
            "omega_ci95": om_ci,
            "phi_circular_mean": ph_mean,
            "phi_circular_std": ph_std,
            "beta_median": be_med,
            "beta_ci95": be_ci,
        },
        "pass_flags": {
            "fit_quality": pass_fit,
            "informative_cohort": pass_info,
            "nonboundary_solution": pass_nonboundary,
            "ci_width": pass_ci,
        },
        "verdict": verdict,
        "selected_rows_head_120": [
            {
                "n": x[0].n,
                "seed": x[0].seed,
                "sign_changes": x[2],
                "osc_strength": x[3],
                "profile": {str(i + 1): float(x[1][i]) for i in range(len(x[1]))},
            }
            for x in selected[:120]
        ],
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1743: OSCILLATORY COHORT DERIVATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Candidates: {len(candidates)}, Selected: {len(selected)}",
        f"- Optimum: omega={th[0]:.6f}, phi={th[1]:.6f}, beta={th[2]:.6f}",
        f"- Median R2={np.median(r2):.4f}, median RMSE={np.median(rmse):.4f}",
        f"- Bootstrap medians: omega={om_med:.6f}, phi_mean={ph_mean:.6f}, beta={be_med:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- fit_quality: {pass_fit}",
        f"- informative_cohort: {pass_info}",
        f"- nonboundary_solution: {pass_nonboundary}",
        f"- ci_width: {pass_ci}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1743] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1743] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
