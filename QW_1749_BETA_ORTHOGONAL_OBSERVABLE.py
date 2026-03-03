#!/usr/bin/env python3
"""
QW-1749: Beta-orthogonal observable from impulse-response envelope.

Intent:
1) Estimate beta_tors from an observable designed to be weakly sensitive to oscillatory phase.
2) Use absolute/time-aggregated tail response to suppress omega/phi nuisance.
3) Quantify stability under phase-scramble surrogates.
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
OUT_JSON = ROOT / "report_qw1749_beta_orthogonal_observable.json"
OUT_MD = ROOT / "RAPORT_QW1749_BETA_ORTHOGONAL_OBSERVABLE.md"


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
    n_steps: int = 900
    tail_start: int = 320
    tail_end: int = 820


def ring_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def smooth_periodic(x: np.ndarray, iters: int = 6) -> np.ndarray:
    y = x.copy()
    for _ in range(iters):
        y = 0.25 * np.roll(y, -1) + 0.5 * y + 0.25 * np.roll(y, 1)
    return y


def build_w(cfg: Cfg, rng: np.random.Generator) -> np.ndarray:
    n = cfg.n
    th = math.pi * np.tanh(smooth_periodic(rng.normal(size=n), iters=8) / 1.2)
    q = rng.integers(-2, 3, size=n).astype(float)

    w = np.zeros((n, n), dtype=float)
    for i in range(n):
        for j in range(i + 1, n):
            d = ring_dist(i, j, n)
            amp = math.exp(-d / cfg.xi)
            if rng.random() < cfg.p_short / (d ** cfg.eta):
                amp += 0.8 * (d ** -0.9) * (1.0 + 0.14 * rng.normal())

            dt = abs(th[i] - th[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.35 * dt

            sgn = 1.0 if (math.cos(th[i] - th[j]) + 0.2 * (q[i] - q[j])) >= 0 else -1.0
            sym = sgn * amp * math.exp(-cfg.tau * tors)
            anti = cfg.rho * math.sin(th[i] - th[j]) / (d ** 1.12)
            w[i, j] = sym + anti
            w[j, i] = sym - anti

    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w /= smax
    return w


def impulse_response_timeseries(cfg: Cfg, w: np.ndarray, src: int) -> np.ndarray:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    x_curr[src] = 1.0

    hist = []
    for _ in range(cfg.n_steps):
        x_next = 2.0 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev)
        x_next = np.tanh(1.35 * x_next)
        x_prev, x_curr = x_curr, x_next
        hist.append(x_curr.copy())
    return np.array(hist, dtype=float)  # [T,N]


def distance_tail_observable(ts: np.ndarray, src: int, tail_start: int, tail_end: int, dmax: int = 12) -> np.ndarray:
    t1 = max(0, tail_start)
    t2 = min(ts.shape[0], tail_end)
    tail = np.abs(ts[t1:t2])  # abs suppresses phase sign
    n = ts.shape[1]
    vals = []
    for d in range(1, min(dmax, n // 2) + 1):
        j1 = (src + d) % n
        j2 = (src - d) % n
        arr = 0.5 * (tail[:, j1] + tail[:, j2])
        vals.append(float(np.median(arr)))
    y = np.array(vals, dtype=float)
    if y[0] > 1e-12:
        y /= y[0]
    return y


def fit_beta_envelope(y: np.ndarray) -> Tuple[float, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    beta_grid = np.linspace(0.0003, 0.25, 3200)
    best_beta = 0.02
    best_rmse = float("inf")
    for b in beta_grid:
        env = 1.0 / (1.0 + b * d)
        s = float(np.dot(y, env) / max(np.dot(env, env), 1e-12))
        pred = s * env
        rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
        if rmse < best_rmse:
            best_rmse = rmse
            best_beta = float(b)
    return best_beta, best_rmse


def phase_scramble_surrogate(y: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    # Preserve envelope scale, randomize alternating sign pattern and local permutations.
    z = y.copy()
    signs = rng.choice([-1.0, 1.0], size=len(z))
    z = np.abs(z) * signs
    perm_idx = np.arange(len(z))
    rng.shuffle(perm_idx)
    z = z[perm_idx]
    z = np.abs(z)  # estimator uses abs/tail; maintain comparable amplitude set
    if z[0] > 1e-12:
        z /= z[0]
    return z


def main() -> None:
    rows: List[Dict[str, object]] = []
    for n in [96, 120]:
        for k in range(18):
            seed = 174900 + 100 * n + k
            cfg = Cfg(
                n=n,
                seed=seed,
                xi=1.55 + 0.22 * ((k % 3) - 1),
                p_short=0.10 + 0.05 * (k % 2),
                eta=1.45 + 0.18 * ((k + 1) % 3),
                tau=0.45 + 0.16 * (k % 3),
                rho=0.18 + 0.10 * (k % 2),
                c2=0.30 + 0.08 * (k % 3),
                damp=0.03 + 0.02 * (k % 2),
            )
            rng = np.random.default_rng(seed)
            w = build_w(cfg, rng)
            src = seed % n
            ts = impulse_response_timeseries(cfg, w, src=src)
            y = distance_tail_observable(ts, src=src, tail_start=cfg.tail_start, tail_end=cfg.tail_end, dmax=12)
            beta_hat, rmse = fit_beta_envelope(y)

            # nuisance test: phase-scramble surrogate should not strongly move beta
            ys = phase_scramble_surrogate(y, rng)
            beta_scr, rmse_scr = fit_beta_envelope(ys)

            rows.append(
                {
                    "n": n,
                    "seed": seed,
                    "beta_hat": beta_hat,
                    "rmse_env": rmse,
                    "beta_scramble_hat": beta_scr,
                    "rmse_scramble": rmse_scr,
                    "delta_beta_scramble": float(abs(beta_hat - beta_scr)),
                    "tail_observable": {str(i + 1): float(y[i]) for i in range(len(y))},
                }
            )

    beta = np.array([r["beta_hat"] for r in rows], dtype=float)
    rmse = np.array([r["rmse_env"] for r in rows], dtype=float)
    delta_scr = np.array([r["delta_beta_scramble"] for r in rows], dtype=float)

    beta_med = float(np.median(beta))
    beta_iqr = float(np.quantile(beta, 0.75) - np.quantile(beta, 0.25))
    beta_ci = [float(np.quantile(beta, 0.025)), float(np.quantile(beta, 0.975))]
    rmse_med = float(np.median(rmse))
    delta_scr_med = float(np.median(delta_scr))
    delta_scr_q90 = float(np.quantile(delta_scr, 0.9))

    pass_fit = rmse_med <= 0.20
    pass_spread = beta_iqr <= 0.09
    pass_phase_orth = delta_scr_q90 <= 0.06
    pass_nonboundary = (beta_ci[0] > 0.0005) and (beta_ci[1] < 0.23)

    if pass_fit and pass_spread and pass_phase_orth and pass_nonboundary:
        verdict = "BETA_ORTHOGONAL_OBSERVABLE_SUPPORTED"
    elif pass_fit and (pass_spread or pass_phase_orth):
        verdict = "BETA_ORTHOGONAL_OBSERVABLE_PARTIAL"
    else:
        verdict = "BETA_ORTHOGONAL_OBSERVABLE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_total_runs": len(rows),
            "observable": "tail median absolute impulse response by distance",
            "beta_fit": "grid fit to A/(1+beta*d)",
            "nuisance_test": "phase-scramble surrogate",
        },
        "summary": {
            "beta_median": beta_med,
            "beta_iqr": beta_iqr,
            "beta_ci95_empirical": beta_ci,
            "median_rmse_env": rmse_med,
            "delta_beta_scramble_median": delta_scr_med,
            "delta_beta_scramble_q90": delta_scr_q90,
        },
        "pass_flags": {
            "fit_quality": bool(pass_fit),
            "spread_control": bool(pass_spread),
            "phase_orthogonality": bool(pass_phase_orth),
            "nonboundary_beta": bool(pass_nonboundary),
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1749: BETA ORTHOGONAL OBSERVABLE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Runs: {len(rows)}",
        f"- beta median/IQR: {beta_med:.6f} / {beta_iqr:.6f}",
        f"- beta empirical CI95: [{beta_ci[0]:.6f}, {beta_ci[1]:.6f}]",
        f"- median envelope RMSE: {rmse_med:.6f}",
        f"- scramble delta beta q90: {delta_scr_q90:.6f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- fit_quality: {pass_fit}",
        f"- spread_control: {pass_spread}",
        f"- phase_orthogonality: {pass_phase_orth}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1749] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1749] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
