#!/usr/bin/env python3
"""
QW-1763: Memory-field microdynamics beta evidence test.

Hypothesis:
- A local memory field coupled to amplitude can generate effective
  envelope attenuation that yields stable, non-null positive beta_tors
  evidence without post-hoc ansatz augmentation at the fit layer.
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
OUT_JSON = ROOT / "report_qw1763_memory_field_microdynamics_beta_test.json"
OUT_MD = ROOT / "RAPORT_QW1763_MEMORY_FIELD_MICRODYNAMICS_BETA_TEST.md"


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
    mem_mu: float
    mem_rate: float
    n_steps: int = 1000
    tail_start: int = 380
    tail_end: int = 920
    dmax: int = 14


def ring_dist(i: int, j: int, n: int) -> int:
    d = abs(i - j)
    return min(d, n - d)


def smooth_periodic(x: np.ndarray, iters: int = 7) -> np.ndarray:
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
                amp += 0.8 * (d ** -0.95) * (1.0 + 0.13 * rng.normal())
            dt = abs(th[i] - th[j]) / math.pi
            dq = abs(q[i] - q[j])
            tors = dq + 0.35 * dt
            sgn = 1.0 if (math.cos(th[i] - th[j]) + 0.2 * (q[i] - q[j])) >= 0 else -1.0
            sym = sgn * amp * math.exp(-cfg.tau * tors)
            anti = cfg.rho * math.sin(th[i] - th[j]) / (d ** 1.1)
            w[i, j] = sym + anti
            w[j, i] = sym - anti
    smax = float(np.linalg.norm(w, 2))
    if smax > 1e-12:
        w /= smax
    return w


def impulse_memory_timeseries(cfg: Cfg, w: np.ndarray, src: int) -> np.ndarray:
    n = cfg.n
    lap = w - np.eye(n)
    x_prev = np.zeros(n, dtype=float)
    x_curr = np.zeros(n, dtype=float)
    x_curr[src] = 1.0
    mem = np.zeros(n, dtype=float)

    hist = []
    for _ in range(cfg.n_steps):
        mem = (1.0 - cfg.mem_rate) * mem + cfg.mem_rate * np.abs(x_curr)
        damping = cfg.mem_mu * mem * x_curr
        x_next = 2.0 * x_curr - x_prev + cfg.c2 * (lap @ x_curr) - cfg.damp * (x_curr - x_prev) - damping
        x_next = np.tanh(1.35 * x_next)
        x_prev, x_curr = x_curr, x_next
        hist.append(x_curr.copy())
    return np.array(hist, dtype=float)


def distance_tail_observable(ts: np.ndarray, src: int, tail_start: int, tail_end: int, dmax: int) -> np.ndarray:
    t1 = max(0, tail_start)
    t2 = min(ts.shape[0], tail_end)
    tail = np.abs(ts[t1:t2])
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


def fit_null_and_beta(y: np.ndarray) -> Dict[str, float]:
    d = np.arange(1, len(y) + 1, dtype=float)

    a0 = float(np.mean(y))
    p0 = np.full_like(y, a0)
    sse0 = float(np.sum((y - p0) ** 2))
    rmse0 = float(np.sqrt(np.mean((y - p0) ** 2)))

    beta_grid = np.logspace(np.log10(1e-5), np.log10(0.25), 2200)
    best = {"beta": 1e-3, "a": a0, "sse": float("inf"), "rmse": float("inf")}
    for b in beta_grid:
        f = 1.0 / (1.0 + b * d)
        den = float(np.dot(f, f))
        if den <= 1e-15:
            continue
        a = float(np.dot(y, f) / den)
        p1 = a * f
        sse = float(np.sum((y - p1) ** 2))
        if sse < best["sse"]:
            best = {
                "beta": float(b),
                "a": a,
                "sse": sse,
                "rmse": float(np.sqrt(np.mean((y - p1) ** 2))),
            }

    n = len(y)
    sse0 = max(sse0, 1e-15)
    sse1 = max(best["sse"], 1e-15)
    bic0 = float(n * np.log(sse0 / n) + 1.0 * np.log(n))
    bic1 = float(n * np.log(sse1 / n) + 2.0 * np.log(n))
    dbic = bic0 - bic1
    return {
        "beta_hat": float(best["beta"]),
        "rmse_null": rmse0,
        "rmse_beta": float(best["rmse"]),
        "delta_bic_null_minus_beta": float(dbic),
    }


def main() -> None:
    rows: List[Dict[str, object]] = []
    for n in [96, 120]:
        for k in range(12):
            seed = 176300 + 100 * n + k
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
                mem_mu=0.04 + 0.06 * (k % 4),
                mem_rate=0.02 + 0.015 * (k % 3),
            )
            rng = np.random.default_rng(seed)
            w = build_w(cfg, rng)
            src = seed % n
            ts = impulse_memory_timeseries(cfg, w, src=src)
            y = distance_tail_observable(ts, src=src, tail_start=cfg.tail_start, tail_end=cfg.tail_end, dmax=cfg.dmax)
            fit = fit_null_and_beta(y)

            rows.append(
                {
                    "n": n,
                    "seed": seed,
                    "mem_mu": cfg.mem_mu,
                    "mem_rate": cfg.mem_rate,
                    "beta_hat": fit["beta_hat"],
                    "rmse_null": fit["rmse_null"],
                    "rmse_beta": fit["rmse_beta"],
                    "delta_bic_null_minus_beta": fit["delta_bic_null_minus_beta"],
                }
            )

    beta = np.array([r["beta_hat"] for r in rows], dtype=float)
    dbic = np.array([r["delta_bic_null_minus_beta"] for r in rows], dtype=float)
    rmse_b = np.array([r["rmse_beta"] for r in rows], dtype=float)
    rmse_gain = np.array([r["rmse_null"] - r["rmse_beta"] for r in rows], dtype=float)

    summary = {
        "n_rows": len(rows),
        "beta_median": float(np.median(beta)),
        "beta_iqr": float(np.quantile(beta, 0.75) - np.quantile(beta, 0.25)),
        "beta_ci95_empirical": [float(np.quantile(beta, 0.025)), float(np.quantile(beta, 0.975))],
        "median_delta_bic_null_minus_beta": float(np.median(dbic)),
        "positive_delta_bic_rate": float(np.mean(dbic > 0.0)),
        "beta_nonboundary_rate": float(np.mean((beta > 5e-4) & (beta < 0.23))),
        "median_rmse_beta": float(np.median(rmse_b)),
        "median_rmse_gain_null_minus_beta": float(np.median(rmse_gain)),
    }

    pass_n = len(rows) >= 20
    pass_evidence = summary["median_delta_bic_null_minus_beta"] >= 1.0 and summary["positive_delta_bic_rate"] >= 0.60
    pass_beta_nonboundary = summary["beta_nonboundary_rate"] >= 0.60
    pass_fit = summary["median_rmse_beta"] <= 0.10 and summary["median_rmse_gain_null_minus_beta"] >= 0.01

    if pass_n and pass_evidence and pass_beta_nonboundary and pass_fit:
        verdict = "MEMORY_FIELD_BETA_EVIDENCE_SUPPORTED"
    elif pass_n and pass_fit and (pass_evidence or pass_beta_nonboundary):
        verdict = "MEMORY_FIELD_BETA_EVIDENCE_PARTIAL"
    else:
        verdict = "MEMORY_FIELD_BETA_EVIDENCE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "protocol": {
            "n_grid": [96, 120],
            "k_per_n": 12,
            "memory_params": {
                "mem_mu": [0.04, 0.10, 0.16, 0.22],
                "mem_rate": [0.02, 0.035, 0.05],
            },
        },
        "summary": summary,
        "pass_flags": {
            "n_runs": bool(pass_n),
            "evidence_strength": bool(pass_evidence),
            "beta_nonboundary": bool(pass_beta_nonboundary),
            "fit_improvement": bool(pass_fit),
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1763: MEMORY FIELD MICRODYNAMICS BETA TEST",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Rows: {len(rows)}",
        f"- beta median/CI95: {summary['beta_median']:.6f} / {summary['beta_ci95_empirical']}",
        f"- median ΔBIC(null-beta): {summary['median_delta_bic_null_minus_beta']:.3f}",
        f"- positive ΔBIC rate: {summary['positive_delta_bic_rate']:.3f}",
        f"- beta nonboundary rate: {summary['beta_nonboundary_rate']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- n_runs: {pass_n}",
        f"- evidence_strength: {pass_evidence}",
        f"- beta_nonboundary: {pass_beta_nonboundary}",
        f"- fit_improvement: {pass_fit}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1763] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1763] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
