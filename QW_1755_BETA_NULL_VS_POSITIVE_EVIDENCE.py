#!/usr/bin/env python3
"""
QW-1755: Beta evidence test against null (beta=0) envelope model.

Goal:
- Verify whether positive beta_tors is statistically required by the
  QW-1749 tail observable, not only numerically fitted near grid edges.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1755_beta_null_vs_positive_evidence.json"
OUT_MD = ROOT / "RAPORT_QW1755_BETA_NULL_VS_POSITIVE_EVIDENCE.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def normalize_weights(w: np.ndarray) -> np.ndarray:
    z = np.nan_to_num(np.array(w, dtype=float), nan=0.0, posinf=0.0, neginf=0.0)
    z = np.clip(z, 0.0, None)
    s = float(np.sum(z))
    if s <= 1e-15:
        return np.ones_like(z) / max(len(z), 1)
    return z / s


def weighted_quantile(values: np.ndarray, weights: np.ndarray, q: np.ndarray) -> np.ndarray:
    x = np.array(values, dtype=float)
    w = normalize_weights(weights)
    qq = np.array(q, dtype=float)
    idx = np.argsort(x)
    xs = x[idx]
    ws = w[idx]
    cw = np.cumsum(ws)
    return np.interp(np.clip(qq, 0.0, 1.0), cw, xs)


def fit_null_and_beta(y: np.ndarray) -> Dict[str, float]:
    n = len(y)
    d = np.arange(1, n + 1, dtype=float)

    # Null model M0: y ~= A
    a0 = float(np.mean(y))
    p0 = np.full_like(y, a0)
    sse0 = float(np.sum((y - p0) ** 2))
    rmse0 = float(np.sqrt(np.mean((y - p0) ** 2)))

    # Positive-beta model M1: y ~= A/(1+beta*d)
    beta_grid = np.logspace(np.log10(1e-5), np.log10(0.3), 4200)
    best = {"beta": 1e-3, "a": a0, "sse": float("inf"), "rmse": float("inf")}
    for b in beta_grid:
        env = 1.0 / (1.0 + b * d)
        den = float(np.dot(env, env))
        if den <= 1e-15:
            continue
        a1 = float(np.dot(y, env) / den)
        p1 = a1 * env
        sse1 = float(np.sum((y - p1) ** 2))
        if sse1 < best["sse"]:
            best = {
                "beta": float(b),
                "a": a1,
                "sse": sse1,
                "rmse": float(np.sqrt(np.mean((y - p1) ** 2))),
            }

    nobs = max(n, 1)
    sse0s = max(sse0, 1e-15)
    sse1s = max(best["sse"], 1e-15)
    bic0 = float(nobs * np.log(sse0s / nobs) + 1.0 * np.log(nobs))
    bic1 = float(nobs * np.log(sse1s / nobs) + 2.0 * np.log(nobs))
    aic0 = float(nobs * np.log(sse0s / nobs) + 2.0 * 1.0)
    aic1 = float(nobs * np.log(sse1s / nobs) + 2.0 * 2.0)

    return {
        "beta_hat": float(best["beta"]),
        "rmse_null": rmse0,
        "rmse_beta": float(best["rmse"]),
        "sse_null": sse0,
        "sse_beta": float(best["sse"]),
        "delta_bic_null_minus_beta": float(bic0 - bic1),
        "delta_aic_null_minus_beta": float(aic0 - aic1),
    }


def main() -> None:
    d1749 = load("report_qw1749_beta_orthogonal_observable.json")
    rows_in = d1749.get("rows", [])
    if len(rows_in) < 12:
        raise RuntimeError("QW-1749 rows too short for evidence analysis.")

    rows: List[Dict[str, object]] = []
    for r in rows_in:
        ymap = r.get("tail_observable", {})
        y = np.array([float(ymap[str(i)]) for i in range(1, len(ymap) + 1)], dtype=float)
        if len(y) < 6:
            continue
        if y[0] > 1e-12:
            y = y / y[0]
        fit = fit_null_and_beta(y)
        rows.append(
            {
                "seed": int(r.get("seed", -1)),
                "n": int(r.get("n", -1)),
                "beta_hat": fit["beta_hat"],
                "rmse_null": fit["rmse_null"],
                "rmse_beta": fit["rmse_beta"],
                "delta_bic_null_minus_beta": fit["delta_bic_null_minus_beta"],
                "delta_aic_null_minus_beta": fit["delta_aic_null_minus_beta"],
            }
        )

    if len(rows) < 12:
        raise RuntimeError("Too few valid rows after envelope extraction.")

    beta = np.array([float(r["beta_hat"]) for r in rows], dtype=float)
    rmse_b = np.array([float(r["rmse_beta"]) for r in rows], dtype=float)
    dbic = np.array([float(r["delta_bic_null_minus_beta"]) for r in rows], dtype=float)
    daic = np.array([float(r["delta_aic_null_minus_beta"]) for r in rows], dtype=float)

    # Row weight favors fit quality and model evidence for positive beta.
    w = (1.0 / np.clip(rmse_b ** 2, 1e-6, None)) * np.exp(np.clip(dbic / 8.0, -4.0, 4.0))
    w = normalize_weights(w)

    bq = weighted_quantile(beta, w, np.array([0.025, 0.5, 0.975], dtype=float))
    p_beta_gt_1e3 = float(np.sum(w * (beta > 1e-3).astype(float)))
    p_beta_gt_3e3 = float(np.sum(w * (beta > 3e-3).astype(float)))

    med_dbic = float(np.median(dbic))
    q75_dbic = float(np.quantile(dbic, 0.75))
    med_daic = float(np.median(daic))
    med_rmse_beta = float(np.median(rmse_b))
    med_rmse_gain = float(np.median(np.array([float(r["rmse_null"]) - float(r["rmse_beta"]) for r in rows], dtype=float)))

    pass_evidence = (med_dbic >= 6.0) and (q75_dbic >= 2.0) and (med_daic >= 4.0)
    pass_positive = (p_beta_gt_1e3 >= 0.80) and (p_beta_gt_3e3 >= 0.55)
    pass_fit = med_rmse_beta <= 0.10 and med_rmse_gain >= 0.015
    pass_nonboundary = bq[0] > 5e-4

    if pass_evidence and pass_positive and pass_fit and pass_nonboundary:
        verdict = "BETA_POSITIVE_EVIDENCE_SUPPORTED"
    elif pass_evidence and (pass_positive or pass_fit):
        verdict = "BETA_POSITIVE_EVIDENCE_PARTIAL"
    else:
        verdict = "BETA_POSITIVE_EVIDENCE_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "input_rows": len(rows),
        "weighted_beta_posterior": {
            "ci95": [float(bq[0]), float(bq[2])],
            "median": float(bq[1]),
            "prob_beta_gt_1e3": p_beta_gt_1e3,
            "prob_beta_gt_3e3": p_beta_gt_3e3,
        },
        "model_evidence": {
            "median_delta_bic_null_minus_beta": med_dbic,
            "q75_delta_bic_null_minus_beta": q75_dbic,
            "median_delta_aic_null_minus_beta": med_daic,
            "median_rmse_beta": med_rmse_beta,
            "median_rmse_gain_null_minus_beta": med_rmse_gain,
        },
        "pass_flags": {
            "evidence_strength": bool(pass_evidence),
            "positive_beta_probability": bool(pass_positive),
            "fit_improvement": bool(pass_fit),
            "nonboundary_beta": bool(pass_nonboundary),
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1755: BETA NULL VS POSITIVE EVIDENCE",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Input rows: {len(rows)}",
        f"- beta median/CI95: {bq[1]:.6f} / [{bq[0]:.6f}, {bq[2]:.6f}]",
        f"- median ΔBIC(null-beta): {med_dbic:.3f}",
        f"- median ΔAIC(null-beta): {med_daic:.3f}",
        f"- P(beta>1e-3): {p_beta_gt_1e3:.3f}",
        f"- P(beta>3e-3): {p_beta_gt_3e3:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- evidence_strength: {pass_evidence}",
        f"- positive_beta_probability: {pass_positive}",
        f"- fit_improvement: {pass_fit}",
        f"- nonboundary_beta: {pass_nonboundary}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1755] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1755] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
