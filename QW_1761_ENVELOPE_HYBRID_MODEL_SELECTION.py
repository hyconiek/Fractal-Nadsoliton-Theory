#!/usr/bin/env python3
"""
QW-1761: Envelope hybrid model selection.

Models:
- M0: constant envelope (null)
- M1: hyperbolic envelope A/(1+beta*d)
- M2: hybrid envelope A*exp(-lambda*d)/(1+beta*d)

Goal:
- Test whether adding a memory-like damping term (lambda) explains
  the persistent weak evidence for positive beta_tors.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1761_envelope_hybrid_model_selection.json"
OUT_MD = ROOT / "RAPORT_QW1761_ENVELOPE_HYBRID_MODEL_SELECTION.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def fit_m0_constant(y: np.ndarray) -> Dict[str, float]:
    a0 = float(np.mean(y))
    pred = np.full_like(y, a0)
    sse = float(np.sum((y - pred) ** 2))
    rmse = float(np.sqrt(np.mean((y - pred) ** 2)))
    return {"a": a0, "sse": sse, "rmse": rmse}


def fit_m1_beta(y: np.ndarray) -> Dict[str, float]:
    d = np.arange(1, len(y) + 1, dtype=float)
    beta_grid = np.logspace(np.log10(1e-5), np.log10(0.3), 1800)
    best = {"beta": 1e-3, "a": float(np.mean(y)), "sse": float("inf"), "rmse": float("inf")}
    for b in beta_grid:
        f = 1.0 / (1.0 + b * d)
        den = float(np.dot(f, f))
        if den <= 1e-15:
            continue
        a = float(np.dot(y, f) / den)
        pred = a * f
        sse = float(np.sum((y - pred) ** 2))
        if sse < best["sse"]:
            best = {
                "beta": float(b),
                "a": a,
                "sse": sse,
                "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
            }
    return best


def fit_m2_beta_lambda(y: np.ndarray) -> Dict[str, float]:
    d = np.arange(1, len(y) + 1, dtype=float)

    beta_coarse = np.logspace(np.log10(1e-5), np.log10(0.3), 280)
    lam_coarse = np.linspace(0.0, 0.25, 180)

    best = {
        "beta": 1e-3,
        "lam": 0.01,
        "a": float(np.mean(y)),
        "sse": float("inf"),
        "rmse": float("inf"),
    }

    for b in beta_coarse:
        inv = 1.0 / (1.0 + b * d)
        for lam in lam_coarse:
            f = np.exp(-lam * d) * inv
            den = float(np.dot(f, f))
            if den <= 1e-15:
                continue
            a = float(np.dot(y, f) / den)
            pred = a * f
            sse = float(np.sum((y - pred) ** 2))
            if sse < best["sse"]:
                best = {
                    "beta": float(b),
                    "lam": float(lam),
                    "a": a,
                    "sse": sse,
                    "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
                }

    b0 = best["beta"]
    l0 = best["lam"]
    b_lo = max(1e-5, b0 * 0.45)
    b_hi = min(0.35, b0 * 1.8 + 1e-6)
    l_lo = max(0.0, l0 - 0.08)
    l_hi = min(0.35, l0 + 0.08)

    beta_ref = np.linspace(b_lo, b_hi, 160)
    lam_ref = np.linspace(l_lo, l_hi, 140)
    for b in beta_ref:
        inv = 1.0 / (1.0 + b * d)
        for lam in lam_ref:
            f = np.exp(-lam * d) * inv
            den = float(np.dot(f, f))
            if den <= 1e-15:
                continue
            a = float(np.dot(y, f) / den)
            pred = a * f
            sse = float(np.sum((y - pred) ** 2))
            if sse < best["sse"]:
                best = {
                    "beta": float(b),
                    "lam": float(lam),
                    "a": a,
                    "sse": sse,
                    "rmse": float(np.sqrt(np.mean((y - pred) ** 2))),
                }

    return best


def bic_from_sse(sse: float, n: int, k: int) -> float:
    s = max(float(sse), 1e-15)
    nn = max(int(n), 1)
    return float(nn * math.log(s / nn) + k * math.log(nn))


def main() -> None:
    d1749 = load("report_qw1749_beta_orthogonal_observable.json")
    rows_in = d1749.get("rows", [])
    if len(rows_in) < 12:
        raise RuntimeError("QW-1749 rows not sufficient for model selection.")

    rows: List[Dict[str, object]] = []
    for r in rows_in:
        ymap = r.get("tail_observable", {})
        y = np.array([float(ymap[str(i)]) for i in range(1, len(ymap) + 1)], dtype=float)
        if len(y) < 8:
            continue
        if y[0] > 1e-12:
            y = y / y[0]

        m0 = fit_m0_constant(y)
        m1 = fit_m1_beta(y)
        m2 = fit_m2_beta_lambda(y)

        bic0 = bic_from_sse(m0["sse"], len(y), 1)
        bic1 = bic_from_sse(m1["sse"], len(y), 2)
        bic2 = bic_from_sse(m2["sse"], len(y), 3)

        bics = {"M0_null": bic0, "M1_beta": bic1, "M2_hybrid": bic2}
        winner = min(bics, key=bics.get)

        rows.append(
            {
                "seed": int(r.get("seed", -1)),
                "n": int(r.get("n", -1)),
                "m0_rmse": float(m0["rmse"]),
                "m1_beta": float(m1["beta"]),
                "m1_rmse": float(m1["rmse"]),
                "m2_beta": float(m2["beta"]),
                "m2_lambda": float(m2["lam"]),
                "m2_rmse": float(m2["rmse"]),
                "bic": bics,
                "winner": winner,
                "delta_bic_m1_minus_m2": float(bic1 - bic2),
                "delta_bic_m0_minus_m2": float(bic0 - bic2),
            }
        )

    if len(rows) < 12:
        raise RuntimeError("Too few valid rows after model fitting.")

    win_counts = {
        "M0_null": int(sum(1 for r in rows if r["winner"] == "M0_null")),
        "M1_beta": int(sum(1 for r in rows if r["winner"] == "M1_beta")),
        "M2_hybrid": int(sum(1 for r in rows if r["winner"] == "M2_hybrid")),
    }
    hybrid_rate = float(win_counts["M2_hybrid"] / len(rows))

    db12 = np.array([float(r["delta_bic_m1_minus_m2"]) for r in rows], dtype=float)
    db02 = np.array([float(r["delta_bic_m0_minus_m2"]) for r in rows], dtype=float)
    b2 = np.array([float(r["m2_beta"]) for r in rows], dtype=float)
    l2 = np.array([float(r["m2_lambda"]) for r in rows], dtype=float)

    b2_ci = [float(np.quantile(b2, 0.025)), float(np.quantile(b2, 0.975))]
    l2_ci = [float(np.quantile(l2, 0.025)), float(np.quantile(l2, 0.975))]
    b2_nonboundary_rate = float(np.mean((b2 > 5e-4) & (b2 < 0.23)))
    l2_positive_rate = float(np.mean(l2 > 0.005))

    med_db12 = float(np.median(db12))
    med_db02 = float(np.median(db02))
    q75_db12 = float(np.quantile(db12, 0.75))

    pass_pref = hybrid_rate >= 0.55
    pass_gain_m1 = med_db12 >= 2.0 and q75_db12 >= 0.5
    pass_gain_m0 = med_db02 >= 2.0
    pass_beta_nonboundary = b2_nonboundary_rate >= 0.60
    pass_lambda = l2_positive_rate >= 0.70

    if pass_pref and pass_gain_m1 and pass_gain_m0 and pass_beta_nonboundary and pass_lambda:
        verdict = "HYBRID_ENVELOPE_MECHANISM_SUPPORTED"
    elif pass_pref and (pass_gain_m1 or pass_gain_m0):
        verdict = "HYBRID_ENVELOPE_MECHANISM_PARTIAL"
    else:
        verdict = "HYBRID_ENVELOPE_MECHANISM_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_rows": len(rows),
        "wins": win_counts,
        "hybrid_preferred_rate": hybrid_rate,
        "summary": {
            "median_delta_bic_m1_minus_m2": med_db12,
            "q75_delta_bic_m1_minus_m2": q75_db12,
            "median_delta_bic_m0_minus_m2": med_db02,
            "m2_beta_median": float(np.median(b2)),
            "m2_beta_ci95_empirical": b2_ci,
            "m2_lambda_median": float(np.median(l2)),
            "m2_lambda_ci95_empirical": l2_ci,
            "m2_beta_nonboundary_rate": b2_nonboundary_rate,
            "m2_lambda_positive_rate": l2_positive_rate,
        },
        "pass_flags": {
            "hybrid_preferred_rate": bool(pass_pref),
            "hybrid_gain_vs_m1": bool(pass_gain_m1),
            "hybrid_gain_vs_m0": bool(pass_gain_m0),
            "m2_beta_nonboundary": bool(pass_beta_nonboundary),
            "m2_lambda_positive": bool(pass_lambda),
        },
        "verdict": verdict,
        "rows": rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1761: ENVELOPE HYBRID MODEL SELECTION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Rows: {len(rows)}",
        f"- Wins: {win_counts}",
        f"- Hybrid preferred rate: {hybrid_rate:.3f}",
        f"- median ΔBIC(M1-M2): {med_db12:.3f}",
        f"- median ΔBIC(M0-M2): {med_db02:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- hybrid_preferred_rate: {pass_pref}",
        f"- hybrid_gain_vs_m1: {pass_gain_m1}",
        f"- hybrid_gain_vs_m0: {pass_gain_m0}",
        f"- m2_beta_nonboundary: {pass_beta_nonboundary}",
        f"- m2_lambda_positive: {pass_lambda}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1761] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1761] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
