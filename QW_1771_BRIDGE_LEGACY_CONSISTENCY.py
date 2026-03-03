#!/usr/bin/env python3
"""
QW-1771: Bridge-to-legacy consistency analysis.

Goal:
- Compare beta_bridge distribution from QW-1770 with legacy envelope beta
  distribution (QW-1755 / QW-1749 lineage).
- Test whether a simple monotonic compression map can align scales without
  destroying ordering.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1771_bridge_legacy_consistency.json"
OUT_MD = ROOT / "RAPORT_QW1771_BRIDGE_LEGACY_CONSISTENCY.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def wasserstein_1d(a: np.ndarray, b: np.ndarray, qn: int = 400) -> float:
    qa = np.quantile(a, np.linspace(0.0, 1.0, qn))
    qb = np.quantile(b, np.linspace(0.0, 1.0, qn))
    return float(np.mean(np.abs(qa - qb)))


def ci95(x: np.ndarray) -> Tuple[float, float]:
    return float(np.quantile(x, 0.025)), float(np.quantile(x, 0.975))


def ci_overlap(ci1: Tuple[float, float], ci2: Tuple[float, float]) -> float:
    l = max(ci1[0], ci2[0])
    r = min(ci1[1], ci2[1])
    if r <= l:
        return 0.0
    den = max(max(ci1[1], ci2[1]) - min(ci1[0], ci2[0]), 1e-12)
    return float((r - l) / den)


def rankdata(x: np.ndarray) -> np.ndarray:
    idx = np.argsort(x)
    out = np.empty_like(idx, dtype=float)
    out[idx] = np.arange(len(x), dtype=float)
    return out


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return float("nan")
    rx = rankdata(x)
    ry = rankdata(y)
    c = np.corrcoef(rx, ry)
    return float(c[0, 1])


def transform_beta(beta: np.ndarray, alpha: float, k: float) -> np.ndarray:
    return alpha * beta / (1.0 + k * beta)


def main() -> None:
    d1770 = load("report_qw1770_kernel_bridge_from_nonenvelope.json")
    d1755 = load("report_qw1755_beta_null_vs_positive_evidence.json")
    d1749 = load("report_qw1749_beta_orthogonal_observable.json")

    rows70 = d1770.get("rows", [])
    pred70 = np.array(d1770.get("predictions", {}).get("beta_bridge_pred", []), dtype=float)
    if len(rows70) == 0 or len(pred70) != len(rows70):
        raise RuntimeError("QW-1770 rows/predictions unavailable or inconsistent.")

    # Use intervention-positive subset for bridge scale.
    btrue = np.array([float(r["beta_bridge_true"]) for r in rows70], dtype=float)
    mask_pos = btrue > 0.02
    bridge = pred70[mask_pos]
    if len(bridge) < 12:
        bridge = pred70.copy()

    # Legacy distribution (preferred: 1755 refit; fallback: 1749)
    r55 = d1755.get("rows", [])
    legacy = np.array([float(r["beta_hat"]) for r in r55], dtype=float) if r55 else np.array([], dtype=float)
    if len(legacy) < 12:
        r49 = d1749.get("rows", [])
        legacy = np.array([float(r["beta_hat"]) for r in r49], dtype=float)
    if len(legacy) < 12:
        raise RuntimeError("Legacy beta distribution too short.")

    # Direct comparison
    wd_direct = wasserstein_1d(bridge, legacy, qn=400)
    med_bridge = float(np.median(bridge))
    med_legacy = float(np.median(legacy))
    ci_bridge = ci95(bridge)
    ci_legacy = ci95(legacy)
    overlap_direct = ci_overlap(ci_bridge, ci_legacy)

    # Monotonic compression map search
    alpha_grid = np.linspace(0.30, 1.20, 120)
    k_grid = np.linspace(0.0, 400.0, 260)
    best = {
        "alpha": 1.0,
        "k": 0.0,
        "wd": float("inf"),
        "med_diff": float("inf"),
        "obj": float("inf"),
    }
    for alpha in alpha_grid:
        for k in k_grid:
            bt = transform_beta(bridge, float(alpha), float(k))
            wd = wasserstein_1d(bt, legacy, qn=250)
            md = abs(float(np.median(bt)) - med_legacy)
            obj = wd + 0.35 * md
            if obj < best["obj"]:
                best = {"alpha": float(alpha), "k": float(k), "wd": float(wd), "med_diff": float(md), "obj": float(obj)}

    bridge_t = transform_beta(bridge, best["alpha"], best["k"])
    wd_t = wasserstein_1d(bridge_t, legacy, qn=500)
    ci_t = ci95(bridge_t)
    overlap_t = ci_overlap(ci_t, ci_legacy)
    med_t = float(np.median(bridge_t))
    red = wd_direct / max(wd_t, 1e-12)
    rho_mon = spearman(bridge, bridge_t)

    # Consistency decisions
    pass_direct = wd_direct <= 0.010 and overlap_direct >= 0.30
    pass_transformed = wd_t <= 0.004 and overlap_t >= 0.45 and red >= 2.0
    pass_map_simple = best["k"] <= 260.0 and 0.35 <= best["alpha"] <= 1.20
    pass_monotonic = np.isfinite(rho_mon) and rho_mon >= 0.995

    if pass_direct:
        verdict = "BRIDGE_LEGACY_CONSISTENCY_DIRECT_SUPPORTED"
    elif pass_transformed and pass_map_simple and pass_monotonic:
        verdict = "BRIDGE_LEGACY_CONSISTENCY_TRANSFORMED_PARTIAL"
    else:
        verdict = "BRIDGE_LEGACY_CONSISTENCY_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_bridge": int(len(bridge)),
        "n_legacy": int(len(legacy)),
        "direct": {
            "wasserstein": wd_direct,
            "median_bridge": med_bridge,
            "median_legacy": med_legacy,
            "ci95_bridge": [ci_bridge[0], ci_bridge[1]],
            "ci95_legacy": [ci_legacy[0], ci_legacy[1]],
            "ci_overlap": overlap_direct,
        },
        "best_monotonic_transform": {
            "alpha": best["alpha"],
            "k": best["k"],
            "wasserstein": wd_t,
            "distance_reduction_factor": red,
            "median_transformed": med_t,
            "ci95_transformed": [ci_t[0], ci_t[1]],
            "ci_overlap_with_legacy": overlap_t,
            "spearman_bridge_vs_transformed": rho_mon,
        },
        "pass_flags": {
            "direct_consistency": bool(pass_direct),
            "transformed_consistency": bool(pass_transformed),
            "simple_map": bool(pass_map_simple),
            "monotonicity": bool(pass_monotonic),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1771: BRIDGE LEGACY CONSISTENCY",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- n_bridge / n_legacy: {len(bridge)} / {len(legacy)}",
        f"- Direct WD: {wd_direct:.6f}, direct CI overlap: {overlap_direct:.3f}",
        f"- Best transform: alpha={best['alpha']:.3f}, k={best['k']:.3f}",
        f"- Transformed WD: {wd_t:.6f}, reduction x{red:.2f}, CI overlap: {overlap_t:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- direct_consistency: {pass_direct}",
        f"- transformed_consistency: {pass_transformed}",
        f"- simple_map: {pass_map_simple}",
        f"- monotonicity: {pass_monotonic}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1771] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1771] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
