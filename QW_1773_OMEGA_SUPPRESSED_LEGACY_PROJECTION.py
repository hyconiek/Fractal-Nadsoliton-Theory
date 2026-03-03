#!/usr/bin/env python3
"""
QW-1773: Omega-suppressed legacy projection.

Projection hypothesis:
    beta_legacy_equiv = alpha * beta_bridge * exp(-gamma / omega^p)

Intent:
- Reparameterize bridge outputs to legacy beta scale using a single
  low-complexity mechanistic suppressor controlled by omega.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1773_omega_suppressed_legacy_projection.json"
OUT_MD = ROOT / "RAPORT_QW1773_OMEGA_SUPPRESSED_LEGACY_PROJECTION.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def wasserstein_1d(a: np.ndarray, b: np.ndarray, qn: int = 500) -> float:
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


def project(beta_bridge: np.ndarray, omega: np.ndarray, alpha: float, gamma: float, p: float) -> np.ndarray:
    om = np.clip(omega, 1e-6, None)
    out = alpha * beta_bridge * np.exp(-gamma / (om ** p))
    return np.clip(out, 0.0, 0.35)


def main() -> None:
    d1770 = load("report_qw1770_kernel_bridge_from_nonenvelope.json")
    d1755 = load("report_qw1755_beta_null_vs_positive_evidence.json")

    rows70 = d1770.get("rows", [])
    pb = np.array(d1770.get("predictions", {}).get("beta_bridge_pred", []), dtype=float)
    po = np.array(d1770.get("predictions", {}).get("omega_pred", []), dtype=float)
    if len(rows70) == 0 or len(pb) != len(rows70) or len(po) != len(rows70):
        raise RuntimeError("QW-1770 inputs inconsistent.")

    btrue = np.array([float(r["beta_bridge_true"]) for r in rows70], dtype=float)
    mask = btrue > 0.02
    if np.sum(mask) < 12:
        mask = np.ones_like(btrue, dtype=bool)
    bridge = pb[mask]
    omega = po[mask]

    r55 = d1755.get("rows", [])
    legacy = np.array([float(r["beta_hat"]) for r in r55], dtype=float)
    if len(legacy) < 12:
        raise RuntimeError("Legacy beta rows too short.")

    wd_direct = wasserstein_1d(bridge, legacy, qn=500)
    ci_leg = ci95(legacy)
    ci_bridge = ci95(bridge)
    overlap_direct = ci_overlap(ci_bridge, ci_leg)

    best = {"obj": float("inf"), "alpha": 1.0, "gamma": 0.0, "p": 1.0, "wd": float("inf"), "med_diff": float("inf")}
    leg_med = float(np.median(legacy))
    rng = np.random.default_rng(1773)

    def eval_candidate(alpha: float, gamma: float, p: float) -> Dict[str, float]:
        bt = project(bridge, omega, float(alpha), float(gamma), float(p))
        wd = wasserstein_1d(bt, legacy, qn=140)
        med_diff = abs(float(np.median(bt)) - leg_med)
        reg = 0.04 * (alpha - 1.0) ** 2 + 0.012 * gamma + 0.02 * abs(p - 1.5)
        obj = wd + 0.30 * med_diff + reg
        return {"obj": float(obj), "alpha": float(alpha), "gamma": float(gamma), "p": float(p), "wd": float(wd), "med_diff": float(med_diff)}

    # Stage 1: random coarse exploration
    for _ in range(7000):
        alpha = float(rng.uniform(0.20, 1.20))
        gamma = float(rng.uniform(0.0, 4.0))
        p = float(rng.uniform(0.8, 2.4))
        cand = eval_candidate(alpha, gamma, p)
        if cand["obj"] < best["obj"]:
            best = cand

    # Stage 2: local refinement around best
    a0, g0, p0 = best["alpha"], best["gamma"], best["p"]
    alpha_grid = np.linspace(max(0.20, a0 - 0.22), min(1.20, a0 + 0.22), 46)
    gamma_grid = np.linspace(max(0.0, g0 - 0.9), min(4.0, g0 + 0.9), 46)
    p_grid = np.linspace(max(0.8, p0 - 0.45), min(2.4, p0 + 0.45), 46)
    for alpha in alpha_grid:
        for gamma in gamma_grid:
            for p in p_grid:
                cand = eval_candidate(float(alpha), float(gamma), float(p))
                if cand["obj"] < best["obj"]:
                    best = cand

    proj = project(bridge, omega, best["alpha"], best["gamma"], best["p"])
    wd_proj = wasserstein_1d(proj, legacy, qn=700)
    red = wd_direct / max(wd_proj, 1e-12)
    ci_proj = ci95(proj)
    overlap_proj = ci_overlap(ci_proj, ci_leg)
    med_proj = float(np.median(proj))

    # Additional checks
    proj_nonboundary_rate = float(np.mean((proj > 5e-4) & (proj < 0.23)))
    proj_upper_tail = float(np.quantile(proj, 0.95))

    pass_reduction = wd_proj <= 0.010 and red >= 8.0
    pass_overlap = overlap_proj >= 0.50
    pass_parsimony = (0.25 <= best["alpha"] <= 1.20) and (best["gamma"] <= 1.5) and (0.9 <= best["p"] <= 2.2)
    pass_projection_shape = proj_nonboundary_rate >= 0.70 and proj_upper_tail <= 0.06

    if pass_reduction and pass_overlap and pass_parsimony and pass_projection_shape:
        verdict = "OMEGA_SUPPRESSED_LEGACY_PROJECTION_SUPPORTED"
    elif pass_reduction and (pass_overlap or pass_parsimony):
        verdict = "OMEGA_SUPPRESSED_LEGACY_PROJECTION_PARTIAL"
    else:
        verdict = "OMEGA_SUPPRESSED_LEGACY_PROJECTION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "n_bridge": int(len(bridge)),
        "n_legacy": int(len(legacy)),
        "direct": {
            "wasserstein": wd_direct,
            "median_bridge": float(np.median(bridge)),
            "median_legacy": leg_med,
            "ci95_bridge": [ci_bridge[0], ci_bridge[1]],
            "ci95_legacy": [ci_leg[0], ci_leg[1]],
            "ci_overlap": overlap_direct,
        },
        "projection": {
            "alpha": best["alpha"],
            "gamma": best["gamma"],
            "p": best["p"],
            "wasserstein": wd_proj,
            "distance_reduction_factor": red,
            "median_projected": med_proj,
            "ci95_projected": [ci_proj[0], ci_proj[1]],
            "ci_overlap_with_legacy": overlap_proj,
            "projected_nonboundary_rate": proj_nonboundary_rate,
            "projected_q95": proj_upper_tail,
        },
        "pass_flags": {
            "distance_reduction": bool(pass_reduction),
            "ci_overlap": bool(pass_overlap),
            "parsimony": bool(pass_parsimony),
            "projection_shape": bool(pass_projection_shape),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1773: OMEGA SUPPRESSED LEGACY PROJECTION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- n_bridge / n_legacy: {len(bridge)} / {len(legacy)}",
        f"- Direct WD: {wd_direct:.6f}, overlap: {overlap_direct:.3f}",
        f"- Projection params: alpha={best['alpha']:.3f}, gamma={best['gamma']:.3f}, p={best['p']:.3f}",
        f"- Projected WD: {wd_proj:.6f}, reduction x{red:.2f}, overlap: {overlap_proj:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- distance_reduction: {pass_reduction}",
        f"- ci_overlap: {pass_overlap}",
        f"- parsimony: {pass_parsimony}",
        f"- projection_shape: {pass_projection_shape}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1773] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1773] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
