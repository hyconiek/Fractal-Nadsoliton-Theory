#!/usr/bin/env python3
"""
QW-1775: OOD validation for omega-suppressed reparameterization.

Validates QW-1773 mapping parameters (alpha, gamma, p) on an independent
cohort with shifted network size and drive ranges.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Tuple

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1775_reparam_ood_validation.json"
OUT_MD = ROOT / "RAPORT_QW1775_REPARAM_OOD_VALIDATION.md"


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
    d1773 = load("report_qw1773_omega_suppressed_legacy_projection.json")
    d1755 = load("report_qw1755_beta_null_vs_positive_evidence.json")

    pr = d1773.get("projection", {})
    alpha = float(pr.get("alpha"))
    gamma = float(pr.get("gamma"))
    pexp = float(pr.get("p"))

    legacy = np.array([float(r["beta_hat"]) for r in d1755.get("rows", [])], dtype=float)
    if len(legacy) < 12:
        raise RuntimeError("Legacy beta rows too short.")
    ci_leg = ci95(legacy)

    # Independent OOD cohort (shifted ranges)
    rng = np.random.default_rng(1775)
    n_cases = 80
    beta_bridge_true = []
    omega_true = []
    for i in range(n_cases):
        # shifted bridge beta support and omega support vs 1770 generation
        b = 0.0 if i % 5 == 0 else float(rng.uniform(0.05, 0.24))
        if i % 7 == 0:
            b = float(rng.uniform(0.0, 0.04))
        om = float(rng.uniform(0.24, 0.48))
        if i % 9 == 0:
            om = float(rng.uniform(0.14, 0.24))
        beta_bridge_true.append(b)
        omega_true.append(om)

    beta_bridge_true = np.array(beta_bridge_true, dtype=float)
    omega_true = np.array(omega_true, dtype=float)

    beta_proj = project(beta_bridge_true, omega_true, alpha, gamma, pexp)
    wd = wasserstein_1d(beta_proj, legacy, qn=500)
    ci_proj = ci95(beta_proj)
    ov = ci_overlap(ci_proj, ci_leg)
    med_proj = float(np.median(beta_proj))
    q95_proj = float(np.quantile(beta_proj, 0.95))
    nonboundary = float(np.mean((beta_proj > 5e-4) & (beta_proj < 0.23)))

    # Parameter sensitivity around projection
    deltas = [-0.10, -0.05, 0.0, 0.05, 0.10]
    wd_sens = []
    ov_sens = []
    for da in deltas:
        for dg in deltas:
            for dp in deltas:
                aa = alpha * (1.0 + da)
                gg = gamma * (1.0 + dg)
                pp = pexp * (1.0 + dp)
                bp = project(beta_bridge_true, omega_true, aa, gg, pp)
                wd_sens.append(wasserstein_1d(bp, legacy, qn=220))
                ov_sens.append(ci_overlap(ci95(bp), ci_leg))
    wd_sens = np.array(wd_sens, dtype=float)
    ov_sens = np.array(ov_sens, dtype=float)
    wd_q90 = float(np.quantile(wd_sens, 0.90))
    ov_q10 = float(np.quantile(ov_sens, 0.10))

    pass_core = wd <= 0.010 and ov >= 0.50
    pass_shape = nonboundary >= 0.70 and q95_proj <= 0.06
    pass_sens = wd_q90 <= 0.016 and ov_q10 >= 0.32

    if pass_core and pass_shape and pass_sens:
        verdict = "REPARAM_OOD_VALIDATION_SUPPORTED"
    elif pass_core and (pass_shape or pass_sens):
        verdict = "REPARAM_OOD_VALIDATION_PARTIAL"
    else:
        verdict = "REPARAM_OOD_VALIDATION_WEAK"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "projection_params_from_1773": {"alpha": alpha, "gamma": gamma, "p": pexp},
        "ood_cohort": {
            "n_cases": n_cases,
            "beta_bridge_true_ci95": [float(np.quantile(beta_bridge_true, 0.025)), float(np.quantile(beta_bridge_true, 0.975))],
            "omega_true_ci95": [float(np.quantile(omega_true, 0.025)), float(np.quantile(omega_true, 0.975))],
        },
        "summary": {
            "wasserstein_vs_legacy": wd,
            "ci95_projected": [ci_proj[0], ci_proj[1]],
            "ci95_legacy": [ci_leg[0], ci_leg[1]],
            "ci_overlap": ov,
            "median_projected": med_proj,
            "q95_projected": q95_proj,
            "nonboundary_rate_projected": nonboundary,
            "sensitivity_wd_q90": wd_q90,
            "sensitivity_overlap_q10": ov_q10,
        },
        "pass_flags": {
            "core_alignment": bool(pass_core),
            "projection_shape": bool(pass_shape),
            "param_sensitivity": bool(pass_sens),
        },
        "verdict": verdict,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1775: REPARAM OOD VALIDATION",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- OOD cases: {n_cases}",
        f"- WD vs legacy: {wd:.6f}, CI overlap: {ov:.3f}",
        f"- Projected median/q95: {med_proj:.6f} / {q95_proj:.6f}",
        f"- Sensitivity WD q90 / overlap q10: {wd_q90:.6f} / {ov_q10:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- core_alignment: {pass_core}",
        f"- projection_shape: {pass_shape}",
        f"- param_sensitivity: {pass_sens}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1775] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1775] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
