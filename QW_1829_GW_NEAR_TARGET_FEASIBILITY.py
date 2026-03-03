#!/usr/bin/env python3
"""
QW-1829: GW near-target feasibility sweep.

Quantifies how much mean shift and variance reduction are required for the GW
shared model to meet critical QW-1828 target:
  P(|H - H_target| <= band) >= 0.7
using Gaussian approximation from QW-1726 row statistics.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np
from scipy.stats import norm


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1829_gw_near_target_feasibility.json"
OUT_MD = ROOT / "RAPORT_QW1829_GW_NEAR_TARGET_FEASIBILITY.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def p_in_band(mu: float, sigma: float, target: float, band: float) -> float:
    lo = (target - band - mu) / sigma
    hi = (target + band - mu) / sigma
    return float(norm.cdf(hi) - norm.cdf(lo))


def sigma_required_for_prob(prob: float, band: float) -> float:
    """
    At mu=target, probability in +/-band is max for given sigma:
      prob = 2*Phi(band/sigma)-1 -> sigma = band / Phi^{-1}((prob+1)/2)
    """
    z = float(norm.ppf((prob + 1.0) / 2.0))
    return float(band / z)


def main() -> None:
    d1726 = load("report_qw1726_gw_fin_projection_retest.json")
    rows = d1726.get("rows", [])
    if not rows:
        raise RuntimeError("No rows in QW-1726.")

    target = 0.31
    band = 0.002
    prob_req = 0.70

    sigma_req_global = sigma_required_for_prob(prob_req, band)

    out_rows: List[Dict[str, float]] = []
    for r in rows:
        snr = float(r["snr_local_over_shared"])
        s = r["shared_background"]

        mu = float(s["mean"])
        sigma = float(s["std"])
        p_now = p_in_band(mu, sigma, target=target, band=band)

        # Required shift with current sigma (best case: move mean exactly to target)
        shift_to_target = target - mu
        abs_shift = abs(shift_to_target)

        # Required sigma shrink if mean can be perfectly centered
        sigma_shrink_factor = sigma / sigma_req_global

        out_rows.append(
            {
                "snr": snr,
                "mu_shared": mu,
                "sigma_shared": sigma,
                "p_now_near_target": p_now,
                "required_mean_shift_to_target": shift_to_target,
                "required_abs_mean_shift": abs_shift,
                "sigma_required_if_centered": sigma_req_global,
                "sigma_shrink_factor_if_centered": sigma_shrink_factor,
            }
        )

    arr_p = np.array([r["p_now_near_target"] for r in out_rows], dtype=float)
    arr_shift = np.array([r["required_abs_mean_shift"] for r in out_rows], dtype=float)
    arr_shrink = np.array([r["sigma_shrink_factor_if_centered"] for r in out_rows], dtype=float)

    summary = {
        "target": target,
        "band": band,
        "required_probability": prob_req,
        "n_snr_rows": len(out_rows),
        "mean_current_p_near_target": float(np.mean(arr_p)),
        "max_current_p_near_target": float(np.max(arr_p)),
        "mean_required_abs_mean_shift": float(np.mean(arr_shift)),
        "min_required_abs_mean_shift": float(np.min(arr_shift)),
        "max_required_abs_mean_shift": float(np.max(arr_shift)),
        "sigma_required_if_centered": float(sigma_req_global),
        "mean_sigma_shrink_factor_if_centered": float(np.mean(arr_shrink)),
        "min_sigma_shrink_factor_if_centered": float(np.min(arr_shrink)),
        "max_sigma_shrink_factor_if_centered": float(np.max(arr_shrink)),
    }

    # Feasibility labels
    hard_mean_shift = summary["mean_required_abs_mean_shift"] >= 0.05
    hard_sigma_shrink = summary["mean_sigma_shrink_factor_if_centered"] >= 5.0

    if hard_mean_shift and hard_sigma_shrink:
        verdict = "GW_NEAR_TARGET_REQUIRES_STRUCTURAL_REDESIGN"
    elif hard_mean_shift or hard_sigma_shrink:
        verdict = "GW_NEAR_TARGET_REQUIRES_MAJOR_REPARAM"
    else:
        verdict = "GW_NEAR_TARGET_TUNABLE"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "flags": {
            "hard_mean_shift": bool(hard_mean_shift),
            "hard_sigma_shrink": bool(hard_sigma_shrink),
        },
        "verdict": verdict,
        "rows": out_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1829: GW NEAR-TARGET FEASIBILITY",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Target H: {target}",
        f"- Band: +/- {band}",
        f"- Required probability: {prob_req}",
        f"- Mean current P(near target): {summary['mean_current_p_near_target']:.6f}",
        f"- Mean required |mean shift|: {summary['mean_required_abs_mean_shift']:.6f}",
        f"- Required sigma if centered: {summary['sigma_required_if_centered']:.6f}",
        f"- Mean sigma shrink factor if centered: {summary['mean_sigma_shrink_factor_if_centered']:.3f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Flags",
        f"- hard_mean_shift: {hard_mean_shift}",
        f"- hard_sigma_shrink: {hard_sigma_shrink}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1829] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1829] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
