#!/usr/bin/env python3
"""
QW-1826: GW shared-vs-control identifiability diagnostic.

Uses QW-1726 simulation table to test whether GW branch has identifiable
signal separation between:
- shared background hypothesis,
- unshared control background.

Focus: robust, distribution-light metrics (effect size, quantile overlap,
near-target prevalence advantage).
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1826_gw_shared_control_identifiability.json"
OUT_MD = ROOT / "RAPORT_QW1826_GW_SHARED_CONTROL_IDENTIFIABILITY.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def pooled_sd(a: float, b: float) -> float:
    return float(np.sqrt(max((a * a + b * b) / 2.0, 1e-12)))


def main() -> None:
    d1726 = load("report_qw1726_gw_fin_projection_retest.json")
    rows = d1726.get("rows", [])
    if not rows:
        raise RuntimeError("No rows in QW-1726 report.")

    out_rows: List[Dict[str, float]] = []
    for r in rows:
        snr = float(r["snr_local_over_shared"])
        s = r["shared_background"]
        u = r["unshared_background_control"]

        mu_s = float(s["mean"])
        sd_s = float(s["std"])
        p10_s = float(s["p10"])
        p90_s = float(s["p90"])
        near_s = float(s["prob_near_031_pm_002"])

        mu_u = float(u["mean"])
        sd_u = float(u["std"])
        p10_u = float(u["p10"])
        p90_u = float(u["p90"])
        near_u = float(u["prob_near_031_pm_002"])

        d_eff = (mu_s - mu_u) / pooled_sd(sd_s, sd_u)

        # Quantile non-overlap fraction proxy (0 = full overlap tendency)
        interval_left = max(p10_s, p10_u)
        interval_right = min(p90_s, p90_u)
        overlap_len = max(0.0, interval_right - interval_left)
        union_left = min(p10_s, p10_u)
        union_right = max(p90_s, p90_u)
        union_len = max(union_right - union_left, 1e-12)
        nonoverlap_frac = 1.0 - overlap_len / union_len

        # Target prevalence advantage of shared over control
        near_adv = near_s - near_u

        out_rows.append(
            {
                "snr": snr,
                "mean_shared": mu_s,
                "mean_control": mu_u,
                "delta_mean": mu_s - mu_u,
                "effect_size_d": float(d_eff),
                "nonoverlap_frac_p10_p90": float(nonoverlap_frac),
                "near_target_advantage": float(near_adv),
            }
        )

    arr_d = np.array([r["effect_size_d"] for r in out_rows], dtype=float)
    arr_abs_d = np.abs(arr_d)
    arr_nonov = np.array([r["nonoverlap_frac_p10_p90"] for r in out_rows], dtype=float)
    arr_near = np.array([r["near_target_advantage"] for r in out_rows], dtype=float)

    summary = {
        "n_snr_rows": len(out_rows),
        "mean_effect_size_d": float(np.mean(arr_d)),
        "mean_abs_effect_size_d": float(np.mean(arr_abs_d)),
        "max_abs_effect_size_d": float(np.max(arr_abs_d)),
        "mean_nonoverlap_frac": float(np.mean(arr_nonov)),
        "min_nonoverlap_frac": float(np.min(arr_nonov)),
        "mean_near_target_advantage": float(np.mean(arr_near)),
        "max_near_target_advantage": float(np.max(arr_near)),
        "prob_shared_beats_control_near_target": float(np.mean(arr_near > 0.0)),
    }

    # Identifiability thresholds: moderate and explicit
    pass_effect = summary["mean_abs_effect_size_d"] >= 0.35
    pass_nonov = summary["mean_nonoverlap_frac"] >= 0.35
    pass_target_adv = summary["prob_shared_beats_control_near_target"] >= 0.70

    if pass_effect and pass_nonov and pass_target_adv:
        verdict = "GW_SHARED_CONTROL_IDENTIFIABLE"
    elif pass_effect and (pass_nonov or pass_target_adv):
        verdict = "GW_SHARED_CONTROL_PARTIAL"
    else:
        verdict = "GW_SHARED_CONTROL_NOT_IDENTIFIABLE"

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "summary": summary,
        "pass_flags": {
            "effect_separation": bool(pass_effect),
            "quantile_nonoverlap": bool(pass_nonov),
            "target_prevalence_advantage": bool(pass_target_adv),
        },
        "verdict": verdict,
        "rows": out_rows,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1826: GW SHARED-CONTROL IDENTIFIABILITY",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- n SNR rows: {summary['n_snr_rows']}",
        f"- Mean |effect size d|: {summary['mean_abs_effect_size_d']:.4f}",
        f"- Mean quantile non-overlap (p10/p90): {summary['mean_nonoverlap_frac']:.4f}",
        f"- P(shared beats control near target): {summary['prob_shared_beats_control_near_target']:.3f}",
        f"- Mean near-target advantage: {summary['mean_near_target_advantage']:.4f}",
        f"- Verdict: **{verdict}**",
        "",
        "## Pass Flags",
        f"- effect_separation: {pass_effect}",
        f"- quantile_nonoverlap: {pass_nonov}",
        f"- target_prevalence_advantage: {pass_target_adv}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1826] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1826] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
