#!/usr/bin/env python3
"""
QW-1821: Likelihood calibration diagnostic for sequence branch.

Assesses metric discordance seen in QW-1817/1818:
- RMSE gain often positive,
- test log-likelihood gain unstable.

Goal:
- quantify whether instability is primarily likelihood-calibration issue
  (variance/tail modeling), not mean-model failure.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1821_likelihood_calibration_diagnostic.json"
OUT_MD = ROOT / "RAPORT_QW1821_LIKELIHOOD_CALIBRATION_DIAGNOSTIC.md"


def load(name: str) -> Dict:
    return json.loads((ROOT / name).read_text(encoding="utf-8"))


def summarize(rep_rows: List[Dict], delta_key: str, rmse_key: str) -> Dict[str, float]:
    delta = np.array([r[delta_key] for r in rep_rows], dtype=float)
    rmse = np.array([r[rmse_key] for r in rep_rows], dtype=float)

    if len(delta) == 0:
        raise RuntimeError("No replication rows found for diagnostic.")

    discord = (rmse > 0.0) & (delta <= 0.0)
    corr = float(np.corrcoef(delta, rmse)[0, 1]) if len(delta) > 2 and np.std(delta) > 1e-12 and np.std(rmse) > 1e-12 else 0.0

    return {
        "n_rep": int(len(delta)),
        "mean_delta_ll": float(np.mean(delta)),
        "std_delta_ll": float(np.std(delta)),
        "prob_delta_ll_positive": float(np.mean(delta > 0.0)),
        "mean_rmse_gain": float(np.mean(rmse)),
        "std_rmse_gain": float(np.std(rmse)),
        "prob_rmse_gain_positive": float(np.mean(rmse > 0.0)),
        "discordance_rate": float(np.mean(discord)),
        "corr_delta_vs_rmse_gain": corr,
    }


def main() -> None:
    d1817 = load("report_qw1817_sequence_oos_validation.json")
    d1818 = load("report_qw1818_robust_sequence_oos.json")

    s17 = summarize(
        d1817["replications"],
        delta_key="delta_ll_m2e_vs_m2",
        rmse_key="rmse_gain_m2_minus_m2e",
    )
    s18 = summarize(
        d1818["replications"],
        delta_key="delta_ll_m2e_vs_m2",
        rmse_key="rmse_gain_m2_minus_m2e",
    )

    # Composite indicator: high RMSE reliability + nontrivial LL discordance.
    mis17 = 0.5 * s17["prob_rmse_gain_positive"] + 0.5 * s17["discordance_rate"]
    mis18 = 0.5 * s18["prob_rmse_gain_positive"] + 0.5 * s18["discordance_rate"]
    mis = float(0.5 * (mis17 + mis18))

    ll_unstable = (s17["std_delta_ll"] > 3.0) or (s18["std_delta_ll"] > 3.0)
    rmse_stable = (s17["prob_rmse_gain_positive"] >= 0.95) and (s18["prob_rmse_gain_positive"] >= 0.95)
    discord_present = (s17["discordance_rate"] >= 0.05) or (s18["discordance_rate"] >= 0.05)

    if ll_unstable and rmse_stable and discord_present:
        verdict = "LIKELIHOOD_MISSPECIFICATION_SIGNAL_STRONG"
    elif (rmse_stable and discord_present) or ll_unstable:
        verdict = "LIKELIHOOD_MISSPECIFICATION_SIGNAL_PARTIAL"
    else:
        verdict = "LIKELIHOOD_MISSPECIFICATION_SIGNAL_WEAK"

    recommendation = (
        "Switch primary gate metric to predictive score robust to variance tails "
        "(e.g., CRPS / quantile loss / Student-t with calibrated df), then retest."
    )

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "qw1817": s17,
        "qw1818": s18,
        "miscalibration_index": mis,
        "flags": {
            "ll_unstable": bool(ll_unstable),
            "rmse_stable_positive": bool(rmse_stable),
            "metric_discordance_present": bool(discord_present),
        },
        "verdict": verdict,
        "recommendation": recommendation,
    }
    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1821: LIKELIHOOD CALIBRATION DIAGNOSTIC",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Miscalibration index: {mis:.3f}",
        f"- Verdict: **{verdict}**",
        f"- Recommendation: **{recommendation}**",
        "",
        "## QW-1817",
        f"- n_rep: {s17['n_rep']}",
        f"- mean delta LL: {s17['mean_delta_ll']:.4f}",
        f"- std delta LL: {s17['std_delta_ll']:.4f}",
        f"- P(delta LL > 0): {s17['prob_delta_ll_positive']:.3f}",
        f"- mean RMSE gain: {s17['mean_rmse_gain']:.6f}",
        f"- P(RMSE gain > 0): {s17['prob_rmse_gain_positive']:.3f}",
        f"- discordance rate: {s17['discordance_rate']:.3f}",
        f"- corr(delta LL, RMSE gain): {s17['corr_delta_vs_rmse_gain']:.3f}",
        "",
        "## QW-1818",
        f"- n_rep: {s18['n_rep']}",
        f"- mean delta LL: {s18['mean_delta_ll']:.4f}",
        f"- std delta LL: {s18['std_delta_ll']:.4f}",
        f"- P(delta LL > 0): {s18['prob_delta_ll_positive']:.3f}",
        f"- mean RMSE gain: {s18['mean_rmse_gain']:.6f}",
        f"- P(RMSE gain > 0): {s18['prob_rmse_gain_positive']:.3f}",
        f"- discordance rate: {s18['discordance_rate']:.3f}",
        f"- corr(delta LL, RMSE gain): {s18['corr_delta_vs_rmse_gain']:.3f}",
        "",
        "## Flags",
        f"- ll_unstable: {ll_unstable}",
        f"- rmse_stable_positive: {rmse_stable}",
        f"- metric_discordance_present: {discord_present}",
        "",
        "## Artifacts",
        f"- JSON: `{OUT_JSON.name}`",
    ]
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1821] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1821] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
