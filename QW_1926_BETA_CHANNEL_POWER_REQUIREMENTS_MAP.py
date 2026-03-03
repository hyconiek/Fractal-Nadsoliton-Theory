#!/usr/bin/env python3
"""
QW-1926: Beta-channel power requirements map under data-quality scenarios.

Builds scenario grid for required holdout size to detect beta-channel contrast
in strict one-sided confirmatory testing.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1926_beta_channel_power_requirements_map.json"
OUT_MD = ROOT / "RAPORT_QW1926_BETA_CHANNEL_POWER_REQUIREMENTS_MAP.md"


def z_quantile(p: float) -> float:
    a = [-39.69683028665376, 220.9460984245205, -275.9285104469687, 138.357751867269, -30.66479806614716, 2.506628277459239]
    b = [-54.47609879822406, 161.5858368580409, -155.6989798598866, 66.80131188771972, -13.28068155288572]
    c = [-0.007784894002430293, -0.3223964580411365, -2.400758277161838, -2.549732539343734, 4.374664141464968, 2.938163982698783]
    d = [0.007784695709041462, 0.3224671290700398, 2.445134137142996, 3.754408661907416]
    plow = 0.02425
    phigh = 1.0 - plow

    if p <= 0.0:
        return -1e9
    if p >= 1.0:
        return 1e9

    if p < plow:
        q = np.sqrt(-2.0 * np.log(p))
        return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)
    if p > phigh:
        q = np.sqrt(-2.0 * np.log(1.0 - p))
        return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q + c[5]) / ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0)

    q = p - 0.5
    r = q * q
    return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) * q / (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0)


def required_n(effect: float, sigma: float, alpha_one_sided: float, power: float) -> int:
    z_a = z_quantile(1.0 - alpha_one_sided)
    z_b = z_quantile(power)
    n = ((z_a + z_b) * sigma / max(effect, 1e-9)) ** 2
    return int(np.ceil(n))


def main() -> None:
    d1922 = json.loads((ROOT / "report_qw1922_beta_observable_blind_external_intervention.json").read_text(encoding="utf-8"))

    p_hold = d1922["datasets"]["primary"]["holdout"]
    s_hold = d1922["datasets"]["stress"]["holdout"]

    c_primary = float(p_hold["contrast"])
    c_stress = float(s_hold["contrast"])

    p_q05 = float(p_hold["contrast_bootstrap"]["q05"])
    p_q95 = float(p_hold["contrast_bootstrap"]["q95"])
    s_q05 = float(s_hold["contrast_bootstrap"]["q05"])
    s_q95 = float(s_hold["contrast_bootstrap"]["q95"])

    sigma_p = float((p_q95 - p_q05) / 3.29)
    sigma_s = float((s_q95 - s_q05) / 3.29)

    base_effect = float(min(c_primary, c_stress))
    base_sigma = float(max(sigma_p, sigma_s, 0.10))

    effect_factors = {
        "optimistic": 0.85,
        "reference": 0.60,
        "conservative": 0.45,
        "stress": 0.30,
    }
    sigma_factors = {
        "optimistic": 0.90,
        "reference": 1.00,
        "conservative": 1.25,
        "stress": 1.60,
    }

    powers = [0.80, 0.90, 0.95]
    alpha = 0.025

    rows = []
    for scen in ["optimistic", "reference", "conservative", "stress"]:
        eff = float(base_effect * effect_factors[scen])
        sig = float(base_sigma * sigma_factors[scen])
        row = {
            "scenario": scen,
            "effect_assumed": eff,
            "sigma_assumed": sig,
            "n_holdout_required_effective_iid": {},
            "n_holdout_required_actual": {},
        }
        for pwr in powers:
            n_eff = required_n(eff, sig, alpha_one_sided=alpha, power=float(pwr))
            # Dependent pair-level structure: enforce strict practical floors.
            n_act = int(np.ceil(1.35 * n_eff))
            if pwr >= 0.90 and scen in {"reference", "optimistic"}:
                n_act = max(500, n_act)
            elif pwr >= 0.90 and scen == "conservative":
                n_act = max(600, n_act)
            elif pwr >= 0.90 and scen == "stress":
                n_act = max(800, n_act)
            elif pwr >= 0.80:
                n_act = max(400, n_act)
            row["n_holdout_required_effective_iid"][f"power_{pwr:.2f}"] = int(n_eff)
            row["n_holdout_required_actual"][f"power_{pwr:.2f}"] = int(n_act)
        rows.append(row)

    # Recommended operational targets.
    ref = next(r for r in rows if r["scenario"] == "reference")
    cons = next(r for r in rows if r["scenario"] == "conservative")

    n_ref_90_adj = int(ref["n_holdout_required_actual"]["power_0.90"])
    n_cons_90_adj = int(cons["n_holdout_required_actual"]["power_0.90"])

    out = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "inputs": {
            "base_effect_min_from_qw1922": base_effect,
            "base_sigma_conservative": base_sigma,
            "alpha_one_sided": alpha,
            "powers": powers,
        },
        "scenario_rows": rows,
        "recommended_targets": {
            "n_holdout_reference_power_0p90": n_ref_90_adj,
            "n_holdout_conservative_power_0p90": n_cons_90_adj,
            "n_total_pairs_reference": int(max(2 * n_ref_90_adj, 1200)),
            "n_total_pairs_conservative": int(max(2 * n_cons_90_adj, 1600)),
        },
        "verdict": "BETA_CHANNEL_POWER_REQUIREMENTS_MAP_READY",
    }

    OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1926: BETA-CHANNEL POWER REQUIREMENTS MAP",
        "",
        f"- Data UTC: {out['generated_utc']}",
        f"- Verdict: **{out['verdict']}**",
        "",
        "## Baseline Inputs",
        f"- base_effect_min_from_qw1922: {base_effect:.4f}",
        f"- base_sigma_conservative: {base_sigma:.4f}",
        f"- alpha_one_sided: {alpha}",
        "",
        "## Scenario Table (n_holdout)",
    ]
    for r in rows:
        lines.append(
            f"- {r['scenario']}: eff={r['effect_assumed']:.4f}, sigma={r['sigma_assumed']:.4f}, "
            f"n80_eff={r['n_holdout_required_effective_iid']['power_0.80']}, "
            f"n90_eff={r['n_holdout_required_effective_iid']['power_0.90']}, "
            f"n95_eff={r['n_holdout_required_effective_iid']['power_0.95']}, "
            f"n80_act={r['n_holdout_required_actual']['power_0.80']}, "
            f"n90_act={r['n_holdout_required_actual']['power_0.90']}, "
            f"n95_act={r['n_holdout_required_actual']['power_0.95']}"
        )

    lines.extend(
        [
            "",
            "## Recommended Targets",
            f"- n_holdout_reference_power_0p90: {n_ref_90_adj}",
            f"- n_holdout_conservative_power_0p90: {n_cons_90_adj}",
            f"- n_total_pairs_reference: {out['recommended_targets']['n_total_pairs_reference']}",
            f"- n_total_pairs_conservative: {out['recommended_targets']['n_total_pairs_conservative']}",
            "",
            "## Artifacts",
            f"- JSON: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1926] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1926] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
