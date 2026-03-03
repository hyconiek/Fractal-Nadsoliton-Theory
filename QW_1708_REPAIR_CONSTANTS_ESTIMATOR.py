#!/usr/bin/env python3
"""
QW-1708: Estimate minimal effective correction terms for FIN 4.5 repaired core.

Purpose:
1) Quantify small correction needed for Weinberg relation.
2) Quantify correction needed for alpha_EM relation.
3) Quantify per-particle mass correction exponent Delta_a for mass formula.
"""

from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT_JSON = ROOT / "report_qw1708_repair_constants_estimator.json"
OUT_MD = ROOT / "RAPORT_QW1708_REPAIR_CONSTANTS_ESTIMATOR.md"


def main() -> None:
    alpha_geo = 4.0 * math.log(2.0)
    beta = 0.01
    gamma = 1.52
    m_top = 173_000.0  # MeV

    # Targets
    sin2_exp = 0.23122
    alpha_inv_exp = 137.035999
    masses_exp = {
        "Top": 173_000.0,
        "Bottom": 4_180.0,
        "Tau": 1_776.9,
        "Charm": 1_270.0,
        "Muon": 105.7,
        "Electron": 0.511,
    }
    q_assign = {
        "Top": 0,
        "Bottom": 7,
        "Tau": 9,
        "Charm": 9,
        "Muon": 14,
        "Electron": 24,
    }

    # 1) Weinberg correction
    sin2_base = alpha_geo / 12.0
    delta_w = sin2_exp / sin2_base - 1.0

    # 2) Fine-structure correction (using repaired baseline with (1-beta))
    alpha_inv_base = (alpha_geo / (2.0 * beta)) * (1.0 - beta)
    delta_vac = alpha_inv_exp / alpha_inv_base - 1.0

    # 3) Mass corrections
    mass_rows = []
    for p, q in q_assign.items():
        pred_base = m_top * (4.0 ** (-(gamma * q / 4.0)))
        exp = masses_exp[p]
        if pred_base > 0 and exp > 0:
            delta_a = math.log(exp / pred_base)
        else:
            delta_a = float("nan")
        mass_rows.append(
            {
                "particle": p,
                "Q": q,
                "pred_base_mev": pred_base,
                "exp_mev": exp,
                "delta_a": delta_a,
                "relative_error_base_pct": abs(pred_base - exp) / exp * 100.0,
            }
        )

    # Summary stats
    deltas_non_top = [r["delta_a"] for r in mass_rows if r["particle"] != "Top"]
    mean_abs_delta = sum(abs(x) for x in deltas_non_top) / len(deltas_non_top)
    max_abs_delta = max(abs(x) for x in deltas_non_top)

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "frozen_params": {
            "alpha_geo": alpha_geo,
            "beta_tors": beta,
            "gamma": gamma,
            "m_top_mev": m_top,
        },
        "effective_corrections": {
            "delta_W_for_weinberg": delta_w,
            "delta_vac_for_alpha_inverse": delta_vac,
        },
        "mass_delta_table": mass_rows,
        "delta_summary": {
            "mean_abs_delta_non_top": mean_abs_delta,
            "max_abs_delta_non_top": max_abs_delta,
        },
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1708: REPAIR CONSTANTS ESTIMATOR",
        "",
        f"- Data UTC: {output['generated_utc']}",
        "",
        "## 1) Poprawka Weinberga",
        f"- sin²θ_W(base)=4ln2/12 = {sin2_base:.9f}",
        f"- sin²θ_W(exp) = {sin2_exp:.9f}",
        f"- wymagane δ_W = {delta_w:+.9e}",
        "",
        "## 2) Poprawka α_EM",
        f"- α⁻¹(base)=4ln2/(2β)*(1-β) = {alpha_inv_base:.9f}",
        f"- α⁻¹(exp) = {alpha_inv_exp:.9f}",
        f"- wymagane δ_vac = {delta_vac:+.9e}",
        "",
        "## 3) Poprawki masowe Δ_a",
    ]
    for r in mass_rows:
        lines.append(
            f"- {r['particle']}: Q={r['Q']}, pred={r['pred_base_mev']:.4f} MeV, exp={r['exp_mev']:.4f} MeV, Δ_a={r['delta_a']:+.6f}, błąd_base={r['relative_error_base_pct']:.2f}%"
        )
    lines.extend(
        [
            "",
            "## 4) Agregat",
            f"- średnia |Δ_a| (bez top): {mean_abs_delta:.6f}",
            f"- max |Δ_a| (bez top): {max_abs_delta:.6f}",
            "",
            "## Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"[QW-1708] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1708] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
