#!/usr/bin/env python3
"""
QW-1703: Claims-vs-computation audit.

Goal:
1) Recompute key README formulas from stated parameters.
2) Compare to experimental targets.
3) Contrast with "0.00% / EXACT" style declarations.
4) Track source document dates used in this audit.
"""

from __future__ import annotations

import json
import math
import re
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parent
README = ROOT / "README.md"
FINAL_SUMMARY = ROOT / "FINAL_SUMMARY_TEORIA_WSZYSTKIEGO.md"

OUT_JSON = ROOT / "report_qw1703_claims_vs_computation_audit.json"
OUT_MD = ROOT / "RAPORT_QW1703_CLAIMS_VS_COMPUTATION_AUDIT.md"


def pct_err(pred: float, target: float) -> float:
    if target == 0:
        return float("nan")
    return abs(pred - target) / abs(target) * 100.0


def file_meta(path: Path) -> dict:
    st = path.stat()
    return {
        "path": str(path.relative_to(ROOT)),
        "mtime_utc": datetime.fromtimestamp(st.st_mtime, tz=timezone.utc).isoformat(),
        "size_bytes": st.st_size,
    }


def main() -> None:
    alpha_geo = 4.0 * math.log(2.0)
    beta_tors = 0.01
    gamma = 1.52
    m_top = 173_000.0  # MeV

    # README table-level targets
    sin2_exp = 0.23122
    alpha_inv_exp = 137.035999

    # 1) Weinberg (as stated in README: sin^2(theta_W) = alpha_geo / 12)
    sin2_pred = alpha_geo / 12.0
    sin2_error = pct_err(sin2_pred, sin2_exp)

    # 2) Fine-structure inverse variants used in repo ecosystem
    alpha_inv_simple = alpha_geo / (2.0 * beta_tors)
    alpha_inv_simple_error = pct_err(alpha_inv_simple, alpha_inv_exp)

    alpha_inv_verify_values = 0.5 * (alpha_geo / beta_tors) * (1 - beta_tors)
    alpha_inv_verify_values_error = pct_err(alpha_inv_verify_values, alpha_inv_exp)

    # 3) Info-geometry identity
    golden_ratio = (1.0 + math.sqrt(5.0)) / 2.0
    info_side = alpha_geo
    geom_side = golden_ratio * math.sqrt(3.0)
    identity_gap_pct = pct_err(info_side, geom_side)

    # 4) Mass formula M(Q) = M_top * 4^(-gamma * Q/4)
    exp_masses = {
        "Top": 173_000.0,
        "Bottom": 4_180.0,
        "Tau": 1_776.9,
        "Charm": 1_270.0,
        "Muon": 105.7,
        "Electron": 0.511,
    }
    q_model = {
        "Top": 0,
        "Bottom": 7,
        "Tau": 9,
        "Charm": 9,
        "Muon": 14,
        "Electron": 24,
    }
    mass_rows = []
    for particle, q in q_model.items():
        pred = m_top * (4.0 ** (-(gamma * q / 4.0)))
        exp = exp_masses[particle]
        err = pct_err(pred, exp)
        mass_rows.append(
            {
                "particle": particle,
                "Q": q,
                "pred_mev": pred,
                "exp_mev": exp,
                "error_pct": err,
            }
        )

    # 5) Read declaration density from key docs
    readme_text = README.read_text(encoding="utf-8", errors="ignore")
    final_text = FINAL_SUMMARY.read_text(encoding="utf-8", errors="ignore")
    exact_mentions = len(re.findall(r"\bexact\b|0\.00%", (readme_text + "\n" + final_text), flags=re.IGNORECASE))
    fitting_mentions = len(re.findall(r"fitting|fit|kalibr|calibr", (readme_text + "\n" + final_text), flags=re.IGNORECASE))

    output = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "sources": [file_meta(README), file_meta(FINAL_SUMMARY)],
        "parameters": {
            "alpha_geo": alpha_geo,
            "beta_tors": beta_tors,
            "gamma": gamma,
            "m_top_mev": m_top,
        },
        "recomputed_observables": {
            "weinberg_from_alpha_over_12": {
                "pred": sin2_pred,
                "exp": sin2_exp,
                "error_pct": sin2_error,
            },
            "fine_structure_inverse_alpha_geo_over_2beta": {
                "pred": alpha_inv_simple,
                "exp": alpha_inv_exp,
                "error_pct": alpha_inv_simple_error,
            },
            "fine_structure_inverse_verify_values_variant": {
                "pred": alpha_inv_verify_values,
                "exp": alpha_inv_exp,
                "error_pct": alpha_inv_verify_values_error,
            },
            "info_geometry_identity": {
                "left_4ln2": info_side,
                "right_phi_sqrt3": geom_side,
                "gap_pct": identity_gap_pct,
            },
        },
        "mass_formula_check": mass_rows,
        "declaration_density": {
            "exact_or_0_00_mentions_in_core_docs": exact_mentions,
            "fitting_or_calibration_mentions_in_core_docs": fitting_mentions,
        },
        "headline_findings": [
            "Formuła Weinberga alpha_geo/12 daje błąd niezerowy (~0.07%), a nie 0.00%.",
            "Formuła alpha_geo/(2*beta) daje ~138.63 (błąd ~1.16%), nie 137.24.",
            "Wariant z czynnikiem (1-beta) daje ~137.24 (błąd ~0.15%), zgodny z częścią raportów.",
            "Tożsamość 4ln2 ≈ phi*sqrt(3) ma ~1% luki, więc jest przybliżeniem, nie równością ścisłą.",
        ],
    }

    OUT_JSON.write_text(json.dumps(output, ensure_ascii=False, indent=2), encoding="utf-8")

    lines = [
        "# RAPORT QW-1703: CLAIMS VS COMPUTATION AUDIT",
        "",
        f"- Data UTC: {output['generated_utc']}",
        f"- Źródła: `{README.name}`, `{FINAL_SUMMARY.name}`",
        f"- Daty źródeł (mtime UTC): {output['sources'][0]['mtime_utc']} ; {output['sources'][1]['mtime_utc']}",
        "",
        "## 1) Recomputacja kluczowych wzorów",
        f"- sin²θ_W = 4ln2/12 = {sin2_pred:.6f}, exp={sin2_exp:.6f}, błąd={sin2_error:.4f}%",
        f"- α⁻¹ (4ln2/(2β), β=0.01) = {alpha_inv_simple:.6f}, exp={alpha_inv_exp:.6f}, błąd={alpha_inv_simple_error:.4f}%",
        f"- α⁻¹ (wariant verify_values) = {alpha_inv_verify_values:.6f}, exp={alpha_inv_exp:.6f}, błąd={alpha_inv_verify_values_error:.4f}%",
        f"- 4ln2 vs φ√3: {info_side:.6f} vs {geom_side:.6f}, luka={identity_gap_pct:.4f}%",
        "",
        "## 2) Test wzoru masowego (Q, γ=1.52)",
    ]
    for row in mass_rows:
        lines.append(
            f"- {row['particle']}: Q={row['Q']}, pred={row['pred_mev']:.4f} MeV, exp={row['exp_mev']:.4f} MeV, błąd={row['error_pct']:.2f}%"
        )
    lines.extend(
        [
            "",
            "## 3) Kontrast deklaracji vs obliczeń",
            f"- Wzmianki `EXACT/0.00%` (README + FINAL_SUMMARY): {exact_mentions}",
            f"- Wzmianki `fitting/calibration` (README + FINAL_SUMMARY): {fitting_mentions}",
            "",
            "## 4) Artefakty",
            f"- JSON szczegółowy: `{OUT_JSON.name}`",
        ]
    )

    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"[QW-1703] Saved JSON: {OUT_JSON.name}")
    print(f"[QW-1703] Saved MD:   {OUT_MD.name}")


if __name__ == "__main__":
    main()
