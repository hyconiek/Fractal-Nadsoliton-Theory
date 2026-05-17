#!/usr/bin/env python3
"""P1913 S863 strict C1/GR common-basis unitarity/background probe."""
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def load(name: str) -> dict:
    p = GEN / name
    if not p.exists():
        return {"_missing": name, "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(p.read_text(encoding="utf-8"))


def cut_row(channel: str, disc_expr: str, cut_expr: str, basis: str) -> dict:
    return {
        "channel": channel,
        "DiscM_common_basis": disc_expr,
        "CutSum_common_basis": cut_expr,
        "common_basis": basis,
        "delta": f"({disc_expr})-({cut_expr})",
        "status": "OPEN_SYMBOLIC_REDUCED_NOT_EVALUATED",
    }


def bg_row(label: str, frw: str, bi: str, basis: str) -> dict:
    return {
        "label": label,
        "Residual_FRW_ren": frw,
        "Residual_BianchiI_ren": bi,
        "common_basis": basis,
        "delta": f"({frw})-({bi})",
        "status": "OPEN_SYMBOLIC_REDUCED_NOT_EVALUATED",
    }


def main() -> None:
    p1912 = load("p1912_s862_strict_c1_gr_first_explicit_A_exports_and_eval_probe.json")

    out = {
        "packet_id": "P1913",
        "stage_id": "S863",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1912_present": "first_explicit_A_exports_v1" in p1912,
            "p1912_a_exports": len(p1912.get("first_explicit_A_exports_v1", [])),
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> A exports -> common-basis Disc/Cut and FRW/Bianchi-I residual rows",
        "common_basis_definition": {
            "loop_basis": ["J_bub", "J_tri", "J_R2", "J_Rmunu2", "J_EH_mix"],
            "kinematic_basis": ["s", "t", "u", "m_f^2", "R"],
            "scheme": "MSbar_candidate",
            "status": "DECLARED",
        },
        "unitarity_common_basis_rows_v1": [
            cut_row("H4_s", "alpha_H4s*J_bub + beta_H4s", "alphaH4s_cut*J_bub + betaH4s_cut", "{J_bub,1}"),
            cut_row("Yf", "alpha_Yf*J_tri + beta_Yf*J_bub + gamma_Yf", "alphaYf_cut*J_tri + betaYf_cut*J_bub + gammaYf_cut", "{J_tri,J_bub,1}"),
            cut_row("grmix", "alpha_gr*J_R2 + beta_gr*J_Rmunu2 + gamma_gr*J_EH_mix", "alpha_gr_cut*J_R2 + beta_gr_cut*J_Rmunu2 + gamma_gr_cut*J_EH_mix", "{J_R2,J_Rmunu2,J_EH_mix}"),
        ],
        "background_common_basis_rows_v1": [
            bg_row(
                "metric_residual_match",
                "r1*J_R2 + r2*J_Rmunu2 + r3",
                "b1*J_R2 + b2*J_Rmunu2 + b3",
                "{J_R2,J_Rmunu2,1}",
            )
        ],
        "first_evaluated_contract_rows": [
            {
                "contract_id": "cutkosky_H4_s",
                "condition": "DiscM_common_basis - CutSum_common_basis == 0",
                "evaluation_state": "OPEN_COEFFICIENT_MATCH_REQUIRED",
            },
            {
                "contract_id": "background_FRW_BI",
                "condition": "Residual_FRW_ren - Residual_BianchiI_ren == 0",
                "evaluation_state": "OPEN_COEFFICIENT_MATCH_REQUIRED",
            },
        ],
        "full_lagrangian_anchor_non_skeleton": {
            "reference": "P1907::full_lagrangian_term_registry_non_skeleton",
            "status": "REQUIRED_AND_ACTIVE",
        },
        "strict_core_closure_statusvector": {
            "renormalization": "OPEN_PARTIAL_A_ONLY",
            "unitarity": "OPEN_COMMON_BASIS_DECLARED",
            "background_independence": "OPEN_COMMON_BASIS_DECLARED",
            "selector_qw2191": "OPEN",
        },
        "false_pass_guard": "Common-basis reduced rows are structural only; closure requires solved coefficient matching and verified zero deltas.",
        "next_honest_step": "Export P1914 with solved coefficient-matching subrows for H4_s and FRW/Bianchi-I residual channel plus explicit PASS/FAIL evaluation stamps.",
        "lay_explanation": "Ujednoliciliśmy obie strony porównań do tej samej bazy matematycznej, więc teraz da się uczciwie sprawdzać równości. Nadal jednak trzeba dopasować współczynniki, żeby uzyskać prawdziwe PASS albo FAIL.",
    }

    path = GEN / "p1913_s863_strict_c1_gr_common_basis_unitarity_background_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
