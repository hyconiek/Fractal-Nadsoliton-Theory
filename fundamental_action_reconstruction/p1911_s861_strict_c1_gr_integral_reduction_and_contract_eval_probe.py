#!/usr/bin/env python3
"""P1911 S861 strict C1/GR integral-reduction and contract-eval probe."""
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


def ibp_row(coeff: str, master_basis: list[str], pole_symbol: str, finite_symbol: str) -> dict:
    return {
        "coefficient_symbol": coeff,
        "master_integral_basis": master_basis,
        "pole_projection": pole_symbol,
        "finite_projection": finite_symbol,
        "status": "OPEN_REDUCTION_PENDING",
    }


def eval_row(contract_id: str, input_needed: list[str], condition: str) -> dict:
    return {
        "contract_id": contract_id,
        "inputs_required": input_needed,
        "evaluation_condition": condition,
        "status": "OPEN_FAIL_BY_MISSING_EVALUATION_INPUTS",
    }


def main() -> None:
    p1910 = load("p1910_s860_strict_c1_gr_coefficient_table_v1_probe.json")

    out = {
        "packet_id": "P1911",
        "stage_id": "S861",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1910_present": "coefficient_table_v1" in p1910,
            "p1910_coeff_count": len(p1910.get("coefficient_table_v1", [])),
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> coefficient_table_v1 -> integral_reduction -> first contract evaluations",
        "integral_reduction_plan_v1": [
            ibp_row("c_scalar4_s", ["J_bub(s)", "J_ct"], "A_s", "F_s"),
            ibp_row("c_scalar4_t", ["J_bub(t)", "J_ct"], "A_t", "F_t"),
            ibp_row("c_scalar4_u", ["J_bub(u)", "J_ct"], "A_u", "F_u"),
            ibp_row("c_yukawa_vertex", ["J_tri(m_f)", "J_bub(m_f)"], "A_yv", "F_yv"),
            ibp_row("c_yukawa_self_energy", ["J_se(m_f)", "J_tad(m_f)"], "A_ys", "F_ys"),
            ibp_row("c_xiH_curvature", ["J_curvH(R)", "J_mix(R,H)"], "A_xi", "F_xi"),
            ibp_row("c_gravity_mixed", ["J_R2", "J_Rmunu2", "J_EH_mix"], "A_gr", "F_gr"),
        ],
        "coefficient_export_template_v2": {
            "required_fields": [
                "coefficient_symbol",
                "A_pole_value",
                "B_log_value",
                "F_finite_value",
                "scheme",
                "scale_mu",
                "derivation_trace_ref",
            ],
            "status": "OPEN_REQUIRED",
        },
        "first_contract_evaluation_queue": [
            eval_row(
                "renorm_scalar4_s",
                ["A_s", "delta_lambda_H_s"],
                "A_s + delta_lambda_H_s == 0",
            ),
            eval_row(
                "renorm_yukawa_vertex",
                ["A_yv", "delta_y_f"],
                "A_yv + delta_y_f == 0",
            ),
            eval_row(
                "st_yukawa_vertex",
                ["Z_yf", "Z_psi", "Z_H", "Z_vertex"],
                "Z_yf - Z_psi**(-1)*Z_H**(-1/2)*Z_vertex == 0",
            ),
            eval_row(
                "cutkosky_H4_s",
                ["DiscM_H4_s", "CutSum_H4_s"],
                "DiscM_H4_s - CutSum_H4_s == 0",
            ),
            eval_row(
                "cutkosky_grmix",
                ["DiscM_grmix", "CutSum_grmix"],
                "DiscM_grmix - CutSum_grmix == 0",
            ),
            eval_row(
                "background_FRW_BI",
                ["Residual_FRW_ren", "Residual_BianchiI_ren", "common_scheme_tag"],
                "Residual_FRW_ren - Residual_BianchiI_ren == 0",
            ),
        ],
        "strict_core_closure_statusvector": {
            "renormalization": "OPEN",
            "unitarity": "OPEN",
            "background_independence": "OPEN",
            "selector_qw2191": "OPEN",
        },
        "false_pass_guard": "Reduction plan and queued evaluations are not solved equalities; strict-core PASS is forbidden until explicit values are exported and checked.",
        "next_honest_step": "Export P1912 with first explicit A_* values from integral reductions and evaluate at least two renorm + one Cutkosky + one FRW/Bianchi-I contract row.",
        "lay_explanation": "To etap przejścia od wzorów ogólnych do konkretnych wyników: wskazujemy jakie całki zredukować i jakie równania sprawdzić. Bez tych liczb teoria wciąż nie jest domknięta.",
    }

    path = GEN / "p1911_s861_strict_c1_gr_integral_reduction_and_contract_eval_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
