#!/usr/bin/env python3
"""P1909 S859 strict C1/GR counterterm-ST-Cutkosky binding probe."""
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


def bind_row(diagram_id: str, ct: str, st: str, cutkosky: str, coeff_tag: str) -> dict:
    return {
        "diagram_id": diagram_id,
        "counterterm_symbol": ct,
        "st_ward_contract": st,
        "cutkosky_equality": cutkosky,
        "coefficient_export": coeff_tag,
        "status": "OPEN_SYMBOLIC_ONLY",
    }


def main() -> None:
    p1908 = load("p1908_s858_strict_c1_gr_divergent_coefficients_table_v1_probe.json")

    out = {
        "packet_id": "P1909",
        "stage_id": "S859",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1908_present": "divergent_coefficients_table_v1" in p1908,
            "p1908_row_count": len(p1908.get("divergent_coefficients_table_v1", [])),
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> divergent coefficients -> ST/Ward + Cutkosky binding",
        "diagram_binding_matrix_v1": [
            bind_row(
                "d_scalar4_bubble_s_channel",
                "delta_lambda_H_s",
                "ST_H4_s: Z_lambda - Z_H^2 consistency",
                "Disc M_H4_s = Sum_n M(HH->n)M*(HH->n)",
                "c_scalar4_s",
            ),
            bind_row(
                "d_scalar4_bubble_t_channel",
                "delta_lambda_H_t",
                "ST_H4_t: crossing-preserved renormalized vertex identity",
                "Disc M_H4_t = Sum_n M(HH->n)M*(HH->n)",
                "c_scalar4_t",
            ),
            bind_row(
                "d_scalar4_bubble_u_channel",
                "delta_lambda_H_u",
                "ST_H4_u: crossing-preserved renormalized vertex identity",
                "Disc M_H4_u = Sum_n M(HH->n)M*(HH->n)",
                "c_scalar4_u",
            ),
            bind_row(
                "d_yukawa_vertex_correction",
                "delta_y_f",
                "ST_Yf: Z_yf = Z_psi^{-1} Z_H^{-1/2} Z_vertex",
                "Disc M_Yf = Sum_n M(fH->n)M*(fH->n)",
                "c_yukawa_vertex",
            ),
            bind_row(
                "d_yukawa_fermion_self_energy",
                "delta_Z_psi_f",
                "Ward_psi: dSigma/dslashp at p=mf gives Z_psi",
                "Disc Sigma_f = Sum_n M(f->n)M*(f->n)",
                "c_yukawa_self_energy",
            ),
            bind_row(
                "d_nonminimal_curvature_scalar_loop",
                "delta_xi_H",
                "Ward_gr_H: improved T_{mu nu} conservation with xi_H term",
                "Disc M_xiH = Sum_n M(HR->n)M*(HR->n)",
                "c_xiH_curvature",
            ),
            bind_row(
                "d_gravity_mixed_counterterm_support",
                "delta_c_R2_bundle",
                "Bianchi+BRST: covariant counterterm closure in R2/Rmunu2 sector",
                "Disc M_grmix = Sum_n M(g_mix->n)M*(g_mix->n)",
                "c_gravity_mixed",
            ),
        ],
        "coefficient_exports_required": {
            "format": "{coefficient_symbol, expression, scheme, scale_mu, epsilon_order, uncertainty_if_numeric}",
            "status": "OPEN_REQUIRED",
            "note": "No coefficient value is asserted here; this packet defines binding contracts only.",
        },
        "pass_fail_contract_v1": [
            {
                "contract_id": "renorm_scalar4_s",
                "condition": "pole_coeff(c_scalar4_s) + delta_lambda_H_s = 0",
                "status": "OPEN_FAIL_BY_MISSING_COEFFICIENT",
            },
            {
                "contract_id": "st_yukawa_vertex",
                "condition": "ST_Yf identity satisfied after renormalization",
                "status": "OPEN_FAIL_BY_MISSING_COEFFICIENT",
            },
            {
                "contract_id": "cutkosky_grmix",
                "condition": "Disc M_grmix equals physical cut sum in common scheme",
                "status": "OPEN_FAIL_BY_MISSING_COEFFICIENT",
            },
        ],
        "closure_blockers_live": {
            "renormalization": "OPEN: explicit c_* coefficients absent",
            "unitarity": "OPEN: channel discontinuity values absent",
            "background_independence": "OPEN: transport witness still not chained to renormalized metric EOM",
            "selector_qw2191": "OPEN: strict-core selector source/premise missing",
        },
        "false_pass_guard": "Symbolic binding contracts are not proof of closure; strict-core closure requires exported coefficients and verified equalities.",
        "next_honest_step": "Export P1910 coefficient_table_v1 with explicit c_* values (symbolic or numeric), then evaluate pass_fail_contract_v1 in the same renormalization and cut scheme.",
        "lay_explanation": "Ten krok łączy każdy diagram z trzema testami: kasowaniem nieskończoności, zgodnością symetrii i warunkiem unitarności. Ale dopóki nie ma konkretnych liczb/wzorów współczynników, wynik pozostaje otwarty.",
    }

    path = GEN / "p1909_s859_strict_c1_gr_counterterm_st_cutkosky_binding_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
