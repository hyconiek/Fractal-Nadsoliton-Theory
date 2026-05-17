#!/usr/bin/env python3
"""P1908 S858 strict C1/GR divergent coefficients table v1 probe."""
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


def row(diagram_id: str, sector: str, integral_basis: str, ct_symbol: str, witness_status: str, note: str) -> dict:
    return {
        "diagram_id": diagram_id,
        "sector": sector,
        "integral_basis": integral_basis,
        "uv_pole_order": "1/epsilon",
        "counterterm_symbol": ct_symbol,
        "scheme_tag": "MSbar_candidate",
        "st_identity_contract": "OPEN_REQUIRED",
        "cutkosky_pairing": "OPEN_REQUIRED",
        "status": witness_status,
        "note": note,
    }


def main() -> None:
    p1906 = load("p1906_s856_strict_c1_gr_diagram_inventory_stub_probe.json")
    p1907 = load("p1907_s857_strict_full_lagrangian_to_eom_witness_matrix_probe.json")

    out = {
        "packet_id": "P1908",
        "stage_id": "S858",
        "status": "OPEN_OBSTRUCTION_WITH_TRACE",
        "route": "strict_only",
        "depends_on": {
            "p1906_present": "diagram_inventory_stub" in p1906,
            "p1907_present": "full_lagrangian_term_registry_non_skeleton" in p1907,
            "p1907_missing_witness_exports_present": "missing_witness_exports" in p1907,
        },
        "strict_chain_step": "K_strict -> coefficients -> full L_SM+L_GR -> EOM -> divergent_coefficients_table_v1",
        "divergent_coefficients_table_v1": [
            row("d_scalar4_bubble_s_channel", "scalar4_sector", "I_scalar4_s(mu,epsilon)", "delta_lambda_H_s", "OPEN_NO_NUMERIC", "Need explicit loop integral reduction and pole coefficient extraction."),
            row("d_scalar4_bubble_t_channel", "scalar4_sector", "I_scalar4_t(mu,epsilon)", "delta_lambda_H_t", "OPEN_NO_NUMERIC", "Crossing-consistent coefficient matching still missing."),
            row("d_scalar4_bubble_u_channel", "scalar4_sector", "I_scalar4_u(mu,epsilon)", "delta_lambda_H_u", "OPEN_NO_NUMERIC", "Joint channel cancellation matrix not exported."),
            row("d_yukawa_vertex_correction", "yukawa_sector", "I_yukawa_v(mu,epsilon)", "delta_y_f", "OPEN_NO_NUMERIC", "Need explicit fermion-flavor block decomposition."),
            row("d_yukawa_fermion_self_energy", "yukawa_sector", "I_yukawa_se(mu,epsilon)", "delta_Z_psi_f", "OPEN_NO_NUMERIC", "Wave-function and mass renormalization split not attached."),
            row("d_nonminimal_curvature_scalar_loop", "nonminimal_gravity_mixed_sector", "I_xiH_R(mu,epsilon)", "delta_xi_H", "OPEN_NO_NUMERIC", "Curvature-mixed divergence basis requires covariant tensor projection."),
            row("d_gravity_mixed_counterterm_support", "nonminimal_gravity_mixed_sector", "I_grmix(mu,epsilon)", "delta_c_R2_bundle", "OPEN_NO_NUMERIC", "Higher-curvature counterterm block is listed but not solved."),
        ],
        "first_pass_fail_rows": [
            {
                "row_id": "trial_scalar4_s",
                "pass_condition": "pole(d_scalar4_bubble_s_channel) + delta_lambda_H_s = 0",
                "current_value": "UNSET",
                "status": "OPEN_FAIL_BY_MISSING_DATA",
            },
            {
                "row_id": "trial_yukawa_vertex",
                "pass_condition": "pole(d_yukawa_vertex_correction) + delta_y_f = 0",
                "current_value": "UNSET",
                "status": "OPEN_FAIL_BY_MISSING_DATA",
            },
            {
                "row_id": "trial_xiH_curvature",
                "pass_condition": "pole(d_nonminimal_curvature_scalar_loop) + delta_xi_H = 0",
                "current_value": "UNSET",
                "status": "OPEN_FAIL_BY_MISSING_DATA",
            },
        ],
        "closure_blockers_live": {
            "renormalization": "OPEN: numeric/symbolic pole coefficients missing",
            "unitarity": "OPEN: diagram rows not yet paired with ImM/Cutkosky equality rows",
            "background_independence": "OPEN: FRW/Bianchi-I transport not yet tied to renormalized EOM closure",
            "selector_qw2191": "OPEN: strict-core selector source/symmetry-breaking export missing",
        },
        "false_pass_guard": "Table seeding with scheme tags is not a renormalization proof; PASS requires explicit coefficients and cancellations.",
        "next_honest_step": "Export P1909 with explicit coefficient values (symbolic or numeric), ST/Ward contracts, and matched ImM/Cutkosky rows for each diagram id.",
        "lay_explanation": "To jak lista zadań do rachunku pętli: wiemy które diagramy i jakie kontrtermy trzeba zestawić, ale bez policzonych współczynników nie da się uczciwie ogłosić domknięcia teorii.",
    }

    path = GEN / "p1908_s858_strict_c1_gr_divergent_coefficients_table_v1_probe.json"
    path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
