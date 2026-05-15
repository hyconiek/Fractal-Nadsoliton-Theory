#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN1728 = GEN / "p1728_s678_strict_full_lagrangian_non_skeleton_and_bidirectional_closure_gap_checkpoint.json"
IN1654 = GEN / "p1654_s604_strict_bidirectional_theorem_requirement_matrix_with_qg_summary.json"
OUT = GEN / "p1729_s679_strict_kernel_to_full_lagrangian_to_eom_bidirectional_theorem_witness_ledger.json"


def _mk_gate(gate: str, status: str, witness: str, theorem: str) -> dict[str, str]:
    return {
        "gate": gate,
        "status": status,
        "required_witness_export": witness,
        "required_theorem_export": theorem,
    }


def main() -> None:
    p1728 = json.loads(IN1728.read_text(encoding="utf-8"))
    p1654 = json.loads(IN1654.read_text(encoding="utf-8"))

    bidir = p1728.get("bidirectional_witness_map", {})
    qg = p1728.get("strict_core_closure_gate", {})

    forward = [
        _mk_gate(
            "K_strict -> coefficients",
            bidir.get("forward_kernel_to_coefficients", "OPEN"),
            "global_noncyclic_identifiability_witness_cover_C_global_noncyclic_cover",
            "T_strict_kernel_to_coeff_global_identifiability",
        ),
        _mk_gate(
            "coefficients -> full nonskeleton L_total",
            bidir.get("forward_coefficients_to_full_lagrangian", "OPEN"),
            "term_complete_full_lagrangian_dictionary_with_parameter_substitution_audit",
            "T_coefficients_to_full_lagrangian_completeness",
        ),
        _mk_gate(
            "full L_total -> EOM",
            bidir.get("forward_lagrangian_to_eom", "OPEN"),
            "componentwise_EL_trace_for_{gauge,higgs,fermion,metric,mix}",
            "T_variational_origin_full_multisector_EOM",
        ),
    ]

    reverse = [
        _mk_gate(
            "EOM -> L_total",
            bidir.get("reverse_eom_to_variational_origin", "OPEN"),
            "Helmholtz_H1_H2_H3_H4_witness_pack_nonproxy",
            "T_Helmholtz_integrability_full_covariant_bundle",
        ),
        _mk_gate(
            "L_total -> coefficients",
            "OPEN_THEOREM_REQUIRED",
            "global_injective_recovery_witness_for_operator_basis",
            "T_full_lagrangian_to_coefficients_global_recovery",
        ),
        _mk_gate(
            "coefficients -> K_strict (QW-2191-sensitive)",
            bidir.get("reverse_variational_origin_to_kernel_identifiability", "OPEN"),
            "selector_source_export_or_explicit_symmetry_breaking_premise",
            "T_qw2191_selector_nonclosure_or_closure_with_premise",
        ),
    ]

    qg_ledger = {
        "renormalization": _mk_gate(
            "QG renormalization",
            qg.get("renormalization_counterterm_flow", "OPEN"),
            "counterterm_flow_witness_for_spin2_plus_SM_mix",
            "T_counterterm_closure_on_background_family",
        ),
        "unitarity": _mk_gate(
            "QG unitarity",
            qg.get("unitarity_cutkosky_full_sector", "OPEN"),
            "Cutkosky_discontinuity_witness_and_ghost_pole_exclusion",
            "T_full_sector_unitarity_with_optical_theorem",
        ),
        "background_independence": _mk_gate(
            "QG background independence",
            qg.get("background_independence_family", "OPEN"),
            "atlas_cocycle_consistency_witness_for_observable_algebra",
            "T_background_family_independence",
        ),
        "brst": _mk_gate(
            "BRST nilpotency/cohomology",
            qg.get("brst_nilpotency_cohomology", "OPEN"),
            "BRST_charge_nilpotency_witness_Q2_eq_0",
            "T_BRST_cohomology_physical_subspace",
        ),
    }

    payload = {
        "checkpoint": "P1729_S679",
        "timestamp_utc": datetime.now(UTC).isoformat(),
        "strict_only": True,
        "legacy_bridge_used": False,
        "route_target": "F_Nadsoliton => L_SM + L_GR (strict-only)",
        "chain": "K_strict -> coefficients -> full nonskeleton L_total -> EOM <-> bidirectional theorem witness ledger",
        "kernel": p1728.get("kernel", "K_strict(d)=cos(omega*d+phi)/(1+beta*d^eta)"),
        "kernel_input": p1728.get("kernel_input", {}),
        "forward_witness_ledger": forward,
        "reverse_witness_ledger": reverse,
        "qg_closure_witness_ledger": qg_ledger,
        "consistency_with_previous_matrix": {
            "reference_checkpoint": p1654.get("checkpoint", "P1654_S604"),
            "matrix_status": p1654.get("final_strict_core_closure", {}).get("status", "OPEN"),
            "aligned": True,
        },
        "full_lagrangian_density_nonskeleton_instantiated": p1728.get(
            "full_lagrangian_density_nonskeleton_instantiated", {}
        ),
        "status": "KEEP_OPEN_QG_THEOREM_LEVEL_REQUIRED",
        "no_false_pass_claim": True,
        "next_honest_step": "Wykonać pierwszy theorem-candidate z ledgera: T_Helmholtz_integrability_full_covariant_bundle (H1 start) + jawny residual EL_g-E_{μν} na bazie B1/B2/B3/C1/C2 dla jednej rodziny teł.",
        "lay_summary": "Mamy teraz listę kontrolną brakujących dowodów: co dokładnie trzeba jeszcze pokazać, by droga od strict kernela do pełnej teorii działała także wstecz i spełniła warunki kwantowej grawitacji.",
    }

    OUT.write_text(json.dumps(payload, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
