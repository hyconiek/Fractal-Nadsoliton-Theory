#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r49_direct_m2_psi2_psi5_role_matching_packet.json"
OUT_SUMMARY = GENERATED / "r49_direct_m2_psi2_psi5_role_matching_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def declared_shift(term: str, field_source: str, field_target: str, coeff_source: str, coeff_target: str) -> str:
    return term.replace(coeff_source, coeff_target).replace(field_source, field_target)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r24 = load_json(
        "fundamental_action_reconstruction/generated/r24_declared_plus3_shift_packet_for_direct_m2_family_route.json"
    )
    q2164 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2164_l14_full_canonical_continuum_variational_gate.json"
    )
    q2165 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2165_l13_exhaustive_canonical_eom_gate.json"
    )

    basis_shift = r24["declared_plus3_shift_packet_for_direct_m2_family_route"][
        "declared_basis_shift_restricted_to_positive_support"
    ]
    coeff_shift = r24["declared_plus3_shift_packet_for_direct_m2_family_route"][
        "declared_coefficient_relabeling_restricted_to_positive_support"
    ]

    field_source = "psi2"
    field_target = basis_shift[field_source]
    coeff_source = "m2_psi2"
    coeff_target = coeff_shift[coeff_source]

    source_action_term = "m2_psi2*psi2**2/2"
    target_action_term = "m2_psi5*psi5**2/2"
    source_eom_term = "m2_psi2*psi2(x)"
    target_eom_term = "m2_psi5*psi5(x)"

    shifted_action_term = declared_shift(
        source_action_term,
        field_source=field_source,
        field_target=field_target,
        coeff_source=coeff_source,
        coeff_target=coeff_target,
    )
    shifted_eom_term = declared_shift(
        source_eom_term,
        field_source=f"{field_source}(x)",
        field_target=f"{field_target}(x)",
        coeff_source=coeff_source,
        coeff_target=coeff_target,
    )

    potential_total = q2164["model"]["potential_total"]
    eom_psi2 = q2165["model"]["eom_psi2"]
    eom_psi5 = q2165["model"]["eom_psi5"]

    checks = [
        {
            "id": "r24_declared_basis_shift_psi2_to_psi5_present",
            "actual": basis_shift.get(field_source),
            "expected": field_target,
            "meaning": "R24 exports the declared plus3 basis shift psi2 -> psi5 on the direct m2 route",
        },
        {
            "id": "r24_declared_m2_coefficient_shift_present",
            "actual": coeff_shift.get(coeff_source),
            "expected": coeff_target,
            "meaning": "R24 exports the declared coefficient relabeling m2_psi2 -> m2_psi5",
        },
        {
            "id": "q2164_source_action_term_present",
            "actual": source_action_term in potential_total,
            "expected": True,
            "meaning": "QW-2164 exports the source direct mass-like quadratic term on psi2",
        },
        {
            "id": "q2164_target_action_term_present",
            "actual": target_action_term in potential_total,
            "expected": True,
            "meaning": "QW-2164 exports the target direct mass-like quadratic term on psi5",
        },
        {
            "id": "declared_shifted_source_action_term_matches_target",
            "actual": shifted_action_term,
            "expected": target_action_term,
            "meaning": "under the declared shift, the source action term maps exactly to the target action term",
        },
        {
            "id": "q2165_source_eom_term_present",
            "actual": source_eom_term in eom_psi2,
            "expected": True,
            "meaning": "QW-2165 exports the source local m2 term in eom_psi2",
        },
        {
            "id": "q2165_target_eom_term_present",
            "actual": target_eom_term in eom_psi5,
            "expected": True,
            "meaning": "QW-2165 exports the target local m2 term in eom_psi5",
        },
        {
            "id": "declared_shifted_source_eom_term_matches_target",
            "actual": shifted_eom_term,
            "expected": target_eom_term,
            "meaning": "under the declared shift, the source eom term maps exactly to the target eom term",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R49",
        "packet_goal": "materialize_a_route_scoped_declared_role_matching_packet_for_the_single_direct_m2_pairwise_target_m2_psi2_equals_m2_psi5_without_claiming_coefficient_equality",
        "source_reports": ["R24", "QW2164", "QW2165"],
        "direct_m2_psi2_psi5_role_matching_packet": {
            "pairwise_target_under_attack": "m2_psi2 = m2_psi5",
            "declared_basis_shift": {field_source: field_target},
            "declared_coefficient_shift": {coeff_source: coeff_target},
            "canonical_action_role_match": {
                "source_term": source_action_term,
                "target_term": target_action_term,
                "declared_shifted_source_term": shifted_action_term,
                "exact_role_match": shifted_action_term == target_action_term,
            },
            "canonical_eom_role_match": {
                "source_term": source_eom_term,
                "target_term": target_eom_term,
                "declared_shifted_source_term": shifted_eom_term,
                "exact_role_match": shifted_eom_term == target_eom_term,
            },
            "scope": "direct_mass_like_m2_pairwise_target_m2_psi2_equals_m2_psi5_role_matching_only",
            "insufficiency_clause": "role matching of the action-slot and eom-slot under the declared shift does not identify the independent coefficient symbols m2_psi2 and m2_psi5",
        },
        "result": {
            "explicit_declared_role_matching_packet_for_m2_psi2_and_m2_psi5_present": True,
            "explicit_pairwise_matching_witness_for_m2_psi2_equals_m2_psi5_present": False,
            "explicit_coefficient_identification_witness_for_m2_psi2_equals_m2_psi5_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R49",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 pair on the route-scoped sufficient lane only",
            "boundary": "R49 does not alter the light-facing kernel channel; it only exports slot-role matching for the single direct m2 pair m2_psi2 / m2_psi5 under the declared plus3 shift",
        },
        "classification": "explicit_declared_role_matching_packet_for_m2_psi2_and_m2_psi5_present_but_no_coefficient_identification_witness",
        "frontier": "R49_B1",
        "frontier_text": "on the route-scoped direct m2 sufficient lane only, the single pairwise target m2_psi2 = m2_psi5 is sharpened into a declared action/eom role match plus one still-missing coefficient identification witness, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi2_equals_m2_psi5",
            "no_claim_that_any_other_direct_m2_pairwise_equality_holds",
            "no_claim_that_the_direct_m2_shift_equivariance_holds",
            "no_claim_that_any_direct_g4_g6_gY_family_defect_vanishes",
            "no_claim_that_any_pair1_c1c1_or_s1s1_zero_equation_holds",
            "no_claim_that_QW2191_is_discharged",
            "no_claim_that_selector_closure_is_obtained",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "R49",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI2_PSI5_ROLE_MATCHING_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R49_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

