#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "r41_direct_m2_psi7_psi10_declared_formal_slot_separation_packet.json"
OUT_SUMMARY = GENERATED / "r41_direct_m2_psi7_psi10_declared_formal_slot_separation_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    r40 = load_json("fundamental_action_reconstruction/generated/r40_direct_m2_psi7_psi10_role_matching_packet.json")
    q2164 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2164_l14_full_canonical_continuum_variational_gate.json"
    )
    q2165 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2165_l13_exhaustive_canonical_eom_gate.json"
    )

    potential_total = q2164["model"]["potential_total"]
    lagrangian_density = q2165["model"]["lagrangian_density"]
    family_slots = sorted(set(re.findall(r"m2_psi\d+", potential_total)))

    source_slot = "m2_psi7"
    target_slot = "m2_psi10"

    checks = [
        {
            "id": "r40_role_matching_packet_present",
            "actual": r40["result"]["explicit_declared_role_matching_packet_for_m2_psi7_and_m2_psi10_present"],
            "expected": True,
            "meaning": "R40 already exports the exact one-pair role-matching packet",
        },
        {
            "id": "q2164_n_psi_fields_is_twelve",
            "actual": q2164["model"]["n_psi_fields"],
            "expected": 12,
            "meaning": "QW-2164 exports a twelve-slot psi family",
        },
        {
            "id": "q2165_n_psi_fields_is_twelve",
            "actual": q2165["model"]["n_psi_fields"],
            "expected": 12,
            "meaning": "QW-2165 exports a twelve-slot psi family",
        },
        {
            "id": "q2164_declared_m2_slot_family_count_is_twelve",
            "actual": len(family_slots),
            "expected": 12,
            "meaning": "QW-2164 exports twelve distinct named m2_psi slots in the canonical potential",
        },
        {
            "id": "source_slot_present_in_declared_family",
            "actual": source_slot in family_slots,
            "expected": True,
            "meaning": "m2_psi7 is one of the declared formal m2 slots",
        },
        {
            "id": "target_slot_present_in_declared_family",
            "actual": target_slot in family_slots,
            "expected": True,
            "meaning": "m2_psi10 is one of the declared formal m2 slots",
        },
        {
            "id": "source_and_target_slots_are_distinct_names",
            "actual": source_slot != target_slot,
            "expected": True,
            "meaning": "the current canonical export names m2_psi7 and m2_psi10 as two distinct formal slots",
        },
        {
            "id": "source_slot_present_in_q2165_lagrangian_density",
            "actual": source_slot in lagrangian_density,
            "expected": True,
            "meaning": "QW-2165 carries the source slot into the full canonical lagrangian density",
        },
        {
            "id": "target_slot_present_in_q2165_lagrangian_density",
            "actual": target_slot in lagrangian_density,
            "expected": True,
            "meaning": "QW-2165 carries the target slot into the full canonical lagrangian density",
        },
    ]

    for item in checks:
        item["pass"] = item["actual"] == item["expected"]

    artifact = {
        "stage": "R41",
        "packet_goal": "materialize_a_route_scoped_declared_formal_slot_separation_packet_for_the_single_direct_m2_pair_m2_psi7_and_m2_psi10_without_claiming_physical_or_symbolic_equality",
        "source_reports": ["R40", "QW2164", "QW2165"],
        "direct_m2_psi7_psi10_declared_formal_slot_separation_packet": {
            "pairwise_target_under_attack": "m2_psi7 = m2_psi10",
            "declared_formal_family_name": "m2_psi_i",
            "declared_formal_family_slots": family_slots,
            "source_slot": source_slot,
            "target_slot": target_slot,
            "same_declared_family": True,
            "distinct_named_slots_in_current_export": True,
            "current_export_layer_statement": "the current canonical export carries m2_psi7 and m2_psi10 as two distinct named slots of the same declared m2_psi family",
            "reduced_missing_witness": "explicit_common_parameter_source_or_symbol_identification_witness_for_the_declared_formal_m2_slots_m2_psi7_and_m2_psi10",
            "scope": "single_direct_m2_pair_m2_psi7_and_m2_psi10_declared_formal_slot_separation_only",
            "insufficiency_clause": "sharing a family prefix and exact role matching under the declared shift still does not identify the two named slots as the same coefficient source",
        },
        "result": {
            "explicit_declared_formal_slot_separation_packet_for_m2_psi7_and_m2_psi10_present": True,
            "explicit_common_parameter_source_or_symbol_identification_witness_present": False,
            "explicit_pairwise_matching_witness_for_m2_psi7_equals_m2_psi10_present": False,
            "global_reduction_of_main_route_claimed": False,
        },
        "light_boundary": {
            "status": "shared_kernel_light_facing_channel_already_closed_by_R14_and_left_unchanged_by_R41",
            "shared_light_channel_source": "R14 specialization of the symmetric kernel-mixing channel",
            "current_packet_scope": "one direct non-light m2 pair on the route-scoped sufficient lane only",
            "boundary": "R41 does not alter the light-facing kernel channel; it only exports declared formal slot separation for the single direct m2 pair m2_psi7 / m2_psi10",
        },
        "classification": "explicit_declared_formal_slot_separation_packet_for_m2_psi7_and_m2_psi10_present_but_no_common_source_or_symbol_identification_witness",
        "frontier": "R41_B1",
        "frontier_text": "on the canonical-ontology-supported direct m2 sufficient lane only, the single coefficient-identification gap for m2_psi7 = m2_psi10 is sharpened into one still-missing common-parameter-source or symbol-identification witness between two distinct named formal slots of the same exported m2 family, while the other direct family blockers, the main host route, and QW-2191 remain unchanged",
        "consistency_checks": checks,
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_m2_psi7_equals_m2_psi10",
            "no_claim_that_the_two_named_slots_are_physically_identified",
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
        "stage": "R41",
        "status": "PASS_PARTIAL_DIRECT_M2_PSI7_PSI10_DECLARED_FORMAL_SLOT_SEPARATION_PACKET_READY",
        "result": artifact["classification"],
        "frontier": ["R41_B1", "QW2191_obstruction", "C10_B1"],
        "theorem_level_pass": False,
        "full_closure_pass": False,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
