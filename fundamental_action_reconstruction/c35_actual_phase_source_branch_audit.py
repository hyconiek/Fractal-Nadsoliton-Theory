#!/usr/bin/env python3
"""C35: actual phase source branch audit without false PASS."""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = ROOT / "generated" / "c35_actual_phase_source_branch_audit_summary.json"

THETA_PAIR_SIGMA_INT_SLOT_FREE = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
THETA_PAIR_CANONICAL_LOCAL_DIAGONAL = GENERATED / "theta_pair_canonical_local_diagonal_strict_derived_v1.json"
R1_POPULATION_SIGMA_INT_SLOT_FREE = (
    GENERATED / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"
)
R1_POPULATION_CANONICAL_LOCAL_DIAGONAL = (
    GENERATED
    / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_canonical_local_diagonal_theta_pair_v1.json"
)


def load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def load_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    return load(path)


def main() -> None:
    q2192 = load(REPO / "report_qw2192_mode_index_selection_axiom_gate.json")
    q2193 = load(REPO / "report_qw2193_selection_axiom_family_robustness_gate.json")

    theta1_num = q2192["numeric_audit"]["theta_pair1_numeric_argmin"]
    theta2_num = q2192["numeric_audit"]["theta_pair2_numeric_argmin"]
    family_rows = q2193["numeric_audit"]["rows"]

    theta_pair_sigma_int_slot_free = load_if_exists(THETA_PAIR_SIGMA_INT_SLOT_FREE)
    theta_pair_canonical_local_diagonal = load_if_exists(THETA_PAIR_CANONICAL_LOCAL_DIAGONAL)
    r1_population_sigma_int_slot_free = load_if_exists(R1_POPULATION_SIGMA_INT_SLOT_FREE)
    r1_population_canonical_local_diagonal = load_if_exists(R1_POPULATION_CANONICAL_LOCAL_DIAGONAL)

    strict_core_theta_source: dict | None = None
    strict_core_theta_source_artifact = None
    if (
        theta_pair_sigma_int_slot_free is not None
        and theta_pair_sigma_int_slot_free.get("object")
        == "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
    ):
        strict_core_theta_source = theta_pair_sigma_int_slot_free
        strict_core_theta_source_artifact = str(THETA_PAIR_SIGMA_INT_SLOT_FREE.relative_to(REPO))
    elif (
        theta_pair_canonical_local_diagonal is not None
        and theta_pair_canonical_local_diagonal.get("object") == "ThetaPair_canonical_local_diagonal_strict_derived_v1"
    ):
        strict_core_theta_source = theta_pair_canonical_local_diagonal
        strict_core_theta_source_artifact = str(THETA_PAIR_CANONICAL_LOCAL_DIAGONAL.relative_to(REPO))

    strict_core_actual_phase_source = "not_shown"
    strict_core_theta_export_status = "not_shown"
    strict_core_theta_pair_values: dict | None = None
    if strict_core_theta_source is not None:
        strict_core_actual_phase_source = "present_strict_core_theta_pair_source"
        strict_core_theta_export_status = "present_strict_core"
        outputs = strict_core_theta_source.get("outputs") or {}
        if isinstance(outputs, dict):
            pair1 = outputs.get("pair1") or {}
            pair2 = outputs.get("pair2") or {}
            if isinstance(pair1, dict) and isinstance(pair2, dict):
                strict_core_theta_pair_values = {
                    "pair1": {
                        k: pair1.get(k)
                        for k in ("m", "theta_1", "theta_1_base_mod_pi", "theta_1_star")
                        if k in pair1
                    },
                    "pair2": {
                        k: pair2.get(k)
                        for k in ("m", "theta_2", "theta_2_base_mod_pi", "theta_2_star")
                        if k in pair2
                    },
                }

    strict_core_orientation_slice_span_present = bool(
        (r1_population_sigma_int_slot_free is not None and r1_population_sigma_int_slot_free.get("object"))
        or (r1_population_canonical_local_diagonal is not None and r1_population_canonical_local_diagonal.get("object"))
    )

    c35_b1 = (
        "no_strict_core_export_of_actual_local_phase_coordinates_theta_1_theta_2_for_the_actual_pair_frames; only an "
        "axiom_augmented_source_branch_theta_star_equals_0_is_currently_available"
    )
    if strict_core_actual_phase_source != "not_shown":
        c35_b1 = (
            "DISCHARGED_ON_CURRENT_REPO_STATE: strict-core theta-pair supply is exported (sigma-int slot-free lane and/or "
            "canonical local-diagonal lane). Residual Z2 sign remains explicit; no implied global selector closure."
        )

    summary = {
        "stage": "C35",
        "status": "C35_EXECUTED_ACTUAL_PHASE_SOURCE_BRANCH_AUDIT_NO_FALSE_PASS",
        "as_of": "2026-03-15",
        "goal": "Reduce C34_B1 by checking whether any packet-ready source branch for actual local phases theta_1, theta_2 already exists in the repo, and whether it belongs to strict core or only to the axiom-augmented branch.",
        "inputs": {
            "strict_admissible": ["C31", "C34", "A10"],
            "contrast_non_strict": ["QW-2192", "QW-2193"]
        },
        "strict_core_state": {
            "local_phase_formula_class_present": "yes_partial",
            "local_reduced_representative_class_present": "yes_partial",
            "actual_theta_1_theta_2_export_for_current_pair_frames": strict_core_theta_export_status,
            "strict_core_theta_source_artifact": strict_core_theta_source_artifact,
            "strict_core_theta_pair_values": strict_core_theta_pair_values,
        },
        "axiom_augmented_branch": {
            "selection_axiom_name": q2192["selection_axiom"]["name"],
            "theta_pair1_numeric_argmin": theta1_num,
            "theta_pair2_numeric_argmin": theta2_num,
            "family_rows_count": len(family_rows),
            "all_family_members_select_theta0_pair1": q2193["flags"]["all_family_members_select_theta0_pair1"],
            "all_family_members_select_theta0_pair2": q2193["flags"]["all_family_members_select_theta0_pair2"],
            "actual_phase_source_branch_present": "yes_non_strict_branch"
        },
        "result": {
            "strict_core_actual_phase_source": strict_core_actual_phase_source,
            "axiom_augmented_actual_phase_source_branch": "present_partial",
            "raw_cross_pair_overlap_scalar_route": "blocked_by_degeneracy",
            "final_basis_level_slice_extraction_present": (
                "present_partial_span_object" if strict_core_orientation_slice_span_present else "not_shown"
            ),
        },
        "residual_blockers": {
            "C35_B1": c35_b1,
            "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
            "C26_B2": "no_explicit_basis_level_embedding_or_extraction_of_the_candidate_two_dimensional_orientation_slice_inside_that_reduced_plane"
        },
        "hard_limits": [
            "no_theorem_level_pass",
            "no_full_closure_pass",
            "no_claim_that_axiom_augmented_branch_discharges_strict_core_blocker",
            "no_claim_that_u_1_u_2_are_materialized_in_strict_core",
            "no_claim_that_final_slice_extraction_is_resolved"
        ],
        "next_step": "C47"
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
