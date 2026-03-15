from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent
generated = root / "generated"

THETA_PAIR_SIGMA_INT_SLOT_FREE = generated / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
THETA_PAIR_CANONICAL_LOCAL_DIAGONAL = generated / "theta_pair_canonical_local_diagonal_strict_derived_v1.json"
R1_POPULATION_SIGMA_INT_SLOT_FREE = (
    generated / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1.json"
)
R1_POPULATION_CANONICAL_LOCAL_DIAGONAL = (
    generated / "r1_residual_orientation_datum_target_slot_population_strict_derived_from_canonical_local_diagonal_theta_pair_v1.json"
)


def load_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def has_object(d: dict | None, obj: str) -> bool:
    return bool(d is not None and d.get("object") == obj)


theta_pair_sigma_int = load_if_exists(THETA_PAIR_SIGMA_INT_SLOT_FREE)
theta_pair_diag = load_if_exists(THETA_PAIR_CANONICAL_LOCAL_DIAGONAL)
r1_pop_sigma_int = load_if_exists(R1_POPULATION_SIGMA_INT_SLOT_FREE)
r1_pop_diag = load_if_exists(R1_POPULATION_CANONICAL_LOCAL_DIAGONAL)

strict_core_theta_export_present = has_object(
    theta_pair_sigma_int, "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1"
) or has_object(theta_pair_diag, "ThetaPair_canonical_local_diagonal_strict_derived_v1")

strict_core_basis_pair_export_present = has_object(
    r1_pop_sigma_int, "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
) or has_object(
    r1_pop_diag, "R1_residual_orientation_datum_target_slot_population_strict_derived_from_canonical_local_diagonal_theta_pair_v1"
)

materialization_dependency = "blocked_by_C35_B1"
if strict_core_theta_export_present:
    materialization_dependency = "no_longer_blocked_by_C35_B1__theta_supply_present"

c47_b1 = (
    "no_explicit_export_of_actual_normalized_basis_pair_u_1_u_2_spanning_the_candidate_two_dimensional_orientation_slice_inside_the_reduced_plane; "
    "materialization_remains_blocked_by_C35_B1"
)
if strict_core_basis_pair_export_present:
    c47_b1 = (
        "DISCHARGED_ON_CURRENT_REPO_STATE: actual basis pair u_1,u_2 and the orientation-slice span object are exported on the "
        "sigma-int slot-free and/or canonical local-diagonal lanes. This does not imply global selector closure."
    )

summary = {
    "step": "C47",
    "status": "C47_EXECUTED_BASIS_LEVEL_ORIENTATION_SLICE_CANDIDATE_AUDIT_NO_FALSE_PASS",
    "goal": "Reduce C26_B2 by checking whether strict core already contains a packet-ready class-level basis candidate for the two-dimensional orientation slice inside the reduced plane.",
    "sources": {
        "C26": "residual restriction blocker split",
        "C28": "local quotient schema",
        "C29": "serialized local projectors",
        "C33": "local phase export formula class",
        "C34": "local reduced representative class",
        "C35": "strict-core theta export status is checked against current generated artifacts",
        "A10": "anti-overclaim boundary",
    },
    "candidate_class": {
        "local_representative": "u_i(theta_i)=cos(theta_i)c_i+sin(theta_i)s_i",
        "candidate_orientation_slice": "S_orient_cand(theta_1,theta_2)=span{u_1(theta_1),u_2(theta_2)}",
        "projector_compatibility": [
            "P_red(theta_i)u_i=u_i",
            "P_tan(theta_i)u_i=0",
        ],
    },
    "findings": {
        "class_level_basis_candidate": "present_partial",
        "actual_theta_export": "present_strict_core" if strict_core_theta_export_present else "not_shown",
        "actual_basis_pair_export": "present_strict_core" if strict_core_basis_pair_export_present else "not_shown",
        "materialization_dependency": materialization_dependency,
    },
    "frontier_after_C47": {
        "C47_B1": c47_b1,
        "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_theta_export_implies_global_selector_closure",
        "no_claim_that_qw_2191_is_discharged",
    ],
    "next_step": "C48",
}

out = root / "generated" / "c47_basis_level_orientation_slice_candidate_audit_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out)
