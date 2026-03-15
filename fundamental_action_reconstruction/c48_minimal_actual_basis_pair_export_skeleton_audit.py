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

population_dependency = "blocked_by_C35_B1"
if strict_core_theta_export_present:
    population_dependency = "no_longer_blocked_by_C35_B1__theta_supply_present"

c48_b1 = (
    "no_explicit_populated_actual_basis_pair_export_instance_even_though_a_minimal_export_skeleton_for_u_1_u_2_is_now_packet_ready; "
    "population_remains_blocked_by_C35_B1"
)
if strict_core_basis_pair_export_present:
    c48_b1 = (
        "DISCHARGED_ON_CURRENT_REPO_STATE: populated u_1,u_2 export instances exist (as strict-core inhabitant instances and/or "
        "theta-pair outputs). This does not promote global selector closure."
    )

summary = {
    "step": "C48",
    "status": "C48_EXECUTED_MINIMAL_ACTUAL_BASIS_PAIR_EXPORT_SKELETON_AUDIT_NO_FALSE_PASS",
    "goal": "Reduce C47_B1 by checking whether strict core already contains a packet-ready minimal export skeleton for the actual basis pair u_1,u_2 spanning the candidate orientation slice.",
    "sources": {
        "C34": "representative class u_i(theta_i)=cos(theta_i)c_i+sin(theta_i)s_i",
        "C35": "strict-core theta export status is checked against current generated artifacts",
        "C40": "minimal field list discipline",
        "C41": "minimal schema artifact discipline",
        "C47": "class-level candidate orientation slice span{u_1(theta_1),u_2(theta_2)}",
        "A10": "anti-overclaim boundary",
    },
    "skeleton": {
        "u_1_formula": "cos(theta_1)c_1 + sin(theta_1)s_1",
        "u_2_formula": "cos(theta_2)c_2 + sin(theta_2)s_2",
        "required_inputs": ["theta_1", "theta_2"],
        "normalization": "class-level ensured",
        "target_role": "basis pair spanning S_orient_cand"
    },
    "findings": {
        "minimal_export_skeleton": "present_partial",
        "actual_theta_export": "present_strict_core" if strict_core_theta_export_present else "not_shown",
        "actual_basis_pair_export": "present_strict_core" if strict_core_basis_pair_export_present else "not_shown",
        "population_dependency": population_dependency,
    },
    "frontier_after_C48": {
        "C48_B1": c48_b1,
        "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12"
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_basis_pair_export_implies_global_selector_closure",
        "no_claim_that_qw_2191_is_discharged"
    ],
    "next_step": "C49"
}

out = root / "generated" / "c48_minimal_actual_basis_pair_export_skeleton_audit_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out)
