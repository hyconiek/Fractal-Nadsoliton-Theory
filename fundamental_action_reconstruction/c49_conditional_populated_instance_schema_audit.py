from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
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

strict_core_populated_instance_present = has_object(
    r1_pop_sigma_int, "R1_residual_orientation_datum_target_slot_population_strict_derived_from_sigma_int_slot_free_theta_pair_v1"
) or has_object(
    r1_pop_diag, "R1_residual_orientation_datum_target_slot_population_strict_derived_from_canonical_local_diagonal_theta_pair_v1"
)

c49_b1 = (
    "no_strict_core_supplied_actual_theta_1_theta_2_values_for_instantiating_the_now_packet_ready_conditional_populated_instance_schema_of_u_1_u_2_and_S_orient_cand"
)
if strict_core_theta_export_present:
    c49_b1 = (
        "DISCHARGED_ON_CURRENT_REPO_STATE: strict-core theta supply is exported and the conditional schema is instantiable in declared "
        "pair1/pair2 scopes (residual Z2 sign remains explicit)."
    )

summary = {
    "step": "C49",
    "status": "C49_EXECUTED_CONDITIONAL_POPULATED_INSTANCE_SCHEMA_AUDIT_NO_FALSE_PASS",
    "goal": "Reduce C48_B1 by checking whether strict core already contains a packet-ready conditional populated-instance schema for u_1,u_2 and the candidate orientation slice, conditional on theta_1,theta_2.",
    "sources": {
        "C34": "class formulas for u_i(theta_i)",
        "C35": "strict-core theta export status is checked against current generated artifacts",
        "C47": "class-level candidate orientation slice",
        "C48": "packet-ready minimal export skeleton",
        "A10": "anti-overclaim boundary"
    },
    "conditional_schema": {
        "inputs": ["theta_1", "theta_2"],
        "u_1": "cos(theta_1)c_1 + sin(theta_1)s_1",
        "u_2": "cos(theta_2)c_2 + sin(theta_2)s_2",
        "orientation_slice_candidate": "span{u_1,u_2}"
    },
    "findings": {
        "conditional_populated_instance_schema": "present_partial",
        "actual_theta_supply": "present_strict_core" if strict_core_theta_export_present else "not_shown",
        "actual_populated_instance": "present_strict_core" if strict_core_populated_instance_present else "not_shown",
    },
    "frontier_after_C49": {
        "C49_B1": c49_b1,
        "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12"
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_theta_supply_implies_global_selector_closure",
        "no_claim_that_qw_2191_is_discharged"
    ],
    "next_step": "C50"
}

out = root / "generated" / "c49_conditional_populated_instance_schema_audit_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out)
