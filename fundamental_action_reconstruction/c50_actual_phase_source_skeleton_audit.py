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

c50_b1: str | None = (
    "no_packet_ready_strict_core_minimal_source_skeleton_for_actual_theta_1_theta_2; only_axiom_augmented_source_branch_is_available"
)
if strict_core_theta_export_present:
    c50_b1 = None

summary = {
    "step": "C50",
    "status": "C50_EXECUTED_ACTUAL_PHASE_SOURCE_SKELETON_AUDIT_NO_FALSE_PASS",
    "goal": "Reduce C49_B1 by checking whether strict core already contains a packet-ready minimal source skeleton for actual theta_1, theta_2, or whether only the axiom-augmented branch remains available.",
    "sources": {
        "C31": "source class alpha_12 = theta_2 - theta_1",
        "C33": "formula class theta_i = atan2(<s_i,u_i>,<c_i,u_i>)",
        "C35": "strict-core actual phase source absent, axiom-augmented branch present",
        "C49": "conditional populated-instance schema depends on actual theta_1, theta_2",
        "A10": "anti-overclaim boundary"
    },
    "findings": {
        "strict_core_minimal_source_skeleton": "present_strict_core" if strict_core_theta_export_present else "not_shown",
        "axiom_augmented_source_branch": "present_non_strict_branch",
        "actual_theta_export": "present_strict_core" if strict_core_theta_export_present else "not_shown",
        "actual_populated_instance": "present_strict_core" if strict_core_populated_instance_present else "not_shown",
    },
    "frontier_after_C50": {
        "C50_B1": c50_b1,
        "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12"
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_theta_supply_implies_global_selector_closure",
        "no_claim_that_qw_2191_is_discharged"
    ],
    "next_step": "C51"
}

out = root / "generated" / "c50_actual_phase_source_skeleton_audit_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out)
