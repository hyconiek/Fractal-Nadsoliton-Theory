from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"

SIGMA_INT_STRICT_DERIVED = GENERATED / "sigma_int_strict_derived_v1.json"
GAUGE_QUOTIENT_WITNESS = GENERATED / "sigma_int_gauge_quotient_safety_witness_v1.json"
THETA_PAIR_SIGMA_INT_SLOT_FREE = GENERATED / "theta_pair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1.json"
ASSIGNMENT_DIAGONAL = GENERATED / "mode_index_assignment_canonical_local_diagonal_strict_derived_v1.json"
ASSIGNMENT_SHANNON = GENERATED / "mode_index_assignment_shannon_element_order_reference_strict_core_v1.json"


def load_if_exists(path: Path) -> dict | None:
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def has_object(d: dict | None, obj: str) -> bool:
    return bool(d is not None and d.get("object") == obj)


sigma_int_strict = load_if_exists(SIGMA_INT_STRICT_DERIVED)
gauge_witness = load_if_exists(GAUGE_QUOTIENT_WITNESS)
theta_pair = load_if_exists(THETA_PAIR_SIGMA_INT_SLOT_FREE)
assignment_diag = load_if_exists(ASSIGNMENT_DIAGONAL)
assignment_shannon = load_if_exists(ASSIGNMENT_SHANNON)

o1_present = has_object(sigma_int_strict, "sigma_int_strict_derived_v1")
o2_present = has_object(gauge_witness, "sigma_int_gauge_quotient_safety_witness_v1")
o3_present = has_object(theta_pair, "ThetaPair_sigma_int_strict_selector_ingredient_o2_cut_slot_free_v1")
o4_supported = has_object(
    assignment_diag, "ModeIndexAssignment_canonical_local_diagonal_strict_derived_v1"
) or has_object(assignment_shannon, "ModeIndexAssignment_shannon_element_order_reference_strict_core_v1")

obligation_matrix = {
    "B3_O1": "strict_datum_exported_premise_based" if o1_present else "candidate_identified",
    "B3_O2": "gauge_quotient_safety_discharged_on_declared_domain" if o2_present else "open",
    "B3_O3": "theta_supply_exported_in_declared_scope" if o3_present else "open",
    "B3_O4": "supported_in_declared_audits" if o4_supported else "open",
    "B3_O5": "executed_no_false_pass",
}

payload = {
    "stage": "B8",
    "status": "B8_UPDATED_NO_FALSE_PASS_SELECTOR_TRACK_RESIDUAL_BLOCKERS_SCOPE_NARROWED_AFTER_STRICT_INTERNAL_SHANNON_THETA_SUPPLY_NO_FALSE_PASS",
    "as_of": "2026-03-15",
    "goal": "Run the anti-overclaim audit for the selector track on current repo state (post strict sigma-int export, gauge-quotient safety, theta supply, and axis-only mode-index assignment exports).",
    "inputs": ["F307/N418", "F308/N419", "F451/N489", "F453/N492", "F454/N496", "A6", "A10"],
    "obligation_matrix": obligation_matrix,
    "residual_blockers": [
        "no_global_discharge_of_QW_2191_kernel_alone_uniqueness",
        "no_strict_core_selector_closure_admissible_S_sel_int",
        "residual_Z2_sign_not_lifted_to_sign_sensitive_physical_orientation_datum_or_no_gauge_irrelevance_proof_for_target_observables",
        "no_theorem_level_derivation_of_FR_sign_map_beyond_declared_strict_side_premise",
        "no_scope_extension_beyond_declared_lanes_and_n12_without_new_strict_objects",
    ],
    "forbidden_claims": [
        "B3_packet_closed",
        "B3_O3_pass",
        "B3_O4_pass",
        "QW_2191_discharged",
        "A6_full_uniqueness_closed",
        "theorem_level_selector_derivation",
        "full_ToE_closure"
    ],
    "next_step": "C1"
}

out = ROOT / "generated" / "b8_selector_track_anti_overclaim_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
