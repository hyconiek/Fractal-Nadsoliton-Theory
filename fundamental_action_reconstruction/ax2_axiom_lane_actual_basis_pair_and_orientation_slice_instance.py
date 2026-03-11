from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent

target = root / "generated" / "axiom_lane_actual_basis_pair_orientation_slice_instance.json"

instance = {
    "lane": "axiom-augmented",
    "axiom": "minimum_harmonic_alignment_with_orientation_convention",
    "source_qw": ["QW-2192", "QW-2193"],
    "theta_1": "0_mod_2pi",
    "theta_2": "0_mod_2pi",
    "u_1": "c_1",
    "u_2": "c_2",
    "orientation_slice": "span{c_1,c_2}",
    "strict_core_status": "not_in_strict_core",
    "forbidden_claims": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_strict_core_discharge",
        "no_qw_2191_discharge",
    ],
}

target.write_text(json.dumps(instance, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX2",
    "status": "AX2_EXECUTED_AXIOM_LANE_ACTUAL_BASIS_PAIR_AND_ORIENTATION_SLICE_INSTANCE_NO_FALSE_PASS",
    "goal": "Materialize the actual basis pair and actual orientation slice on the explicit axiom-augmented selector lane opened by AX1.",
    "created_file": {
        "relative_path": "generated/axiom_lane_actual_basis_pair_orientation_slice_instance.json",
        "exists_after_step": target.exists(),
        "content_keys": list(instance.keys()),
    },
    "result": {
        "actual_theta_values_available": "yes_axiom_lane_only",
        "actual_basis_pair_available": "yes_axiom_lane_only",
        "actual_orientation_slice_available": "yes_axiom_lane_only",
        "strict_core_changed": False,
    },
    "residual_frontier": {
        "T12_B1": "the typing judgment with totality and uniqueness is specified but not discharged for the current selector track",
        "T2_B1": "the bridge theorem is specified but not discharged; strict-core target slot and sign-only export-map object exist, but target-slot population (theta_1,theta_2) remains absent",
        "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_axiom_lane_equals_strict_core",
        "no_claim_that_qw_2191_is_discharged",
    ],
    "next_step": "AX3",
}

out = root / "generated" / "ax2_axiom_lane_actual_basis_pair_and_orientation_slice_instance_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out)
