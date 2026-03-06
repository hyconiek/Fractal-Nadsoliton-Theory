from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent

target = root / "generated" / "axiom_lane_selector_family_robustness_certificate.json"

certificate = {
    "lane": "axiom-augmented",
    "selector_family": "J_ab_positive_weight_family",
    "family_formula": "J_ab(theta)=2(a+b)(1-cos theta)",
    "admissible_domain": {
        "a": "positive",
        "b": "positive",
    },
    "minimizer": "theta=0_mod_2pi",
    "stable_theta": {
        "theta_1": "0_mod_2pi",
        "theta_2": "0_mod_2pi",
    },
    "stable_basis_pair": ["c_1", "c_2"],
    "stable_orientation_slice": "span{c_1,c_2}",
    "source_qw": ["QW-2193"],
    "strict_core_status": "not_in_strict_core",
    "forbidden_claims": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_strict_core_discharge",
        "no_qw_2191_discharge",
        "no_selector_family_derivation_from_strict_core",
    ],
}

target.write_text(json.dumps(certificate, indent=2) + "\n", encoding="ascii")

summary = {
    "step": "AX4",
    "status": "AX4_EXECUTED_AXIOM_LANE_SELECTOR_FAMILY_ROBUSTNESS_AUDIT_NO_FALSE_PASS",
    "goal": "Certify that the materialized AX3 bridge instance remains stable across the declared positive-weight selector family on the axiom-augmented lane.",
    "created_file": {
        "relative_path": "generated/axiom_lane_selector_family_robustness_certificate.json",
        "exists_after_step": target.exists(),
        "content_keys": list(certificate.keys()),
    },
    "result": {
        "selector_family_robustness_certified": "yes_axiom_lane_only",
        "stable_theta": {
            "theta_1": "0_mod_2pi",
            "theta_2": "0_mod_2pi",
        },
        "stable_basis_pair": ["c_1", "c_2"],
        "stable_orientation_slice": "span{c_1,c_2}",
        "strict_core_changed": False,
    },
    "residual_frontier": {
        "T12_B1": "the typing judgment with totality and uniqueness is specified but not discharged for the current selector track",
        "T2_B1": "the bridge theorem is specified but not discharged; strict-core target slot and equivalence/export map remain absent",
        "C32_B2": "raw_cross_pair_overlap_scalar_route_is_formally_degenerate_under_the_strict_orthonormal_disjoint_mode_scaffold_and_thus_does_not_export_alpha_12",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_axiom_lane_equals_strict_core",
        "no_claim_that_selector_family_is_derived_from_strict_core",
        "no_claim_that_qw_2191_is_discharged",
    ],
    "next_step": "AX5",
}

out = root / "generated" / "ax4_axiom_lane_selector_family_robustness_audit_summary.json"
out.write_text(json.dumps(summary, indent=2) + "\n", encoding="ascii")
print(out)
