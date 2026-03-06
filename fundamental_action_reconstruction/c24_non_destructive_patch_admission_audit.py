from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent

q2165_py = repo / "QW_2165_L13_EXHAUSTIVE_CANONICAL_EOM_GATE.py"
src_2165 = q2165_py.read_text(encoding="utf-8")

admission_conditions = {
    "existing_export_carrier_present": 'OUT_JSON = ROOT / "report_qw2165_l13_exhaustive_canonical_eom_gate.json"'
    in src_2165,
    "existing_model_clause_present": '"model": {' in src_2165,
    "finite_family_size_present": "N = 12" in src_2165,
    "all_rows_family_present": "eom_psi = [euler_lagrange(l_density, psi[i]) for i in range(N)]"
    in src_2165,
    "patch_candidate_is_serialization_only": True,
    "patch_candidate_does_not_change_lagrangian_or_eom_definitions": True,
    "anti_overclaim_boundary_retained": True,
}

admission_allowed = all(admission_conditions.values())

payload = {
    "stage": "C24",
    "status": "C24_EXECUTED_NON_DESTRUCTIVE_PATCH_ADMISSION_AUDIT_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Check whether a minimal serialization-only patch-candidate for the full 12-row model clause is admissible as a non-destructive next step, without claiming that the patch has already been applied or rerun.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "C23",
            "A10",
        ]
    },
    "admission_conditions": admission_conditions,
    "result": {
        "non_destructive_patch_candidate_admission_allowed": admission_allowed,
        "patch_applied": "not_shown",
        "rerun_with_full_12_row_output": "not_shown",
    },
    "residual_blockers": {
        "C24_B1": "no_applied_patch_candidate_and_no_rerun_validated_report_for_the_full_12_row_model_clause_even_though_non_destructive_patch_admission_is_allowed",
        "C24_B2": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_the_patch_has_been_applied",
        "no_claim_that_the_rerun_has_been_executed",
        "no_orientation_slice_restriction_claim",
    ],
    "next_step": "C25",
}

out = root / "generated" / "c24_non_destructive_patch_admission_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
