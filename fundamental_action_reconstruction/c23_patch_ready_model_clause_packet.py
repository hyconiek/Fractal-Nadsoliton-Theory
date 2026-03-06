from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent

q2165_py = repo / "QW_2165_L13_EXHAUSTIVE_CANONICAL_EOM_GATE.py"
src_2165 = q2165_py.read_text(encoding="utf-8")

patch_ready_schema = '\n'.join(
    [
        '"model": {',
        '    "n_psi_fields": N,',
        '    "lagrangian_density": str(l_density),',
        '    "eom_phi": str(eom_phi),',
        '    **{f"eom_psi{i}": str(eom_psi[i]) for i in range(N)},',
        '},',
    ]
)

payload = {
    "stage": "C23",
    "status": "C23_EXECUTED_PATCH_READY_MODEL_CLAUSE_PACKET_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce the missing model-clause-schema blocker by constructing a minimal patch-ready all-12-rows model clause for QW-2165, without claiming that the patch has already been applied or rerun.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "QW-2166",
            "QW-2180",
            "C22",
            "A10",
        ]
    },
    "patch_ready_packet": {
        "existing_model_clause_present": '"model": {' in src_2165,
        "finite_family_size_present": "N = 12" in src_2165,
        "eom_family_present": "eom_psi = [euler_lagrange(l_density, psi[i]) for i in range(N)]"
        in src_2165,
        "patch_ready_schema_present_in_packet": True,
        "patch_ready_schema": patch_ready_schema,
    },
    "result": {
        "minimal_patch_ready_schema_packet_present": "yes",
        "patch_applied_inside_qw2165": "not_shown",
        "rerun_and_persisted_full_12_row_export": "not_shown",
    },
    "residual_blockers": {
        "C23_B1": "no_applied_and_rerun_materialization_of_the_patch_ready_all_12_rows_model_clause_inside_qw2165_export_carrier",
        "C23_B2": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_the_patch_has_been_applied",
        "no_claim_that_qw2165_has_been_rerun_with_full_12_row_export",
        "no_orientation_slice_restriction_claim",
    ],
    "next_step": "C24",
}

out = root / "generated" / "c23_patch_ready_model_clause_packet_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
