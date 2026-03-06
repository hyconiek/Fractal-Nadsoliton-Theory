from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent

q2165_py = repo / "QW_2165_L13_EXHAUSTIVE_CANONICAL_EOM_GATE.py"
q2165_report = repo / "report_qw2165_l13_exhaustive_canonical_eom_gate.json"

src_2165 = q2165_py.read_text(encoding="utf-8")
report_2165 = json.loads(q2165_report.read_text(encoding="utf-8"))
model_2165 = report_2165["model"]
flags_2165 = report_2165["flags"]

sample_keys_2165 = sorted(k for k in model_2165 if k.startswith("sample_eom_psi"))

payload = {
    "stage": "C20",
    "status": "C20_EXECUTED_FINITE_MATERIALIZATION_RECIPE_AUDIT_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce the persisted-serialization blocker by showing that the strict core already contains a finite, persisted materialization recipe for all 12 Psi rows, even if no executed 12-row persisted serialization run has been stored yet.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "QW-2166",
            "QW-2180",
            "C19",
            "A10",
        ]
    },
    "finite_materialization_recipe": {
        "finite_family_size_n_present": "N = 12" in src_2165,
        "psi_family_constructor_present": 'psi = [sp.Function(f"psi{i}")(x) for i in range(N)]'
        in src_2165,
        "persisted_lagrangian_density_present": "lagrangian_density" in model_2165,
        "euler_lagrange_constructor_present": "def euler_lagrange(lagr, f):" in src_2165,
        "row_family_materialization_comprehension_present": "eom_psi = [euler_lagrange(l_density, psi[i]) for i in range(N)]"
        in src_2165,
        "all_fields_execution_flag_present": flags_2165["euler_lagrange_executed_for_all_13_fields"],
    },
    "persisted_outputs": {
        "sample_keys": sample_keys_2165,
        "persisted_12_row_serialization_run_present": len(sample_keys_2165) >= 12,
        "persisted_full_12_row_table_present": False,
        "orientation_slice_restriction_present": False,
    },
    "result": {
        "finite_persisted_materialization_recipe_present": "yes",
        "executed_persisted_12_row_serialization_run_present": "not_shown",
        "frontier_can_be_reduced_from_missing_recipe_to_missing_executed_serialization_run": "yes",
    },
    "residual_blockers": {
        "C20_B1": "no_explicit_executed_and_persisted_12_row_serialization_run_from_the_already_present_finite_materialization_recipe",
        "C20_B2": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_the_12_row_serialization_run_has_already_been_executed_and_persisted",
        "no_claim_that_the_full_12x12_table_is_already_present",
        "no_orientation_slice_restriction_claim",
    ],
    "next_step": "C21",
}

out = root / "generated" / "c20_finite_materialization_recipe_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
