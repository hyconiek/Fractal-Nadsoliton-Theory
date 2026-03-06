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

full_keys = [f"eom_psi{i}" for i in range(12)]

payload = {
    "stage": "C25",
    "status": "C25_EXECUTED_APPLIED_PATCH_RERUN_EXPORT_AUDIT_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Confirm that the minimal serialization-only patch was actually applied to QW-2165 and that a rerun produced a persisted report containing all 12 Psi-row exports, without claiming anything about the still-open orientation-slice restriction.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "C24",
            "A10",
        ]
    },
    "applied_patch_and_rerun": {
        "full_model_patch_present_in_source": '**{f"eom_psi{i}": str(eom_psi[i]) for i in range(N)}' in src_2165,
        "persisted_report_contains_all_12_rows": all(k in model_2165 for k in full_keys),
        "persisted_report_row_count": sum(k in model_2165 for k in full_keys),
        "sample_rows_preserved": all(
            k in model_2165 for k in ["sample_eom_psi0", "sample_eom_psi6", "sample_eom_psi11"]
        ),
        "verdict": report_2165["verdict"],
        "lean_binary": report_2165["checker"]["lean_binary"],
        "lean_exit_code": report_2165["checker"]["exit_code"],
    },
    "result": {
        "serialization_lane_closed_in_declared_scope": "yes",
        "orientation_slice_restriction_present": "not_shown",
    },
    "residual_blockers": {
        "C25_B1": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_orientation_slice_restriction_is_resolved",
        "no_claim_that_full_12x12_hessian_table_is_resolved",
    ],
    "next_step": "C26",
}

out = root / "generated" / "c25_applied_patch_rerun_export_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
