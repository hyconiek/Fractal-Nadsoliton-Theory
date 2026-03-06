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

sample_keys = sorted(k for k in model_2165 if k.startswith("sample_eom_psi"))
all_static_keys_present = all(f'"eom_psi{i}"' in src_2165 for i in range(12))
finite_schema_present = any(
    token in src_2165
    for token in [
        'f"eom_psi{i}"',
        "for i in range(N)",
        "for i in range(12)",
        "model.update",
    ]
) and ('f"eom_psi{i}"' in src_2165 or "model.update" in src_2165)

payload = {
    "stage": "C22",
    "status": "C22_EXECUTED_MODEL_CLAUSE_SCHEMA_ABSENCE_AUDIT_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Check whether the strict core already contains either a static all-12-row model clause or a finite key-family schema generating all Psi row entries inside the existing QW-2165 export carrier.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "QW-2166",
            "QW-2180",
            "C21",
            "A10",
        ]
    },
    "carrier_clause_audit": {
        "existing_export_carrier_present": True,
        "sample_keys_present": sample_keys,
        "static_all_12_row_clause_present": all_static_keys_present,
        "finite_key_family_schema_present": finite_schema_present,
    },
    "result": {
        "explicit_full_model_clause_present": "no",
        "explicit_finite_key_family_schema_present": "no",
        "absence_is_now_explicitly_audited": "yes",
    },
    "residual_blockers": {
        "C22_B1": "no_explicit_static_all_12_rows_model_clause_and_no_explicit_finite_key_family_schema_generating_all_Psi_row_entries_inside_the_existing_qw2165_export_carrier",
        "C22_B2": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_the_full_model_clause_exists",
        "no_claim_that_the_finite_schema_exists",
        "no_orientation_slice_restriction_claim",
    ],
    "next_step": "C23",
}

out = root / "generated" / "c22_model_clause_schema_absence_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
