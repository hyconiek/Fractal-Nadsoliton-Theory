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

sample_keys_2165 = sorted(k for k in model_2165 if k.startswith("sample_eom_psi"))

payload = {
    "stage": "C21",
    "status": "C21_EXECUTED_EXISTING_EXPORT_CARRIER_AUDIT_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce the missing executed serialization-run blocker by showing that the strict core already contains an existing persisted export carrier for the Psi family, even if the model payload still serializes only three sample rows instead of all 12 rows.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "QW-2166",
            "QW-2180",
            "C20",
            "A10",
        ]
    },
    "existing_export_carrier": {
        "out_json_path_present": 'OUT_JSON = ROOT / "report_qw2165_l13_exhaustive_canonical_eom_gate.json"'
        in src_2165,
        "json_write_present": "OUT_JSON.write_text(json.dumps(out, ensure_ascii=False, indent=2), encoding=\"utf-8\")"
        in src_2165,
        "existing_model_payload_present": '"model": {' in src_2165,
        "sample_row_serialization_keys_present": all(
            key in src_2165
            for key in ['"sample_eom_psi0"', '"sample_eom_psi6"', '"sample_eom_psi11"']
        ),
        "sample_key_count_in_persisted_report": len(sample_keys_2165),
    },
    "result": {
        "persisted_export_carrier_present": "yes",
        "all_rows_model_serialization_clause_present": "not_shown",
        "frontier_can_be_reduced_from_missing_executed_run_to_missing_full_model_clause_inside_existing_export_carrier": "yes",
    },
    "residual_blockers": {
        "C21_B1": "no_explicit_all_12_rows_model_serialization_clause_inside_the_already_existing_qw2165_persisted_export_carrier",
        "C21_B2": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_the_existing_export_carrier_already_serializes_all_12_rows",
        "no_claim_that_the_full_12x12_table_is_already_present",
        "no_orientation_slice_restriction_claim",
    ],
    "next_step": "C22",
}

out = root / "generated" / "c21_existing_export_carrier_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
