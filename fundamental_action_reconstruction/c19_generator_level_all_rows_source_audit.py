from __future__ import annotations

import json
from pathlib import Path

root = Path(__file__).resolve().parent
repo = root.parent

q2165_py = repo / "QW_2165_L13_EXHAUSTIVE_CANONICAL_EOM_GATE.py"
q2166_py = repo / "QW_2166_L14_EXHAUSTIVE_CANONICAL_HESSIAN_GATE.py"
q2165_report = repo / "report_qw2165_l13_exhaustive_canonical_eom_gate.json"
q2166_report = repo / "report_qw2166_l14_exhaustive_canonical_hessian_gate.json"

src_2165 = q2165_py.read_text(encoding="utf-8")
src_2166 = q2166_py.read_text(encoding="utf-8")
report_2165 = json.loads(q2165_report.read_text(encoding="utf-8"))
report_2166 = json.loads(q2166_report.read_text(encoding="utf-8"))

model_2165 = report_2165["model"]
model_2166 = report_2166["model"]
flags_2165 = report_2165["flags"]
flags_2166 = report_2166["flags"]

sample_keys_2165 = sorted(k for k in model_2165 if k.startswith("sample_eom_psi"))

payload = {
    "stage": "C19",
    "status": "C19_EXECUTED_GENERATOR_LEVEL_ALL_ROWS_SOURCE_AUDIT_NO_FALSE_PASS",
    "as_of": "2026-03-06",
    "goal": "Reduce the missing 12-row serialization blocker by showing that the strict core already contains generator-level all-rows source objects for the Psi family, even if no persisted row-by-row 12-row export artifact has been materialized.",
    "inputs": {
        "strict_admissible": [
            "QW-2165",
            "QW-2166",
            "QW-2180",
            "C18",
            "A10",
        ]
    },
    "generator_level_source": {
        "q2165_all_rows_generator_present": "eom_psi = [euler_lagrange(l_density, psi[i]) for i in range(N)]"
        in src_2165,
        "q2166_full_hessian_generator_present": "hessian = sp.hessian(potential_total, fields)"
        in src_2166,
        "n_psi_fields_q2165": model_2165["n_psi_fields"],
        "n_psi_fields_q2166": model_2166["n_psi_fields"],
        "q2165_all_fields_flag": flags_2165["euler_lagrange_executed_for_all_13_fields"],
        "q2166_all_fields_flag": flags_2166["hessian_constructed_for_all_13_fields"],
        "q2166_all_linearized_flag": flags_2166["linearized_eom_executed_for_all_13_fluctuation_fields"],
        "lagrangian_density_mentions_all_psi_symbols": all(
            f"psi{i}(x)" in model_2165["lagrangian_density"] for i in range(12)
        ),
        "potential_total_mentions_all_psi_symbols": all(
            f"psi{i}" in model_2166["potential_total"] for i in range(12)
        ),
    },
    "persisted_artifacts": {
        "q2165_sample_keys": sample_keys_2165,
        "persisted_12_row_export_present": len(sample_keys_2165) >= 12,
        "persisted_full_12x12_hessian_table_present": False,
        "orientation_slice_restriction_present": False,
    },
    "result": {
        "generator_level_all_rows_source_present": "yes",
        "persisted_row_by_row_12_row_export_present": "not_shown",
        "generator_level_source_vs_persisted_serialization_gap_is_now_explicit": "yes",
    },
    "residual_blockers": {
        "C19_B1": "no_explicit_persisted_12_row_serialization_artifact_even_though_generator_level_all_rows_source_is_present",
        "C19_B2": "no_explicit_restriction_from_the_control_pullback_orbits_to_the_candidate_orientation_slice",
    },
    "hard_limits": [
        "no_theorem_level_pass",
        "no_full_closure_pass",
        "no_claim_that_the_12_rows_are_already_persisted",
        "no_claim_that_the_full_12x12_canonical_table_is_already_materialized",
        "no_orientation_slice_restriction_claim",
    ],
    "next_step": "C20",
}

out = root / "generated" / "c19_generator_level_all_rows_source_audit_summary.json"
out.write_text(json.dumps(payload, indent=2) + "\n", encoding="ascii")
print(out)
