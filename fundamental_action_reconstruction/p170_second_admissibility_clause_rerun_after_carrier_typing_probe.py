#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p170_second_admissibility_clause_rerun_after_carrier_typing_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n188 = load_json(
        "fundamental_action_reconstruction/generated/n188_current_first_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    f82 = load_json(
        "fundamental_action_reconstruction/generated/f82_first_exported_preobserver_source_object_second_clause_typing_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "first_clause_already_discharged",
            "actual": n188["theorem_result"]["discharged"],
            "expected": True,
        },
        {
            "id": "typed_carrier_decomposition_present",
            "actual": f82["carrier"]["typed_sum"],
            "expected": ["V_topo", "L_int", "M_int"],
        },
        {
            "id": "topological_seed_slot_present",
            "actual": f82["carrier"]["topological_seed_slot"],
            "expected": "u_T",
        },
        {
            "id": "light_transport_slot_present",
            "actual": f82["carrier"]["light_transport_slot"],
            "expected": "u_L",
        },
        {
            "id": "matter_encoding_slot_present",
            "actual": f82["carrier"]["matter_encoding_slot"],
            "expected": "u_M",
        },
        {
            "id": "nonzero_light_support",
            "actual": f82["state_support"]["nonzero_light_support"],
            "expected": True,
        },
        {
            "id": "future_orientation_export_target_frame",
            "actual": f82["future_orientation_export_target_frame"],
            "expected": ["u_T", "u_L"],
        },
        {
            "id": "E_orient_not_yet_exported",
            "actual": f82["E_orient_exported"],
            "expected": False,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P170",
            "lane": "second_admissibility_clause_rerun_after_carrier_typing_only",
            "status": "P170_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SECOND_CLAUSE_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P170",
            "lane": "second_admissibility_clause_rerun_after_carrier_typing_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_SECOND_ADMISSIBILITY_CLAUSE_AFTER_P170",
            "exported_object": "S_preLM_strict_core_source_object_v1",
            "second_clause": "carrier_typed_enough_for_later_export",
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "source_seed_only",
                "non_substitutive",
                "selector_acceptance_independent",
                "future_bridge_compatible",
            ],
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
