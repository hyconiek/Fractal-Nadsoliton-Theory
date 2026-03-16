#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = (
    GENERATED
    / "p649_current_second_admissibility_clause_rerun_after_carrier_typing_for_s_sel_int_strict_core_source_object_v1_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p649_current_second_admissibility_clause_rerun_after_carrier_typing_for_s_sel_int_strict_core_source_object_v1_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    n540 = load_json(
        "fundamental_action_reconstruction/generated/n540_current_first_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    f649 = load_json(
        "fundamental_action_reconstruction/generated/f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet.json"
    )
    f649_summary = load_json(
        "fundamental_action_reconstruction/generated/f649_first_exported_s_sel_int_strict_core_source_object_second_clause_typing_packet_summary.json"
    )

    exported_object = str(
        f649_summary.get("exported_object")
        or f649.get("exported_object")
        or "S_sel_int_strict_core_source_object_v1"
    )

    carrier = f649.get("carrier") or {}
    basis = (carrier.get("basis_functions_by_x") or {}) if isinstance(carrier, dict) else {}
    c1 = basis.get("c1") if isinstance(basis, dict) else None
    s1 = basis.get("s1") if isinstance(basis, dict) else None

    second_clause = "carrier_typed_enough_for_later_E_orient_export_required"

    checks_spec = [
        {
            "id": "f34_second_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"][second_clause],
            "expected": True,
            "meaning": "F34 keeps the carrier-typing clause active for admissible S_sel_int",
        },
        {
            "id": "first_clause_already_discharged",
            "actual": bool(n540["theorem_result"]["discharged"]),
            "expected": True,
            "meaning": "N540 discharges the first admissibility clause for S_sel_int on current repo state",
        },
        {
            "id": "f649_typing_packet_status",
            "actual": f649_summary.get("status"),
            "expected": "F649_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SECOND_CLAUSE_TYPING_PACKET_NO_FALSE_PASS",
            "meaning": "F649 typing packet summary status matches the executed no-false-pass marker",
        },
        {
            "id": "exported_object_identity",
            "actual": exported_object,
            "expected": "S_sel_int_strict_core_source_object_v1",
            "meaning": "the typing packet targets the strict-core exported S_sel_int source object v1 identity",
        },
        {
            "id": "typed_carrier_decomposition_present",
            "actual": f649_summary["carrier"]["typed_sum"],
            "expected": ["L2_Z12_real_v1"],
            "meaning": "a minimal explicit typed carrier is frozen (L2(Z_12) as a real carrier)",
        },
        {
            "id": "index_set_present",
            "actual": f649_summary["carrier"]["index_set"],
            "expected": "Z_12",
            "meaning": "the typed carrier index set is explicit",
        },
        {
            "id": "basis_frame_present",
            "actual": f649_summary["carrier"]["basis_frame"],
            "expected": ["c1", "s1"],
            "meaning": "the future orientation-export target frame is expressed in an explicit basis frame",
        },
        {
            "id": "nonzero_s1_support",
            "actual": bool(f649_summary["state_support"]["nonzero_s1_support"]),
            "expected": True,
            "meaning": "the designated sine axis has nonzero support in the exported state support marker",
        },
        {
            "id": "basis_function_c1_length_12",
            "actual": None if not isinstance(c1, list) else len(c1),
            "expected": 12,
            "meaning": "the carrier exports an explicit c1(x) array on Z_12",
        },
        {
            "id": "basis_function_s1_length_12",
            "actual": None if not isinstance(s1, list) else len(s1),
            "expected": 12,
            "meaning": "the carrier exports an explicit s1(x) array on Z_12",
        },
        {
            "id": "future_orientation_export_target_frame",
            "actual": f649_summary["future_orientation_export_target_frame"],
            "expected": ["c1", "s1"],
            "meaning": "a future E_orient target frame is explicit (only as a target frame)",
        },
        {
            "id": "E_orient_not_yet_exported",
            "actual": bool(f649_summary.get("E_orient_exported")),
            "expected": False,
            "meaning": "E_orient is not falsely claimed as already exported",
        },
        {
            "id": "f649_no_false_pass",
            "actual": bool(f649_summary.get("no_false_pass")),
            "expected": True,
            "meaning": "F649 summary self-reports no-false-pass discipline",
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
                "meaning": item["meaning"],
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        status = "P649_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_SECOND_CLAUSE_STATE_FOR_S_SEL_INT_AFTER_CARRIER_TYPING"
        artifact: dict[str, Any] = {
            "stage": "P649",
            "lane": "current_second_admissibility_clause_rerun_after_carrier_typing_only",
            "status": status,
            "exported_object": exported_object,
            "second_clause": second_clause,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
        summary: dict[str, Any] = {
            "stage": "P649",
            "lane": artifact["lane"],
            "status": status,
            "checks": checks,
            "blocking_mismatches": mismatches,
            "no_false_pass": True,
        }
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_SECOND_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P649"
        artifact = {
            "stage": "P649",
            "lane": "current_second_admissibility_clause_rerun_after_carrier_typing_only",
            "status": status,
            "exported_object": exported_object,
            "second_clause": second_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": [
                "source_seed_only_not_counted_as_E_orient_or_bridge",
                "strict_core_only_required (still not full admissibility)",
                "selector_acceptance_outside_strict_core_may_not_count_as_source_construction",
                "future_bridge_compatible_required",
            ],
            "no_false_pass": True,
        }
        summary = {
            "stage": "P649",
            "lane": artifact["lane"],
            "status": status,
            "exported_object": exported_object,
            "second_clause": second_clause,
            "checks": checks,
            "remaining_admissibility_clauses_unresolved": artifact["remaining_admissibility_clauses_unresolved"],
            "no_false_pass": True,
        }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

