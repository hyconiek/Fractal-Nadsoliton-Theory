#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P653_SUMMARY = (
    GENERATED
    / "p653_current_sixth_admissibility_clause_rerun_after_future_bridge_compatibility_packet_for_s_sel_int_strict_core_source_object_v1_probe_summary.json"
)

OUT_JSON = (
    GENERATED
    / "p676_current_first_admissible_s_sel_int_source_object_discharge_probe.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p676_current_first_admissible_s_sel_int_source_object_discharge_probe_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not P653_SUMMARY.exists():
        status = "P676_NOT_COMPUTABLE_MISSING_P653_SUMMARY"
        artifact = {
            "stage": "P676",
            "lane": "current_first_admissible_s_sel_int_source_object_discharge_probe_only",
            "status": status,
            "missing": str(P653_SUMMARY.relative_to(REPO)),
            "no_false_pass": True,
        }
        summary = {
            "stage": "P676",
            "lane": artifact["lane"],
            "status": status,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p653 = load_json(P653_SUMMARY)
    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_STRICT_CORE_SOURCE_OBJECT_SATISFYING_THE_SIXTH_ADMISSIBILITY_CLAUSE_FOR_S_SEL_INT_AFTER_P653"
    )

    status_ok = p653.get("status") == expected_status
    remaining = p653.get("remaining_admissibility_clauses_unresolved")
    remaining_ok = isinstance(remaining, list) and len(remaining) == 0

    checks = [
        {
            "id": "sixth_clause_rerun_positive",
            "actual": p653.get("status"),
            "expected": expected_status,
            "pass": status_ok,
        },
        {
            "id": "no_remaining_admissibility_clauses_unresolved",
            "actual": remaining,
            "expected": [],
            "pass": remaining_ok,
        },
    ]

    ok = status_ok and remaining_ok

    status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_STRICT_CORE_SOURCE_OBJECT_FOR_S_SEL_INT_IN_THE_SENSE_OF_F34_AFTER_P676"
        if ok
        else "P676_REQUIRES_REVIEW_SOURCE_OBJECT_ADMISSIBILITY_CONTRACT_NOT_FULLY_SATISFIED_ON_CURRENT_REPO_STATE"
    )

    exported_object = p653.get("exported_object")
    artifact = {
        "stage": "P676",
        "lane": "current_first_admissible_s_sel_int_source_object_discharge_probe_only",
        "status": status,
        "inputs": {"p653_summary": str(P653_SUMMARY.relative_to(REPO))},
        "checks": checks,
        "result": {
            "admissible_S_sel_int_source_object_in_F34_sense": ok,
            "exported_object": exported_object,
            "counts_as": (
                "admissible_source_object_for_S_sel_int (F34 contract satisfied)"
                if ok
                else "not_yet_verified_as_admissible_source_object_for_S_sel_int"
            ),
        },
        "hard_limits": [
            "no_strict_core_selector_closure",
            "no_global_kernel_alone_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P676",
        "lane": artifact["lane"],
        "status": status,
        "admissible_S_sel_int_source_object_in_F34_sense": ok,
        "exported_object": exported_object,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

