#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "f50_genuinely_additive_source_object_construction_contract_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f50_genuinely_additive_source_object_construction_contract_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n126 = load_json(
        "fundamental_action_reconstruction/generated/n126_current_repo_exports_no_admissible_strict_core_internal_selector_source_object_theorem_summary.json"
    )
    n136 = load_json(
        "fundamental_action_reconstruction/generated/n136_current_first_clause_obstruction_theorem_for_s_sel_int_candidate_seed_v0_summary.json"
    )
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n126_no_admissible_object_present",
            "actual": n126["theorem_result"][
                "admissible_strict_core_internal_selector_source_object_present"
            ],
            "expected": False,
            "meaning": "N126 already excludes currently exported objects from admissible source-object status",
        },
        {
            "id": "n136_first_clause_obstruction_discharged",
            "actual": n136["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N136 already excludes packaged reuse from counting as a genuinely new source object",
        },
        {
            "id": "n149_constructive_frontier_exhausted",
            "actual": n149["theorem_result"][
                "constructive_selector_frontier_exhausted_on_current_repo_state"
            ],
            "expected": True,
            "meaning": "N149 already closes the constructive selector frontier inside present exports",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "F50",
        "lane": "future_genuinely_additive_source_object_construction_only",
        "goal": "freeze_the_minimal_contract_for_any_future_positive_move_that_would_count_as_a_genuinely_additive_new_strict_core_source_object_construction",
        "status": "F50_EXECUTED_GENUINELY_ADDITIVE_SOURCE_OBJECT_CONSTRUCTION_CONTRACT_PACKET_NO_FALSE_PASS",
        "minimal_additive_construction_contract": {
            "new_exported_object_identity": True,
            "not_just_packaged_reuse": True,
            "strict_core_only": True,
            "no_external_selector_import": True,
            "kernel_split_safe": True,
            "future_admissibility_compatible": True,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F50",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "minimal_additive_construction_contract": artifact[
            "minimal_additive_construction_contract"
        ],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
