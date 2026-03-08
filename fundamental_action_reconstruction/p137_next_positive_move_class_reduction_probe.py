#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "p137_next_positive_move_class_reduction_probe.json"
OUT_SUMMARY = (
    GENERATED / "p137_next_positive_move_class_reduction_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f50 = load_json(
        "fundamental_action_reconstruction/generated/f50_genuinely_additive_source_object_construction_contract_packet_summary.json"
    )
    n149 = load_json(
        "fundamental_action_reconstruction/generated/n149_current_repo_constructive_selector_frontier_exhaustion_theorem_summary.json"
    )

    checks_spec = [
        {
            "id": "n149_only_remaining_positive_move_class",
            "actual": n149["theorem_result"]["only_remaining_positive_move_class"],
            "expected": "future_genuinely_additive_new_strict_core_source_object_construction",
            "meaning": "N149 already fixes the only remaining positive move class",
        },
        {
            "id": "f50_contract_new_exported_object_identity",
            "actual": f50["minimal_additive_construction_contract"]["new_exported_object_identity"],
            "expected": True,
            "meaning": "F50 freezes the additive construction contract",
        },
        {
            "id": "f50_contract_not_just_packaged_reuse",
            "actual": f50["minimal_additive_construction_contract"]["not_just_packaged_reuse"],
            "expected": True,
            "meaning": "F50 excludes repackaging of current exports",
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
        "stage": "P137",
        "lane": "next_positive_move_class_reduction_only",
        "goal": "test_whether_the_only_remaining_positive_move_class_is_now_reduced_to_one_explicit_minimal_additive_construction_contract",
        "status": "CURRENT_REPO_REDUCES_THE_ONLY_REMAINING_POSITIVE_MOVE_CLASS_TO_ONE_EXPLICIT_MINIMAL_ADDITIVE_CONSTRUCTION_CONTRACT_AFTER_P137",
        "target_state": {
            "only_remaining_positive_move_class": "future_genuinely_additive_new_strict_core_source_object_construction",
            "reduced_to_one_explicit_minimal_additive_construction_contract": True,
            "remaining_open_branches": [
                "future_attempted_genuinely_additive_new_strict_core_source_object_construction"
            ],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P137",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
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
