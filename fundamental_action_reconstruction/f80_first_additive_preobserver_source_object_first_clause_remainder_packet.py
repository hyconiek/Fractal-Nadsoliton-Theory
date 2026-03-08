#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f80_first_additive_preobserver_source_object_first_clause_remainder_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p166 = load_json(
        "fundamental_action_reconstruction/generated/p166_genuinely_new_strict_core_source_object_clause_probe_for_s_prelm_additive_candidate_v1_summary.json"
    )
    n186 = load_json(
        "fundamental_action_reconstruction/generated/n186_current_additive_preobserver_source_object_nonreduction_theorem_summary.json"
    )

    blocking = [
        item
        for item in p166["blocking_obstructions"]
        if item
        in {
            "future_only_attempt",
            "constructed_source_object_exported",
            "n182_not_nonpromoted_only",
            "n183_not_upgrade_target_only",
        }
    ]

    summary = {
        "stage": "F80",
        "lane": "first_additive_preobserver_source_object_first_clause_remainder_only",
        "goal": "remove_the_simplest_packaging_reduction_reading_and_freeze_the_remaining_first_clause_obstruction_package",
        "status": "F80_EXECUTED_FIRST_ADDITIVE_PREOBSERVER_SOURCE_OBJECT_FIRST_CLAUSE_REMAINDER_PACKET_NO_FALSE_PASS",
        "construction_attempt": "S_preLM_additive_candidate_v1",
        "first_clause": "genuinely_new_strict_core_source_object_required",
        "nonreduction_witness_present": n186["theorem_result"]["delta_norm"] > 0.0,
        "simple_packaging_reduction_removed": n186["theorem_result"]["delta_norm"] > 0.0,
        "remaining_first_clause_obstructions": blocking,
        "remaining_first_clause_package": "realized_constructed_source_object_export_package",
        "package_meaning": [
            "future_only_attempt_scope_still_active",
            "constructed_source_object_export_still_absent",
            "construction_attempt_scope_not_yet_promoted",
            "admissibility_upgrade_target_scope_not_yet_realized",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
