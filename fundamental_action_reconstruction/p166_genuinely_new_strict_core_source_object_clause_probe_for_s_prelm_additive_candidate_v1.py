#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p166_genuinely_new_strict_core_source_object_clause_probe_for_s_prelm_additive_candidate_v1_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f76 = load_json(
        "fundamental_action_reconstruction/generated/f76_first_additive_preobserver_source_object_construction_attempt_packet_summary.json"
    )
    n182 = load_json(
        "fundamental_action_reconstruction/generated/n182_current_first_additive_preobserver_source_object_construction_attempt_theorem_summary.json"
    )
    n183 = load_json(
        "fundamental_action_reconstruction/generated/n183_current_first_additive_preobserver_source_object_admissibility_upgrade_target_theorem_summary.json"
    )
    n184 = load_json(
        "fundamental_action_reconstruction/generated/n184_current_first_additive_preobserver_source_object_first_clause_admissibility_target_theorem_summary.json"
    )

    checks = [
        {
            "id": "future_only_attempt",
            "actual": f76["guardrails"]["future_only"],
            "expected": False,
            "pass": f76["guardrails"]["future_only"] is False,
            "meaning": "A genuinely new strict-core source object may not remain only future-only attempt packaging",
        },
        {
            "id": "constructed_source_object_exported",
            "actual": f76["constructed_source_object"],
            "expected": True,
            "pass": f76["constructed_source_object"] is True,
            "meaning": "The object would need to be exported as a constructed source object",
        },
        {
            "id": "n182_not_nonpromoted_only",
            "actual": n182["theorem_result"]["full_closure_pass"],
            "expected": True,
            "pass": n182["theorem_result"]["full_closure_pass"] is True,
            "meaning": "Current theorem state would need promotion beyond construction-attempt scope",
        },
        {
            "id": "n183_not_upgrade_target_only",
            "actual": n183["theorem_result"]["admissibility_upgrade_target"],
            "expected": "constructed_source_object_export",
            "pass": n183["theorem_result"]["admissibility_upgrade_target"] == "constructed_source_object_export",
            "meaning": "Current state would need a realized object, not only an upgrade target",
        },
        {
            "id": "n184_first_clause_fixed",
            "actual": n184["theorem_result"]["first_clause"],
            "expected": "genuinely_new_strict_core_source_object_required",
            "pass": n184["theorem_result"]["first_clause"] == "genuinely_new_strict_core_source_object_required",
            "meaning": "The first clause is indeed the clause under test",
        },
    ]

    blocking = [item["id"] for item in checks if not item["pass"] and item["id"] != "n184_first_clause_fixed"]

    if blocking:
        summary = {
            "stage": "P166",
            "lane": "first_additive_preobserver_source_object_first_clause_probe_only",
            "status": "CURRENT_REPO_DOES_NOT_YET_SHOW_THAT_S_PRELM_ADDITIVE_CANDIDATE_V1_SATISFIES_THE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_CLAUSE_AFTER_P166",
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "checks": checks,
            "blocking_obstructions": blocking,
            "no_false_pass": True,
        }
    else:
        summary = {
            "stage": "P166",
            "lane": "first_additive_preobserver_source_object_first_clause_probe_only",
            "status": "P166_REQUIRES_REVIEW_UNEXPECTED_POSITIVE_FIRST_CLAUSE_STATE",
            "construction_attempt": "S_preLM_additive_candidate_v1",
            "first_clause": "genuinely_new_strict_core_source_object_required",
            "checks": checks,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
