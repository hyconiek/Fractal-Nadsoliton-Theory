#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f38_first_source_seed_first_clause_refinement_packet.json"
OUT_SUMMARY = GENERATED / "f38_first_source_seed_first_clause_refinement_packet_summary.json"


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    p123 = load_json(
        "fundamental_action_reconstruction/generated/p123_first_source_seed_admissibility_upgrade_target_probe_summary.json"
    )

    checks_spec = [
        {
            "id": "p123_target_fixed",
            "actual": p123["target_state"]["next_constructive_move_reduced_to_one_first_admissibility_upgrade_target"],
            "expected": True,
            "meaning": "P123 already fixes one and only one first admissibility-upgrade target",
        },
        {
            "id": "f34_first_clause_present",
            "actual": f34["minimal_source_seed_construction_contract"]["genuinely_new_strict_core_source_object_required"],
            "expected": True,
            "meaning": "F34 explicitly keeps the genuinely-new strict-core object clause active",
        },
        {
            "id": "p123_target_uses_sigma_candidate",
            "actual": p123["target_state"]["admissibility_upgrade_target"]["candidate_seed_instance"]["construction_route"]["internal_binary_candidate"],
            "expected": "sigma_int_candidate",
            "meaning": "the current target still reuses a current internal binary candidate",
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

    first_clause_target = {
        "candidate_seed_name": p123["target_state"]["admissibility_upgrade_target"]["candidate_seed_name"],
        "target_admissible_source_name": p123["target_state"]["admissibility_upgrade_target"]["target_admissible_source_name"],
        "first_clause_name": "genuinely_new_strict_core_source_object_required",
        "first_clause_expected_value": True,
        "why_first": [
            "candidate_seed_still_reuses_current_route_material",
            "minimal_contract_requires_genuinely_new_source_object",
            "no_earlier_clause_ordering_ambiguity_remains",
        ],
    }

    artifact = {
        "stage": "F38",
        "lane": "first_source_seed_first_clause_only",
        "goal": "freeze_the_first_clause_level_admissibility_target_for_S_sel_int_candidate_seed_v0",
        "status": "F38_EXECUTED_FIRST_SOURCE_SEED_FIRST_CLAUSE_REFINEMENT_PACKET_NO_FALSE_PASS",
        "first_clause_target": first_clause_target,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F38",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "first_clause_target": artifact["first_clause_target"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
