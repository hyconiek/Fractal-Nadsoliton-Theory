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
    / "p634_genuinely_new_strict_core_source_object_clause_probe_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p634_genuinely_new_strict_core_source_object_clause_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f34 = load_json(
        "fundamental_action_reconstruction/generated/f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
    )
    f318 = load_json(
        "fundamental_action_reconstruction/generated/f318_first_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_packet_summary.json"
    )
    p392 = load_json(
        "fundamental_action_reconstruction/generated/p392_current_strict_side_source_seed_strict_sigma_int_upgrade_candidate_instance_probe_summary.json"
    )
    sigma_int = load_json(
        "fundamental_action_reconstruction/generated/sigma_int_strict_derived_v1.json"
    )

    # Conservative: F318/P392 upgrade only the sigma-int provenance slot of the seed.
    # It does not, by itself, export an admissible strict-core internal selector source object S_sel_int.
    clause_satisfied = False

    checks_spec = [
        {
            "id": "p392_seed_v1_exported",
            "actual": bool(p392.get("overall_pass")),
            "expected": True,
            "meaning": "P392 confirms the strict sigma-int upgraded seed candidate instance exists",
        },
        {
            "id": "f34_first_clause_required",
            "actual": f34["minimal_source_seed_construction_contract"][
                "genuinely_new_strict_core_source_object_required"
            ],
            "expected": True,
            "meaning": "F34 keeps the genuinely-new strict-core source-object clause active",
        },
        {
            "id": "f318_candidate_seed_name",
            "actual": f318["candidate_seed_instance"]["candidate_seed_name"],
            "expected": "S_sel_int_candidate_seed_v1",
            "meaning": "F318 exports the strict sigma-int upgraded seed candidate instance name",
        },
        {
            "id": "f318_candidate_built_from_strict_sigma_int",
            "actual": f318["candidate_seed_instance"]["construction_route"][
                "internal_binary_candidate"
            ],
            "expected": "sigma_int_strict_derived_v1",
            "meaning": "F318 candidate seed uses the strict sigma-int source-upgrade datum (no hybrid FR reuse)",
        },
        {
            "id": "sigma_int_is_actual_strict_source_upgrade_value",
            "actual": str(sigma_int.get("status") or ""),
            "expected": "actual_exported_strict_source_upgrade_value__premise_based",
            "meaning": "sigma_int_strict_derived_v1 is exported as an actual strict-side source-upgrade value object (premise provenance kept explicit)",
        },
        {
            "id": "first_clause_currently_satisfied",
            "actual": clause_satisfied,
            "expected": False,
            "meaning": "the seed-v1 packet remains below an exported admissible strict-core internal selector source object (no false pass)",
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
        "stage": "P634",
        "lane": "first_source_seed_first_clause_test_only",
        "goal": "test_whether_S_sel_int_candidate_seed_v1_already_satisfies_the_genuinely_new_strict_core_source_object_clause",
        "status": "CURRENT_REPO_DOES_NOT_YET_SHOW_THAT_S_SEL_INT_CANDIDATE_SEED_V1_SATISFIES_THE_GENUINELY_NEW_STRICT_CORE_SOURCE_OBJECT_CLAUSE_AFTER_P634",
        "clause_test_result": {
            "candidate_seed_name": "S_sel_int_candidate_seed_v1",
            "target_admissible_source_name": "S_sel_int",
            "first_clause_name": "genuinely_new_strict_core_source_object_required",
            "currently_satisfied": clause_satisfied,
            "strict_sigma_int_slot_upgraded": True,
            "current_reason": (
                "seed_v1_upgrades_only_the_sigma_int_provenance_slot_(sigma_int_candidate -> sigma_int_strict_derived_v1)_"
                "but_does_not_export_a_genuinely_new_strict_core_internal_selector_source_object_identity_admissible_as_S_sel_int"
            ),
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P634",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "clause_test_result": artifact["clause_test_result"],
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

