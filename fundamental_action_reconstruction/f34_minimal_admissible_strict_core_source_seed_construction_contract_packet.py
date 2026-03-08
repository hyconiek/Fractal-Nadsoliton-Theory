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
    / "f34_minimal_admissible_strict_core_source_seed_construction_contract_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f34_minimal_admissible_strict_core_source_seed_construction_contract_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f33 = load_json(
        "fundamental_action_reconstruction/generated/f33_future_strict_core_source_seed_construction_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f33_source_seed_name",
            "actual": f33["source_seed_construction_target"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "F33 already freezes S_sel_int as the first source-seed target",
        },
        {
            "id": "f33_requires_strict_core_only",
            "actual": f33["source_seed_construction_target"]["strict_core_only_required"],
            "expected": True,
            "meaning": "the first source target already requires strict-core scope",
        },
        {
            "id": "f33_requires_later_orientation_export",
            "actual": f33["source_seed_construction_target"][
                "source_carrying_enough_for_later_E_orient_export_required"
            ],
            "expected": True,
            "meaning": "the first source target already carries a later E_orient export requirement",
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

    minimal_source_seed_construction_contract = {
        "strict_core_source_object": "S_sel_int",
        "genuinely_new_strict_core_source_object_required": True,
        "carrier_typed_enough_for_later_E_orient_export_required": True,
        "source_seed_only_not_counted_as_E_orient_or_bridge": True,
        "strict_core_only_required": True,
        "silent_legacy_to_strict_substitution_forbidden": True,
        "selector_acceptance_outside_strict_core_may_not_count_as_source_construction": True,
        "future_bridge_compatible_required": True,
    }

    artifact = {
        "stage": "F34",
        "lane": "minimal_admissible_strict_core_source_seed_construction_contract_only",
        "goal": "freeze_the_minimal_admissible_construction_contract_for_S_sel_int",
        "status": "F34_EXECUTED_MINIMAL_ADMISSIBLE_STRICT_CORE_SOURCE_SEED_CONSTRUCTION_CONTRACT_PACKET_NO_FALSE_PASS",
        "minimal_source_seed_construction_contract": minimal_source_seed_construction_contract,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F34",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "minimal_source_seed_construction_contract": artifact["minimal_source_seed_construction_contract"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n",
        encoding="ascii",
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
