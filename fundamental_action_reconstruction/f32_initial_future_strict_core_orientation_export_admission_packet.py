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
    / "f32_initial_future_strict_core_orientation_export_admission_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f32_initial_future_strict_core_orientation_export_admission_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f29 = load_json(
        "fundamental_action_reconstruction/generated/f29_genuine_strict_core_internal_selector_source_admission_packet_summary.json"
    )
    f31 = load_json(
        "fundamental_action_reconstruction/generated/f31_initial_future_strict_core_selector_source_seed_target_packet_summary.json"
    )

    checks_spec = [
        {
            "id": "f29_requires_internal_orientation_discharge",
            "actual": f29["admission_contract"]["internal_orientation_discharge_required"],
            "expected": True,
            "meaning": "the future source admission contract already requires an internal orientation discharge",
        },
        {
            "id": "f31_seed_starts_with_future_source",
            "actual": f31["initial_seed_target"]["strict_core_source_object"],
            "expected": "S_sel_int",
            "meaning": "the current seed target begins with the future source object",
        },
        {
            "id": "f31_seed_second_node_is_orientation_export",
            "actual": f31["initial_seed_target"]["internal_orientation_export"],
            "expected": "E_orient",
            "meaning": "the current seed target already names the required orientation export node",
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

    orientation_export_contract = {
        "derived_from_future_source_object_required": True,
        "strict_core_only_export_required": True,
        "internal_orientation_datum_or_equivalent_required": True,
        "selector_bearing_without_external_anchor_required": True,
        "quotient_or_gauge_safe_required": True,
        "bridge_ready_for_B_sel_required": True,
        "silent_legacy_to_strict_substitution_forbidden": True,
        "axiom_augmented_selector_acceptance_may_not_count_as_export": True,
    }

    artifact = {
        "stage": "F32",
        "lane": "initial_future_strict_core_orientation_export_admission_only",
        "goal": "freeze_the_admissible_contract_for_the_orientation_export_node_of_the_future_seed",
        "status": "F32_EXECUTED_INITIAL_FUTURE_STRICT_CORE_ORIENTATION_EXPORT_ADMISSION_PACKET_NO_FALSE_PASS",
        "orientation_export_name": "E_orient",
        "upstream_source_name": "S_sel_int",
        "orientation_export_contract": orientation_export_contract,
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F32",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "orientation_export_name": artifact["orientation_export_name"],
        "upstream_source_name": artifact["upstream_source_name"],
        "orientation_export_contract": artifact["orientation_export_contract"],
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
