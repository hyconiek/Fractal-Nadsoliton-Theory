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
    / "f11_legacy_weinberg_object_successor_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f11_legacy_weinberg_object_successor_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def registry_groups_have_id(payload: dict[str, Any], entry_id: str) -> bool:
    for items in payload.get("groups", {}).values():
        for item in items:
            if item.get("id") == entry_id:
                return True
    return False


def entries_have_id(payload: dict[str, Any], entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("entries", []))


def updates_have_id(payload: dict[str, Any], entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("updates", []))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p72 = load_json(
        "fundamental_action_reconstruction/generated/p72_legacy_weinberg_replaced_successor_subbranch_probe_summary.json"
    )
    release_text = load_text("RELEASE_4_9_TEXTBOOK_EN_PL.md")
    qw2068 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2068_sm_gr_parameter_registry.json"
    )
    qw2069 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
    )
    qw2098 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json"
    )
    qw2094 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2094_strict_rigor_defect_sweep.json"
    )

    object_chain_present = all(
        [
            registry_groups_have_id(qw2068, "sin2_theta_w_mz"),
            entries_have_id(qw2069, "sin2_theta_w_mz"),
            updates_have_id(qw2098, "sin2_theta_w_mz"),
            "sin2_theta_w_mz" in release_text,
        ]
    )
    object_chain_consistency_present = (
        qw2094.get("checks", {}).get("QW-2069_sin2_theta_w_mz_matches_qw2098_status")
        is True
        and qw2094.get("checks", {}).get("QW-2069_sin2_theta_w_mz_matches_qw2098_method")
        is True
    )

    textual_object_successor_verdict_present = False
    object_lineage_upgrade_verdict_present = False

    checks_spec = [
        {
            "id": "p72_object_replaced_verdict_absent",
            "actual": p72["object_replaced_verdict_present"],
            "expected": False,
            "meaning": "P72 already keeps the object-successor replaced sub-branch absent",
        },
        {
            "id": "object_chain_present",
            "actual": object_chain_present,
            "expected": True,
            "meaning": "the repo exports the sin2_theta_w_mz object chain across QW-2068/QW-2069/QW-2098/Release 4.9",
        },
        {
            "id": "object_chain_consistency_present",
            "actual": object_chain_consistency_present,
            "expected": True,
            "meaning": "QW-2094 confirms consistency of the object chain across QW-2069 and QW-2098",
        },
        {
            "id": "textual_object_successor_verdict_present",
            "actual": textual_object_successor_verdict_present,
            "expected": False,
            "meaning": "no explicit textual object-successor verdict is currently exported",
        },
        {
            "id": "object_lineage_upgrade_verdict_present",
            "actual": object_lineage_upgrade_verdict_present,
            "expected": False,
            "meaning": "no explicit object-lineage-upgrade verdict is currently exported",
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
        "stage": "F11",
        "lane": "legacy_weinberg_object_successor_refinement_current_repo_state_only",
        "goal": "refine_the_missing_object_successor_verdict_into_textual_vs_lineage_upgrade_subbranches",
        "status": "F11_EXECUTED_LEGACY_WEINBERG_OBJECT_SUCCESSOR_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P72 keeps the object-successor sub-branch absent, while the repo already exports a real sin2_theta_w_mz object chain and QW-2094 object-side consistency markers; therefore the narrowest honest refinement is textual object-successor verdict vs object-lineage-upgrade verdict",
        "candidate_state": {
            "object_chain_present": object_chain_present,
            "object_chain_consistency_present": object_chain_consistency_present,
            "textual_object_successor_verdict_present": textual_object_successor_verdict_present,
            "object_lineage_upgrade_verdict_present": object_lineage_upgrade_verdict_present,
        },
        "remaining_missing_objects": [
            "explicit_textual_object_successor_verdict_identifying_sin2_theta_w_mz_as_the_strict_side_successor_object_replacing_the_legacy_weinberg_angle_role",
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_sin2_theta_w_mz_candidate_chain_into_replacement_semantics_for_the_legacy_weinberg_angle_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F11",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "candidate_state": artifact["candidate_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
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
