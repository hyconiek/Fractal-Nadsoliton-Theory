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
    / "f16_legacy_fine_structure_semantic_transfer_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f16_legacy_fine_structure_semantic_transfer_refinement_packet_summary.json"
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

    p82 = load_json(
        "fundamental_action_reconstruction/generated/p82_strict_side_role_equivalence_probe_for_legacy_fine_structure_role_summary.json"
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
            registry_groups_have_id(qw2068, "alpha_em_inv_mz"),
            entries_have_id(qw2069, "alpha_em_inv_mz"),
            updates_have_id(qw2098, "alpha_em_inv_mz"),
            "alpha_em_inv_mz" in release_text,
        ]
    )
    object_chain_consistency_present = (
        qw2094.get("checks", {}).get("QW-2069_alpha_em_inv_mz_matches_qw2098_status")
        is True
        and qw2094.get("checks", {}).get("QW-2069_alpha_em_inv_mz_matches_qw2098_method")
        is True
    )

    checks_spec = [
        {
            "id": "p82_candidate_present",
            "actual": p82["strict_side_candidate_object"]["present"],
            "expected": True,
            "meaning": "P82 already confirms that a real strict-side fine-structure candidate object is exported",
        },
        {
            "id": "object_chain_present",
            "actual": object_chain_present,
            "expected": True,
            "meaning": "the repo exports the alpha_em_inv_mz candidate chain across QW-2068/QW-2069/QW-2098/Release 4.9",
        },
        {
            "id": "object_chain_consistency_present",
            "actual": object_chain_consistency_present,
            "expected": True,
            "meaning": "QW-2094 confirms consistency of the alpha_em_inv_mz candidate chain across QW-2069 and QW-2098",
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
        "stage": "F16",
        "lane": "legacy_fine_structure_semantic_transfer_refinement_current_repo_state_only",
        "goal": "refine_the_missing_retained_semantic_transfer_verdict_into_textual_vs_object_lineage_upgrade_subbranches",
        "status": "F16_EXECUTED_LEGACY_FINE_STRUCTURE_SEMANTIC_TRANSFER_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P82 keeps the retained semantic-transfer verdict absent, while the repo already exports a real alpha_em_inv_mz candidate object, a real strict-side candidate chain, and QW-2094 consistency markers; therefore the narrowest honest refinement is textual retained-successor verdict vs object-lineage-upgrade verdict",
        "candidate_state": {
            "strict_side_candidate_object_present": p82["strict_side_candidate_object"]["present"],
            "strict_side_candidate_object": "alpha_em_inv_mz",
            "object_chain_present": object_chain_present,
            "object_chain_consistency_present": object_chain_consistency_present,
        },
        "remaining_missing_objects": [
            "explicit_textual_retained_successor_verdict_identifying_alpha_em_inv_mz_as_the_retained_strict_side_successor_of_the_legacy_fine_structure_role",
            "explicit_object_lineage_upgrade_verdict_elevating_the_existing_alpha_em_inv_mz_candidate_chain_into_retained_strict_side_fine_structure_role_transfer",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F16",
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
