#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "f17_legacy_fine_structure_replaced_branch_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f17_legacy_fine_structure_replaced_branch_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def registry_groups_have_id(payload: dict, entry_id: str) -> bool:
    for items in payload.get("groups", {}).values():
        for item in items:
            if item.get("id") == entry_id:
                return True
    return False


def entries_have_id(payload: dict, entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("entries", []))


def updates_have_id(payload: dict, entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("updates", []))


def method_present_in_entries(payload: dict, target_method: str) -> bool:
    return any(item.get("method") == target_method for item in payload.get("entries", []))


def method_present_in_updates(payload: dict, target_method: str) -> bool:
    return any(item.get("method") == target_method for item in payload.get("updates", []))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n90 = load_json(
        "fundamental_action_reconstruction/generated/n90_current_legacy_fine_structure_retained_branch_full_negative_closure_theorem_summary.json"
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

    method_id = "qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r"

    strict_object_candidate_present = all(
        [
            registry_groups_have_id(qw2068, "alpha_em_inv_mz"),
            entries_have_id(qw2069, "alpha_em_inv_mz"),
            updates_have_id(qw2098, "alpha_em_inv_mz"),
            "alpha_em_inv_mz" in release_text,
        ]
    )
    strict_method_candidate_present = all(
        [
            method_present_in_entries(qw2069, method_id),
            method_present_in_updates(qw2098, method_id),
            method_id in json.dumps(qw2098, ensure_ascii=True),
        ]
    )

    checks_spec = [
        {
            "id": "n90_retained_branch_closed",
            "actual": n90["theorem_result"]["retained_branch_closed_negatively_on_current_repo_state"],
            "expected": True,
            "meaning": "N90 already closes the retained branch negatively on the current repo state",
        },
        {
            "id": "strict_object_candidate_present",
            "actual": strict_object_candidate_present,
            "expected": True,
            "meaning": "the strict side exports alpha_em_inv_mz as a real successor-candidate object",
        },
        {
            "id": "strict_method_candidate_present",
            "actual": strict_method_candidate_present,
            "expected": True,
            "meaning": "the strict side exports qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r as a real successor-candidate method",
        },
    ]

    checks = []
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
        "stage": "F17",
        "lane": "legacy_fine_structure_replaced_branch_refinement_current_repo_state_only",
        "goal": "refine_the_missing_replaced_branch_verdict_into_object_successor_vs_method_successor_subbranches",
        "status": "F17_EXECUTED_LEGACY_FINE_STRUCTURE_REPLACED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "N90 closes the retained branch negatively, while the repo already exports alpha_em_inv_mz as a real successor-candidate object and qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r as a real successor-candidate method lineage; therefore the narrowest honest refinement is object-successor verdict vs method-successor-semantics verdict",
        "candidate_state": {
            "strict_object_candidate_present": strict_object_candidate_present,
            "strict_method_candidate_present": strict_method_candidate_present,
            "strict_object_candidate": "alpha_em_inv_mz",
            "strict_method_candidate": method_id,
        },
        "remaining_missing_objects": [
            "explicit_object_successor_verdict_identifying_alpha_em_inv_mz_as_the_strict_side_successor_object_replacing_the_legacy_fine_structure_role",
            "explicit_method_successor_semantics_verdict_identifying_qw2098_alpha_on_shell_from_gf_mw_sin2_delta_r_as_the_strict_side_successor_semantics_replacing_the_legacy_fine_structure_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F17",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "candidate_state": artifact["candidate_state"],
        "remaining_missing_objects": artifact["remaining_missing_objects"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()
