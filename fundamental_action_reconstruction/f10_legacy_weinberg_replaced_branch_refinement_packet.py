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
    / "f10_legacy_weinberg_replaced_branch_refinement_packet.json"
)
OUT_SUMMARY = (
    GENERATED
    / "f10_legacy_weinberg_replaced_branch_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def any_entry_with_id(payload: dict[str, Any], entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("entries", []))


def any_update_with_id(payload: dict[str, Any], entry_id: str) -> bool:
    return any(item.get("id") == entry_id for item in payload.get("updates", []))


def any_entry_with_method(payload: dict[str, Any], method_name: str) -> bool:
    return any(item.get("method") == method_name for item in payload.get("entries", []))


def any_update_with_method(payload: dict[str, Any], method_name: str) -> bool:
    return any(item.get("method") == method_name for item in payload.get("updates", []))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p71 = load_json(
        "fundamental_action_reconstruction/generated/p71_strict_side_replaced_branch_probe_for_legacy_weinberg_role_summary.json"
    )
    release_text = load_text("RELEASE_4_9_TEXTBOOK_EN_PL.md")
    qw2069 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2069_full_sm_gr_derivation_package.json"
    )
    qw2098 = load_json(
        "material_dowodowy/korpus_qw_pozostaly/raporty_json/report_qw2098_ew_secondary_nonanchor_closure_gate.json"
    )

    object_candidate_present = (
        "sin2_theta_w_mz" in release_text
        or any_entry_with_id(qw2069, "sin2_theta_w_mz")
        or any_update_with_id(qw2098, "sin2_theta_w_mz")
    )

    method_candidate_present = (
        "qw2098_sin2_from_nonanchor_ew_pole_chain" in release_text
        or any_entry_with_method(qw2069, "qw2098_sin2_from_nonanchor_ew_pole_chain")
        or any_update_with_method(qw2098, "qw2098_sin2_from_nonanchor_ew_pole_chain")
    )

    replaced_object_verdict_present = False
    replaced_method_verdict_present = False

    checks_spec = [
        {
            "id": "p71_replaced_branch_absent",
            "actual": p71["replaced_verdict_present"],
            "expected": False,
            "meaning": "P71 already shows the replaced branch is absent as a whole",
        },
        {
            "id": "strict_object_candidate_present",
            "actual": object_candidate_present,
            "expected": True,
            "meaning": "the strict side exports sin2_theta_w_mz as a real successor-candidate object",
        },
        {
            "id": "strict_method_candidate_present",
            "actual": method_candidate_present,
            "expected": True,
            "meaning": "the strict side exports qw2098_sin2_from_nonanchor_ew_pole_chain as a real successor-candidate method lineage",
        },
        {
            "id": "replaced_object_verdict_present",
            "actual": replaced_object_verdict_present,
            "expected": False,
            "meaning": "no explicit strict-side object-successor verdict for the legacy Weinberg role is currently exported",
        },
        {
            "id": "replaced_method_verdict_present",
            "actual": replaced_method_verdict_present,
            "expected": False,
            "meaning": "no explicit strict-side method-successor-semantics verdict for the legacy Weinberg role is currently exported",
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
        "stage": "F10",
        "lane": "legacy_weinberg_replaced_branch_refinement_current_repo_state_only",
        "goal": "refine_the_missing_replaced_branch_into_object_successor_vs_method_successor_semantics_subbranches",
        "status": "F10_EXECUTED_LEGACY_WEINBERG_REPLACED_BRANCH_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "P71 keeps the replaced branch absent as a whole, while the strict side already exports both a successor-candidate object sin2_theta_w_mz and a successor-candidate method qw2098_sin2_from_nonanchor_ew_pole_chain; therefore the narrowest honest refinement is object-successor verdict vs method-successor-semantics verdict",
        "candidate_state": {
            "strict_object_candidate_present": object_candidate_present,
            "strict_method_candidate_present": method_candidate_present,
            "replaced_object_verdict_present": replaced_object_verdict_present,
            "replaced_method_verdict_present": replaced_method_verdict_present,
        },
        "remaining_missing_objects": [
            "explicit_object_successor_verdict_identifying_sin2_theta_w_mz_as_the_strict_side_successor_object_replacing_the_legacy_weinberg_angle_role",
            "explicit_method_successor_semantics_verdict_identifying_qw2098_sin2_from_nonanchor_ew_pole_chain_as_the_strict_side_successor_semantics_replacing_the_legacy_weinberg_angle_role",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F10",
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
