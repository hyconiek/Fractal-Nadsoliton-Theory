#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = GENERATED / "f9_legacy_weinberg_semantic_transfer_refinement_packet.json"
OUT_SUMMARY = (
    GENERATED / "f9_legacy_weinberg_semantic_transfer_refinement_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def load_text(repo_relative_path: str) -> str:
    return (REPO / repo_relative_path).read_text(encoding="utf-8")


def any_entry_with_id(obj: Any, target_id: str) -> bool:
    if isinstance(obj, dict):
        if obj.get("id") == target_id:
            return True
        return any(any_entry_with_id(value, target_id) for value in obj.values())
    if isinstance(obj, list):
        return any(any_entry_with_id(value, target_id) for value in obj)
    return False


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    p67 = load_json(
        "fundamental_action_reconstruction/generated/p67_strict_side_role_equivalence_probe_for_legacy_weinberg_role_summary.json"
    )
    qw2093_output = load_json(
        "t1_nonanchor_observables_input_qw2085_2086.json"
    )

    candidate_present = p67["strict_side_candidate_object"]["present"]
    qw2093_text = json.dumps(qw2093_output, ensure_ascii=True)
    lineage_touchpoint_present = (
        "sin2_theta_w_eff" in qw2093_text
        and "alpha_geo" in qw2093_text
        and "micro-radiative shift relation" in qw2093_text
    )

    checks_spec = [
        {
            "id": "p67_candidate_present",
            "actual": candidate_present,
            "expected": True,
            "meaning": "P67 already confirms presence of the strict-side Weinberg candidate object",
        },
        {
            "id": "qw2093_alpha_geo_lineage_touchpoint_present",
            "actual": lineage_touchpoint_present,
            "expected": True,
            "meaning": "the QW-2093 output artifact exports a real alpha_geo-bearing lineage touchpoint for the strict-side Weinberg chain",
        },
        {
            "id": "p67_explicit_transfer_verdict_absent",
            "actual": p67["explicit_role_equivalence_verdict_present"],
            "expected": False,
            "meaning": "P67 already confirms that no explicit retained semantic-transfer verdict is currently exported",
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
        "stage": "F9",
        "lane": "legacy_weinberg_semantic_transfer_refinement_current_repo_state_only",
        "goal": "refine_the_remaining_weinberg_semantic_transfer_blocker_into_textual_successor_versus_lineage_upgrade_subbranches",
        "status": "F9_EXECUTED_LEGACY_WEINBERG_SEMANTIC_TRANSFER_REFINEMENT_PACKET_NO_FALSE_PASS",
        "reason": "the strict side already exports sin2_theta_w_mz and QW-2093 already exports an alpha_geo-bearing lineage touchpoint, but neither fact yet upgrades to an explicit retained semantic-transfer verdict",
        "semantic_transfer_state": {
            "strict_side_candidate_object_present": candidate_present,
            "qw2093_alpha_geo_lineage_touchpoint_present": lineage_touchpoint_present,
            "explicit_transfer_verdict_present": p67["explicit_role_equivalence_verdict_present"],
        },
        "remaining_missing_objects": [
            "explicit_textual_retained_successor_verdict_identifying_sin2_theta_w_mz_as_the_retained_strict_side_successor_of_the_legacy_weinberg_angle_role",
            "explicit_lineage_upgrade_verdict_elevating_the_qw2093_alpha_geo_touchpoint_into_retained_strict_side_weinberg_role_transfer",
        ],
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "F9",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "semantic_transfer_state": artifact["semantic_transfer_state"],
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
