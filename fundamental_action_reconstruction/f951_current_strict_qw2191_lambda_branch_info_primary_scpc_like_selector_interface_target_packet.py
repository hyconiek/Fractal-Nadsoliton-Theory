#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1012 = GENERATED / "p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe.json"
IN_F950 = GENERATED / "f950_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_packet_summary.json"

OUT = GENERATED / "f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet.json"
OUT_SUMMARY = GENERATED / "f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1012, IN_F950]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F951",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1012 = load_json(IN_P1012)
    f950 = load_json(IN_F950)
    theorem = p1012.get("theorem_result") or {}

    status = (
        "F951_EXECUTED_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
        if p1012.get("status")
        == "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDITED"
        and theorem.get("exact_selector_interface_exported_on_current_repo_state") is False
        and theorem.get("next_honest_move_is_freeze_exact_selector_interface_target") is True
        and f950.get("status")
        == "F950_EXECUTED_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
        else "F951_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F951",
        "packet_name": "CurrentStrictQW2191LambdaBranchInfoPrimarySCPCLikeSelectorInterfaceTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1012_selector_interface_nonexport_audit_probe": rel(IN_P1012),
            "f950_candidate_reference_packet_summary": rel(IN_F950),
        },
        "exported_object": {
            "object_id": p1012.get("candidate_selector_interface_target_id"),
            "goal": "Freeze one exact strict selector-interface target for the admitted info-primary SCPC-like candidate lane.",
            "reference_lane_packet_ref": f950.get("exported_object_id"),
            "candidate_reference_lane_grade": f950.get("candidate_reference_lane_grade"),
            "strict_selector_interface_status": "future_only_target_exported",
            "strict_selector_source_status": "blocked_nonexport",
            "provider_class_shift_realization_status": "not_realized",
        },
        "current_honest_reading": [
            "The candidate lane still has no exact selector interface on the current repo state.",
            "So the strongest honest positive move is now one future-only selector-interface target.",
            "This remains below selector-source export and below QW-2191 discharge.",
        ],
        "recommended_next_move": "Audit whether the new selector-interface target has any actual realization or whether one exact actual interface attempt must be frozen next.",
        "hard_limits": [
            "Does not claim strict selector interface realization.",
            "Does not claim strict selector source export.",
            "Does not claim T176 export.",
            "Does not claim QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F951",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "strict_selector_interface_status": artifact["exported_object"]["strict_selector_interface_status"],
        "provider_class_shift_realization_status": artifact["exported_object"][
            "provider_class_shift_realization_status"
        ],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
