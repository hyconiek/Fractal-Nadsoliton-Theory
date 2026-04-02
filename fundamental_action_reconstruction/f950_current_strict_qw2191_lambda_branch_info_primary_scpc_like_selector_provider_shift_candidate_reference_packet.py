#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1011 = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe.json"
IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"

OUT = GENERATED / "f950_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_packet.json"
OUT_SUMMARY = GENERATED / "f950_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1011, IN_F949]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F950",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1011 = load_json(IN_P1011)
    f949 = load_json(IN_F949)

    theorem_result = p1011.get("theorem_result") or {}
    lane = p1011.get("candidate_reference_lane") or {}

    status = (
        "F950_EXECUTED_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
        if p1011.get("status")
        == "P1011_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_SELECTOR_INTERFACE_BLOCKED"
        and theorem_result.get("info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted")
        is True
        and theorem_result.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
        and theorem_result.get("strict_selector_interface_exported") is False
        and theorem_result.get("strict_selector_source_exported") is False
        else "F950_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "F950",
        "packet_name": "CurrentStrictQW2191LambdaBranchInfoPrimarySCPCLikeSelectorProviderShiftCandidateReference_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1011_candidate_reference_lane_admission_probe": rel(IN_P1011),
            "f949_noncyclic_pivot_packet_summary": rel(IN_F949),
        },
        "exported_object": {
            "object_id": lane.get("candidate_id"),
            "goal": "Export the active Lambda branch only as an info-primary SCPC-like selector-provider shift candidate reference lane for QW-2191.",
            "pivot_family_ref": f949.get("preferred_noncyclic_pivot_family"),
            "pivot_branch_ref": f949.get("preferred_first_pivot_branch"),
            "candidate_reference_lane_grade": theorem_result.get("candidate_reference_lane_grade"),
            "reading_contract": lane.get("reading_contract"),
            "information_primary_ontology_constraint_active": theorem_result.get(
                "information_primary_ontology_constraint_active"
            ),
            "strict_selector_interface_status": "blocked_nonexport",
            "strict_selector_source_status": "blocked_nonexport",
            "provider_class_shift_realization_status": "not_realized",
            "identity_contract": {
                "predictive_coding_like_reference_context_only": True,
                "neural_network_identity_exported": False,
                "energy_based_identity_exported": False,
            },
        },
        "current_honest_reading": [
            "The active Lambda branch is now exportable as a narrower selector-provider shift candidate reference lane.",
            "That reading is explicitly information-primary and predictive-coding-like only by reference-context grade.",
            "The strict selector interface remains blocked, so no selector-source or QW-2191 closure claim is admitted.",
        ],
        "recommended_next_move": "Audit whether the candidate lane already exports any exact strict selector interface before freezing one new exact interface target.",
        "hard_limits": [
            "Does not claim strict selector interface export.",
            "Does not claim strict selector source export.",
            "Does not claim theorem-level neural-network identity.",
            "Does not claim energy-based selector identity.",
            "Does not claim T176 export.",
            "Does not claim QW-2191 discharge.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F950",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "candidate_reference_lane_grade": artifact["exported_object"]["candidate_reference_lane_grade"],
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
