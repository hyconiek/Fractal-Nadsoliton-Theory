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
IN_F950 = GENERATED / "f950_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_packet_summary.json"

OUT = GENERATED / "p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe_summary.json"

SELECTOR_INTERFACE_OBJECT_ID = (
    "lambda_nad12_sigma_residual_shannon_information_primary_"
    "sparse_competitive_predictive_coding_like_selector_interface_target_v1"
)
SELECTOR_INTERFACE_STEM = "information_primary_sparse_competitive_predictive_coding_like_selector_interface"

EXCLUDED_PATHS = {
    "fundamental_action_reconstruction/P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDIT_PROBE.md",
    "fundamental_action_reconstruction/N845_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDIT_THEOREM.md",
    "fundamental_action_reconstruction/F951_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_TARGET_PACKET.md",
    "fundamental_action_reconstruction/p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe.py",
    "fundamental_action_reconstruction/f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet.py",
    "fundamental_action_reconstruction/generated/p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe.json",
    "fundamental_action_reconstruction/generated/p1012_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_nonexport_audit_probe_summary.json",
    "fundamental_action_reconstruction/generated/f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet.json",
    "fundamental_action_reconstruction/generated/f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet_summary.json",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def search_repo_for_selector_interface() -> list[dict[str, str]]:
    hits: list[dict[str, str]] = []
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        repo_rel = rel(path)
        if repo_rel in EXCLUDED_PATHS:
            continue
        if path.suffix.lower() not in {".md", ".py", ".json", ".tex"}:
            continue
        text = path.read_text(encoding="utf-8", errors="ignore")
        if SELECTOR_INTERFACE_OBJECT_ID in text or SELECTOR_INTERFACE_STEM in text:
            hits.append({"path": repo_rel, "match": "candidate_specific_selector_interface_marker"})
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1011, IN_F950]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1012",
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
    f950 = load_json(IN_F950)
    theorem_result = p1011.get("theorem_result") or {}
    selector_interface_hits = search_repo_for_selector_interface()

    checks = [
        {
            "id": "p1011_already_admits_candidate_reference_lane_only",
            "pass": (
                p1011.get("status")
                == "P1011_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_SELECTOR_INTERFACE_BLOCKED"
                and theorem_result.get("info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted")
                is True
                and theorem_result.get("candidate_reference_lane_grade") == "reference_context_candidate_only"
            ),
            "details": "P1011 already admits the active Lambda branch only as a reference-context candidate lane.",
        },
        {
            "id": "f950_already_exports_candidate_reference_packet_only",
            "pass": (
                f950.get("status")
                == "F950_EXECUTED_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
                and f950.get("strict_selector_interface_status") == "blocked_nonexport"
                and f950.get("provider_class_shift_realization_status") == "not_realized"
            ),
            "details": "F950 already exports the candidate lane only as reference-context grade.",
        },
        {
            "id": "current_artifacts_still_keep_selector_interface_unexported",
            "pass": theorem_result.get("strict_selector_interface_exported") is False,
            "details": "The candidate-lane admission artifact still keeps strict selector interface explicitly unexported.",
        },
        {
            "id": "repo_has_no_exact_candidate_specific_selector_interface_marker_yet",
            "pass": len(selector_interface_hits) == 0,
            "details": "Repo scan finds no exact candidate-specific selector-interface target marker outside the files being introduced now.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem = {
        "exact_selector_interface_exported_on_current_repo_state": False if all_pass else None,
        "exact_selector_interface_target_already_exists_on_current_repo_state": False if all_pass else None,
        "next_honest_move_is_freeze_exact_selector_interface_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDITED"
        if all_pass
        else "P1012_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1012",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1011_candidate_reference_lane_admission_probe": rel(IN_P1011),
            "f950_candidate_reference_packet_summary": rel(IN_F950),
        },
        "theorem_result": theorem,
        "candidate_selector_interface_target_id": SELECTOR_INTERFACE_OBJECT_ID,
        "selector_interface_hits_outside_current_introduction": selector_interface_hits,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The candidate lane now exists only at reference-context grade.",
            "No exact strict selector interface is yet exported for that lane on the current repo state.",
            "The honest next positive move is therefore one exact selector-interface target, not a source or closure claim.",
        ],
        "recommended_next_packet": {
            "id": "F951_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_TARGET_PACKET",
            "goal": "Freeze one exact strict selector-interface target for the admitted candidate lane.",
            "export_object_id": SELECTOR_INTERFACE_OBJECT_ID,
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1012",
        "status": status,
        "as_of": AS_OF,
        "exact_selector_interface_exported_on_current_repo_state": theorem[
            "exact_selector_interface_exported_on_current_repo_state"
        ],
        "exact_selector_interface_target_already_exists_on_current_repo_state": theorem[
            "exact_selector_interface_target_already_exists_on_current_repo_state"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
