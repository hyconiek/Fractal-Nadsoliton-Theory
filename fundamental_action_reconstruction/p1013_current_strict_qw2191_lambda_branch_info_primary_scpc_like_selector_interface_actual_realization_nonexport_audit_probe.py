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
IN_F951 = GENERATED / "f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet_summary.json"

OUT = GENERATED / "p1013_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1013_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_nonexport_audit_probe_summary.json"

ACTUAL_REALIZATION_OBJECT_ID = (
    "lambda_nad12_sigma_residual_shannon_information_primary_"
    "sparse_competitive_predictive_coding_like_selector_interface_"
    "actual_realization_v1"
)
ACTUAL_REALIZATION_STEM = (
    "information_primary_sparse_competitive_predictive_coding_like_"
    "selector_interface_actual_realization"
)

EXCLUDED_PATHS = {
    "fundamental_action_reconstruction/P1013_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
    "fundamental_action_reconstruction/N846_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
    "fundamental_action_reconstruction/T287_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
    "fundamental_action_reconstruction/N847_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
    "fundamental_action_reconstruction/p1013_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_nonexport_audit_probe.py",
    "fundamental_action_reconstruction/p1014_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_probe.py",
    "fundamental_action_reconstruction/generated/p1013_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_nonexport_audit_probe.json",
    "fundamental_action_reconstruction/generated/p1013_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_nonexport_audit_probe_summary.json",
    "fundamental_action_reconstruction/generated/p1014_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_probe.json",
    "fundamental_action_reconstruction/generated/p1014_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_probe_summary.json",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def search_repo_for_actual_realization() -> list[dict[str, str]]:
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
        if ACTUAL_REALIZATION_OBJECT_ID in text or ACTUAL_REALIZATION_STEM in text:
            hits.append({"path": repo_rel, "match": "candidate_specific_selector_interface_actual_realization_marker"})
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1012, IN_F951]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1013",
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
    f951 = load_json(IN_F951)
    theorem = p1012.get("theorem_result") or {}
    actual_realization_hits = search_repo_for_actual_realization()

    checks = [
        {
            "id": "p1012_keeps_selector_interface_export_false",
            "pass": (
                p1012.get("status")
                == "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDITED"
                and theorem.get("exact_selector_interface_exported_on_current_repo_state") is False
                and theorem.get("exact_selector_interface_target_already_exists_on_current_repo_state") is False
            ),
            "details": "P1012 already keeps exact selector-interface export false on the current repo state.",
        },
        {
            "id": "f951_exports_future_only_selector_interface_target",
            "pass": (
                f951.get("status")
                == "F951_EXECUTED_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and f951.get("strict_selector_interface_status") == "future_only_target_exported"
                and f951.get("provider_class_shift_realization_status") == "not_realized"
            ),
            "details": "F951 already exports one future-only exact selector-interface target.",
        },
        {
            "id": "repo_has_no_exact_candidate_specific_selector_interface_actual_realization_marker_yet",
            "pass": len(actual_realization_hits) == 0,
            "details": "Repo scan finds no exact candidate-specific actual-realization marker for that interface outside the files being introduced now.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "current_repo_has_exported_actual_realization_of_selector_interface_target": False if all_pass else None,
        "selector_interface_target_still_remains_future_only_not_actual_export": True if all_pass else None,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_selector_interface_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1013_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if all_pass
        else "P1013_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1013",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1012_selector_interface_nonexport_audit_probe": rel(IN_P1012),
            "f951_selector_interface_target_packet_summary": rel(IN_F951),
        },
        "theorem_result": theorem_result,
        "candidate_selector_interface_actual_realization_id": ACTUAL_REALIZATION_OBJECT_ID,
        "actual_realization_hits_outside_current_introduction": actual_realization_hits,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The frozen selector-interface target still remains future-only on the current repo state.",
            "No exact actual realization of that target is exported yet.",
            "The honest next positive move is therefore one exact actual-realization attempt below the same target.",
        ],
        "recommended_next_probe": {
            "id": "P1014_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_PROBE",
            "goal": "Export one exact actual-realization attempt below the same selector-interface target.",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1013",
        "status": status,
        "as_of": AS_OF,
        "current_repo_has_exported_actual_realization_of_selector_interface_target": theorem_result[
            "current_repo_has_exported_actual_realization_of_selector_interface_target"
        ],
        "selector_interface_target_still_remains_future_only_not_actual_export": theorem_result[
            "selector_interface_target_still_remains_future_only_not_actual_export"
        ],
        "recommended_next_probe": artifact["recommended_next_probe"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
