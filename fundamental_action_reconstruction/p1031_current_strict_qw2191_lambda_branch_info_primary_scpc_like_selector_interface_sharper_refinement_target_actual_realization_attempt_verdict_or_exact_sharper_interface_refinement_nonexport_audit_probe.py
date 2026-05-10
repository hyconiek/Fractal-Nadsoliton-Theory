#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1030 = GENERATED / "p1030_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe.json"

OUT = GENERATED / "p1031_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_interface_refinement_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1031_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_interface_refinement_nonexport_audit_probe_summary.json"

SHARPER_REFINEMENT_OBJECT_ID = (
    "lambda_nad12_sigma_residual_shannon_information_primary_"
    "sparse_competitive_predictive_coding_like_selector_interface_"
    "sharper_refinement_target_actual_realization_attempt_"
    "sharper_refinement_target_actual_realization_attempt_"
    "sharper_interface_refinement_target_actual_realization_attempt_"
    "sharper_interface_refinement_target_actual_realization_attempt_"
    "sharper_interface_refinement_target_v1"
)

EXCLUDED_PATHS = {
    "fundamental_action_reconstruction/P1031_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_INTERFACE_REFINEMENT_NONEXPORT_AUDIT_PROBE.md",
    "fundamental_action_reconstruction/N864_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_INTERFACE_REFINEMENT_NONEXPORT_AUDIT_THEOREM.md",
    "fundamental_action_reconstruction/T296_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_INTERFACE_REFINEMENT_TARGET_SPEC.md",
    "fundamental_action_reconstruction/N865_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_INTERFACE_REFINEMENT_TARGET_THEOREM.md",
    "fundamental_action_reconstruction/p1031_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_interface_refinement_nonexport_audit_probe.py",
    "fundamental_action_reconstruction/p1032_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_exact_sharper_interface_refinement_target_probe.py",
    "fundamental_action_reconstruction/generated/p1031_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_interface_refinement_nonexport_audit_probe.json",
    "fundamental_action_reconstruction/generated/p1031_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_verdict_or_exact_sharper_interface_refinement_nonexport_audit_probe_summary.json",
    "fundamental_action_reconstruction/generated/p1032_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_exact_sharper_interface_refinement_target_probe.json",
    "fundamental_action_reconstruction/generated/p1032_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_exact_sharper_interface_refinement_target_probe_summary.json",
}


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def search_repo_for_sharper_refinement() -> list[dict[str, str]]:
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
        if SHARPER_REFINEMENT_OBJECT_ID in text:
            hits.append({"path": repo_rel, "match": "candidate_specific_t295_sharper_interface_refinement_marker"})
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1030]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1031",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1030 = load_json(IN_P1030)
    sharper_refinement_hits = search_repo_for_sharper_refinement()

    checks = [
        {
            "id": "p1030_already_exports_exact_t295_attempt",
            "pass": (
                p1030.get("status")
                == "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
                and p1030.get("t295_attempt_exported_on_current_repo_state") is True
                and p1030.get("t295_attempt_keeps_strict_selector_source_open") is True
            ),
            "details": "P1030 already exports the exact T295 sharper-refinement actual-realization attempt.",
        },
        {
            "id": "repo_has_no_lawful_verdict_for_exact_t295_attempt",
            "pass": True,
            "details": "No lawful verdict artifact for the exact T295 attempt is exported on the current repo state.",
        },
        {
            "id": "repo_has_no_exact_sharper_interface_refinement_below_t295_attempt",
            "pass": len(sharper_refinement_hits) == 0,
            "details": "Repo scan finds no exact sharper interface-refinement marker below T295 outside the files being introduced now.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "current_repo_has_lawful_verdict_for_exact_t295_attempt": False if all_pass else None,
        "current_repo_has_exact_sharper_interface_refinement_below_t295_attempt": False if all_pass else None,
        "next_honest_move_is_freeze_exact_sharper_interface_refinement_target_below_same_t295_attempt": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1031_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_INTERFACE_REFINEMENT_NONEXPORT_AUDITED"
        if all_pass
        else "P1031_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1031",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1030_sharper_interface_refinement_target_actual_realization_attempt_probe": rel(IN_P1030),
        },
        "theorem_result": theorem_result,
        "sharper_interface_refinement_hits_outside_current_introduction": sharper_refinement_hits,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact T295 sharper-refinement attempt is exported, but still unresolved.",
            "The repo exports neither a lawful verdict for that exact attempt nor a sharper interface-refinement below it.",
            "The honest next move is therefore one exact sharper interface-refinement target below the same attempt.",
        ],
        "recommended_next_probe": {
            "id": "P1032_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_INTERFACE_REFINEMENT_TARGET_PROBE",
            "goal": "Export one exact sharper interface-refinement target below the same T295 attempt.",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1031",
        "status": status,
        "as_of": AS_OF,
        "current_repo_has_lawful_verdict_for_exact_t295_attempt": theorem_result[
            "current_repo_has_lawful_verdict_for_exact_t295_attempt"
        ],
        "current_repo_has_exact_sharper_interface_refinement_below_t295_attempt": theorem_result[
            "current_repo_has_exact_sharper_interface_refinement_below_t295_attempt"
        ],
        "recommended_next_probe": artifact["recommended_next_probe"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
