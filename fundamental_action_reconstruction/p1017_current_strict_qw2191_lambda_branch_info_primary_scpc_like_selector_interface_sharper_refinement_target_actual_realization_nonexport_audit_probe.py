#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1016 = GENERATED / "p1016_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_exact_sharper_interface_refinement_target_probe.json"

OUT = GENERATED / "p1017_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1017_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_nonexport_audit_probe_summary.json"

ACTUAL_REALIZATION_OBJECT_ID = (
    "lambda_nad12_sigma_residual_shannon_information_primary_"
    "sparse_competitive_predictive_coding_like_selector_interface_"
    "sharper_refinement_target_actual_realization_v1"
)
ACTUAL_REALIZATION_STEM = (
    "information_primary_sparse_competitive_predictive_coding_like_"
    "selector_interface_sharper_refinement_target_actual_realization"
)

EXCLUDED_PATHS = {
    "fundamental_action_reconstruction/P1017_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
    "fundamental_action_reconstruction/N850_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_THEOREM.md",
    "fundamental_action_reconstruction/T289_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
    "fundamental_action_reconstruction/N851_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_THEOREM.md",
    "fundamental_action_reconstruction/p1017_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_nonexport_audit_probe.py",
    "fundamental_action_reconstruction/p1018_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe.py",
    "fundamental_action_reconstruction/generated/p1017_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_nonexport_audit_probe.json",
    "fundamental_action_reconstruction/generated/p1017_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_nonexport_audit_probe_summary.json",
    "fundamental_action_reconstruction/generated/p1018_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe.json",
    "fundamental_action_reconstruction/generated/p1018_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_sharper_refinement_target_actual_realization_attempt_probe_summary.json",
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
            hits.append({"path": repo_rel, "match": "candidate_specific_sharper_interface_actual_realization_marker"})
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1016]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1017",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1016 = load_json(IN_P1016)
    actual_realization_hits = search_repo_for_actual_realization()

    checks = [
        {
            "id": "p1016_already_exports_exact_t288_target",
            "pass": (
                p1016.get("status")
                == "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_INTERFACE_REFINEMENT_TARGET_EXPORTED"
                and p1016.get("t288_target_exported_on_current_repo_state") is True
                and p1016.get("t288_target_keeps_strict_selector_source_open") is True
            ),
            "details": "P1016 already exports the exact T288 sharper interface-refinement target.",
        },
        {
            "id": "repo_has_no_exact_actual_realization_of_t288_target",
            "pass": len(actual_realization_hits) == 0,
            "details": "Repo scan finds no exact actual-realization marker for the T288 sharper target outside the files being introduced now.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "current_repo_has_exported_actual_realization_of_t288_target": False if all_pass else None,
        "t288_target_still_remains_future_only_not_actual_export": True if all_pass else None,
        "next_honest_move_is_exact_actual_realization_attempt_of_same_t288_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1017_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if all_pass
        else "P1017_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1017",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1016_sharper_interface_refinement_target_probe": rel(IN_P1016),
        },
        "theorem_result": theorem_result,
        "actual_realization_hits_outside_current_introduction": actual_realization_hits,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The exact T288 sharper target is exported, but still future-only.",
            "No exact actual realization of that target is exported on the current repo state.",
            "The honest next move is therefore one exact actual-realization attempt below that same sharper target.",
        ],
        "recommended_next_probe": {
            "id": "P1018_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_SHARPER_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_PROBE",
            "goal": "Export one exact actual-realization attempt below the same T288 sharper target.",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1017",
        "status": status,
        "as_of": AS_OF,
        "current_repo_has_exported_actual_realization_of_t288_target": theorem_result[
            "current_repo_has_exported_actual_realization_of_t288_target"
        ],
        "t288_target_still_remains_future_only_not_actual_export": theorem_result[
            "t288_target_still_remains_future_only_not_actual_export"
        ],
        "recommended_next_probe": artifact["recommended_next_probe"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
