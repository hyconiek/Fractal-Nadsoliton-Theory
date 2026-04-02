#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1013 = GENERATED / "p1013_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_nonexport_audit_probe.json"
IN_F951 = GENERATED / "f951_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_target_packet_summary.json"

OUT = GENERATED / "p1014_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_probe.json"
OUT_SUMMARY = GENERATED / "p1014_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_interface_actual_realization_attempt_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P1013, IN_F951]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1014",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1013 = load_json(IN_P1013)
    f951 = load_json(IN_F951)
    theorem = p1013.get("theorem_result") or {}

    status = (
        "PASS_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
        if p1013.get("status")
        == "P1013_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        and theorem.get("current_repo_has_exported_actual_realization_of_selector_interface_target") is False
        and theorem.get("next_honest_move_is_exact_actual_realization_attempt_of_same_selector_interface_target")
        is True
        and f951.get("status")
        == "F951_EXECUTED_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
        else "P1014_REQUIRES_REVIEW"
    )

    attempt_name = (
        "lambda_nad12_sigma_residual_shannon_information_primary_"
        "sparse_competitive_predictive_coding_like_selector_interface_"
        "actual_realization_attempt_v1"
    )

    artifact = {
        "stage": "P1014",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p1013_selector_interface_actual_realization_nonexport_audit_probe": rel(IN_P1013),
            "f951_selector_interface_target_packet_summary": rel(IN_F951),
        },
        "t287_attempt_name": attempt_name,
        "t287_attempt_shape": {
            "reference_lane_packet_ref": "lambda_nad12_sigma_residual_shannon_information_primary_sparse_competitive_predictive_coding_like_selector_provider_shift_candidate_reference_lane_v1",
            "selector_interface_target_ref": "lambda_nad12_sigma_residual_shannon_information_primary_sparse_competitive_predictive_coding_like_selector_interface_target_v1",
            "pivot_family_ref": "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1",
            "pivot_branch_ref": "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1",
        },
        "t287_attempt_exported_on_current_repo_state": status.startswith("PASS_"),
        "t287_attempt_keeps_strict_selector_source_open": True,
        "hard_limits": [
            "No actual selector-interface realization is claimed.",
            "No actual strict selector source export is claimed.",
            "No T176 export is claimed.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1014",
        "status": status,
        "as_of": AS_OF,
        "t287_attempt_exported_on_current_repo_state": artifact["t287_attempt_exported_on_current_repo_state"],
        "t287_attempt_keeps_strict_selector_source_open": artifact[
            "t287_attempt_keeps_strict_selector_source_open"
        ],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
