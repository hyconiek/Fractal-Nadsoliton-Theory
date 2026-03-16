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
    / "p645_current_downstream_completion_branch_discharge_probe_for_s_sel_int_candidate_seed_v1.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p645_current_downstream_completion_branch_discharge_probe_for_s_sel_int_candidate_seed_v1_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f645 = load_json(
        "fundamental_action_reconstruction/generated/f645_last_downstream_completion_branch_packet_for_s_sel_int_candidate_seed_v1_summary.json"
    )
    n535 = load_json(
        "fundamental_action_reconstruction/generated/n535_current_orientation_export_branch_obstruction_theorem_for_s_sel_int_candidate_seed_v1_summary.json"
    )

    explicit_downstream_completion_branch_discharge_exported = False
    try:
        n546 = load_json(
            "fundamental_action_reconstruction/generated/n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
        )
        n547 = load_json(
            "fundamental_action_reconstruction/generated/n547_current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_theorem_summary.json"
        )
        n548 = load_json(
            "fundamental_action_reconstruction/generated/n548_current_exported_s_sel_int_strict_core_source_object_selector_reduction_operator_theorem_summary.json"
        )
        n549 = load_json(
            "fundamental_action_reconstruction/generated/n549_current_exported_s_sel_int_strict_core_source_object_selector_output_operator_theorem_summary.json"
        )
        explicit_downstream_completion_branch_discharge_exported = bool(
            n546["theorem_result"]["admissible_E_orient"]
            and n547["theorem_result"]["admissible_B_sel"]
            and n548["theorem_result"]["admissible_R_sel"]
            and n549["theorem_result"]["admissible_O_sel"]
        )
    except Exception:
        explicit_downstream_completion_branch_discharge_exported = False

    checks_spec = [
        {
            "id": "f645_downstream_branch_last",
            "actual": f645["last_lower_branch_to_attack"],
            "expected": "future_completion_of_B_sel_R_sel_O_sel_after_new_source_object_construction",
            "meaning": "F645 freezes downstream completion as the last remaining lower branch (seed v1)",
        },
        {
            "id": "n535_orientation_export_obstruction_discharged",
            "actual": n535["theorem_result"]["discharged"],
            "expected": True,
            "meaning": "N535 packages the seed-v1 orientation-export obstruction",
        },
        {
            "id": "explicit_downstream_completion_branch_discharge_exported",
            "actual": explicit_downstream_completion_branch_discharge_exported,
            "expected": True,
            "meaning": "the current repo exports an explicit downstream-completion discharge in the seed-v1 lane (B_sel/R_sel/O_sel chain)",
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

    any_fail = any(not c["pass"] for c in checks)
    status = (
        "P645_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_DOWNSTREAM_COMPLETION_BRANCH_STATE_FOR_SEED_V1"
        if any_fail
        else "CURRENT_REPO_EXPORTS_AN_EXPLICIT_DOWNSTREAM_COMPLETION_BRANCH_DISCHARGE_FOR_SEED_V1_AFTER_P645"
    )

    artifact = {
        "stage": "P645",
        "lane": "current_downstream_completion_branch_discharge_only",
        "goal": "test_whether_the_current_repo_already_exports_an_explicit_downstream_completion_branch_discharge_for_seed_v1_last_remaining_lower_branch",
        "status": status,
        "target_state": {
            "downstream_completion_branch_selected_as_next_attack": True,
            "explicit_downstream_completion_branch_discharge_exported": explicit_downstream_completion_branch_discharge_exported,
            "remaining_open_branches": [],
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P645",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(
        json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    OUT_SUMMARY.write_text(
        json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii"
    )
    print(OUT_JSON)


if __name__ == "__main__":
    main()
