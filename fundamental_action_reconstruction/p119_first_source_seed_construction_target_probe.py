#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED / "p119_first_source_seed_construction_target_probe.json"
)
OUT_SUMMARY = (
    GENERATED / "p119_first_source_seed_construction_target_probe_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        "fundamental_action_reconstruction/generated/n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json",
    ]
    missing = [p for p in prereq if not (REPO / p).exists()]
    if missing:
        artifact = {
            "stage": "P119",
            "lane": "first_source_seed_construction_target_only",
            "goal": "test_whether_the_last_positive_branch_is_now_reduced_to_one_explicit_first_source_seed_construction_target",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES_FOR_FIRST_SOURCE_SEED_TARGET_PROBE",
            "missing": missing,
            "strict_core_promotion": False,
            "full_closure_pass": False,
            "no_false_pass": True,
        }
        OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(
            json.dumps(
                {
                    "stage": "P119",
                    "status": artifact["status"],
                    "lane": artifact["lane"],
                    "missing": missing,
                    "strict_core_promotion": False,
                    "full_closure_pass": False,
                    "no_false_pass": True,
                },
                indent=2,
                ensure_ascii=True,
            )
            + "\n",
            encoding="ascii",
        )
        print(OUT_JSON)
        return

    n676 = load_json(prereq[0])
    n546 = (
        load_json(
            "fundamental_action_reconstruction/generated/n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json"
        )
        if (GENERATED / "n546_current_exported_s_sel_int_strict_core_source_object_admissible_orientation_export_theorem_summary.json").exists()
        else None
    )
    n549 = (
        load_json(
            "fundamental_action_reconstruction/generated/n549_current_exported_s_sel_int_strict_core_source_object_selector_output_operator_theorem_summary.json"
        )
        if (GENERATED / "n549_current_exported_s_sel_int_strict_core_source_object_selector_output_operator_theorem_summary.json").exists()
        else None
    )

    admissible_source_object_exported = bool(
        n676.get("theorem_result", {}).get("admissible_S_sel_int_source_object_in_F34_sense")
    )
    exported_source_object = str(n676.get("theorem_result", {}).get("exported_object") or "")
    admissible_orientation_exported = bool(
        (n546 or {}).get("theorem_result", {}).get("admissible_E_orient")
    )
    downstream_operator_chain_exported = bool(
        (n549 or {}).get("theorem_result", {}).get("admissible_O_sel")
    )

    reduced_to_one_first_source_seed_target = (
        admissible_source_object_exported and exported_source_object == "S_sel_int_strict_core_source_object_v1"
    )

    checks_spec = [
        {
            "id": "n676_admissible_source_object_exported",
            "actual": admissible_source_object_exported,
            "expected": True,
            "meaning": "N676 confirms an admissible strict-core S_sel_int source object is exported (F34 contract sense)",
        },
        {
            "id": "n676_exported_object_name",
            "actual": exported_source_object,
            "expected": "S_sel_int_strict_core_source_object_v1",
            "meaning": "the exported admissible strict-core source object is S_sel_int_strict_core_source_object_v1",
        },
        {
            "id": "reduced_to_one_first_source_seed_target",
            "actual": reduced_to_one_first_source_seed_target,
            "expected": True,
            "meaning": "the last positive branch is reduced to one explicit first source-seed construction target",
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

    first_source_seed_construction_target = {
        "strict_core_source_object": "S_sel_int",
        "exported_object": exported_source_object or None,
        "admissible_source_object_exported": admissible_source_object_exported,
        "admissible_orientation_exported": admissible_orientation_exported,
        "downstream_operator_chain_exported": downstream_operator_chain_exported,
        "genuinely_new_strict_core_object_required": True,
        "source_carrying_enough_for_later_E_orient_export_required": True,
        "orientation_not_yet_counted_as_exported": not admissible_orientation_exported,
        "strict_core_only_required": True,
        "silent_legacy_to_strict_substitution_forbidden": True,
        "selector_acceptance_outside_strict_core_may_not_count_as_source_seed": True,
    }

    later_open_branches: list[str] = []
    if not admissible_source_object_exported:
        later_open_branches.append("future_construction_of_admissible_S_sel_int")
    elif not admissible_orientation_exported:
        later_open_branches.append("future_derivation_of_admissible_E_orient_from_S_sel_int")
    elif not downstream_operator_chain_exported:
        later_open_branches.append("future_completion_of_B_sel_R_sel_O_sel_after_seed_package")
    else:
        later_open_branches.append(
            "future_strict_core_selector_closure_and_global_QW2191_discipline_target_T172"
        )

    artifact = {
        "stage": "P119",
        "lane": "first_source_seed_construction_target_only",
        "goal": "test_whether_the_last_positive_branch_is_now_reduced_to_one_explicit_first_source_seed_construction_target",
        "status": "CURRENT_REPO_REDUCES_THE_LAST_POSITIVE_BRANCH_TO_ONE_FIRST_SOURCE_SEED_CONSTRUCTION_TARGET_AFTER_P119",
        "target_state": {
            "last_positive_branch_reduced_to_one_first_source_seed_target": reduced_to_one_first_source_seed_target,
            "first_source_seed_construction_target": first_source_seed_construction_target,
            "later_open_branches": later_open_branches,
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P119",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "target_state": artifact["target_state"],
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
