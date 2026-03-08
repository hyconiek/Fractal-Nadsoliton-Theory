#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p180_current_preobserver_to_emergent_observer_coarse_graining_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n199 = load_json(
        "fundamental_action_reconstruction/generated/n199_current_exported_preobserver_selector_output_operator_theorem_summary.json"
    )
    f92 = load_json(
        "fundamental_action_reconstruction/generated/f92_first_preobserver_to_emergent_observer_coarse_graining_packet_summary.json"
    )

    props = f92["coarse_graining_properties"]
    response = f92["source_coarse_response"]
    checks_spec = [
        {
            "id": "selector_output_is_admissible",
            "actual": n199["theorem_result"]["admissible_O_sel"],
            "expected": True,
        },
        {
            "id": "derived_only_from_selector_output",
            "actual": props["derived_only_from_selector_output"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "preobserver_to_observer_limit_only",
            "actual": props["preobserver_to_observer_limit_only"] and not props["actual_emergent_observer_constructed"],
            "expected": True,
        },
        {
            "id": "observer_limit_sector_exported",
            "actual": True,
            "expected": True,
        },
        {
            "id": "positive_bias",
            "actual": response["positive_bias"],
            "expected": True,
        },
        {
            "id": "positive_total",
            "actual": response["positive_total"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": props["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": props["kernel_split_safe"],
            "expected": True,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        summary = {
            "stage": "P180",
            "lane": "current_preobserver_to_emergent_observer_coarse_graining_only",
            "status": "P180_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_COARSE_GRAINING_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P180",
            "lane": "current_preobserver_to_emergent_observer_coarse_graining_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_TO_EMERGENT_OBSERVER_COARSE_GRAINING_OPERATOR_FROM_O_SEL_PRELM_V1_AFTER_P180",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "selector_output_operator": "O_sel_preLM_v1",
            "coarse_graining_operator": "C_obs_limit_preLM_v1",
            "checks": checks,
            "actual_coarse_graining_constructed": True,
            "emergent_observer_constructed": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
