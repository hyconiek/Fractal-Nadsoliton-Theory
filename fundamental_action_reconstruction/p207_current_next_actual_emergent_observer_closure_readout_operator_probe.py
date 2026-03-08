#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p207_current_next_actual_emergent_observer_closure_readout_operator_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n226 = load_json(
        "fundamental_action_reconstruction/generated/n226_current_next_actual_emergent_observer_closure_support_object_theorem_summary.json"
    )
    f119 = load_json(
        "fundamental_action_reconstruction/generated/f119_next_actual_emergent_observer_closure_readout_operator_packet_summary.json"
    )

    props = f119["observer_actual_closure_readout_properties"]
    response = f119["source_observer_actual_closure_readout"]
    checks_spec = [
        {
            "id": "actual_closure_support_is_admissible",
            "actual": n226["theorem_result"]["admissible_AI_obs_actual_closure_support_object"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_support_state",
            "actual": props["derived_only_from_actual_closure_support_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_readout_only",
            "actual": props["downstream_actual_closure_readout_only"]
            and not props["actual_emergent_observer_closure"],
            "expected": True,
        },
        {
            "id": "positive_commit_amplitude",
            "actual": response["positive_commit_amplitude"],
            "expected": True,
        },
        {
            "id": "zero_gap_channel",
            "actual": response["zero_gap_channel"],
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
            "stage": "P207",
            "lane": "current_next_actual_emergent_observer_closure_readout_only",
            "status": "P207_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P207",
            "lane": "current_next_actual_emergent_observer_closure_readout_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_FROM_AI_OBS_ACTUAL_CLOSURE_SUPPORT_PRELM_V1_AFTER_P207",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_actual_closure_support_object": "AI_obs_actual_closure_support_preLM_v1",
            "observer_actual_closure_readout_operator": "AJ_obs_actual_closure_readout_preLM_v1",
            "checks": checks,
            "actual_closure_readout_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
