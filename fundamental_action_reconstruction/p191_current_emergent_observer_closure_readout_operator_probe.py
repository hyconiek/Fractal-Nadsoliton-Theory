#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p191_current_emergent_observer_closure_readout_operator_probe_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    n210 = load_json(
        "fundamental_action_reconstruction/generated/n210_current_emergent_observer_closure_support_object_theorem_summary.json"
    )
    f103 = load_json(
        "fundamental_action_reconstruction/generated/f103_first_emergent_observer_closure_readout_operator_packet_summary.json"
    )

    props = f103["observer_closure_readout_properties"]
    response = f103["source_observer_closure_readout"]
    checks_spec = [
        {
            "id": "observer_closure_support_is_admissible",
            "actual": n210["theorem_result"]["admissible_S_obs_closure_support"],
            "expected": True,
        },
        {
            "id": "derived_only_from_closure_support_state",
            "actual": props["derived_only_from_closure_support_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": props["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_closure_readout_only",
            "actual": props["downstream_closure_readout_only"]
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
            "stage": "P191",
            "lane": "current_emergent_observer_closure_readout_only",
            "status": "P191_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_EMERGENT_OBSERVER_CLOSURE_READOUT_STATE",
            "checks": checks,
            "blocking_mismatches": mismatches,
        }
    else:
        summary = {
            "stage": "P191",
            "lane": "current_emergent_observer_closure_readout_only",
            "status": "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_EMERGENT_OBSERVER_CLOSURE_READOUT_OPERATOR_FROM_S_OBS_CLOSURE_SUPPORT_PRELM_V1_AFTER_P191",
            "source_object": "S_preLM_strict_core_source_object_v1",
            "observer_closure_support_map": "S_obs_closure_support_preLM_v1",
            "observer_closure_readout_operator": "T_obs_closure_readout_preLM_v1",
            "checks": checks,
            "closure_readout_exported": True,
            "actual_emergent_observer_closure": False,
            "strict_core_selector_closure": False,
            "no_false_pass": True,
        }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
