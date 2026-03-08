#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p203_current_actual_emergent_observer_closure_commit_map_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n222 = load_json(GENERATED / "n222_current_actual_emergent_observer_closure_object_candidate_theorem_summary.json")
    f115 = load_json(GENERATED / "f115_first_actual_emergent_observer_closure_commit_map_packet_summary.json")

    checks = [
        {
            "id": "actual_closure_object_candidate_is_admissible",
            "actual": n222["theorem_result"]["admissible_AE_obs_actual_closure_object_candidate"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_object_candidate_state",
            "actual": f115["observer_actual_closure_commit_properties"]["derived_only_from_actual_closure_object_candidate_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": f115["observer_actual_closure_commit_properties"]["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_commit_only",
            "actual": f115["observer_actual_closure_commit_properties"]["downstream_actual_closure_commit_only"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_commit_amplitude",
            "actual": f115["source_observer_actual_closure_commit"]["positive_actual_closure_commit_amplitude"],
            "expected": True,
        },
        {
            "id": "zero_actual_closure_residual",
            "actual": f115["source_observer_actual_closure_commit"]["zero_actual_closure_residual"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": f115["observer_actual_closure_commit_properties"]["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": f115["observer_actual_closure_commit_properties"]["kernel_split_safe"],
            "expected": True,
        },
    ]

    for check in checks:
        check["pass"] = check["actual"] == check["expected"]

    status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_COMMIT_MAP_FROM_AE_OBS_ACTUAL_CLOSURE_OBJECT_CANDIDATE_PRELM_V1_AFTER_P203"
        if all(check["pass"] for check in checks)
        else "CURRENT_REPO_DOES_NOT_EXPORT_AN_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_COMMIT_MAP_FROM_AE_OBS_ACTUAL_CLOSURE_OBJECT_CANDIDATE_PRELM_V1_AFTER_P203"
    )

    summary = {
        "stage": "P203",
        "lane": "current_actual_emergent_observer_closure_commit_only",
        "status": status,
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_object_candidate_operator": "AE_obs_actual_closure_object_candidate_preLM_v1",
        "observer_actual_closure_commit_map": "AF_obs_actual_closure_commit_preLM_v1",
        "checks": checks,
        "actual_closure_commit_exported": all(check["pass"] for check in checks),
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
