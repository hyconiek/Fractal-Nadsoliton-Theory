#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p204_current_actual_emergent_observer_closure_realization_map_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n223 = load_json(GENERATED / "n223_current_actual_emergent_observer_closure_commit_map_theorem_summary.json")
    f116 = load_json(GENERATED / "f116_first_actual_emergent_observer_closure_realization_map_packet_summary.json")

    checks = [
        {
            "id": "actual_closure_commit_is_admissible",
            "actual": n223["theorem_result"]["admissible_AF_obs_actual_closure_commit"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_commit_state",
            "actual": f116["observer_actual_closure_realization_properties"]["derived_only_from_actual_closure_commit_state"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": f116["observer_actual_closure_realization_properties"]["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_realization_only",
            "actual": f116["observer_actual_closure_realization_properties"]["downstream_actual_closure_realization_only"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_realization_amplitude",
            "actual": f116["source_observer_actual_closure_realization"]["positive_actual_closure_realization_amplitude"],
            "expected": True,
        },
        {
            "id": "actual_residual_annihilated",
            "actual": f116["source_observer_actual_closure_realization"]["actual_residual_annihilated"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": f116["observer_actual_closure_realization_properties"]["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": f116["observer_actual_closure_realization_properties"]["kernel_split_safe"],
            "expected": True,
        },
    ]

    for check in checks:
        check["pass"] = check["actual"] == check["expected"]

    status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_FROM_AF_OBS_ACTUAL_CLOSURE_COMMIT_PRELM_V1_AFTER_P204"
        if all(check["pass"] for check in checks)
        else "CURRENT_REPO_DOES_NOT_EXPORT_AN_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_REALIZATION_MAP_FROM_AF_OBS_ACTUAL_CLOSURE_COMMIT_PRELM_V1_AFTER_P204"
    )

    summary = {
        "stage": "P204",
        "lane": "current_actual_emergent_observer_closure_realization_only",
        "status": status,
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_commit_map": "AF_obs_actual_closure_commit_preLM_v1",
        "observer_actual_closure_realization_map": "AG_obs_actual_closure_realization_preLM_v1",
        "checks": checks,
        "actual_closure_realization_exported": all(check["pass"] for check in checks),
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
