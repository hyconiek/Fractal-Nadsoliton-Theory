#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "p202_current_actual_emergent_observer_closure_object_candidate_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    n221 = load_json(GENERATED / "n221_current_actual_emergent_observer_closure_readout_operator_theorem_summary.json")
    f114 = load_json(GENERATED / "f114_first_actual_emergent_observer_closure_object_candidate_packet_summary.json")

    checks = [
        {
            "id": "observer_actual_closure_readout_is_admissible",
            "actual": n221["theorem_result"]["admissible_AD_obs_actual_closure_readout"],
            "expected": True,
        },
        {
            "id": "derived_only_from_actual_closure_readout",
            "actual": f114["observer_actual_closure_object_candidate_properties"]["derived_only_from_actual_closure_readout"],
            "expected": True,
        },
        {
            "id": "strict_core_only",
            "actual": f114["observer_actual_closure_object_candidate_properties"]["strict_core_only"],
            "expected": True,
        },
        {
            "id": "downstream_actual_closure_object_candidate_only",
            "actual": f114["observer_actual_closure_object_candidate_properties"]["downstream_actual_closure_object_candidate_only"],
            "expected": True,
        },
        {
            "id": "positive_actual_closure_object_candidate_amplitude",
            "actual": f114["source_observer_actual_closure_object_candidate"]["positive_actual_closure_object_candidate_amplitude"],
            "expected": True,
        },
        {
            "id": "actual_gap_annihilated",
            "actual": f114["source_observer_actual_closure_object_candidate"]["actual_gap_annihilated"],
            "expected": True,
        },
        {
            "id": "observer_information_deficit_downstream",
            "actual": f114["observer_actual_closure_object_candidate_properties"]["observer_information_deficit_is_downstream_symptom"],
            "expected": True,
        },
        {
            "id": "kernel_split_safe",
            "actual": f114["observer_actual_closure_object_candidate_properties"]["kernel_split_safe"],
            "expected": True,
        },
    ]

    for check in checks:
        check["pass"] = check["actual"] == check["expected"]

    status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_CANDIDATE_FROM_AD_OBS_ACTUAL_CLOSURE_READOUT_PRELM_V1_AFTER_P202"
        if all(check["pass"] for check in checks)
        else "CURRENT_REPO_DOES_NOT_EXPORT_AN_ADMISSIBLE_ACTUAL_EMERGENT_OBSERVER_CLOSURE_OBJECT_CANDIDATE_FROM_AD_OBS_ACTUAL_CLOSURE_READOUT_PRELM_V1_AFTER_P202"
    )

    summary = {
        "stage": "P202",
        "lane": "current_actual_emergent_observer_closure_object_candidate_only",
        "status": status,
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_readout_operator": "AD_obs_actual_closure_readout_preLM_v1",
        "observer_actual_closure_object_candidate_operator": "AE_obs_actual_closure_object_candidate_preLM_v1",
        "checks": checks,
        "actual_closure_object_candidate_exported": all(check["pass"] for check in checks),
        "actual_emergent_observer_closure": False,
        "strict_core_selector_closure": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
