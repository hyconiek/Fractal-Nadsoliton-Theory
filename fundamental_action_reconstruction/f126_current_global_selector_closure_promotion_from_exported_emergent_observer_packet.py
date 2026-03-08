#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f126_current_global_selector_closure_promotion_from_exported_emergent_observer_packet_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    p213 = load_json(
        GENERATED / "p213_current_next_actual_emergent_observer_closure_readout_operator_probe_summary.json"
    )
    n118 = load_json(
        GENERATED
        / "n118_current_selector_or_symmetry_breaking_requirement_theorem_for_qw2191_summary.json"
    )
    n163 = load_json(
        GENERATED
        / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n233 = load_json(
        GENERATED
        / "n233_current_next_actual_emergent_observer_closure_readout_operator_theorem_summary.json"
    )

    summary = {
        "stage": "F126",
        "status": "F126_EXECUTED_CURRENT_GLOBAL_SELECTOR_CLOSURE_PROMOTION_FROM_EXPORTED_EMERGENT_OBSERVER_PACKET_NO_FALSE_PASS",
        "scope": "current_global_selector_closure_and_qw2191_discharge_promotion_question_only",
        "fixed_exported_composite": [
            "S_preLM_strict_core_source_object_v1",
            "E_orient_preLM_v1",
            "B_sel_preLM_v1",
            "R_sel_preLM_v1",
            "O_sel_preLM_v1",
            "AO_obs_actual_closure_support_preLM_v1",
            "AP_obs_actual_closure_readout_preLM_v1",
        ],
        "promotion_targets": [
            "global_selector_closure_from_exported_emergent_observer_composite",
            "global_QW2191_discharge_from_exported_emergent_observer_composite",
        ],
        "required_promotion_criteria": [
            {
                "id": "admissible_local_positive_composite_with_stable_asymmetry",
                "currently_exported": p213["actual_closure_readout_exported"],
            },
            {
                "id": "theorem_level_actual_emergent_observer_closure",
                "currently_exported": False,
            },
            {
                "id": "observer_not_downstream_only",
                "currently_exported": False,
            },
            {
                "id": "theorem_level_strict_core_selector_closure",
                "currently_exported": False,
            },
            {
                "id": "basis_independent_global_promotion_witness",
                "currently_exported": False,
            },
            {
                "id": "theorem_level_global_qw2191_discharge",
                "currently_exported": False,
            },
            {
                "id": "kernel_split_safe_promotion_semantics",
                "currently_exported": True,
            },
        ],
        "anchors": {
            "p213_status": p213["status"],
            "n233_status": n233["status"],
            "n118_status": n118["status"],
            "n163_status": n163["status"],
        },
        "global_promotion_not_yet_decided": True,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
