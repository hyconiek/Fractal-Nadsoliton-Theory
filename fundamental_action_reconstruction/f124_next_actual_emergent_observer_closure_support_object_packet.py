#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f124_next_actual_emergent_observer_closure_support_object_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f123 = load_json(
        GENERATED / "f123_next_actual_emergent_observer_closure_fixed_point_test_packet_summary.json"
    )

    aa_actual_closure_fix_v2 = f123["source_observer_actual_closure_fixed_point"][
        "aa_actual_closure_fix_v2"
    ]
    summary = {
        "stage": "F124",
        "status": "F124_EXECUTED_NEXT_ACTUAL_EMERGENT_OBSERVER_CLOSURE_SUPPORT_OBJECT_PACKET_NO_FALSE_PASS",
        "input_object": "AN_obs_actual_closure_fixed_point_test_preLM_v1",
        "source_object": "S_preLM_strict_core_source_object_v1",
        "observer_actual_closure_support_operator": {
            "operator": "AO_obs_actual_closure_support_preLM_v1",
            "domain_basis": ["aa_actual_closure_fix"],
            "codomain_basis": ["ab_actual_closure_support"],
            "matrix": [[1.0]],
            "actual_closure_support_dimension": 1,
        },
        "source_observer_actual_closure_support": {
            "input_basis": ["aa_actual_closure_fix"],
            "input_vector": [aa_actual_closure_fix_v2],
            "output_basis": ["ab_actual_closure_support"],
            "output_vector": [aa_actual_closure_fix_v2],
            "ab_actual_closure_support_v2": aa_actual_closure_fix_v2,
            "positive_support_amplitude": aa_actual_closure_fix_v2 > 0.0,
            "one_dimensional_actual_closure_support": True,
        },
        "observer_actual_closure_support_properties": {
            "derived_only_from_actual_closure_fixed_point_test": True,
            "strict_core_only": True,
            "downstream_actual_closure_support_only": True,
            "kernel_split_safe": True,
            "actual_emergent_observer_closure": False,
            "observer_information_deficit_is_downstream_symptom": True,
        },
        "hard_limits": [
            "not_actual_emergent_observer_closure",
            "no_QW2191_discharge",
            "no_strict_core_selector_closure",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
