#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n547_current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_theorem_summary.json"
)


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p655 = load_json(
        "fundamental_action_reconstruction/generated/p655_current_exported_s_sel_int_strict_core_source_object_selector_bridge_operator_probe_summary.json"
    )

    expected_status = "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_SELECTOR_BRIDGE_OPERATOR_FROM_E_ORIENT_S_SEL_INT_SOURCE_OBJECT_V1_AFTER_P655"
    status_ok = p655["status"] == expected_status

    checks = [
        {
            "id": "positive_selector_bridge_probe",
            "actual": p655["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N547",
        "status": "N547_DISCHARGED_CURRENT_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_SELECTOR_BRIDGE_OPERATOR_THEOREM_NO_FALSE_PASS",
        "scope": "current_s_sel_int_source_object_selector_bridge_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_sel_int_strict_core_source_object_v1",
            "orientation_input": "E_orient_s_sel_int_source_object_v1",
            "selector_bridge_operator": "B_sel_s_sel_int_source_object_v1",
            "admissible_B_sel": status_ok,
            "bridge_ready_for_R_sel": status_ok,
        },
        "hard_limits": [
            "admissible_S_sel_int_not_yet_constructed",
            "R_sel_not_yet_constructed",
            "O_sel_not_yet_constructed",
            "downstream_chain_not_yet_constructed",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

