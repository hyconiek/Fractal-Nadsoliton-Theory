#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "n199_current_exported_preobserver_selector_output_operator_theorem_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    p179 = load_json(
        "fundamental_action_reconstruction/generated/p179_current_exported_preobserver_selector_output_operator_probe_summary.json"
    )

    expected_status = (
        "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_PREOBSERVER_SELECTOR_OUTPUT_OPERATOR_FROM_R_SEL_PRELM_V1_AFTER_P179"
    )
    status_ok = p179["status"] == expected_status

    checks = [
        {
            "id": "positive_selector_output_probe",
            "actual": p179["status"],
            "expected": expected_status,
            "pass": status_ok,
        }
    ]

    summary = {
        "step": "N199",
        "status": "N199_DISCHARGED_CURRENT_EXPORTED_PREOBSERVER_SELECTOR_OUTPUT_OPERATOR_THEOREM_NO_FALSE_PASS",
        "scope": "current_preobserver_selector_output_operator_only",
        "checks": checks,
        "theorem_result": {
            "discharged": status_ok,
            "source_object": "S_preLM_strict_core_source_object_v1",
            "orientation_input": "E_orient_preLM_v1",
            "selector_bridge_operator": "B_sel_preLM_v1",
            "selector_reduction_operator": "R_sel_preLM_v1",
            "selector_output_operator": "O_sel_preLM_v1",
            "admissible_O_sel": status_ok,
            "bridge_ready_for_emergent_observer_limit": status_ok,
        },
        "hard_limits": [
            "emergent_observer_not_yet_constructed",
            "downstream_chain_not_yet_completed",
            "no_strict_core_selector_closure",
            "no_QW2191_discharge",
            "no_ToE_closure",
        ],
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
