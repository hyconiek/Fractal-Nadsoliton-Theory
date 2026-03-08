#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f87_first_exported_preobserver_source_object_future_bridge_compatibility_packet_summary.json"


def load_json(repo_relative_path: str) -> dict:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    f77 = load_json(
        "fundamental_action_reconstruction/generated/f77_first_additive_preobserver_source_object_admissibility_upgrade_target_packet_summary.json"
    )
    n189 = load_json(
        "fundamental_action_reconstruction/generated/n189_current_second_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n190 = load_json(
        "fundamental_action_reconstruction/generated/n190_current_third_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n192 = load_json(
        "fundamental_action_reconstruction/generated/n192_current_fifth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )
    n193 = load_json(
        "fundamental_action_reconstruction/generated/n193_current_sixth_admissibility_clause_discharge_theorem_for_s_prelm_strict_core_source_object_v1_summary.json"
    )

    guardrails = f77["guardrails"]

    summary = {
        "stage": "F87",
        "lane": "first_exported_preobserver_source_object_future_bridge_compatibility_only",
        "status": "F87_EXECUTED_FIRST_EXPORTED_PREOBSERVER_SOURCE_OBJECT_FUTURE_BRIDGE_COMPATIBILITY_PACKET_NO_FALSE_PASS",
        "exported_object": "S_preLM_strict_core_source_object_v1",
        "second_clause_discharged": n189["theorem_result"]["discharged"],
        "third_clause_discharged": n190["theorem_result"]["discharged"],
        "fifth_clause_discharged": n192["theorem_result"]["discharged"],
        "sixth_clause_discharged": n193["theorem_result"]["discharged"],
        "source_object_first": guardrails["source_object_first"],
        "upstream_of_observer": guardrails["upstream_of_observer"],
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
