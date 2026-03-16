#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "f653_first_exported_s_sel_int_strict_core_source_object_future_bridge_compatibility_packet_summary.json"
)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f647 = load_json(
        "fundamental_action_reconstruction/generated/f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt_summary.json"
    )
    n541 = load_json(
        "fundamental_action_reconstruction/generated/n541_current_second_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    n542 = load_json(
        "fundamental_action_reconstruction/generated/n542_current_third_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )
    n544 = load_json(
        "fundamental_action_reconstruction/generated/n544_current_fifth_admissibility_clause_discharge_theorem_for_s_sel_int_strict_core_source_object_v1_summary.json"
    )

    props = f647["strict_core_export_properties"]

    summary = {
        "stage": "F653",
        "lane": "first_exported_s_sel_int_strict_core_source_object_future_bridge_compatibility_only",
        "status": "F653_EXECUTED_FIRST_EXPORTED_S_SEL_INT_STRICT_CORE_SOURCE_OBJECT_FUTURE_BRIDGE_COMPATIBILITY_PACKET_NO_FALSE_PASS",
        "exported_object": "S_sel_int_strict_core_source_object_v1",
        "second_clause_discharged": bool(n541["theorem_result"]["discharged"]),
        "third_clause_discharged": bool(n542["theorem_result"]["discharged"]),
        "selector_acceptance_independence_clause_discharged": bool(n544["theorem_result"]["discharged"]),
        "kernel_split_safe": bool(props["kernel_split_safe"]),
        "no_external_selector_import": bool(props["no_external_selector_import"]),
        "upstream_of_observer": bool(props["upstream_of_observer"]),
        "full_admissibility_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

