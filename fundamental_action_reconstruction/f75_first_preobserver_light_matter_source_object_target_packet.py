#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f75_first_preobserver_light_matter_source_object_target_packet_summary.json"


def main() -> None:
    phi = 0.16250
    kernel_core_amplitude = math.cos(phi)

    summary = {
        "stage": "F75",
        "lane": "future_preobserver_light_matter_source_object_target_only",
        "goal": "reduce_the_first_guardrail_consistent_preobserver_provider_packet_to_one_explicit_upstream_source_object_target",
        "status": "F75_EXECUTED_FIRST_PREOBSERVER_LIGHT_MATTER_SOURCE_OBJECT_TARGET_PACKET_NO_FALSE_PASS",
        "provider_packet": "preobserver_light_matter_source_provider_packet_v1",
        "source_object_target": "preobserver_light_matter_source_object_target_v1",
        "carrier": {
            "used": ["V_topo", "L_int", "M_int"],
            "observer_excluded": True,
        },
        "profile": {
            "T_flow": "cos(phi) * e_topo",
            "L_seed": "P_NL^(0) T_flow^(0)",
            "M_seed": "P_LM(d) P_NL^(0) T_flow^(0)",
            "assembled_target": "T_flow^(0) ⊕ L_seed^(0) ⊕ M_seed(d)",
            "kernel_core_amplitude": kernel_core_amplitude,
        },
        "guardrails": {
            "future_only": True,
            "genuinely_additive_target_only": True,
            "upstream_of_observer": True,
            "light_before_observer": True,
            "matter_as_encoding_not_primary_selector_source": True,
            "kernel_split_safe": True,
            "no_external_selector_import": True,
            "source_object_first": True,
        },
        "theorem_level_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
