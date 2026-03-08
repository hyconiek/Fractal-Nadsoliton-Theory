#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUT = ROOT / "generated" / "f128_first_source_topology_selector_promotion_target_packet_summary.json"


def main() -> None:
    summary = {
        "packet_id": "F128",
        "status": "F128_EXECUTED_FIRST_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": "tau_src_candidate_v1",
        "promotion_target": "Pi_sel_src_target_v1",
        "codomain_packet": {
            "name": "Sigma_sel_src_target_v1",
            "components": [
                "selector_axis_class_v1",
                "selector_signed_split_class_v1",
                "preobserver_scope_tag_v1",
            ],
        },
        "downstream_chart_realization_candidates": [
            "E_orient_preLM_v1",
            "B_sel_preLM_v1",
            "R_sel_preLM_v1",
            "O_sel_preLM_v1",
        ],
        "observer_role": "downstream_only",
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
