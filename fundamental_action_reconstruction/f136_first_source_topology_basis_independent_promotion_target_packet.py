#!/usr/bin/env python3
import json
from pathlib import Path


def main() -> None:
    out = {
        "status": "F136_EXECUTED_FIRST_SOURCE_TOPOLOGY_BASIS_INDEPENDENT_PROMOTION_TARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "promotion_target_name": "Upsilon_sel_basis_target_v1",
        "domain_targets": [
            "Theta_src_nontriv_discharge_target_v1",
            "Pi_sel_src_target_v1",
        ],
        "codomain_target": "Sigma_sel_basis_free_target_v1",
        "codomain_packet": {
            "name": "Sigma_sel_basis_free_target_v1",
            "components": [
                "selector_axis_basis_free_class_v1",
                "selector_signed_split_basis_free_class_v1",
                "preobserver_basis_free_scope_tag_v1",
            ],
        },
        "scope": {
            "source_side_only": True,
            "future_only": True,
            "below_actual_basis_independent_selector_promotion": True,
            "below_qw2191_quotient_safe_resolution": True,
            "no_false_pass": True,
        },
    }
    path = Path(__file__).resolve().parent / "generated" / "f136_first_source_topology_basis_independent_promotion_target_packet_summary.json"
    path.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
