#!/usr/bin/env python3
import json
from pathlib import Path


def main() -> None:
    out = {
        "status": "F135_EXECUTED_FIRST_SOURCE_TOPOLOGY_FULL_NONTRIVIALITY_DISCHARGE_TARGET_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "discharge_target_name": "Theta_src_nontriv_discharge_target_v1",
        "domain_target": "Mu_src_nontriv_assembly_target_v1",
        "codomain_target": "actual_full_source_topology_nontriviality_discharge_target_v1",
        "scope": {
            "source_side_only": True,
            "future_only": True,
            "below_actual_full_nontriviality_discharge": True,
            "below_selector_promotion": True,
            "below_qw2191_resolution": True,
            "no_false_pass": True,
        },
    }
    path = Path(__file__).resolve().parent / "generated" / "f135_first_source_topology_full_nontriviality_discharge_target_packet_summary.json"
    path.write_text(json.dumps(out, indent=2) + "\n")


if __name__ == "__main__":
    main()
