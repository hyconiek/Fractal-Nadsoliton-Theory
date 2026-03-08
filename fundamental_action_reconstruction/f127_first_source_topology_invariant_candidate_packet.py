#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f74 = load_json(
        GENERATED / "f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
    )
    n163 = load_json(
        GENERATED / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
    )
    n234 = load_json(
        GENERATED
        / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
    )

    summary = {
        "stage": "F127",
        "lane": "future_source_topology_selector_route_only",
        "status": "F127_EXECUTED_FIRST_SOURCE_TOPOLOGY_INVARIANT_CANDIDATE_PACKET_FUTURE_ROUTE_ONLY_NO_FALSE_PASS",
        "candidate_packet": "tau_src_candidate_v1",
        "components": {
            "source_limit_tag_v1": "d -> 0",
            "phi_barrier_tag_v1": "fixed_nonzero_phi_kernel_core_barrier",
            "T_flow_0": f74["topological_flow_vector"]["definition"],
        },
        "anchors": {
            "kernel_core_amplitude": f74["topological_flow_vector"]["kernel_core_amplitude"],
            "observer_downstream_only": n163["theorem_result"][
                "observer_information_deficit_downstream_symptom_on_current_repo_state"
            ],
            "global_selector_closure_justified_on_current_repo_state": n234[
                "theorem_result"
            ]["global_selector_closure_justified_on_current_repo_state"],
            "global_qw2191_discharge_justified_on_current_repo_state": n234[
                "theorem_result"
            ]["global_qw2191_discharge_justified_on_current_repo_state"],
        },
        "properties": {
            "future_route_only": True,
            "upstream_of_observer": True,
            "external_selector_import": False,
            "basis_independent_selector_promotion_exported": False,
            "qw2191_quotient_safe_promotion_exported": False,
            "current_selector_closure": False,
            "current_qw2191_discharge": False,
            "kernel_split_safe": True,
        },
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
