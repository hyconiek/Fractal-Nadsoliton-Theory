#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
IN_F74 = GENERATED / "f74_first_preobserver_light_matter_asymmetric_provider_packet_summary.json"
IN_F127 = GENERATED / "f127_first_source_topology_invariant_candidate_packet_summary.json"
IN_F141 = GENERATED / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
IN_F143 = GENERATED / "f143_first_actual_source_topology_nonzero_flow_witness_packet_summary.json"
IN_N163 = GENERATED / "n163_current_observer_information_deficit_downstream_symptom_theorem_summary.json"
IN_N234 = GENERATED / "n234_current_global_selector_closure_and_qw2191_discharge_promotion_obstruction_theorem_summary.json"
OUT = GENERATED / "f142_first_actual_source_topology_observer_free_scope_witness_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    f74 = load_json(IN_F74)
    f127 = load_json(IN_F127)
    f141 = load_json(IN_F141)
    f143 = load_json(IN_F143)
    n163 = load_json(IN_N163)
    n234 = load_json(IN_N234)

    actual_observer_free_scope_witness_exported = (
        f127["properties"]["upstream_of_observer"] is True
        and f74["observer_nonparticipation"]["dP_NL_dO_zero"] is True
        and f74["observer_nonparticipation"]["dP_LM_dO_zero"] is True
        and f74["observer_nonparticipation"]["observer_to_upstream_blocks_zero"] is True
        and n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"] is True
        and n163["theorem_result"]["primary_missing_selector_source_gap_upstream_of_observer"] is True
        and n234["theorem_result"]["observer_downstream_only"] is True
    )

    support_packet = {
        "source_limit_tag_v1": f127["components"]["source_limit_tag_v1"],
        "phi_barrier_tag_v1": f127["components"]["phi_barrier_tag_v1"],
        "T_flow_0": f127["components"]["T_flow_0"],
        "ordering": f74["cascade"]["ordering"],
        "dP_NL_dO_zero": f74["observer_nonparticipation"]["dP_NL_dO_zero"],
        "dP_LM_dO_zero": f74["observer_nonparticipation"]["dP_LM_dO_zero"],
        "observer_to_upstream_blocks_zero": f74["observer_nonparticipation"]["observer_to_upstream_blocks_zero"],
        "observer_information_deficit_downstream_symptom_on_current_repo_state": n163["theorem_result"]["observer_information_deficit_downstream_symptom_on_current_repo_state"],
        "primary_missing_selector_source_gap_upstream_of_observer": n163["theorem_result"]["primary_missing_selector_source_gap_upstream_of_observer"],
        "observer_downstream_only": n234["theorem_result"]["observer_downstream_only"],
    }

    summary = {
        "packet_id": "F142",
        "status": "F142_EXECUTED_FIRST_ACTUAL_SOURCE_TOPOLOGY_OBSERVER_FREE_SCOPE_WITNESS_PACKET_NO_FALSE_PASS",
        "as_of": "2026-03-08",
        "input_packet": f127["candidate_packet"],
        "witness": "Omega_src_observer_free_scope_actual_witness_v1",
        "codomain": "observer_free_scope_tag_v1",
        "support_packet_id": "W_src_observer_free_scope_support_packet_v1",
        "support_packet": support_packet,
        "observer_role": "downstream_only",
        "actual_observer_free_scope_witness_exported": actual_observer_free_scope_witness_exported,
        "actual_observer_free_scope_discharged": actual_observer_free_scope_witness_exported,
        "barrier_protected_sign_discharged": f141["barrier_protected_sign_discharged"],
        "actual_nonzero_flow_discharged": f143["actual_nonzero_flow_discharged"],
        "full_source_topology_nontriviality_discharged": False,
        "basis_independence_discharged": False,
        "qw2191_quotient_safe_discharged": False,
        "current_selector_closure": False,
        "current_global_qw2191_discharge": False,
        "kernel_split_safe": True,
        "legacy_to_strict_bridge_claimed": False,
        "no_false_pass": True,
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
