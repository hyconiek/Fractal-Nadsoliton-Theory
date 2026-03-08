#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F129 = ROOT / "generated" / "f129_first_source_topology_invariant_nontriviality_target_packet_summary.json"
OUT = ROOT / "generated" / "p217_current_source_topology_invariant_nontriviality_target_probe_summary.json"


def main() -> None:
    f129 = json.loads(IN_F129.read_text())

    passed = (
        f129["input_packet"] == "tau_src_candidate_v1"
        and f129["nontriviality_target"] == "Nu_src_nontriv_target_v1"
        and f129["codomain_packet"]["name"] == "Lambda_src_nontriv_target_v1"
        and f129["observer_role"] == "downstream_only"
        and f129["future_only"] is True
        and f129["actual_nontriviality_discharged"] is False
        and f129["basis_independence_discharged"] is False
        and f129["qw2191_quotient_safe_discharged"] is False
        and f129["current_selector_closure"] is False
        and f129["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P217",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_INVARIANT_NONTRIVIALITY_TARGET_BELOW_ACTUAL_NONTRIVIALITY_DISCHARGE_AFTER_P217"
            if passed
            else "P217_FAIL"
        ),
        "input_packet": f129["input_packet"],
        "nontriviality_target": f129["nontriviality_target"],
        "codomain_packet": f129["codomain_packet"]["name"],
        "observer_role": f129["observer_role"],
        "future_only": f129["future_only"],
        "actual_nontriviality_discharged": f129["actual_nontriviality_discharged"],
        "basis_independence_discharged": f129["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f129["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f129["current_selector_closure"],
        "current_global_qw2191_discharge": f129["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
