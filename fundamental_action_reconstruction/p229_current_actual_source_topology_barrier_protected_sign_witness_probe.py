#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F141 = ROOT / "generated" / "f141_first_actual_source_topology_barrier_protected_sign_witness_packet_summary.json"
OUT = ROOT / "generated" / "p229_current_actual_source_topology_barrier_protected_sign_witness_probe_summary.json"


def main() -> None:
    f141 = json.loads(IN_F141.read_text(encoding="utf-8"))

    passed = (
        f141["input_packet"] == "tau_src_candidate_v1"
        and f141["witness"] == "Psi_src_barrier_sign_actual_witness_v1"
        and f141["codomain"] == "barrier_protected_sign_class_v1"
        and f141["actual_barrier_protected_sign_witness_exported"] is True
        and f141["barrier_protected_sign_discharged"] is True
        and f141["observer_role"] == "downstream_only"
        and f141["full_source_topology_nontriviality_discharged"] is False
        and f141["basis_independence_discharged"] is False
        and f141["qw2191_quotient_safe_discharged"] is False
        and f141["current_selector_closure"] is False
        and f141["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P229",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_ACTUAL_SOURCE_TOPOLOGY_BARRIER_PROTECTED_SIGN_WITNESS_BELOW_FULL_SOURCE_TOPOLOGY_NONTRIVIALITY_AFTER_P229"
            if passed
            else "P229_FAIL"
        ),
        "input_packet": f141["input_packet"],
        "witness": f141["witness"],
        "codomain": f141["codomain"],
        "support_packet_id": f141["support_packet_id"],
        "observer_role": f141["observer_role"],
        "actual_barrier_protected_sign_witness_exported": f141["actual_barrier_protected_sign_witness_exported"],
        "barrier_protected_sign_discharged": f141["barrier_protected_sign_discharged"],
        "full_source_topology_nontriviality_discharged": f141["full_source_topology_nontriviality_discharged"],
        "basis_independence_discharged": f141["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f141["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f141["current_selector_closure"],
        "current_global_qw2191_discharge": f141["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()
