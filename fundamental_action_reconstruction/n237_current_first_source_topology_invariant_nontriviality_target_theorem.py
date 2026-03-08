#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P217 = ROOT / "generated" / "p217_current_source_topology_invariant_nontriviality_target_probe_summary.json"
OUT = ROOT / "generated" / "n237_current_first_source_topology_invariant_nontriviality_target_theorem_summary.json"


def main() -> None:
    p217 = json.loads(IN_P217.read_text())

    passed = (
        p217["status"]
        == "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_INVARIANT_NONTRIVIALITY_TARGET_BELOW_ACTUAL_NONTRIVIALITY_DISCHARGE_AFTER_P217"
        and p217["future_only"] is True
        and p217["actual_nontriviality_discharged"] is False
        and p217["basis_independence_discharged"] is False
        and p217["qw2191_quotient_safe_discharged"] is False
        and p217["current_selector_closure"] is False
        and p217["current_global_qw2191_discharge"] is False
    )

    summary = {
        "theorem_id": "N237",
        "status": (
            "N237_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_INVARIANT_NONTRIVIALITY_TARGET_THEOREM_NO_FALSE_PASS"
            if passed
            else "N237_FAIL"
        ),
        "input_packet": p217["input_packet"],
        "nontriviality_target": p217["nontriviality_target"],
        "codomain_packet": p217["codomain_packet"],
        "observer_role": p217["observer_role"],
        "future_only": p217["future_only"],
        "actual_nontriviality_discharged": p217["actual_nontriviality_discharged"],
        "basis_independence_discharged": p217["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p217["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p217["current_selector_closure"],
        "current_global_qw2191_discharge": p217["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
