#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_P216 = ROOT / "generated" / "p216_current_source_topology_selector_promotion_target_probe_summary.json"
OUT = ROOT / "generated" / "n236_current_first_source_topology_selector_promotion_target_theorem_summary.json"


def main() -> None:
    p216 = json.loads(IN_P216.read_text())

    discharged = p216["status"] == (
        "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_BELOW_BASIS_INDEPENDENCE_AND_QW2191_QUOTIENT_SAFETY_AFTER_P216"
    )

    summary = {
        "theorem_id": "N236",
        "status": (
            "N236_DISCHARGED_CURRENT_FIRST_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_THEOREM_NO_FALSE_PASS"
            if discharged
            else "N236_FAIL"
        ),
        "source_packet": p216["input_packet"],
        "promotion_target": p216["promotion_target"],
        "codomain_packet": p216["codomain_packet"],
        "observer_role": p216["observer_role"],
        "future_only": p216["future_only"],
        "basis_independence_discharged": p216["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": p216["qw2191_quotient_safe_discharged"],
        "current_selector_closure": p216["current_selector_closure"],
        "current_global_qw2191_discharge": p216["current_global_qw2191_discharge"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
