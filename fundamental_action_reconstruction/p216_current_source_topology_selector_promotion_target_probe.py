#!/usr/bin/env python3

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_F128 = ROOT / "generated" / "f128_first_source_topology_selector_promotion_target_packet_summary.json"
OUT = ROOT / "generated" / "p216_current_source_topology_selector_promotion_target_probe_summary.json"


def main() -> None:
    f128 = json.loads(IN_F128.read_text())

    passed = (
        f128["promotion_target"] == "Pi_sel_src_target_v1"
        and f128["input_packet"] == "tau_src_candidate_v1"
        and f128["observer_role"] == "downstream_only"
        and f128["basis_independence_discharged"] is False
        and f128["qw2191_quotient_safe_discharged"] is False
        and f128["current_selector_closure"] is False
        and f128["current_global_qw2191_discharge"] is False
    )

    summary = {
        "probe_id": "P216",
        "status": (
            "CURRENT_REPO_EXPORTS_ONE_FUTURE_SOURCE_TOPOLOGY_SELECTOR_PROMOTION_TARGET_BELOW_BASIS_INDEPENDENCE_AND_QW2191_QUOTIENT_SAFETY_AFTER_P216"
            if passed
            else "P216_FAIL"
        ),
        "input_packet": f128["input_packet"],
        "promotion_target": f128["promotion_target"],
        "codomain_packet": f128["codomain_packet"]["name"],
        "observer_role": f128["observer_role"],
        "future_only": True,
        "basis_independence_discharged": f128["basis_independence_discharged"],
        "qw2191_quotient_safe_discharged": f128["qw2191_quotient_safe_discharged"],
        "current_selector_closure": f128["current_selector_closure"],
        "current_global_qw2191_discharge": f128["current_global_qw2191_discharge"],
        "downstream_chart_realization_candidates": f128["downstream_chart_realization_candidates"],
    }
    OUT.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
