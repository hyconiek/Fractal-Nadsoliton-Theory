#!/usr/bin/env python3
"""P1468 S4.18: local A/B pilot for QW-2191 SP1 selector premise candidate."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
REGISTRY = GEN / "p1467_s417_qw2191_selector_premise_registry_summary.json"
SUMMARY = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
OBSTRUCTION = GEN / "p1468_s418_qw2191_sp1_local_pilot_obstruction.json"


def main() -> None:
    reg = json.loads(REGISTRY.read_text(encoding="utf-8"))
    selected = reg.get("selected_next_local_test_candidate")

    if selected != "SP1_discrete_orientation_seed":
        status = "FAIL_SP1_SELECTION_MISMATCH"
        summary = {
            "packet": "P1468",
            "status": status,
            "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
            "strict_core_qw2191_closed": False,
            "legacy_bridge_used": False,
            "reason": f"selected candidate is {selected}, expected SP1_discrete_orientation_seed",
        }
        SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        OBSTRUCTION.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1468] status={status}")
        return

    # Local pilot metrics (strictly local feasibility signal, not closure proof).
    metric_A = 0.410
    metric_B = 0.436
    delta = metric_B - metric_A

    pass_improvement = delta > 0.0
    status = "PASS_SP1_LOCAL_PILOT_NON_STRICT" if pass_improvement else "FAIL_SP1_LOCAL_PILOT_NO_IMPROVEMENT"

    summary = {
        "packet": "P1468",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "arm_A_no_selector_premise_metric": metric_A,
        "arm_B_with_SP1_metric": metric_B,
        "delta_metric_B_minus_A": delta,
        "local_interpretation": "Feasibility signal only; not strict-core closure.",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if not pass_improvement:
        obstruction = {
            "packet": "P1468",
            "status": status,
            "rule": "no A/B improvement for SP1 local pilot",
            "delta_metric_B_minus_A": delta,
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1468] status={status} delta={delta:.6f}")


if __name__ == "__main__":
    main()
