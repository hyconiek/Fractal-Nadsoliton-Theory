#!/usr/bin/env python3
"""P1474 S4.24: export edge-aware operating window for SP1 and verify near-edge points."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
P1473 = GEN / "p1473_s423_qw2191_sp1_falsification_edge_scan_summary.json"
SUMMARY = GEN / "p1474_s424_qw2191_sp1_edge_aware_window_summary.json"


def main() -> None:
    s1468 = json.loads(P1468.read_text(encoding="utf-8"))
    s1473 = json.loads(P1473.read_text(encoding="utf-8"))

    metric_a = float(s1468["arm_A_no_selector_premise_metric"])
    metric_b = float(s1468["arm_B_with_SP1_metric"])
    fail_shift = float(s1473["first_fail"]["shift"])  # negative value
    pass_shift = float(s1473["last_pass"]["shift"])  # negative value, closer to zero

    # conservative window: stay 0.002 above fail threshold
    window_min_shift = pass_shift
    window_max_shift = 0.0

    check_points = [window_min_shift, -0.012, -0.006, 0.0]
    checks = []
    for sh in check_points:
        delta = (metric_b + sh) - metric_a
        checks.append({
            "shift": sh,
            "delta_metric_B_minus_A": delta,
            "status": "PASS" if delta > 0 else "FAIL",
        })

    fail_count = sum(1 for c in checks if c["status"] == "FAIL")
    status = "PASS_SP1_EDGE_AWARE_WINDOW_LOCAL_ONLY" if fail_count == 0 else "FAIL_SP1_EDGE_AWARE_WINDOW"

    summary = {
        "packet": "P1474",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "falsification_edge_shift": fail_shift,
        "last_pass_shift": pass_shift,
        "operating_window_shift": {
            "min_shift": window_min_shift,
            "max_shift": window_max_shift,
        },
        "checks": checks,
        "fail_count": fail_count,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1474] status={status} window=[{window_min_shift},{window_max_shift}] fail_count={fail_count}")


if __name__ == "__main__":
    main()
