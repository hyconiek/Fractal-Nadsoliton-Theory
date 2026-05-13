#!/usr/bin/env python3
"""P1476 S4.26: transfer-check SP1 window under slight baseline-context variation."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
P1474 = GEN / "p1474_s424_qw2191_sp1_edge_aware_window_summary.json"
SUMMARY = GEN / "p1476_s426_qw2191_sp1_window_transfer_check_summary.json"
OBSTRUCTION = GEN / "p1476_s426_qw2191_sp1_window_transfer_check_obstruction.json"


def main() -> None:
    s1468 = json.loads(P1468.read_text(encoding="utf-8"))
    s1474 = json.loads(P1474.read_text(encoding="utf-8"))

    metric_a = float(s1468["arm_A_no_selector_premise_metric"])
    metric_b = float(s1468["arm_B_with_SP1_metric"])
    min_shift = float(s1474["operating_window_shift"]["min_shift"])

    # simulate slight context change in baseline A
    context_offsets = [0.0, 0.001, -0.001]
    test_shifts = [min_shift, -0.018, -0.012, -0.006, 0.0]

    rows = []
    first_fail = None
    for off in context_offsets:
        a_ctx = metric_a + off
        for sh in test_shifts:
            delta = (metric_b + sh) - a_ctx
            ok = delta > 0.0
            row = {
                "context_offset_A": off,
                "shift": sh,
                "delta_metric_B_minus_Actx": delta,
                "status": "PASS" if ok else "FAIL",
            }
            rows.append(row)
            if (not ok) and first_fail is None:
                first_fail = row

    fail_count = sum(1 for r in rows if r["status"] == "FAIL")
    status = "PASS_SP1_WINDOW_TRANSFER_LOCAL_ONLY" if fail_count == 0 else "FAIL_SP1_WINDOW_TRANSFER"

    summary = {
        "packet": "P1476",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "context_offsets_A": context_offsets,
        "tested_shifts": test_shifts,
        "rows": rows,
        "fail_count": fail_count,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1476",
            "status": "FAIL_SP1_WINDOW_TRANSFER",
            "first_fail": first_fail,
            "rule": "window does not transfer under mild baseline-context change",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1476] status={status} fail_count={fail_count}")


if __name__ == "__main__":
    main()
