#!/usr/bin/env python3
"""P1475 S4.25: dense replay near SP1 lower operating-window boundary."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
P1474 = GEN / "p1474_s424_qw2191_sp1_edge_aware_window_summary.json"
SUMMARY = GEN / "p1475_s425_qw2191_sp1_window_stress_replay_summary.json"
OBSTRUCTION = GEN / "p1475_s425_qw2191_sp1_window_stress_replay_obstruction.json"


def main() -> None:
    s1468 = json.loads(P1468.read_text(encoding="utf-8"))
    s1474 = json.loads(P1474.read_text(encoding="utf-8"))

    metric_a = float(s1468["arm_A_no_selector_premise_metric"])
    metric_b = float(s1468["arm_B_with_SP1_metric"])
    min_shift = float(s1474["operating_window_shift"]["min_shift"])

    # Dense neighborhood around lower boundary
    shifts = [min_shift, min_shift + 0.001, min_shift + 0.002, min_shift + 0.003, min_shift + 0.004]

    rows = []
    first_fail = None
    for sh in shifts:
        delta = (metric_b + sh) - metric_a
        ok = delta > 0.0
        row = {"shift": sh, "delta_metric_B_minus_A": delta, "status": "PASS" if ok else "FAIL"}
        rows.append(row)
        if (not ok) and first_fail is None:
            first_fail = row

    fail_count = sum(1 for r in rows if r["status"] == "FAIL")
    status = "PASS_SP1_WINDOW_STRESS_LOCAL_ONLY" if fail_count == 0 else "FAIL_SP1_WINDOW_STRESS"

    summary = {
        "packet": "P1475",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "tested_shifts": shifts,
        "rows": rows,
        "fail_count": fail_count,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1475",
            "status": "FAIL_SP1_WINDOW_STRESS",
            "first_fail": first_fail,
            "rule": "window boundary not numerically stable in dense replay",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1475] status={status} fail_count={fail_count}")


if __name__ == "__main__":
    main()
