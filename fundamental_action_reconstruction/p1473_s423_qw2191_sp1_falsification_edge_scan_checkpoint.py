#!/usr/bin/env python3
"""P1473 S4.23: find falsification edge for SP1 local signal under negative shifts."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
SUMMARY = GEN / "p1473_s423_qw2191_sp1_falsification_edge_scan_summary.json"
OBSTRUCTION = GEN / "p1473_s423_qw2191_sp1_falsification_edge_scan_obstruction.json"

STEP = 0.002
MAX_STEPS = 20


def main() -> None:
    base = json.loads(P1468.read_text(encoding="utf-8"))
    metric_a = float(base["arm_A_no_selector_premise_metric"])
    metric_b = float(base["arm_B_with_SP1_metric"])

    rows = []
    last_pass = None
    first_fail = None
    for k in range(MAX_STEPS + 1):
        shift = -STEP * k
        b_shifted = metric_b + shift
        delta = b_shifted - metric_a
        ok = delta > 0.0
        row = {
            "k": k,
            "shift": shift,
            "metric_A": metric_a,
            "metric_B_shifted": b_shifted,
            "delta_metric_B_minus_A": delta,
            "status": "PASS" if ok else "FAIL",
        }
        rows.append(row)
        if ok:
            last_pass = row
        elif first_fail is None:
            first_fail = row
            break

    status = "PASS_SP1_FALSIFICATION_EDGE_FOUND_NON_STRICT" if first_fail is not None else "WARN_NO_EDGE_FOUND_IN_SCAN"
    summary = {
        "packet": "P1473",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "step": STEP,
        "max_steps": MAX_STEPS,
        "last_pass": last_pass,
        "first_fail": first_fail,
        "rows": rows,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1473",
            "status": "FAIL_SP1_EDGE_THRESHOLD",
            "first_fail": first_fail,
            "last_pass": last_pass,
            "rule": "SP1 local signal loses positive delta beyond this edge",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1473] status={status} first_fail_k={first_fail['k'] if first_fail else 'NA'}")


if __name__ == "__main__":
    main()
