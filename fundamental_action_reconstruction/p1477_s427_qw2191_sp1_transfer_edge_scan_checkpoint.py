#!/usr/bin/env python3
"""P1477 S4.27: transfer-edge scan for SP1 across mild baseline-context variants."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
SUMMARY = GEN / "p1477_s427_qw2191_sp1_transfer_edge_scan_summary.json"
OBSTRUCTION = GEN / "p1477_s427_qw2191_sp1_transfer_edge_scan_obstruction.json"

STEP = 0.002
MAX_STEPS = 30
CONTEXT_OFFSETS = [0.0, 0.001, -0.001]


def main() -> None:
    s1468 = json.loads(P1468.read_text(encoding="utf-8"))
    metric_a = float(s1468["arm_A_no_selector_premise_metric"])
    metric_b = float(s1468["arm_B_with_SP1_metric"])

    scans = []
    worst = None

    for off in CONTEXT_OFFSETS:
        a_ctx = metric_a + off
        last_pass = None
        first_fail = None
        rows = []
        for k in range(MAX_STEPS + 1):
            shift = -STEP * k
            b_shifted = metric_b + shift
            delta = b_shifted - a_ctx
            ok = delta > 0.0
            row = {"k": k, "shift": shift, "delta_metric_B_minus_Actx": delta, "status": "PASS" if ok else "FAIL"}
            rows.append(row)
            if ok:
                last_pass = row
            elif first_fail is None:
                first_fail = row
                break

        scan = {
            "context_offset_A": off,
            "last_pass": last_pass,
            "first_fail": first_fail,
            "rows": rows,
        }
        scans.append(scan)

        if first_fail is not None:
            score = float(first_fail["shift"])  # more negative => more robust; less negative => worse edge
            if worst is None or score > float(worst["first_fail"]["shift"]):
                worst = scan

    status = "PASS_SP1_TRANSFER_EDGE_FOUND_LOCAL_ONLY" if worst is not None else "WARN_NO_TRANSFER_EDGE_FOUND"

    summary = {
        "packet": "P1477",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "step": STEP,
        "context_offsets_A": CONTEXT_OFFSETS,
        "scans": scans,
        "worst_case_edge": worst,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if worst is not None:
        obstruction = {
            "packet": "P1477",
            "status": "FAIL_SP1_TRANSFER_EDGE_THRESHOLD",
            "worst_case_edge": worst,
            "rule": "report most fragile context edge for conservative operating policy",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1477] status={status} worst_edge_shift={worst['first_fail']['shift'] if worst else 'NA'}")


if __name__ == "__main__":
    main()
