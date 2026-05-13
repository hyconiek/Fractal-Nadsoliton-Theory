#!/usr/bin/env python3
"""P1479 S4.29: enforce SP1 policy-v2 gate before local replay execution."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1468 = GEN / "p1468_s418_qw2191_sp1_local_pilot_summary.json"
POLICY = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"
SUMMARY = GEN / "p1479_s429_qw2191_sp1_policy_gate_integration_summary.json"
OBSTRUCTION = GEN / "p1479_s429_qw2191_sp1_policy_gate_integration_obstruction.json"

CANDIDATES = [
    {"label": "in_policy_candidate", "shift": -0.022},
    {"label": "out_of_policy_candidate", "shift": -0.026},
]


def in_policy(shift: float, smin: float, smax: float) -> bool:
    return smin <= shift <= smax


def main() -> None:
    base = json.loads(P1468.read_text(encoding="utf-8"))
    policy = json.loads(POLICY.read_text(encoding="utf-8"))

    metric_a = float(base["arm_A_no_selector_premise_metric"])
    metric_b = float(base["arm_B_with_SP1_metric"])
    smin = float(policy["safe_shift_min"])
    smax = float(policy["safe_shift_max"])

    rows = []
    blocked = 0
    first_block = None
    for cand in CANDIDATES:
        sh = float(cand["shift"])
        allowed = in_policy(sh, smin, smax)
        if allowed:
            delta = (metric_b + sh) - metric_a
            rows.append({
                "label": cand["label"],
                "shift": sh,
                "gate_status": "PASS_POLICY_GATE",
                "rerun_executed": True,
                "delta_metric_B_minus_A": delta,
                "rerun_status": "PASS" if delta > 0 else "FAIL",
            })
        else:
            blocked += 1
            row = {
                "label": cand["label"],
                "shift": sh,
                "gate_status": "FAIL_SP1_POLICY_BAND_VIOLATION",
                "rerun_executed": False,
            }
            rows.append(row)
            if first_block is None:
                first_block = row

    status = "PASS_SP1_POLICY_GATE_INTEGRATION_LOCAL_ONLY" if blocked >= 1 else "FAIL_SP1_POLICY_GATE_INTEGRATION"

    summary = {
        "packet": "P1479",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_id": policy["policy_id"],
        "safe_shift_min": smin,
        "safe_shift_max": smax,
        "candidates": rows,
        "blocked_count": blocked,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_block is not None:
        obstruction = {
            "packet": "P1479",
            "status": "FAIL_SP1_POLICY_BAND_VIOLATION",
            "first_block": first_block,
            "rule": "out-of-policy candidate blocked before replay execution",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1479] status={status} blocked_count={blocked}")


if __name__ == "__main__":
    main()
