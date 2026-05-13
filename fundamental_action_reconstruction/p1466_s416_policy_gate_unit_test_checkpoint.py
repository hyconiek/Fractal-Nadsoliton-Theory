#!/usr/bin/env python3
"""P1466 S4.16: point-based unit tests for shared policy gate semantics."""

from __future__ import annotations

import json
from pathlib import Path

from p1465_policy_gate_core import PolicyBand, gate_decision

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
SUMMARY = GEN / "p1466_s416_policy_gate_unit_test_summary.json"
OBSTRUCTION = GEN / "p1466_s416_policy_gate_unit_test_obstruction.json"


TESTS = [
    ("below_band", 0.0009, False),
    ("edge_min", 0.0010, True),
    ("inside", 0.0015, True),
    ("edge_max", 0.0020, True),
    ("above_band", 0.0021, False),
]


def main() -> None:
    band = PolicyBand(delta_min=0.001, delta_max=0.002)
    rows = []
    first_fail = None

    for label, delta, expected_allow in TESTS:
        decision = gate_decision(delta, band)
        allow = bool(decision["allow_rerun"])
        pass_case = allow == expected_allow
        row = {
            "label": label,
            "delta": delta,
            "expected_allow": expected_allow,
            "observed_allow": allow,
            "gate_status": decision["gate_status"],
            "status": "PASS" if pass_case else "FAIL",
        }
        rows.append(row)
        if (not pass_case) and first_fail is None:
            first_fail = row

    ok = first_fail is None
    summary = {
        "packet": "P1466",
        "status": "PASS_POLICY_GATE_UNIT_TEST_LOCAL_ONLY" if ok else "FAIL_POLICY_GATE_UNIT_TEST",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "band": {"delta_min": band.delta_min, "delta_max": band.delta_max},
        "cases": rows,
        "fail_count": sum(1 for r in rows if r["status"] == "FAIL"),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if not ok:
        obstruction = {
            "packet": "P1466",
            "status": "FAIL_POLICY_GATE_UNIT_TEST",
            "first_fail": first_fail,
            "rule": "shared policy gate semantics mismatch",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1466] status={summary['status']} fail_count={summary['fail_count']}")


if __name__ == "__main__":
    main()
