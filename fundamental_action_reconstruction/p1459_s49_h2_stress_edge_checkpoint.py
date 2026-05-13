#!/usr/bin/env python3
"""P1459 S4.9: local stress-edge replay to find first failing h2 remediation delta."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
HOLDOUT = GEN / "p1452_s43_holdout_input.json"
PROP = GEN / "p1453_s44_h2_remediation_proposal.json"
SUMMARY = GEN / "p1459_s49_h2_stress_edge_summary.json"
OBSTRUCTION = GEN / "p1459_s49_h2_stress_edge_obstruction.json"

STEP = 0.0001
MAX_STEPS = 30


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def eval_h2(case: dict, min_gain: float, replay_tol: float, delta: float) -> dict:
    margin_before = float(case["margin_before"])
    margin_after = float(case["margin_after"]) + delta
    gain = margin_after - margin_before
    gain_ok = gain >= min_gain
    replay_gap = float(case["replay_gap"])
    replay_ok = replay_gap <= replay_tol
    return {
        "id": case["id"],
        "delta": round(delta, 6),
        "margin_before": margin_before,
        "margin_after": margin_after,
        "gain": gain,
        "min_gain": min_gain,
        "replay_gap": replay_gap,
        "replay_tol": replay_tol,
        "gain_ok": gain_ok,
        "replay_ok": replay_ok,
        "status": "PASS" if gain_ok and replay_ok else "FAIL",
    }


def main() -> None:
    holdout = readj(HOLDOUT)
    prop = readj(PROP)

    min_gain = float(holdout["min_gain"])
    replay_tol = float(holdout["replay_tol"])
    base_delta = float(prop["delta_margin_boost_h2"])
    h2 = next(c for c in holdout["cases"] if c["id"] == "h2")

    trajectory = []
    first_fail = None
    last_pass = None
    for k in range(MAX_STEPS + 1):
        delta = max(0.0, base_delta - STEP * k)
        row = eval_h2(h2, min_gain, replay_tol, delta)
        trajectory.append(row)
        if row["status"] == "PASS":
            last_pass = row
        elif first_fail is None:
            first_fail = row
            break
        if delta == 0.0:
            break

    status = "PASS_STRESS_EDGE_FOUND" if first_fail is not None else "WARN_NO_FAIL_IN_SCAN"

    summary = {
        "packet": "P1459",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "base_delta_margin_boost_h2": base_delta,
        "step": STEP,
        "max_steps": MAX_STEPS,
        "last_pass": last_pass,
        "first_fail": first_fail,
        "trajectory": trajectory,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1459",
            "status": "FAIL_H2_STRESS_EDGE_THRESHOLD",
            "first_fail": first_fail,
            "last_pass": last_pass,
            "rule": "immediate obstruction export at first fail threshold",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1459] status=PASS_STRESS_EDGE_FOUND first_fail_delta={first_fail['delta']}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print("[P1459] status=WARN_NO_FAIL_IN_SCAN")


if __name__ == "__main__":
    main()
