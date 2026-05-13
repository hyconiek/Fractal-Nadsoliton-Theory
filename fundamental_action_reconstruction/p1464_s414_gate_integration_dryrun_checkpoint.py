#!/usr/bin/env python3
"""P1464 S4.14: dry-run integration showing policy gate blocks out-of-band before replay."""

from __future__ import annotations

import json
from pathlib import Path

from p1465_policy_gate_core import PolicyBand, gate_decision

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
POLICY = GEN / "p1462_s412_h2_safety_band_policy.json"
HOLDOUT = GEN / "p1452_s43_holdout_input.json"
SUMMARY = GEN / "p1464_s414_gate_integration_dryrun_summary.json"
OBSTRUCTION = GEN / "p1464_s414_gate_integration_dryrun_obstruction.json"

CANDIDATES = [
    {"label": "in_band_candidate", "delta": 0.0015},
    {"label": "out_of_band_candidate", "delta": 0.0009},
]


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def replay_h2_status(h2_case: dict, delta: float, min_gain: float, replay_tol: float) -> dict:
    before = float(h2_case["margin_before"])
    after = float(h2_case["margin_after"]) + delta
    gain = after - before
    gain_ok = gain >= min_gain
    replay_gap = float(h2_case["replay_gap"])
    replay_ok = replay_gap <= replay_tol
    return {
        "gain": gain,
        "gain_ok": gain_ok,
        "replay_gap": replay_gap,
        "replay_ok": replay_ok,
        "status": "PASS" if gain_ok and replay_ok else "FAIL",
    }


def main() -> None:
    policy = readj(POLICY)
    holdout = readj(HOLDOUT)

    dmin = float(policy["delta_min"])
    dmax = float(policy["delta_max"])
    min_gain = float(holdout["min_gain"])
    replay_tol = float(holdout["replay_tol"])
    h2 = next(c for c in holdout["cases"] if c["id"] == "h2")

    rows = []
    blocked = 0
    for cand in CANDIDATES:
        delta = float(cand["delta"])
        decision = gate_decision(delta, PolicyBand(dmin, dmax))
        allowed = bool(decision["allow_rerun"])
        if allowed:
            replay = replay_h2_status(h2, delta, min_gain, replay_tol)
            rows.append(
                {
                    "label": cand["label"],
                    "delta": decision["delta"],
                    "gate_status": decision["gate_status"],
                    "rerun_executed": True,
                    "rerun_result": replay,
                }
            )
        else:
            blocked += 1
            rows.append(
                {
                    "label": cand["label"],
                    "delta": decision["delta"],
                    "gate_status": decision["gate_status"],
                    "rerun_executed": False,
                    "rerun_result": None,
                }
            )

    status = "PASS_GATE_INTEGRATION_DRYRUN" if blocked >= 1 else "FAIL_GATE_INTEGRATION_DRYRUN"
    summary = {
        "packet": "P1464",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_band": {"delta_min": dmin, "delta_max": dmax},
        "candidates": rows,
        "blocked_count": blocked,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if blocked >= 1:
        obstruction = {
            "packet": "P1464",
            "status": "FAIL_POLICY_BAND_VIOLATION",
            "first_blocked": next(r for r in rows if not r["rerun_executed"]),
            "rule": "out-of-band candidate blocked prior to replay execution",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1464] status={status} blocked={blocked}")


if __name__ == "__main__":
    main()
