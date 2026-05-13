#!/usr/bin/env python3
"""P1463 S4.13: apply policy-band gate before running local rerun candidates."""

from __future__ import annotations

import json
from pathlib import Path

from p1465_policy_gate_core import PolicyBand, gate_decision

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
POLICY = GEN / "p1462_s412_h2_safety_band_policy.json"
PROPOSAL = GEN / "p1453_s44_h2_remediation_proposal.json"
SUMMARY = GEN / "p1463_s413_policy_aware_rerun_gate_summary.json"
OBSTRUCTION = GEN / "p1463_s413_policy_aware_rerun_gate_obstruction.json"

# candidate queue emulates future rerun requests
QUEUE_DELTAS = [0.0015, 0.0009]


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    policy = readj(POLICY)
    proposal = readj(PROPOSAL)

    dmin = float(policy["delta_min"])
    dmax = float(policy["delta_max"])

    candidate_deltas = [float(proposal["delta_margin_boost_h2"])] + QUEUE_DELTAS
    seen = set()
    ordered = []
    for d in candidate_deltas:
        r = round(d, 6)
        if r not in seen:
            seen.add(r)
            ordered.append(r)

    decisions = []
    first_block = None
    for delta in ordered:
        decision = gate_decision(delta, PolicyBand(dmin, dmax))
        allow = bool(decision["allow_rerun"])
        decisions.append(decision)
        if (not allow) and first_block is None:
            first_block = decision

    summary = {
        "packet": "P1463",
        "status": "PASS_POLICY_GATE_ACTIVE",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_id": policy["policy_id"],
        "band": {"delta_min": dmin, "delta_max": dmax},
        "decisions": decisions,
        "blocked_count": sum(1 for d in decisions if not d["allow_rerun"]),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_block is not None:
        obstruction = {
            "packet": "P1463",
            "status": "FAIL_POLICY_BAND_VIOLATION",
            "first_blocked_candidate": first_block,
            "rule": "rerun blocked before execution when candidate delta is outside policy band",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1463] status=PASS_POLICY_GATE_ACTIVE blocked={summary['blocked_count']}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print("[P1463] status=PASS_POLICY_GATE_ACTIVE blocked=0")


if __name__ == "__main__":
    main()
