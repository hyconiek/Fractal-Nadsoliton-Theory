#!/usr/bin/env python3
"""P1462 S4.12: export safety-band policy and validate candidate deltas against it."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
BAND = GEN / "p1460_s410_h2_delta_safety_band_summary.json"
PROPOSAL = GEN / "p1453_s44_h2_remediation_proposal.json"
POLICY = GEN / "p1462_s412_h2_safety_band_policy.json"
VALIDATION = GEN / "p1462_s412_h2_policy_validation_summary.json"
OBSTRUCTION = GEN / "p1462_s412_h2_policy_validation_obstruction.json"


EXTRA_TEST_DELTAS = [0.0009, 0.0010, 0.0015, 0.0020, 0.0021]


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def in_band(delta: float, dmin: float, dmax: float) -> bool:
    return dmin <= delta <= dmax


def main() -> None:
    band = readj(BAND)["h2_safety_band"]
    proposal = readj(PROPOSAL)

    dmin = float(band["delta_min"])
    dmax = float(band["delta_max"])

    policy = {
        "packet": "P1462",
        "policy_id": "h2_delta_safety_band_v1",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "delta_min": dmin,
        "delta_max": dmax,
        "source": "P1460",
    }
    POLICY.write_text(json.dumps(policy, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    candidate_deltas = [float(proposal["delta_margin_boost_h2"])] + EXTRA_TEST_DELTAS
    seen = set()
    ordered = []
    for d in candidate_deltas:
        r = round(d, 6)
        if r not in seen:
            ordered.append(r)
            seen.add(r)

    rows = []
    first_violation = None
    for delta in ordered:
        ok = in_band(delta, dmin, dmax)
        row = {
            "delta": delta,
            "in_policy_band": ok,
            "status": "PASS_POLICY_BAND" if ok else "FAIL_POLICY_BAND_VIOLATION",
        }
        rows.append(row)
        if (not ok) and first_violation is None:
            first_violation = row

    summary = {
        "packet": "P1462",
        "status": "PASS_POLICY_GUARD_ACTIVE",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_id": "h2_delta_safety_band_v1",
        "delta_min": dmin,
        "delta_max": dmax,
        "checked_candidates": rows,
        "violation_count": sum(1 for r in rows if not r["in_policy_band"]),
    }
    VALIDATION.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_violation is not None:
        obstruction = {
            "packet": "P1462",
            "status": "FAIL_POLICY_BAND_VIOLATION",
            "first_violation": first_violation,
            "rule": "outside-band rerun requires new policy update step",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1462] status=PASS_POLICY_GUARD_ACTIVE violations={summary['violation_count']}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print("[P1462] status=PASS_POLICY_GUARD_ACTIVE violations=0")


if __name__ == "__main__":
    main()
