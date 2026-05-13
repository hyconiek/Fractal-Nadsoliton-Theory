#!/usr/bin/env python3
"""P1478 S4.28: conservative SP1 operating policy update from transfer-edge evidence."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
P1477 = GEN / "p1477_s427_qw2191_sp1_transfer_edge_scan_summary.json"
POLICY = GEN / "p1478_s428_qw2191_sp1_operating_policy_v2.json"
SUMMARY = GEN / "p1478_s428_qw2191_sp1_policy_update_summary.json"
OBSTRUCTION = GEN / "p1478_s428_qw2191_sp1_policy_update_obstruction.json"

MARGIN = 0.002


def main() -> None:
    s1477 = json.loads(P1477.read_text(encoding="utf-8"))
    worst = s1477["worst_case_edge"]
    edge_worst = float(worst["first_fail"]["shift"])

    safe_min_shift = edge_worst + MARGIN
    safe_max_shift = 0.0

    policy = {
        "packet": "P1478",
        "policy_id": "sp1_operating_policy_v2",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "premise_id": "SP1_discrete_orientation_seed",
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
        "source_packet": "P1477",
        "worst_case_edge_shift": edge_worst,
        "safety_margin": MARGIN,
        "safe_shift_min": safe_min_shift,
        "safe_shift_max": safe_max_shift,
    }
    POLICY.write_text(json.dumps(policy, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    # quick sanity check candidates
    candidates = [-0.028, -0.024, -0.023, -0.022, -0.01, 0.0]
    checks = []
    first_block = None
    for sh in candidates:
        ok = safe_min_shift <= sh <= safe_max_shift
        row = {"shift": sh, "in_policy": ok, "status": "PASS" if ok else "BLOCK"}
        checks.append(row)
        if (not ok) and first_block is None:
            first_block = row

    summary = {
        "packet": "P1478",
        "status": "PASS_SP1_POLICY_UPDATE_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "policy_id": policy["policy_id"],
        "checks": checks,
        "blocked_count": sum(1 for c in checks if c["status"] == "BLOCK"),
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_block is not None:
        obstruction = {
            "packet": "P1478",
            "status": "FAIL_SP1_POLICY_BAND_VIOLATION",
            "first_block": first_block,
            "rule": "candidate outside conservative SP1 policy window",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1478] status={summary['status']} safe_min_shift={safe_min_shift}")


if __name__ == "__main__":
    main()
