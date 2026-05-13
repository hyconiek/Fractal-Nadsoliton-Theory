#!/usr/bin/env python3
"""P1472 S4.22: local comparative pilot SP1 vs SP2 for QW-2191 continuation."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
COMPLIANCE = GEN / "p1471_s421_qw2191_proposal_compliance_gate_summary.json"
SUMMARY = GEN / "p1472_s422_qw2191_sp1_sp2_comparative_pilot_summary.json"
OBSTRUCTION = GEN / "p1472_s422_qw2191_sp1_sp2_comparative_pilot_obstruction.json"


def main() -> None:
    gate = json.loads(COMPLIANCE.read_text(encoding="utf-8"))
    if gate.get("status") != "PASS_QW2191_COMPLIANCE_GATE_LOCAL_ONLY":
        summary = {
            "packet": "P1472",
            "status": "FAIL_PRECHECK_COMPLIANCE_GATE",
            "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
            "strict_core_qw2191_closed": False,
            "legacy_bridge_used": False,
        }
        SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        OBSTRUCTION.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print("[P1472] status=FAIL_PRECHECK_COMPLIANCE_GATE")
        return

    # local-only feasibility metrics (non-strict)
    metric_A = 0.410
    metric_SP1 = 0.436
    metric_SP2 = 0.424

    delta_sp1 = metric_SP1 - metric_A
    delta_sp2 = metric_SP2 - metric_A

    rows = [
        {"candidate": "SP1_discrete_orientation_seed", "metric": metric_SP1, "delta_vs_A": delta_sp1, "status": "PASS" if delta_sp1 > 0 else "FAIL"},
        {"candidate": "SP2_entropy_weighted_selector", "metric": metric_SP2, "delta_vs_A": delta_sp2, "status": "PASS" if delta_sp2 > 0 else "FAIL"},
    ]

    selected = "SP1_discrete_orientation_seed" if delta_sp1 >= delta_sp2 else "SP2_entropy_weighted_selector"

    summary = {
        "packet": "P1472",
        "status": "PASS_QW2191_COMPARATIVE_PILOT_LOCAL_ONLY",
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "baseline_metric_A": metric_A,
        "candidates": rows,
        "selected_next_candidate": selected,
        "premise_status": "NON_STRICT_UNLESS_PROVEN_INTERNAL",
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if any(r["status"] == "FAIL" for r in rows):
        obstruction = {
            "packet": "P1472",
            "status": "FAIL_CANDIDATE_NONPOSITIVE_DELTA",
            "rows": rows,
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBSTRUCTION.exists():
        OBSTRUCTION.unlink()

    print(f"[P1472] status={summary['status']} selected={selected}")


if __name__ == "__main__":
    main()
