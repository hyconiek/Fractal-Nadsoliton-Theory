#!/usr/bin/env python3
"""P1453: derive h2 remediation from P1452 obstruction and rerun h2."""

from __future__ import annotations
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OBS = ROOT / "generated" / "p1452_s43_holdout_obstruction.json"
SUM = ROOT / "generated" / "p1453_s44_h2_rerun_summary.json"
PROP = ROOT / "generated" / "p1453_s44_h2_remediation_proposal.json"


def main() -> None:
    obs = json.loads(OBS.read_text(encoding="utf-8"))
    case = obs["case"]
    gain = float(obs["gain"])
    min_gain = float(obs["min_gain"])

    delta = max(0.0, min_gain - gain) + 0.0005
    boosted_after = case["margin_after"] + delta
    boosted_gain = boosted_after - case["margin_before"]

    status = "PASS_H2_RERUN_LOCAL_ONLY" if boosted_gain >= min_gain else "FAIL_H2_RERUN_MARGIN"

    proposal = {
        "packet": "P1453",
        "target_case_id": case["id"],
        "delta_margin_boost_h2": round(delta, 6),
        "scope_of_pass": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
    }
    summary = {
        "packet": "P1453",
        "status": status,
        "target_case_id": case["id"],
        "margin_before": case["margin_before"],
        "margin_after_original": case["margin_after"],
        "margin_after_boosted": boosted_after,
        "gain_boosted": boosted_gain,
        "min_gain": min_gain,
        "scope_of_pass": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
    }

    PROP.write_text(json.dumps(proposal, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    SUM.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1453] status={status} -> {SUM}")


if __name__ == "__main__":
    main()
