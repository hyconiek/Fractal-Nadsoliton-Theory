#!/usr/bin/env python3
"""P1451 checkpoint: targeted rerun on first obstruction using P1450 proposal."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OBSTRUCTION = ROOT / "generated" / "p1449_s4_first_obstruction.json"
INPUT = ROOT / "generated" / "p1449_s4_stress_scenarios_input.json"
PROPOSAL = ROOT / "generated" / "p1450_s41_noncyclic_remediation_proposal.json"
OUT = ROOT / "generated" / "p1451_s42_targeted_rerun_summary.json"


def main() -> None:
    obs = json.loads(OBSTRUCTION.read_text(encoding="utf-8"))
    base = json.loads(INPUT.read_text(encoding="utf-8"))
    proposal = json.loads(PROPOSAL.read_text(encoding="utf-8"))

    target_id = obs["scenario"]["id"]
    replay_tol = float(base["replay_tol"])
    min_gain = float(base["min_gain"])
    boost = float(proposal["delta_margin_boost"])

    target = next(s for s in base["scenarios"] if s["id"] == target_id)
    boosted_after = target["margin_after"] + boost
    gain = boosted_after - target["margin_before"]

    if target["replay_gap"] > replay_tol:
        status = "FAIL_TARGETED_RERUN_REPLAY"
    elif gain < min_gain:
        status = "FAIL_TARGETED_RERUN_MARGIN"
    else:
        status = "PASS_TARGETED_RERUN_LOCAL_ONLY"

    out = {
        "packet": "P1451",
        "status": status,
        "target_scenario_id": target_id,
        "replay_tol": replay_tol,
        "min_gain": min_gain,
        "margin_before": target["margin_before"],
        "margin_after_original": target["margin_after"],
        "delta_margin_boost": boost,
        "margin_after_boosted": boosted_after,
        "gain_boosted": gain,
        "replay_gap": target["replay_gap"],
        "scope_of_pass": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
    }

    OUT.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1451] status={status} -> {OUT}")


if __name__ == "__main__":
    main()
