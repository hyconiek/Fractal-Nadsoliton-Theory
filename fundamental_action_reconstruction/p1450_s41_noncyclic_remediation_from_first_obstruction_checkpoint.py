#!/usr/bin/env python3
"""P1450 checkpoint: export noncyclic remediation proposal from first obstruction."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OBS = ROOT / "generated" / "p1449_s4_first_obstruction.json"
OUT = ROOT / "generated" / "p1450_s41_noncyclic_remediation_proposal.json"
SUM = ROOT / "generated" / "p1450_s41_noncyclic_remediation_summary.json"


def main() -> None:
    if not OBS.exists():
        summary = {"packet": "P1450", "status": "NO_OBSTRUCTION_INPUT", "legacy_bridge_used": False}
        SUM.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1450] status=NO_OBSTRUCTION_INPUT -> {SUM}")
        return

    obs = json.loads(OBS.read_text(encoding="utf-8"))
    scenario = obs.get("scenario", {})
    gain = float(obs.get("gain", 0.0))
    min_gain = float(obs.get("min_gain", 0.0))

    needed = max(0.0, min_gain - gain)
    proposal = {
        "packet": "P1450",
        "proposal_type": "NONCYCLIC_PROVIDER_LOCAL_MARGIN_BOOST_V1",
        "target_scenario_id": scenario.get("id", "unknown"),
        "delta_margin_boost": round(needed + 0.0005, 6),
        "scope_of_pass": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "note": "Proposal only; requires replay validation in rerun of P1449.",
    }

    summary = {
        "packet": "P1450",
        "status": "PROPOSAL_EXPORTED",
        "source_obstruction_reason": obs.get("reason"),
        "target_scenario_id": proposal["target_scenario_id"],
        "legacy_bridge_used": False,
    }

    OUT.write_text(json.dumps(proposal, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    SUM.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1450] status=PROPOSAL_EXPORTED -> {OUT}")


if __name__ == "__main__":
    main()
