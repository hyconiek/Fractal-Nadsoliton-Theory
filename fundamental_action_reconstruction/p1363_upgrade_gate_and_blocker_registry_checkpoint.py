#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {}


def main() -> None:
    b = _load(GEN / "p1362_strict_candidate_residual_benchmark_summary.json")
    items = b.get("strict_candidates_benchmark", [])

    upgrades = []
    blockers = []
    for it in items:
        cid = it.get("claim_id")
        if it.get("strict_upgrade_ready") is True:
            upgrades.append({"claim_id": cid, "new_status": "strict_verified"})
            continue

        reasons = []
        if it.get("residual_status") != "PASS":
            reasons.append("residual_not_pass")
        ub = it.get("uncertainty_budget", {})
        if ub.get("declared") is not True:
            reasons.append("uncertainty_budget_missing")
        if cid == "C3_fine_structure_successor":
            reasons.append("successor_role_equivalence_theorem_missing")

        blockers.append({"claim_id": cid, "blockers": reasons})

    out = {
        "packet": "P1363",
        "as_of": "2026-05-12",
        "upgrade_actions": upgrades,
        "blocker_registry": blockers,
        "upgrade_count": len(upgrades),
        "blocked_count": len(blockers),
        "next_priority": "P1364_BLOCKER_TARGETED_CORRECTIVE_RUNS",
    }

    path = GEN / "p1363_upgrade_gate_and_blocker_registry_summary.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1363] wrote {path}; upgrades={len(upgrades)} blocked={len(blockers)}")


if __name__ == "__main__":
    main()
