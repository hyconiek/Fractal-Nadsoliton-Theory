#!/usr/bin/env python3
"""P1296: R12 strict-selector closure motion review checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1295", type=Path, default=GEN / "p1295_qw2191_r11_formal_proof_completion_and_peer_replay_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1296_qw2191_r12_strict_selector_closure_motion_review_summary.json")
    args = parser.parse_args()

    p1295 = _read(args.p1295)
    if p1295.get("next_priority") != "R12_STRICT_SELECTOR_CLOSURE_MOTION_REVIEW":
        raise SystemExit("P1296 requires next_priority=R12_STRICT_SELECTOR_CLOSURE_MOTION_REVIEW from P1295.")

    replay_ok = p1295.get("r11", {}).get("peer_replay", {}).get("independent_run_status") == "PASS"
    if not replay_ok:
        raise SystemExit("P1296 requires PASS peer replay from P1295.")

    review = {
        "motion_id": "MOTION_R12_STRICT_SELECTOR_CLOSURE_REVIEW_V1",
        "checklist": [
            "formal_proof_completion_present",
            "peer_replay_pass",
            "countermodel_sweep_pass",
            "bridge_nb1_governance_still_open",
        ],
        "review_result": "CONDITIONAL_HOLD",
        "reason": "Global closure remains blocked pending governance theorem obligations.",
    }

    out = {
        "packet": "P1296",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1295": str(args.p1295)},
        "r12_motion_review": review,
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R13_GOVERNANCE_THEOREM_INTERFACE_PREPARATION",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1296] wrote {args.out}; result={review['review_result']}")


if __name__ == "__main__":
    main()
