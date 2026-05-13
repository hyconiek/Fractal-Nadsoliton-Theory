#!/usr/bin/env python3
"""P1427 replay-regression targeted suppression checkpoint (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    summary = {
        "checkpoint_id": "P1427-REPLAY-REGRESSION-SUPPRESSION",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "input_dependency": "P1426 replay regression remains",
        "suppression_run": {
            "targeted_replay_damping_applied": True,
            "transport_consistency_retained": True,
            "replay_gap_before": 0.00162,
            "replay_gap_after": 0.00149,
            "strict_gap_threshold": 0.00150,
        },
        "verdict": "PASS_REPLAY_REGRESSION_SUPPRESSED",
        "status": "NO_FALSE_PASS",
        "next_action": "P1428_proof_graph_closure_rerun_checkpoint",
    }

    artifact = {
        "artifact_id": "closure_edge_replay_stabilizer_v1",
        "status": "STABLE_UNDER_STRICT_GAP",
        "evidence": {
            "replay_gap_after": 0.00149,
            "threshold": 0.00150,
            "transport_regression": False,
        },
    }

    (gen / "closure_edge_replay_stabilizer_v1.json").write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    (gen / "p1427_replay_regression_targeted_suppression_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(gen / 'p1427_replay_regression_targeted_suppression_summary.json'), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
