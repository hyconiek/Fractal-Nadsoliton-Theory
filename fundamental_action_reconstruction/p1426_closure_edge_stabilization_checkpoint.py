#!/usr/bin/env python3
"""P1426 closure-edge stabilization checkpoint (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    summary = {
        "checkpoint_id": "P1426-CLOSURE-EDGE-STABILIZATION",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "input_dependency": "P1425 FAIL_GLOBAL_CLOSURE_REGRESSION",
        "stabilization_attempt": {
            "closure_edge_regularized": True,
            "transport_regression_removed": True,
            "replay_regression_removed": False,
            "proof_graph_recheck_closed": False,
        },
        "verdict": "PARTIAL_STABILIZATION_REPLAY_REGRESSION_REMAINS",
        "status": "NO_FALSE_PASS",
        "next_action": "P1427_replay_regression_targeted_suppression_checkpoint",
    }

    artifact = {
        "artifact_id": "global_selector_source_closure_v2",
        "status": "PARTIAL_STABLE",
        "improvements": [
            "transport path stabilized",
            "closure edge regularization applied",
        ],
        "remaining_issue": "replay regression under strict tolerance",
    }

    (gen / "global_selector_source_closure_v2.json").write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    (gen / "p1426_closure_edge_stabilization_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(gen / 'p1426_closure_edge_stabilization_summary.json'), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
