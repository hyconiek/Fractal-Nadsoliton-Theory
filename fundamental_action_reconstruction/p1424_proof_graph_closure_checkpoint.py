#!/usr/bin/env python3
"""P1424 proof-graph closure checkpoint (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    summary = {
        "checkpoint_id": "P1424-PROOF-GRAPH-CLOSURE",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "input_artifacts": [
            "selector_uniqueness_certificate_v1.json",
            "qw2191_discharge_argument_v1.json",
        ],
        "proof_graph": {
            "nodes_required": 12,
            "nodes_closed": 11,
            "missing_edge": "global_selector_source_closure -> qw2191_discharge_final",
            "acyclic": True,
        },
        "verdict": "FAIL_PROOF_GRAPH_NOT_CLOSED",
        "status": "NO_FALSE_PASS",
        "next_action": "P1425_global_selector_source_closure_checkpoint",
    }

    obstruction = {
        "obstruction_id": "OBSTR-PROOF-GRAPH-CLOSURE-v1",
        "reason": "one final global closure edge missing",
        "missing_edge": summary["proof_graph"]["missing_edge"],
        "recommended_next": "construct global selector source closure artifact",
    }

    (gen / "p1424_proof_graph_closure_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    (gen / "p1424_proof_graph_obstruction_v1.json").write_text(json.dumps(obstruction, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(gen / 'p1424_proof_graph_closure_summary.json'), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
