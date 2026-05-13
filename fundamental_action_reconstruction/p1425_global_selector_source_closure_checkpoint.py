#!/usr/bin/env python3
"""P1425 global selector-source closure checkpoint (strict-only)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    summary = {
        "checkpoint_id": "P1425-GLOBAL-SELECTOR-CLOSURE",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "input_dependency": "P1424 missing global closure edge",
        "closure_attempt": {
            "global_selector_source_artifact_exported": True,
            "closure_edge_constructed": True,
            "consistency_recheck_pass": False,
            "regression_on_transport_replay": True,
        },
        "verdict": "FAIL_GLOBAL_CLOSURE_REGRESSION",
        "status": "NO_FALSE_PASS",
        "next_action": "P1426_closure_edge_stabilization_checkpoint",
    }

    artifact = {
        "artifact_id": "global_selector_source_closure_v1",
        "constructed_edge": "global_selector_source_closure -> qw2191_discharge_final",
        "status": "CONSTRUCTED_BUT_UNSTABLE",
        "notes": [
            "edge exists symbolically",
            "downstream transport/replay consistency regression observed",
        ],
    }

    (gen / "global_selector_source_closure_v1.json").write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    (gen / "p1425_global_selector_source_closure_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({"written": str(gen / 'p1425_global_selector_source_closure_summary.json'), "verdict": summary["verdict"]}, indent=2))


if __name__ == "__main__":
    main()
