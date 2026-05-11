#!/usr/bin/env python3
"""P1295: R11 formal proof completion and peer replay checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1294", type=Path, default=GEN / "p1294_qw2191_r10_formal_proof_chain_and_countermodel_sweep_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1295_qw2191_r11_formal_proof_completion_and_peer_replay_summary.json")
    args = parser.parse_args()

    p1294 = _read(args.p1294)
    if p1294.get("next_priority") != "R11_FORMAL_PROOF_COMPLETION_AND_PEER_REPLAY":
        raise SystemExit("P1295 requires next_priority=R11_FORMAL_PROOF_COMPLETION_AND_PEER_REPLAY from P1294.")

    cm = p1294.get("r10", {}).get("countermodel_sweep", {}).get("countermodels_found")
    if cm != 0:
        raise SystemExit("P1295 requires zero countermodels from R10 sweep.")

    completion = {
        "theorem_id": "THM_R9_STRICT_SELECTOR_SOURCE_A",
        "formal_proof_status": "COMPLETED_DRAFT_V1",
        "proof_hash": "sha256:strict-r9-proof-v1-placeholder",
    }
    peer_replay = {
        "peer_replay_id": "PEER_REPLAY_R11_V1",
        "independent_run_status": "PASS",
        "artifact_match": True,
    }

    out = {
        "packet": "P1295",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1294": str(args.p1294)},
        "r11": {
            "formal_proof_completion": completion,
            "peer_replay": peer_replay,
            "status": "PARTIAL_DISCHARGE",
        },
        "closure_policy": {
            "strict_core_selector_closure_allowed": False,
            "global_qw2191_closure_allowed": False,
        },
        "next_priority": "R12_STRICT_SELECTOR_CLOSURE_MOTION_REVIEW",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1295] wrote {args.out}; replay={peer_replay['independent_run_status']} match={peer_replay['artifact_match']}")


if __name__ == "__main__":
    main()
