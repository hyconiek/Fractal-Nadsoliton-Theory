#!/usr/bin/env python3
"""P1256: formal specification checkpoint for legacy->strict bridge vs non-bridge theorem.

This stage defines the minimal theorem-spec surface and proof-obligation map
required before any strict-core closure claim can be considered.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1255", type=Path, default=GEN / "p1255_bridge_or_nonbridge_decision_gate_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1256_bridge_nonbridge_theorem_spec_summary.json")
    args = parser.parse_args()

    p1255 = _read(args.p1255)
    if p1255.get("decision") != "PROCEED_TO_BRIDGE_OR_NONBRIDGE_FORMALIZATION_ONLY":
        raise SystemExit("P1256 requires P1255 decision gate in formalization-only mode.")

    theorem_spec = {
        "bridge_theorem": {
            "statement_id": "B1",
            "goal": "Establish explicit map from legacy ontological kernel to strict gate kernel under declared hypotheses.",
            "required_objects": [
                "domain-alignment map",
                "parameter-role transfer map",
                "error/obstruction bounds",
            ],
        },
        "nonbridge_theorem": {
            "statement_id": "NB1",
            "goal": "Prove impossibility (or inconsistency) of legacy->strict identification under canonical assumptions.",
            "required_objects": [
                "obstruction witness",
                "invariance violation or no-go lemma",
                "selector-consistent countermodel class",
            ],
        },
    }

    proof_obligations = [
        {
            "id": "O1",
            "lane": "strict",
            "text": "No strict-core closure claim without resolved B1 or NB1.",
            "status": "OPEN",
        },
        {
            "id": "O2",
            "lane": "strict",
            "text": "QW-2191 compatibility must be explicit in selected theorem lane.",
            "status": "OPEN",
        },
        {
            "id": "O3",
            "lane": "cross-lane",
            "text": "No implicit legacy physical-role transfer onto strict kernel.",
            "status": "OPEN",
        },
    ]

    out = {
        "packet": "P1256",
        "as_of": "2026-05-11",
        "input": {"p1255": str(args.p1255)},
        "entry_gate": p1255.get("gate_pass", False),
        "theorem_spec": theorem_spec,
        "proof_obligations": proof_obligations,
        "closure_policy": "STRICT_CLOSURE_FORBIDDEN_UNTIL_B1_OR_NB1_DISCHARGED",
        "next_action": "Select B1 (bridge) or NB1 (non-bridge) and instantiate formal assumptions.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1256] wrote {args.out}")


if __name__ == "__main__":
    main()
