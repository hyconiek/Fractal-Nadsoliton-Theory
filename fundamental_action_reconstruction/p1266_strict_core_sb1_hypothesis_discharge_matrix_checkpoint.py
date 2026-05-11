#!/usr/bin/env python3
"""P1266: strict-core SB1 hypothesis discharge matrix checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1265", type=Path, default=GEN / "p1265_strict_core_symmetry_breaking_source_theorem_scaffold_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1266_strict_core_sb1_hypothesis_discharge_matrix_summary.json")
    args = parser.parse_args()

    p1265 = _read(args.p1265)
    if p1265.get("lane") != "STRICT_CORE_ONLY":
        raise SystemExit("P1266 requires STRICT_CORE_ONLY lane from P1265.")

    matrix = [
        {
            "hypothesis": "H1",
            "lemma_id": "L_SB1_1",
            "required_artifact": "generated/p1246_w1_strict_step2_chain_summary.json",
            "status": "PARTIAL",
            "note": "Locality envelope partially evidenced by strict-step chain continuity.",
        },
        {
            "hypothesis": "H2",
            "lemma_id": "L_SB1_2",
            "required_artifact": "generated/p1231_w1_nonstrict_axiom_scenario_summary.json",
            "status": "OPEN",
            "note": "Need strict exclusion proof that selector-breaking step imports no non-strict axioms.",
        },
        {
            "hypothesis": "H3",
            "lemma_id": "L_SB1_3",
            "required_artifact": "generated/p1255_bridge_or_nonbridge_decision_gate_summary.json",
            "status": "PARTIAL",
            "note": "Kernel-split policy is active; formal non-transfer lemma still open.",
        },
        {
            "hypothesis": "H4",
            "lemma_id": "L_SB1_4",
            "required_artifact": "generated/p1249_w1_step5_obstruction_interface_summary.json",
            "status": "OPEN",
            "note": "Bounded obstruction interface must be converted into strict theorem form.",
        },
        {
            "hypothesis": "QW2191_COMPAT",
            "lemma_id": "L_SB1_QW2191",
            "required_artifact": "generated/p1250_w1_strict_uniqueness_witness_candidate_summary.json",
            "status": "OPEN",
            "note": "Explicit coexistence/resolution theorem with QW-2191 not yet exported.",
        },
    ]

    open_count = sum(1 for r in matrix if r["status"] == "OPEN")
    out = {
        "packet": "P1266",
        "as_of": "2026-05-11",
        "lane": "STRICT_CORE_ONLY",
        "input": {"p1265": str(args.p1265)},
        "sb1_hypothesis_matrix": matrix,
        "open_count": open_count,
        "sb1_discharge_ready": open_count == 0,
        "closure_policy": "STRICT_KERNEL_CLOSURE_FORBIDDEN_UNTIL_SB1_MATRIX_OPEN_COUNT_ZERO_AND_B1_OR_NB1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1266] wrote {args.out}; open_count={open_count}")


if __name__ == "__main__":
    main()
