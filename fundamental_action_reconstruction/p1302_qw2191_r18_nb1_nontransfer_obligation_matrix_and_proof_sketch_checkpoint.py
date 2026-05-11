#!/usr/bin/env python3
"""P1302: R18 NB1 non-transfer obligation matrix and proof sketch checkpoint."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _read(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1301", type=Path, default=GEN / "p1301_qw2191_r17_nb1_formal_nontransfer_theorem_draft_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1302_qw2191_r18_nb1_nontransfer_obligation_matrix_and_proof_sketch_summary.json")
    args = parser.parse_args()

    p1301 = _read(args.p1301)
    if p1301.get("next_priority") != "R18_NB1_NONTRANSFER_OBLIGATION_MATRIX_AND_PROOF_SKETCH":
        raise SystemExit("P1302 requires next_priority=R18_NB1_NONTRANSFER_OBLIGATION_MATRIX_AND_PROOF_SKETCH from P1301.")

    r17 = p1301.get("r17_nb1_nontransfer_theorem", {})
    if r17.get("status") != "DRAFT_WITH_OBLIGATIONS":
        raise SystemExit("P1302 requires DRAFT_WITH_OBLIGATIONS from P1301.")

    matrix = [
        {
            "id": "O1_DOMAIN_MISMATCH",
            "lemma": "Legacy parameter-role codomain is not equal to strict theorem-role codomain.",
            "evidence_type": "typed_role_map_counterexample_or_nonisomorphism_proof",
            "pass_criterion": "No total role-preserving map legacy->strict exists.",
        },
        {
            "id": "O2_TRANSPORT_FUNCTOR_ABSENCE",
            "lemma": "No admissible transport functor preserves strict closure predicates.",
            "evidence_type": "functorial_obstruction_certificate",
            "pass_criterion": "Any candidate functor violates at least one strict closure invariant.",
        },
        {
            "id": "O3_QW2191_PERSISTENCE",
            "lemma": "QW-2191 obstruction persists without new internal selector source.",
            "evidence_type": "obstruction_replay_with_provider_null_result",
            "pass_criterion": "Obstruction class remains non-empty under NB assumptions.",
        },
    ]

    out = {
        "packet": "P1302",
        "as_of": "2026-05-11",
        "lane": "NB1_NONBRIDGE_TRACK",
        "input": {"p1301": str(args.p1301)},
        "r18_obligation_matrix": matrix,
        "proof_sketch_status": "MATRIX_DRAFTED",
        "next_priority": "R19_NB1_NONTRANSFER_LEMMA_PACK_V1",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1302] wrote {args.out}; obligations={len(matrix)}")


if __name__ == "__main__":
    main()
