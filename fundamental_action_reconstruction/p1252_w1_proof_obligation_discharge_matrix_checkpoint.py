#!/usr/bin/env python3
"""P1252: build theorem-draft proof-obligation discharge matrix."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1251", type=Path, default=GEN / "p1251_w1_candidate_to_theorem_promotion_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1252_w1_proof_obligation_discharge_matrix_summary.json")
    args = parser.parse_args()

    p1251 = json.loads(args.p1251.read_text(encoding="utf-8"))
    draft = p1251.get("theorem_draft", {}) if isinstance(p1251.get("theorem_draft"), dict) else {}
    lemmas = draft.get("lemmas", []) if isinstance(draft.get("lemmas", []), list) else []

    matrix = []
    for lemma in lemmas:
        if "L1" in lemma:
            source = "generated/p1245_w1_strict_step1_inference_summary.json"
        elif "L2" in lemma:
            source = "generated/p1248_w1_step4_formal_chain_extension_summary.json"
        else:
            source = "generated/p1249_w1_step5_obstruction_interface_summary.json"

        matrix.append(
            {
                "lemma": lemma,
                "status": "OPEN",
                "pass_condition": "export theorem-grade derivation text with explicit symbolic rule chain",
                "current_source_artifact": source,
            }
        )

    open_count = sum(1 for row in matrix if row["status"] == "OPEN")

    out = {
        "packet": "P1252",
        "as_of": "2026-05-11",
        "theorem_id": draft.get("theorem_id", "UNKNOWN"),
        "proof_obligation_matrix": matrix,
        "open_obligation_count": open_count,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Proof-obligation discharge matrix for theorem-grade promotion workflow.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1252] open_obligation_count={open_count} wrote {args.out}")


if __name__ == "__main__":
    main()
