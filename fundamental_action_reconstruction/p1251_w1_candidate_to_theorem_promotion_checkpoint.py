#!/usr/bin/env python3
"""P1251: promote strict uniqueness candidate into theorem-grade draft skeleton."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1250", type=Path, default=GEN / "p1250_w1_strict_uniqueness_witness_candidate_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1251_w1_candidate_to_theorem_promotion_summary.json")
    args = parser.parse_args()

    p1250 = json.loads(args.p1250.read_text(encoding="utf-8"))
    candidate = p1250.get("candidate", {}) if isinstance(p1250.get("candidate"), dict) else {}

    candidate_ready = bool(p1250.get("all_gaps_candidate_resolved", False))

    theorem_draft = {
        "packet": "P1251_THEOREM_DRAFT",
        "as_of": "2026-05-11",
        "theorem_id": "T_W1_STRICT_UNIQUENESS_V1_DRAFT",
        "hypotheses": [
            "H1: strict selector-source reference exists and is bound to strict lane",
            "H2: strict chain C1.1 -> C1.2 -> C1.3 is supported",
            "H3: non-strict lane is not imported as dependency in strict chain",
        ],
        "lemmas": [
            "L1: local strict-lane admissibility (from C1.1)",
            "L2: strict chain consistency extension (from C1.2, C1.3)",
            "L3: obstruction-interface candidate closure at candidate-level only",
        ],
        "conclusion": "Under H1-H3 and L1-L3, W1 strict uniqueness candidate is theorem-grade draft ready for formal proof obligations.",
        "source_hooks": candidate.get("evidence_hooks", []),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
    }

    out = {
        "packet": "P1251",
        "as_of": "2026-05-11",
        "candidate_ready_for_promotion": candidate_ready,
        "theorem_draft_emitted": candidate_ready,
        "theorem_draft": theorem_draft if candidate_ready else None,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Candidate-to-theorem promotion checkpoint; proof obligations remain to be discharged.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1251] emitted={out['theorem_draft_emitted']} wrote {args.out}")


if __name__ == "__main__":
    main()
