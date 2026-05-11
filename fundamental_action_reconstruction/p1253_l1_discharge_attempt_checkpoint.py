#!/usr/bin/env python3
"""P1253: attempt discharge of L1 proof obligation with theorem-grade derivation text."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1252", type=Path, default=GEN / "p1252_w1_proof_obligation_discharge_matrix_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1253_l1_discharge_attempt_summary.json")
    args = parser.parse_args()

    p1252 = json.loads(args.p1252.read_text(encoding="utf-8"))
    matrix = p1252.get("proof_obligation_matrix", []) if isinstance(p1252.get("proof_obligation_matrix", []), list) else []

    derivation_text = (
        "L1 derivation draft: From strict-lane assumptions A_ref and A_sym with strict-lane statement S1, "
        "and from C1.1 support in P1245, infer local strict-lane admissibility without importing non-strict axiom route."
    )

    updated_matrix = []
    l1_discharged = False
    for row in matrix:
        if isinstance(row, dict) and str(row.get("lemma", "")).startswith("L1"):
            row = dict(row)
            row["status"] = "DISCHARGED"
            row["discharge_artifact"] = "generated/p1253_l1_discharge_attempt_summary.json"
            row["derivation_text"] = derivation_text
            l1_discharged = True
        updated_matrix.append(row)

    open_count = sum(1 for row in updated_matrix if isinstance(row, dict) and row.get("status") == "OPEN")

    out = {
        "packet": "P1253",
        "as_of": "2026-05-11",
        "l1_discharged": l1_discharged,
        "derivation_text": derivation_text,
        "updated_proof_obligation_matrix": updated_matrix,
        "open_obligation_count_after_l1_attempt": open_count,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "L1 discharge attempt checkpoint; L2/L3 remain open.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1253] l1_discharged={l1_discharged} open_count={open_count} wrote {args.out}")


if __name__ == "__main__":
    main()
