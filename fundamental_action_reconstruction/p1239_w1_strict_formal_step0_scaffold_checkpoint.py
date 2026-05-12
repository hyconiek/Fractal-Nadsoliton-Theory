#!/usr/bin/env python3
"""P1239: create strict formal step-0 scaffold gated by P1238 pass."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1238", type=Path, default=GEN / "p1238_w1_strict_formalization_entry_gate_summary.json")
    parser.add_argument("--scaffold", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold_summary.json")
    args = parser.parse_args()

    gate = json.loads(args.p1238.read_text(encoding="utf-8"))
    gate_pass = bool(gate.get("strict_formalization_entry_gate_pass", False))

    scaffold_written = False
    if gate_pass:
        scaffold = {
            "packet": "P1239_SCAFFOLD",
            "as_of": "2026-05-11",
            "w1_strict_formal_statement": "TODO_FORMAL_STATEMENT_V1",
            "w1_strict_formal_assumptions": ["TODO_ASSUMPTION_1", "TODO_ASSUMPTION_2"],
            "w1_scope_boundary": "strict_lane_only_no_global_closure_claim",
            "verification_hooks": ["hook_symbolic_consistency", "hook_selector_uniqueness_trace"],
            "strict_closure_claim_allowed": False,
            "theory_closure_status": "OPEN",
        }
        args.scaffold.write_text(json.dumps(scaffold, indent=2, sort_keys=True) + "\n", encoding="utf-8")
        scaffold_written = True

    out = {
        "packet": "P1239",
        "as_of": "2026-05-11",
        "entry_gate_pass": gate_pass,
        "scaffold_written": scaffold_written,
        "scaffold_path": str(args.scaffold) if scaffold_written else None,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Strict formal step-0 scaffold is created only when P1238 gate passes.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1239] gate_pass={gate_pass} scaffold_written={scaffold_written} wrote {args.out}")


if __name__ == "__main__":
    main()
