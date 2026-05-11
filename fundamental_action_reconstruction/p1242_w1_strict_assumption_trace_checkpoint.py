#!/usr/bin/env python3
"""P1242: trace strict step-0 assumptions to concrete checkpoint artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scaffold", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1242_w1_strict_assumption_trace_summary.json")
    args = parser.parse_args()

    scaffold = json.loads(args.scaffold.read_text(encoding="utf-8"))
    assumptions = scaffold.get("w1_strict_formal_assumptions", []) if isinstance(scaffold.get("w1_strict_formal_assumptions", []), list) else []

    trace_map = {
        "strict_selector_source_theorem_exported=true for the active candidate reference": "generated/p1227_w1_execution_summary.json",
        "symmetry_state is UNBROKEN_IN_STRICT_CORE on the strict lane": "generated/p1233_w1_strict_lane_readiness_gate_summary.json",
        "theory_closure_status remains OPEN and no global closure claim is made": "generated/p1238_w1_strict_formalization_entry_gate_summary.json",
    }

    traces = []
    all_traced = True
    for a in assumptions:
        target = trace_map.get(a)
        exists = bool(target) and (ROOT / target).exists()
        traced = bool(target) and exists
        all_traced = all_traced and traced
        traces.append({
            "assumption": a,
            "artifact": target,
            "artifact_exists": exists,
            "traced": traced,
        })

    out = {
        "packet": "P1242",
        "as_of": "2026-05-11",
        "assumption_count": len(assumptions),
        "all_assumptions_traced": all_traced,
        "traces": traces,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Assumption-to-artifact trace checkpoint for strict step-0 draft.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1242] all_assumptions_traced={all_traced} wrote {args.out}")


if __name__ == "__main__":
    main()
