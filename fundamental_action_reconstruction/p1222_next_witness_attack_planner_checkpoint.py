#!/usr/bin/env python3
"""P1222: build execution-ready plan for the next open witness obligation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1221", type=Path, default=GEN / "p1221_post_w4_discharge_continuation_summary.json")
    parser.add_argument("--p1192", type=Path, default=GEN / "p1192_selector_premise_witness_map_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1222_next_witness_attack_planner_summary.json")
    args = parser.parse_args()

    p1221 = json.loads(args.p1221.read_text(encoding="utf-8"))
    p1192 = json.loads(args.p1192.read_text(encoding="utf-8"))

    next_target = p1221.get("next_target")
    witness_map = p1192.get("witness_obligations", []) if isinstance(p1192.get("witness_obligations"), list) else []

    selected = None
    for w in witness_map:
        if isinstance(w, dict) and w.get("status") == "OPEN":
            selected = w
            break

    out = {
        "packet": "P1222",
        "as_of": "2026-05-11",
        "continuation_signal": next_target,
        "selected_open_witness": selected,
        "execution_plan": {
            "input_artifact": "next_open_witness_symbolic_or_numeric_payload.json",
            "expected_output": "witness_check_summary.json",
            "pass_criterion": "witness_status == DISCHARGED and strict_closure_claim_allowed == false",
            "fail_criterion": "witness_status != DISCHARGED",
        },
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Planner checkpoint only; no witness proof generated here.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1222] selected_open_witness={bool(selected)} wrote {args.out}")


if __name__ == "__main__":
    main()
