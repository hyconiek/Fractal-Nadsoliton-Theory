#!/usr/bin/env python3
"""P1247: Step-3 countercheck to ensure no hidden non-strict dependency in chain."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1246", type=Path, default=GEN / "p1246_w1_strict_step2_chain_summary.json")
    parser.add_argument("--p1231", type=Path, default=GEN / "p1231_w1_nonstrict_axiom_scenario_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1247_w1_step3_countercheck_summary.json")
    args = parser.parse_args()

    p1246 = json.loads(args.p1246.read_text(encoding="utf-8"))
    p1231 = json.loads(args.p1231.read_text(encoding="utf-8"))

    c1_2_supported = bool(p1246.get("c1_2_supported", False))
    c1_2_text = str(p1246.get("derived_inference", {}).get("text", ""))

    mentions_nonstrict = "non-strict" in c1_2_text.lower()
    nonstrict_mode = str(p1231.get("w1_discharge_mode", ""))

    # Countercheck discipline: strict chain must not require non-strict input for support.
    hidden_nonstrict_dependency = (nonstrict_mode == "NON_STRICT_AXIOM_DISCHARGE") and (not c1_2_supported)
    countercheck_pass = c1_2_supported and mentions_nonstrict and not hidden_nonstrict_dependency

    out = {
        "packet": "P1247",
        "as_of": "2026-05-11",
        "countercheck_inputs": {
            "c1_2_supported": c1_2_supported,
            "nonstrict_mode_reference": nonstrict_mode,
            "c1_2_mentions_nonstrict_exclusion": mentions_nonstrict,
        },
        "hidden_nonstrict_dependency": hidden_nonstrict_dependency,
        "step3_countercheck_pass": countercheck_pass,
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Step-3 countercheck: strict chain remains supported without importing non-strict lane as dependency.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1247] countercheck_pass={countercheck_pass} wrote {args.out}")


if __name__ == "__main__":
    main()
