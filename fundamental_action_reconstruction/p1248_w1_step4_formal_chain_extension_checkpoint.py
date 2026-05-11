#!/usr/bin/env python3
"""P1248: Step-4 formal chain extension deriving C1.3."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1246", type=Path, default=GEN / "p1246_w1_strict_step2_chain_summary.json")
    parser.add_argument("--p1247", type=Path, default=GEN / "p1247_w1_step3_countercheck_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1248_w1_step4_formal_chain_extension_summary.json")
    args = parser.parse_args()

    p1246 = json.loads(args.p1246.read_text(encoding="utf-8"))
    p1247 = json.loads(args.p1247.read_text(encoding="utf-8"))

    c1_2_supported = bool(p1246.get("c1_2_supported", False))
    countercheck_pass = bool(p1247.get("step3_countercheck_pass", False))

    c1_3_supported = c1_2_supported and countercheck_pass
    c1_3_text = (
        "C1.3: From C1.2 and Step-3 countercheck pass, strict-lane local witness chain remains admissible "
        "through Step-4 without non-strict dependence."
    ) if c1_3_supported else "C1.3 not derivable from current Step-2/Step-3 state"

    out = {
        "packet": "P1248",
        "as_of": "2026-05-11",
        "chain_rule": "C1.2 + Countercheck_pass -> C1.3",
        "c1_2_supported": c1_2_supported,
        "step3_countercheck_pass": countercheck_pass,
        "c1_3_supported": c1_3_supported,
        "derived_inference": {
            "id": "C1.3",
            "text": c1_3_text,
        },
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Step-4 formal chain extension checkpoint.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1248] C1.3_supported={c1_3_supported} wrote {args.out}")


if __name__ == "__main__":
    main()
