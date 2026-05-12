#!/usr/bin/env python3
"""P1241: fill strict step-0 scaffold with first draft content and rerun fill-check."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--scaffold", type=Path, default=GEN / "p1239_w1_strict_formal_step0_scaffold.json")
    parser.add_argument("--p1240-out", type=Path, default=GEN / "p1240_w1_strict_statement_fill_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1241_w1_strict_step0_content_draft_summary.json")
    args = parser.parse_args()

    scaffold = json.loads(args.scaffold.read_text(encoding="utf-8"))
    scaffold["w1_strict_formal_statement"] = (
        "On the strict lane only, W1 is locally discharged when a strict selector-source reference is present "
        "and strict-path discharge is realized without invoking non-strict symmetry-breaking axioms."
    )
    scaffold["w1_strict_formal_assumptions"] = [
        "strict_selector_source_theorem_exported=true for the active candidate reference",
        "symmetry_state is UNBROKEN_IN_STRICT_CORE on the strict lane",
        "theory_closure_status remains OPEN and no global closure claim is made",
    ]
    args.scaffold.write_text(json.dumps(scaffold, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    root = Path(__file__).resolve().parent
    subprocess.run(
        [
            "python3",
            str(root / "p1240_w1_strict_statement_fill_checkpoint.py"),
            "--scaffold",
            str(args.scaffold),
            "--out",
            str(args.p1240_out),
        ],
        check=True,
    )

    p1240 = json.loads(args.p1240_out.read_text(encoding="utf-8"))
    out = {
        "packet": "P1241",
        "as_of": "2026-05-11",
        "draft_content_written": True,
        "post_fill_strict_formal_step0_ready": p1240.get("strict_formal_step0_ready"),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "First strict step-0 content draft applied and validated through P1240 fill-check.",
    }
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1241] ready={out['post_fill_strict_formal_step0_ready']} wrote {args.out}")


if __name__ == "__main__":
    main()
