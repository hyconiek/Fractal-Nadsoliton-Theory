#!/usr/bin/env python3
"""P1232: strict-vs-nonstrict comparative summary for W1."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--strict", type=Path, default=GEN / "p1227_w1_real_candidate_run_summary.json")
    parser.add_argument("--nonstrict", type=Path, default=GEN / "p1231_w1_nonstrict_axiom_scenario_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1232_w1_strict_vs_nonstrict_comparative_summary.json")
    args = parser.parse_args()

    strict = _load(args.strict)
    nonstrict = _load(args.nonstrict)

    strict_row = {
        "lane": "STRICT_CANDIDATE",
        "discharge_mode": strict.get("w1_discharge_mode"),
        "witness_status": strict.get("w1_witness_status"),
        "strict_core_closure_eligibility": True if strict.get("w1_discharge_mode") == "STRICT_PATH_DISCHARGE" else False,
        "symmetry_state": "UNBROKEN_IN_STRICT_CORE",
    }
    nonstrict_row = {
        "lane": "NONSTRICT_AXIOM",
        "discharge_mode": nonstrict.get("w1_discharge_mode"),
        "witness_status": nonstrict.get("w1_witness_status"),
        "strict_core_closure_eligibility": bool(nonstrict.get("strict_core_closure_eligibility", False)),
        "symmetry_state": nonstrict.get("symmetry_state", "AXIOM_DECLARED_NON_STRICT"),
    }

    out = {
        "packet": "P1232",
        "as_of": "2026-05-11",
        "comparison": [strict_row, nonstrict_row],
        "delta": {
            "discharge_mode_differs": strict_row["discharge_mode"] != nonstrict_row["discharge_mode"],
            "strict_core_eligibility_differs": strict_row["strict_core_closure_eligibility"] != nonstrict_row["strict_core_closure_eligibility"],
            "symmetry_state_differs": strict_row["symmetry_state"] != nonstrict_row["symmetry_state"],
        },
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Comparative checkpoint only; keeps strict and non-strict lanes explicitly separated.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1232] wrote {args.out}")


if __name__ == "__main__":
    main()
