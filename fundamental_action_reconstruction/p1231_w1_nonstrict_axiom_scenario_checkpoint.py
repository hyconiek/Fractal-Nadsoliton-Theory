#!/usr/bin/env python3
"""P1231: explicit non-strict axiom scenario note without mixing strict-core claim."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _dump(path: Path, data: dict) -> None:
    path.write_text(json.dumps(data, indent=2, sort_keys=True) + "\n", encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--base-input", type=Path, default=GEN / "p1223_w1_input_payload.json")
    parser.add_argument("--scenario-input", type=Path, default=GEN / "p1231_w1_input_payload_nonstrict_axiom_scenario.json")
    parser.add_argument("--p1223-out", type=Path, default=GEN / "p1231_w1_execution_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1231_w1_nonstrict_axiom_scenario_summary.json")
    args = parser.parse_args()

    payload = _load(args.base_input)
    avail = payload.setdefault("available_inputs", {})
    avail["strict_selector_source_theorem_exported"] = False
    avail["explicit_symmetry_breaking_axiom_declared"] = True
    avail["scoped_non_strict_tag_if_axiom_used"] = True
    _dump(args.scenario_input, payload)

    root = Path(__file__).resolve().parent
    subprocess.run([
        "python3", str(root / "p1223_w1_execution_checkpoint.py"),
        "--input", str(args.scenario_input), "--out", str(args.p1223_out)
    ], check=True)

    p1223 = _load(args.p1223_out)
    out = {
        "packet": "P1231",
        "as_of": "2026-05-11",
        "scenario": "NONSTRICT_AXIOM_ONLY",
        "w1_witness_status": p1223.get("witness_status"),
        "w1_discharge_mode": p1223.get("evaluation", {}).get("discharge_mode"),
        "strict_core_closure_eligibility": p1223.get("evaluation", {}).get("strict_core_closure_eligibility"),
        "symmetry_state": "AXIOM_DECLARED_NON_STRICT",
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Scenario recorded explicitly as non-strict; not merged into strict-core closure lane.",
    }
    _dump(args.out, out)
    print(f"[P1231] mode={out['w1_discharge_mode']} strict_elig={out['strict_core_closure_eligibility']} wrote {args.out}")


if __name__ == "__main__":
    main()
