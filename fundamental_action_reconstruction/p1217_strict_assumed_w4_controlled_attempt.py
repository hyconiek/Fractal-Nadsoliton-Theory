#!/usr/bin/env python3
"""P1217: run controlled W4 attempt on STRICT_ASSUMED gate mode."""
from __future__ import annotations

import argparse
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run_cmd(args: list[str]) -> None:
    subprocess.run(args, check=True, cwd=ROOT.parent)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=GEN / "p1217_strict_assumed_w4_controlled_attempt_summary.json")
    args = parser.parse_args()

    gate_path = GEN / "p1217_p1209_strict_assumed_summary.json"
    attempt_path = GEN / "p1217_p1210_strict_assumed_summary.json"

    run_cmd(["python3", str(ROOT / "p1209_symbolic_attempt_gate_open.py"), "--assume-strict-artifact", "--out", str(gate_path)])
    run_cmd(["python3", str(ROOT / "p1210_controlled_w4_symbolic_attempt.py"), "--p1209", str(gate_path), "--out", str(attempt_path)])

    gate = read_json(gate_path)
    attempt = read_json(attempt_path)

    out = {
        "packet": "P1217",
        "as_of": "2026-05-11",
        "gate_mode": gate.get("gate_mode"),
        "strict_assumption_applied": gate.get("strict_assumption_applied"),
        "symbolic_attempt_gate_open": gate.get("symbolic_attempt_gate_open"),
        "attempt_executed": attempt.get("attempt_executed"),
        "w4_status": attempt.get("w4_status"),
        "w4_discharge_pass": attempt.get("w4_discharge_pass"),
        "theory_closure_status": "OPEN",
        "strict_closure_claim_allowed": False,
        "note": "Operational continuation checkpoint under STRICT_ASSUMED mode.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1217] wrote {args.out}")


if __name__ == "__main__":
    main()
