#!/usr/bin/env python3
"""P1215: compare STRICT_CERTIFIED vs PROVISIONAL_UNCERTIFIED W4 gate modes."""
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
    parser.add_argument("--out", type=Path, default=GEN / "p1215_w4_gate_mode_comparison_summary.json")
    args = parser.parse_args()

    strict_gate = GEN / "p1215_p1209_strict_summary.json"
    provisional_gate = GEN / "p1215_p1209_provisional_summary.json"
    strict_attempt = GEN / "p1215_p1210_strict_summary.json"
    provisional_attempt = GEN / "p1215_p1210_provisional_summary.json"

    run_cmd(["python3", str(ROOT / "p1209_symbolic_attempt_gate_open.py"), "--out", str(strict_gate)])
    run_cmd(["python3", str(ROOT / "p1210_controlled_w4_symbolic_attempt.py"), "--p1209", str(strict_gate), "--out", str(strict_attempt)])

    run_cmd([
        "python3", str(ROOT / "p1209_symbolic_attempt_gate_open.py"),
        "--allow-uncertified-artifact", "--out", str(provisional_gate),
    ])
    run_cmd(["python3", str(ROOT / "p1210_controlled_w4_symbolic_attempt.py"), "--p1209", str(provisional_gate), "--out", str(provisional_attempt)])

    s_gate = read_json(strict_gate)
    p_gate = read_json(provisional_gate)
    s_attempt = read_json(strict_attempt)
    p_attempt = read_json(provisional_attempt)

    out = {
        "packet": "P1215",
        "as_of": "2026-05-11",
        "strict_mode": {
            "gate_mode": s_gate.get("gate_mode"),
            "symbolic_attempt_gate_open": s_gate.get("symbolic_attempt_gate_open"),
            "attempt_executed": s_attempt.get("attempt_executed"),
            "w4_status": s_attempt.get("w4_status"),
        },
        "provisional_mode": {
            "gate_mode": p_gate.get("gate_mode"),
            "symbolic_attempt_gate_open": p_gate.get("symbolic_attempt_gate_open"),
            "attempt_executed": p_attempt.get("attempt_executed"),
            "w4_status": p_attempt.get("w4_status"),
        },
        "theory_closure_status": "OPEN",
        "strict_closure_claim_allowed": False,
        "note": "Methodology comparison only: gate behavior vs attempt behavior.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1215] wrote {args.out}")


if __name__ == "__main__":
    main()
