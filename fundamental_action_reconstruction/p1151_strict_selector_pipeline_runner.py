#!/usr/bin/env python3
"""P1151 strict selector pipeline runner.

Runs mandatory sequence:
1) P1150 admissibility gate on provided candidate premise JSON,
2) P1146 probe,
3) P1147 probe,
4) P1148 probe,
5) P1149 reproducibility audit.

Exports unified run summary JSON and fails fast on any nonzero exit.
"""
from __future__ import annotations
import json
import subprocess
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def run(cmd: list[str]) -> dict:
    proc = subprocess.run(cmd, capture_output=True, text=True)
    return {
        "cmd": cmd,
        "returncode": proc.returncode,
        "stdout": proc.stdout.strip(),
        "stderr": proc.stderr.strip(),
    }


def main() -> None:
    premise = (
        Path(sys.argv[1]).resolve()
        if len(sys.argv) > 1
        else (GEN / "p1150_candidate_input_example.json").resolve()
    )

    payload = {}
    try:
        payload = json.loads(premise.read_text(encoding="utf-8"))
    except Exception:
        payload = {}
    force_fail_stage = payload.get("force_fail_stage")

    steps = [
        [sys.executable, str(ROOT / "p1150_strict_premise_admissibility_gate.py"), str(premise)],
        [sys.executable, str(ROOT / "p1146_strict_shannon_selector_candidate_grid_nonclosure_probe.py")],
        [sys.executable, str(ROOT / "p1147_strict_selector_obstruction_location_probe.py")],
        [sys.executable, str(ROOT / "p1148_strict_phase_shifted_selector_family_probe.py")],
        [sys.executable, str(ROOT / "p1149_strict_selector_probe_reproducibility_audit.py")],
    ]

    logs = []
    overall_pass = True
    for idx, s in enumerate(steps):
        rec = run(s)
        logs.append(rec)
        if isinstance(force_fail_stage, int) and idx == force_fail_stage:
            rec["returncode"] = 1
            rec["stderr"] = (rec.get("stderr") or "") + " | forced_fail_stage"
        if rec["returncode"] != 0:
            overall_pass = False
            break

    failed_step = None
    for idx, rec in enumerate(logs):
        if rec["returncode"] != 0:
            failed_step = {"index": idx, "cmd": rec["cmd"]}
            break

    out = {
        "packet": "P1151",
        "as_of": "2026-05-10",
        "premise_input": str(premise),
        "steps": logs,
        "overall_pass": overall_pass,
        "failed_step": failed_step,
    }

    out_path = GEN / "p1151_strict_selector_pipeline_runner_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1151] overall_pass={overall_pass} wrote {out_path}")
    if not overall_pass:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
