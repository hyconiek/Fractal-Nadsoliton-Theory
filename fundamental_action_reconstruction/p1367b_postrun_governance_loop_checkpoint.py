#!/usr/bin/env python3
from __future__ import annotations

import json
import subprocess
from pathlib import Path

FAR = Path(__file__).resolve().parent
GEN = FAR / "generated"


def run(cmd: list[str]) -> dict:
    p = subprocess.run(cmd, capture_output=True, text=True)
    return {"cmd": " ".join(cmd), "returncode": p.returncode, "stdout": p.stdout.strip(), "stderr": p.stderr.strip()}


def main() -> None:
    runs = [
        run(["python3", str(FAR / "p1362_strict_candidate_residual_benchmark_checkpoint.py")]),
        run(["python3", str(FAR / "p1363_upgrade_gate_and_blocker_registry_checkpoint.py")]),
        run(["python3", str(FAR / "p1365_verify_artifact_tracking_checkpoint.py")]),
    ]
    ok = all(r["returncode"] == 0 for r in runs)
    out = {
        "packet": "P1367b",
        "as_of": "2026-05-12",
        "runs": runs,
        "loop_status": "PASS" if ok else "FAIL",
    }
    path = GEN / "p1367b_postrun_governance_loop_summary.json"
    path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1367b] wrote {path}; loop_status={out['loop_status']}")


if __name__ == "__main__":
    main()
