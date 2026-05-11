#!/usr/bin/env python3
"""P1219: execute a concrete controlled W4 run using a minimal ready P1205 artifact."""
from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def canonical_trace_bytes(trace_payload: dict) -> bytes:
    return json.dumps(trace_payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def run_cmd(args: list[str]) -> None:
    subprocess.run(args, check=True, cwd=ROOT.parent)


def read_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--out", type=Path, default=GEN / "p1219_w4_execution_checkpoint_summary.json")
    args = parser.parse_args()

    trace_payload = {
        "backend": "external_cas_mock_minimal",
        "reduced_zero_ok": False,
        "witness": "W4",
        "note": "Execution checkpoint artifact; not a discharge certificate.",
    }
    trace_hash = hashlib.sha256(canonical_trace_bytes(trace_payload)).hexdigest()

    p1205_path = GEN / "p1219_p1205_minimal_ready_artifact.json"
    p1205 = {
        "packet": "P1205_EXTERNAL_MINIMAL_READY",
        "as_of": "2026-05-11",
        "sympy_available": True,
        "trace_payload": trace_payload,
        "trace_hash_sha256": trace_hash,
        "strict_closure_claim_allowed": False,
    }
    p1205_path.write_text(json.dumps(p1205, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    gate_path = GEN / "p1219_p1209_gate_summary.json"
    seal_path = GEN / "p1219_p1212_seal_summary.json"
    attempt_path = GEN / "p1219_p1210_attempt_summary.json"

    run_cmd(["python3", str(ROOT / "p1209_symbolic_attempt_gate_open.py"), "--p1205-only", "--p1205", str(p1205_path), "--out", str(gate_path)])
    run_cmd(["python3", str(ROOT / "p1212_p1205_trace_provenance_seal.py"), "--p1205", str(p1205_path), "--out", str(seal_path)])
    run_cmd(["python3", str(ROOT / "p1210_controlled_w4_symbolic_attempt.py"), "--p1205", str(p1205_path), "--p1209", str(gate_path), "--p1212", str(seal_path), "--out", str(attempt_path)])

    gate = read_json(gate_path)
    attempt = read_json(attempt_path)

    out = {
        "packet": "P1219",
        "as_of": "2026-05-11",
        "artifact_path": str(p1205_path),
        "gate_mode": gate.get("gate_mode"),
        "symbolic_attempt_gate_open": gate.get("symbolic_attempt_gate_open"),
        "attempt_executed": attempt.get("attempt_executed"),
        "w4_status": attempt.get("w4_status"),
        "w4_discharge_pass": attempt.get("w4_discharge_pass"),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Controlled runtime checkpoint for W4 execution path.",
    }
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1219] wrote {args.out}")


if __name__ == "__main__":
    main()
