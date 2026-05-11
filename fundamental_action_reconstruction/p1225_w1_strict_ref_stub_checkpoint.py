#!/usr/bin/env python3
"""P1225: minimal strict-reference stub run for W1 continuity."""
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
    parser.add_argument("--stub-input", type=Path, default=GEN / "p1225_w1_input_payload_with_strict_ref_stub.json")
    parser.add_argument("--p1223-out", type=Path, default=GEN / "p1225_w1_execution_summary.json")
    parser.add_argument("--p1224-out", type=Path, default=GEN / "p1225_w1_binding_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1225_w1_strict_ref_stub_checkpoint_summary.json")
    args = parser.parse_args()

    payload = _load(args.base_input)
    available = payload.setdefault("available_inputs", {})
    available["strict_selector_source_theorem_exported"] = True
    payload["strict_selector_source_theorem_ref"] = "internal://w1/strict_selector_source_theorem_stub_v1"
    _dump(args.stub_input, payload)

    root = Path(__file__).resolve().parent
    subprocess.run([
        "python3", str(root / "p1223_w1_execution_checkpoint.py"),
        "--input", str(args.stub_input), "--out", str(args.p1223_out)
    ], check=True)
    subprocess.run([
        "python3", str(root / "p1224_w1_artifact_binding_checkpoint.py"),
        "--input", str(args.stub_input), "--p1223", str(args.p1223_out), "--out", str(args.p1224_out)
    ], check=True)

    p1223 = _load(args.p1223_out)
    p1224 = _load(args.p1224_out)

    out = {
        "packet": "P1225",
        "as_of": "2026-05-11",
        "strict_ref_stub_applied": True,
        "w1_witness_status_after_stub": p1223.get("witness_status"),
        "w1_discharge_mode_after_stub": p1223.get("evaluation", {}).get("discharge_mode"),
        "strict_binding_ok_after_stub": p1224.get("strict_binding_ok"),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Pipeline continuity run with minimal strict theorem reference stub only.",
    }
    _dump(args.out, out)
    print(f"[P1225] witness={out['w1_witness_status_after_stub']} binding={out['strict_binding_ok_after_stub']} wrote {args.out}")


if __name__ == "__main__":
    main()
