#!/usr/bin/env python3
"""P1227: run W1 pipeline with a non-stub strict reference candidate."""
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
    parser.add_argument(
        "--candidate-ref",
        default="internal://w1/strict_selector_source_theorem_candidate_v1",
    )
    parser.add_argument("--input", type=Path, default=GEN / "p1227_w1_input_payload_real_candidate_ref.json")
    parser.add_argument("--p1223-out", type=Path, default=GEN / "p1227_w1_execution_summary.json")
    parser.add_argument("--p1224-out", type=Path, default=GEN / "p1227_w1_binding_summary.json")
    parser.add_argument("--p1226-out", type=Path, default=GEN / "p1227_w1_reference_tier_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1227_w1_real_candidate_run_summary.json")
    args = parser.parse_args()

    payload = _load(args.base_input)
    payload.setdefault("available_inputs", {})["strict_selector_source_theorem_exported"] = True
    payload["strict_selector_source_theorem_ref"] = args.candidate_ref
    _dump(args.input, payload)

    root = Path(__file__).resolve().parent
    subprocess.run(["python3", str(root / "p1223_w1_execution_checkpoint.py"), "--input", str(args.input), "--out", str(args.p1223_out)], check=True)
    subprocess.run(["python3", str(root / "p1224_w1_artifact_binding_checkpoint.py"), "--input", str(args.input), "--p1223", str(args.p1223_out), "--out", str(args.p1224_out)], check=True)
    subprocess.run(["python3", str(root / "p1226_w1_reference_tier_checkpoint.py"), "--input", str(args.input), "--p1223-out", str(args.p1223_out), "--p1224-out", str(args.p1224_out), "--out", str(args.p1226_out)], check=True)

    p1223 = _load(args.p1223_out)
    p1224 = _load(args.p1224_out)
    p1226 = _load(args.p1226_out)

    out = {
        "packet": "P1227",
        "as_of": "2026-05-11",
        "candidate_ref": args.candidate_ref,
        "w1_witness_status": p1223.get("witness_status"),
        "w1_discharge_mode": p1223.get("evaluation", {}).get("discharge_mode"),
        "strict_binding_ok": p1224.get("strict_binding_ok"),
        "strict_reference_tier": p1226.get("strict_reference_tier"),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Real-candidate continuity run; no global closure claim.",
    }
    _dump(args.out, out)
    print(f"[P1227] tier={out['strict_reference_tier']} witness={out['w1_witness_status']} wrote {args.out}")


if __name__ == "__main__":
    main()
