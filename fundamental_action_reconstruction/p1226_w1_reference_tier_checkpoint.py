#!/usr/bin/env python3
"""P1226: classify W1 strict-reference tier and rerun continuity checks."""
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


def _tier_from_ref(ref: str) -> str:
    if "_stub_" in ref or ref.endswith("_stub_v1"):
        return "STUB"
    return "REAL_CANDIDATE"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=GEN / "p1225_w1_input_payload_with_strict_ref_stub.json")
    parser.add_argument("--p1223-out", type=Path, default=GEN / "p1226_w1_execution_summary.json")
    parser.add_argument("--p1224-out", type=Path, default=GEN / "p1226_w1_binding_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1226_w1_reference_tier_summary.json")
    args = parser.parse_args()

    payload = _load(args.input)
    ref = str(payload.get("strict_selector_source_theorem_ref", ""))
    tier = _tier_from_ref(ref) if ref else "MISSING"

    root = Path(__file__).resolve().parent
    subprocess.run(["python3", str(root / "p1223_w1_execution_checkpoint.py"), "--input", str(args.input), "--out", str(args.p1223_out)], check=True)
    subprocess.run(["python3", str(root / "p1224_w1_artifact_binding_checkpoint.py"), "--input", str(args.input), "--p1223", str(args.p1223_out), "--out", str(args.p1224_out)], check=True)

    p1223 = _load(args.p1223_out)
    p1224 = _load(args.p1224_out)

    out = {
        "packet": "P1226",
        "as_of": "2026-05-11",
        "strict_reference_tier": tier,
        "strict_reference_value": ref,
        "w1_witness_status": p1223.get("witness_status"),
        "w1_discharge_mode": p1223.get("evaluation", {}).get("discharge_mode"),
        "strict_binding_ok": p1224.get("strict_binding_ok"),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Tier-classification checkpoint; does not claim global closure.",
    }
    _dump(args.out, out)
    print(f"[P1226] tier={tier} witness={out['w1_witness_status']} wrote {args.out}")


if __name__ == "__main__":
    main()
