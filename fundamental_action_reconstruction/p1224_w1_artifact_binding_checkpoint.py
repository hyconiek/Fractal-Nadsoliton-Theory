#!/usr/bin/env python3
"""P1224: minimal artifact-binding check for W1 discharge continuity."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=GEN / "p1223_w1_input_payload.json")
    parser.add_argument("--p1223", type=Path, default=GEN / "p1223_w1_execution_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1224_w1_artifact_binding_summary.json")
    args = parser.parse_args()

    p1223_input = _load(args.input)
    p1223 = _load(args.p1223)

    available = p1223_input.get("available_inputs", {})
    strict_theorem_flag = bool(available.get("strict_selector_source_theorem_exported", False))

    theorem_ref = p1223_input.get("strict_selector_source_theorem_ref")
    theorem_hash = p1223_input.get("strict_selector_source_theorem_sha256")

    # Single-researcher mode: require only explicit theorem reference when strict flag is raised.
    strict_binding_ok = (not strict_theorem_flag) or bool(theorem_ref)

    out = {
        "packet": "P1224",
        "as_of": "2026-05-11",
        "witness_id": p1223_input.get("witness_id", "W1_selector_uniqueness_bridge"),
        "continuation_from": "P1223",
        "p1223_witness_status": p1223.get("witness_status"),
        "strict_theorem_flag": strict_theorem_flag,
        "strict_binding_ok": strict_binding_ok,
        "binding_rule": "if strict_selector_source_theorem_exported=true then theorem_ref is required (sha256 optional in single-researcher mode)",
        "strict_selector_source_theorem_ref_present": bool(theorem_ref),
        "strict_selector_source_theorem_sha256_present": bool(theorem_hash),
        "strict_closure_claim_allowed": False,
        "theory_closure_status": "OPEN",
        "note": "Continuity checkpoint only; does not alter prior P1223 verdict.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1224] strict_binding_ok={strict_binding_ok} wrote {args.out}")


if __name__ == "__main__":
    main()
