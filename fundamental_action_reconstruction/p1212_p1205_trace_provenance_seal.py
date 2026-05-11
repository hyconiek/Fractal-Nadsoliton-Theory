#!/usr/bin/env python3
"""P1212: verify P1205 trace payload hash integrity (provenance seal)."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"


def canonical_trace_bytes(trace_payload: dict) -> bytes:
    return json.dumps(trace_payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode("utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1205", type=Path, default=GEN / "p1205_w4_sympy_cas_runner_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1212_p1205_trace_provenance_seal_summary.json")
    args = parser.parse_args()

    p1205 = json.loads(args.p1205.read_text(encoding="utf-8"))
    trace = p1205.get("trace_payload") if isinstance(p1205.get("trace_payload"), dict) else None
    declared_hash = p1205.get("trace_hash_sha256") if isinstance(p1205.get("trace_hash_sha256"), str) else ""

    trace_present = trace is not None
    recomputed_hash = hashlib.sha256(canonical_trace_bytes(trace)).hexdigest() if trace_present else ""
    hash_format_ok = len(declared_hash) == 64 and all(ch in "0123456789abcdef" for ch in declared_hash.lower())
    hash_match = bool(trace_present and hash_format_ok and declared_hash.lower() == recomputed_hash)

    out = {
        "packet": "P1212",
        "as_of": "2026-05-11",
        "input_path": str(args.p1205),
        "trace_payload_present": trace_present,
        "declared_trace_hash_sha256": declared_hash,
        "recomputed_trace_hash_sha256": recomputed_hash,
        "declared_hash_format_ok": hash_format_ok,
        "trace_hash_match": hash_match,
        "provenance_seal_pass": hash_match,
        "strict_closure_claim_allowed": False,
        "note": "Integrity-only seal for external symbolic trace payload.",
    }

    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1212] provenance_seal_pass={hash_match} wrote {args.out}")


if __name__ == "__main__":
    main()
