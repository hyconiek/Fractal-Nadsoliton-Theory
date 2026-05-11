#!/usr/bin/env python3
"""P1201: trace-hash probe for W4 adapter auditability."""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"


def main() -> None:
    p1200 = json.loads((GEN / "p1200_w4_backend_adapter_minimal_impl_summary.json").read_text(encoding="utf-8"))

    trace_payload = {
        "step": "expand",
        "expr": "a*psi4**2 - a*psi4**2",
        "ok": False,
        "reason": "backend_not_integrated",
    }
    trace_json = json.dumps(trace_payload, sort_keys=True)
    trace_hash = hashlib.sha256(trace_json.encode("utf-8")).hexdigest()

    out = {
        "packet": "P1201",
        "as_of": "2026-05-10",
        "input_methods_present": bool(p1200.get("methods_present", False)),
        "trace_payload": trace_payload,
        "trace_hash_sha256": trace_hash,
        "auditability_ready": True,
        "symbolic_execution_ready": False,
        "strict_closure_claim_allowed": False,
        "note": "Trace hashing ready even though symbolic execution remains unavailable.",
    }

    out_path = GEN / "p1201_w4_trace_hash_probe_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1201] auditability_ready=True symbolic_execution_ready=False wrote {out_path}")


if __name__ == "__main__":
    main()
