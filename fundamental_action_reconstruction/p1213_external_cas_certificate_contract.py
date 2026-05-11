#!/usr/bin/env python3
"""P1213: enforce minimal external CAS certificate contract for P1205 artifacts."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

GEN = Path(__file__).resolve().parent / "generated"
REQUIRED_TEXT = ["cas_engine", "cas_engine_version", "deterministic_replay_seed", "trace_schema_version"]


def _is_nonempty_str(v: object) -> bool:
    return isinstance(v, str) and bool(v.strip())


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--p1205", type=Path, default=GEN / "p1205_w4_sympy_cas_runner_summary.json")
    parser.add_argument("--out", type=Path, default=GEN / "p1213_external_cas_certificate_contract_summary.json")
    args = parser.parse_args()

    data = json.loads(args.p1205.read_text(encoding="utf-8"))

    missing_text = [k for k in REQUIRED_TEXT if not _is_nonempty_str(data.get(k))]
    hash_value = data.get("trace_hash_sha256")
    hash_format_ok = isinstance(hash_value, str) and len(hash_value) == 64 and all(c in "0123456789abcdef" for c in hash_value.lower())
    as_of_ok = data.get("as_of") == "2026-05-11"
    strict_guard_ok = data.get("strict_closure_claim_allowed") is False

    contract_pass = bool(not missing_text and hash_format_ok and as_of_ok and strict_guard_ok)

    out = {
        "packet": "P1213",
        "as_of": "2026-05-11",
        "input_path": str(args.p1205),
        "missing_or_empty_text_fields": missing_text,
        "trace_hash_format_ok": hash_format_ok,
        "as_of_ok": as_of_ok,
        "strict_guard_ok": strict_guard_ok,
        "external_certificate_contract_pass": contract_pass,
        "strict_closure_claim_allowed": False,
        "note": "Minimal external CAS provenance/replay contract gate.",
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1213] external_certificate_contract_pass={contract_pass} wrote {args.out}")


if __name__ == "__main__":
    main()
