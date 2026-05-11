#!/usr/bin/env python3
"""P1150 strict premise admissibility gate.

Single gate for candidate premises before running selector probes.
Checks:
- strict-side provenance declaration
- explicit noncyclic anchor declaration
- explicit no-legacy-role-transfer declaration
- explicit non-closure / no-QW2191-discharge declaration

Input: JSON file path with fields listed in REQUIRED_KEYS.
Output: generated/p1150_strict_premise_admissibility_gate_summary.json
Exit code 0 only if all checks pass.
"""
from __future__ import annotations
import json
from pathlib import Path
import sys

REQUIRED_KEYS = {
    "premise_id": str,
    "strict_side_provenance": bool,
    "noncyclic_anchor_declared": bool,
    "no_legacy_role_transfer": bool,
    "no_closure_claim": bool,
    "no_qw2191_discharge_claim": bool,
}


def validate(payload: dict) -> dict:
    checks = {}
    for k, t in REQUIRED_KEYS.items():
        checks[f"has_{k}"] = k in payload
        checks[f"type_{k}_ok"] = isinstance(payload.get(k), t) if k in payload else False

    semantic = {
        "strict_side_provenance_true": payload.get("strict_side_provenance") is True,
        "noncyclic_anchor_true": payload.get("noncyclic_anchor_declared") is True,
        "no_legacy_role_transfer_true": payload.get("no_legacy_role_transfer") is True,
        "no_closure_claim_true": payload.get("no_closure_claim") is True,
        "no_qw2191_discharge_claim_true": payload.get("no_qw2191_discharge_claim") is True,
    }

    overall = all(checks.values()) and all(semantic.values())
    return {"schema_checks": checks, "semantic_checks": semantic, "overall_pass": overall}


def main() -> None:
    if len(sys.argv) != 2:
        print("usage: p1150_strict_premise_admissibility_gate.py <premise_json>")
        raise SystemExit(2)

    in_path = Path(sys.argv[1]).resolve()
    payload = json.loads(in_path.read_text(encoding="utf-8"))
    result = validate(payload)
    out = {
        "packet": "P1150",
        "as_of": "2026-05-10",
        "input_file": str(in_path),
        "premise_id": payload.get("premise_id"),
        **result,
    }

    out_path = Path(__file__).resolve().parent / "generated" / "p1150_strict_premise_admissibility_gate_summary.json"
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1150] overall_pass={out['overall_pass']} wrote {out_path}")
    if not out["overall_pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
