#!/usr/bin/env python3
"""P1793 S743 strict SV1->SV5 evidence intake validator."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
SCHEMA_PATH = GENERATED / "p1792_s742_strict_execution_evidence_intake_schema_checkpoint.json"
DEFAULT_INPUT = GENERATED / "p1793_s743_sv1_sv5_evidence_intake_input.json"
DEFAULT_OUTPUT = GENERATED / "p1793_s743_sv1_sv5_evidence_intake_validation_result.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def validate(schema: dict[str, Any], payload: dict[str, Any]) -> dict[str, Any]:
    required_fields = schema.get("required_fields", [])
    allowed_result_kind = set(schema.get("validation_rules", {}).get("allowed_result_kind", []))

    ledgers = payload.get("ledgers", [])
    tags: list[str] = []
    issues: list[dict[str, Any]] = []

    freeze_triplets = set()
    result_kinds: list[str] = []

    for i, ledger in enumerate(ledgers):
        missing = [f for f in required_fields if f not in ledger]
        if missing:
            tags.append("INVALID_LEDGER_SCHEMA")
            issues.append({"ledger_index": i, "type": "missing_fields", "fields": missing})
            continue

        rk = ledger.get("result_kind")
        result_kinds.append(rk)
        if rk not in allowed_result_kind:
            tags.append("INVALID_RESULT_KIND")
            issues.append({"ledger_index": i, "type": "invalid_result_kind", "value": rk})

        if rk == "PASS_ZERO" and "residual_vector" not in ledger:
            tags.append("INVALID_PASS_CLAIM")
            issues.append({"ledger_index": i, "type": "pass_without_residual_vector"})

        freeze_triplets.add(
            (
                ledger.get("background_family_id"),
                ledger.get("index_convention_id"),
                ledger.get("boundary_clause_id"),
            )
        )

    if len(freeze_triplets) > 1:
        tags.append("FREEZE_MISMATCH")
        issues.append({"type": "freeze_mismatch", "distinct_freezes": len(freeze_triplets)})

    unique_tags = sorted(set(tags))

    if not ledgers:
        unique_tags = sorted(set(unique_tags + ["EMPTY_LEDGER_SET"]))

    all_pass_zero = bool(result_kinds) and all(rk == "PASS_ZERO" for rk in result_kinds)
    if not all_pass_zero:
        unique_tags = sorted(set(unique_tags + ["NON_PASS_LEDGER_PRESENT"]))

    verdict = "PASS_ZERO" if (not unique_tags and all_pass_zero) else "OPEN_OBSTRUCTION_WITH_TRACE"

    return {
        "verdict": verdict,
        "obstruction_tags": unique_tags,
        "issues": issues,
        "ledger_count": len(ledgers),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    schema = load_json(SCHEMA_PATH)
    payload = load_json(args.input)
    result = validate(schema, payload)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(args.output)


if __name__ == "__main__":
    main()
