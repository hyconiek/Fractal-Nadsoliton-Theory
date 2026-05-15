#!/usr/bin/env python3
"""P1792 S742 strict execution evidence intake schema checkpoint."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"


def load(name: str) -> dict:
    p = GENERATED / name
    if not p.exists():
        return {"_missing": name}
    return json.loads(p.read_text(encoding="utf-8"))


def main() -> None:
    p1790 = load("p1790_s740_strict_sv1_to_sv5_coordinated_run_pack_contract_checkpoint.json")
    p1791 = load("p1791_s741_strict_sv1_to_sv5_ledger_completeness_and_verdict_protocol_checkpoint.json")

    required_fields = [
        "ledger_id",
        "produced_by",
        "background_family_id",
        "index_convention_id",
        "boundary_clause_id",
        "component_basis",
        "result_kind",
        "residual_vector",
        "obstruction_tags",
        "timestamp_utc",
    ]

    out = {
        "checkpoint_id": "P1792_S742",
        "title": "STRICT_EXECUTION_EVIDENCE_INTAKE_SCHEMA",
        "input_reuse": [
            "p1790_s740_strict_sv1_to_sv5_coordinated_run_pack_contract_checkpoint.json",
            "p1791_s741_strict_sv1_to_sv5_ledger_completeness_and_verdict_protocol_checkpoint.json",
        ],
        "applies_to_ledgers": p1791.get("required_ledgers", []),
        "required_fields": required_fields,
        "validation_rules": {
            "missing_required_field": "INVALID_LEDGER_SCHEMA",
            "freeze_mismatch_background_or_index_or_boundary": "FREEZE_MISMATCH",
            "pass_zero_without_explicit_residual_vector": "INVALID_PASS_CLAIM",
            "allowed_result_kind": ["PASS_ZERO", "OPEN_OBSTRUCTION_WITH_TRACE"],
        },
        "intake_outcome_mapping": {
            "schema_invalid": "OPEN_OBSTRUCTION_WITH_TRACE",
            "freeze_mismatch": "OPEN_OBSTRUCTION_WITH_TRACE",
            "nonzero_or_obstructed_residual": "OPEN_OBSTRUCTION_WITH_TRACE",
            "all_valid_and_zero": "PASS_ZERO",
        },
        "dependency_snapshot": {
            "p1790_status": p1790.get("status", "MISSING"),
            "p1791_status": p1791.get("status", "MISSING"),
        },
        "next_honest_step": "Implement intake validator that checks required fields and freeze coherence before passing evidence into P1791 verdict protocol.",
        "status": "INTAKE_SCHEMA_DEFINED_VALIDATOR_PENDING",
    }

    GENERATED.mkdir(parents=True, exist_ok=True)
    out_path = GENERATED / "p1792_s742_strict_execution_evidence_intake_schema_checkpoint.json"
    out_path.write_text(json.dumps(out, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(out_path)


if __name__ == "__main__":
    main()
