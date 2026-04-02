#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-27"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1062 = GENERATED / "p1062_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_admission_probe_summary.json"
OUT_JSON = GENERATED / "f961_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_packet.json"
OUT_SUMMARY = GENERATED / "f961_current_strict_t173_t176_existing_t183_negative_3_cycle_sign_holonomy_obstruction_reference_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1062.exists():
        artifact = {
            "stage": "F961",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1062.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1062 = load_json(IN_P1062)
    admitted = bool(p1062.get("admissible_as_reference_only") is True)
    status = (
        "F961_EXECUTED_CURRENT_STRICT_T173_T176_EXISTING_T183_NEGATIVE_3_CYCLE_SIGN_HOLONOMY_OBSTRUCTION_REFERENCE_PACKET_NO_FALSE_PASS"
        if admitted
        else "F961_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_REFERENCE_ADMISSION_STATE"
    )

    artifact = {
        "stage": "F961",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1062_reference_admission_probe_summary": str(IN_P1062.relative_to(REPO)),
        },
        "exported_reference": {
            "object_id": "Negative3CycleSignHolonomyObstructionReference_global_C_v1_strict_v1",
            "grade": "obstruction_reference_only",
            "supports_search_for": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
            "typed_meaning": (
                "one strict reference-only obstruction object freezing the currently exported negative 3-cycle witness as "
                "a search clue for a future inversion-sensitive source-side provider class on the active T183 frontier"
            ),
        },
        "admission_properties": {
            "triangle_witness_present": True,
            "counts_as_strict_physical_orientation_datum": False,
            "counts_as_strict_berry_holonomy_primitive": False,
            "counts_as_lawful_supplier": False,
            "counts_as_reference_only": admitted,
        },
        "hard_limits": [
            "Does not discharge T183.",
            "Does not discharge T176.",
            "Does not export a strict physical orientation datum.",
            "Does not export a strict Berry/holonomy primitive.",
            "Does not discharge QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "reference_object_id": artifact["exported_reference"]["object_id"],
        "reference_grade": artifact["exported_reference"]["grade"],
        "supports_search_for": artifact["exported_reference"]["supports_search_for"],
        "counts_as_strict_physical_orientation_datum": False,
        "counts_as_lawful_supplier": False,
        "counts_as_reference_only": admitted,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
