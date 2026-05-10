#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-28"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1090 = GENERATED / "p1090_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_family_failure_map_audit_probe_summary.json"
OUT_JSON = GENERATED / "f965_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_failure_map_packet.json"
OUT_SUMMARY = GENERATED / "f965_current_strict_t173_t176_nadsoliton_information_process_to_orientation_supplier_failure_map_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1090.exists():
        artifact = {
            "stage": "F965",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1090.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1090 = load_json(IN_P1090)
    admitted = bool(p1090.get("process_reading_of_nadsoliton_is_repo_real"))
    status = (
        "F965_EXECUTED_CURRENT_STRICT_T173_T176_NADSOLITON_INFORMATION_PROCESS_TO_ORIENTATION_SUPPLIER_FAILURE_MAP_PACKET_NO_FALSE_PASS"
        if admitted
        else "F965_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_FAILURE_MAP_STATE"
    )

    artifact = {
        "stage": "F965",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1090_failure_map_audit_probe_summary": str(IN_P1090.relative_to(REPO)),
        },
        "exported_failure_map": {
            "object_id": "NadsolitonInformationProcessToOrientationSupplierFailureMap_global_C_v1_strict_v1",
            "grade": "frontier_failure_map_only",
            "supports_search_for": "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1",
            "typed_meaning": (
                "one strict frontier-only map freezing the exact repo-tested families that already tried to move from "
                "nadsoliton as informational process to lawful orientation supplier, together with the exact failure boundary of each family"
            ),
        },
        "admission_properties": {
            "process_reading_of_nadsoliton_is_repo_real": admitted,
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
            "tested_family_count": int(p1090.get("tested_family_count", 0)),
        },
        "hard_limits": [
            "Does not export a lawful process-to-orientation supplier.",
            "Does not discharge T183.",
            "Does not discharge T176.",
            "Does not discharge QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "failure_map_object_id": artifact["exported_failure_map"]["object_id"],
        "failure_map_grade": artifact["exported_failure_map"]["grade"],
        "supports_search_for": artifact["exported_failure_map"]["supports_search_for"],
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "process_reading_of_nadsoliton_is_repo_real": admitted,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
