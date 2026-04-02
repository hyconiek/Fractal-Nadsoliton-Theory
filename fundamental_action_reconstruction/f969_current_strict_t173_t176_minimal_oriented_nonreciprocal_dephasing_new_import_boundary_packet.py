#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1094 = GENERATED / "p1094_current_strict_t173_t176_minimal_oriented_nonreciprocal_dephasing_new_import_boundary_admission_probe_summary.json"
OUT_JSON = GENERATED / "f969_current_strict_t173_t176_minimal_oriented_nonreciprocal_dephasing_new_import_boundary_packet.json"
OUT_SUMMARY = GENERATED / "f969_current_strict_t173_t176_minimal_oriented_nonreciprocal_dephasing_new_import_boundary_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1094.exists():
        artifact = {
            "stage": "F969",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1094.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1094 = load_json(IN_P1094)
    admitted = p1094.get("status") == "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ORIENTED_NONRECIPROCAL_DEPHASING_NEW_IMPORT_BOUNDARY_ADMITTED"
    status = (
        "F969_EXECUTED_CURRENT_STRICT_T173_T176_MINIMAL_ORIENTED_NONRECIPROCAL_DEPHASING_NEW_IMPORT_BOUNDARY_PACKET_NO_FALSE_PASS"
        if admitted
        else "F969_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_BOUNDARY_STATE"
    )

    artifact = {
        "stage": "F969",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1094_boundary_admission_summary": str(IN_P1094.relative_to(REPO)),
        },
        "exported_boundary_packet": {
            "object_id": "MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1",
            "grade": p1094.get("boundary_grade"),
            "minimal_new_import_boundary": p1094.get("minimal_new_import_boundary"),
            "frontier_context_bridge": p1094.get("frontier_context_bridge"),
            "exact_reduction_to_frontier_context_bridge_exported": False,
            "typed_meaning": "one explicit irreducible new-import boundary for the oriented/nonreciprocal dephasing lane, frozen only as a candidate provider-class seed",
        },
        "admission_properties": {
            "admissible_as_candidate_provider_class_seed_only": bool(p1094.get("admissible_as_candidate_provider_class_seed_only")),
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "hard_limits": [
            "Does not export the full operator.",
            "Does not export a lawful supplier.",
            "Does not export a solution.",
            "Does not export a strict physical orientation datum.",
            "Does not export exact reduction to the frontier context bridge.",
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
        "boundary_object_id": artifact["exported_boundary_packet"]["object_id"],
        "boundary_grade": artifact["exported_boundary_packet"]["grade"],
        "minimal_new_import_boundary": artifact["exported_boundary_packet"]["minimal_new_import_boundary"],
        "frontier_context_bridge": artifact["exported_boundary_packet"]["frontier_context_bridge"],
        "exact_reduction_to_frontier_context_bridge_exported": False,
        "admissible_as_candidate_provider_class_seed_only": artifact["admission_properties"]["admissible_as_candidate_provider_class_seed_only"],
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
