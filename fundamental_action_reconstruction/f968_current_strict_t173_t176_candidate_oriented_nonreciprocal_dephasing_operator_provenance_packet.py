#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1093 = GENERATED / "p1093_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_audit_probe_summary.json"
OUT_JSON = GENERATED / "f968_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_packet.json"
OUT_SUMMARY = GENERATED / "f968_current_strict_t173_t176_candidate_oriented_nonreciprocal_dephasing_operator_provenance_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1093.exists():
        artifact = {
            "stage": "F968",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1093.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1093 = load_json(IN_P1093)
    admitted = p1093.get("status") == "PASS_CURRENT_STRICT_T173_T176_CANDIDATE_ORIENTED_NONRECIPROCAL_DEPHASING_OPERATOR_PROVENANCE_AUDITED"
    status = (
        "F968_EXECUTED_CURRENT_STRICT_T173_T176_CANDIDATE_ORIENTED_NONRECIPROCAL_DEPHASING_OPERATOR_PROVENANCE_PACKET_NO_FALSE_PASS"
        if admitted
        else "F968_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_PROVENANCE_STATE"
    )

    artifact = {
        "stage": "F968",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1093_provenance_audit_summary": str(IN_P1093.relative_to(REPO)),
        },
        "exported_provenance_packet": {
            "object_id": "CandidateOrientedNonreciprocalDephasingOperatorProvenanceMap_v1",
            "repo_rooted_component_count": p1093.get("repo_rooted_component_count"),
            "genuine_new_import_component_count": p1093.get("genuine_new_import_component_count"),
            "complete_candidate_operator_exported": p1093.get("complete_candidate_operator_exported"),
            "minimal_new_import_boundary": p1093.get("minimal_new_import_boundary"),
            "active_missing_bridge": p1093.get("active_missing_bridge"),
            "typed_meaning": "one provenance packet freezing which parts of the candidate oriented/nonreciprocal dephasing operator are already repo-rooted motifs and which parts are still genuine new imports",
        },
        "admission_properties": {
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "hard_limits": [
            "Does not export the full operator.",
            "Does not export a lawful supplier.",
            "Does not export a strict physical orientation datum.",
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
        "provenance_object_id": artifact["exported_provenance_packet"]["object_id"],
        "repo_rooted_component_count": artifact["exported_provenance_packet"]["repo_rooted_component_count"],
        "genuine_new_import_component_count": artifact["exported_provenance_packet"]["genuine_new_import_component_count"],
        "complete_candidate_operator_exported": artifact["exported_provenance_packet"]["complete_candidate_operator_exported"],
        "minimal_new_import_boundary": artifact["exported_provenance_packet"]["minimal_new_import_boundary"],
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
