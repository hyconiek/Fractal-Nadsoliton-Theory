#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path

AS_OF = "2026-03-29"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1092 = GENERATED / "p1092_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_audit_probe_summary.json"
OUT_JSON = GENERATED / "f967_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_packet.json"
OUT_SUMMARY = GENERATED / "f967_current_strict_t173_t176_time_arrow_light_resonance_oriented_dephasing_competing_extension_packet_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1092.exists():
        artifact = {
            "stage": "F967",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1092.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1092 = load_json(IN_P1092)
    admitted = p1092.get("status") == "PASS_CURRENT_STRICT_T173_T176_TIME_ARROW_LIGHT_RESONANCE_ORIENTED_DEPHASING_COMPETING_EXTENSION_AUDITED"
    status = (
        "F967_EXECUTED_CURRENT_STRICT_T173_T176_TIME_ARROW_LIGHT_RESONANCE_ORIENTED_DEPHASING_COMPETING_EXTENSION_PACKET_NO_FALSE_PASS"
        if admitted
        else "F967_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_HYPOTHESIS_STATE"
    )

    artifact = {
        "stage": "F967",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1092_hypothesis_audit_summary": str(IN_P1092.relative_to(REPO)),
        },
        "exported_hypothesis_packet": {
            "object_id": "TimeArrowLightResonanceOrientedDephasingCompetingExtensionHypothesis_v1",
            "hypothesis_grade": p1092.get("hypothesis_grade"),
            "closest_existing_repo_lane": p1092.get("closest_existing_repo_lane"),
            "preferred_candidate_form": p1092.get("preferred_candidate_form"),
            "active_missing_bridge": p1092.get("active_missing_bridge"),
            "typed_meaning": "one frontier-only packet freezing that time-arrow plus light-resonance plus dephasing is currently admissible only as a competing extension hypothesis, with the strongest honest future form being an oriented/nonreciprocal dephasing-or-retardation operator",
        },
        "admission_properties": {
            "plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker": bool(p1092.get("plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker")),
            "requires_additional_oriented_transport_or_selector_slot": bool(p1092.get("requires_additional_oriented_transport_or_selector_slot")),
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "hard_limits": [
            "Does not export an actual oriented dephasing operator.",
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
        "hypothesis_object_id": artifact["exported_hypothesis_packet"]["object_id"],
        "hypothesis_grade": artifact["exported_hypothesis_packet"]["hypothesis_grade"],
        "preferred_candidate_form": artifact["exported_hypothesis_packet"]["preferred_candidate_form"],
        "plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker": artifact["admission_properties"]["plain_time_constant_dephasing_scalar_counts_as_symmetry_breaker"],
        "requires_additional_oriented_transport_or_selector_slot": artifact["admission_properties"]["requires_additional_oriented_transport_or_selector_slot"],
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
