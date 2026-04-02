#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P808 = GENERATED / "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.json"
IN_F807 = GENERATED / "f807_current_strict_alpha_s_provider_class_shift_requirement_packet.json"

OUT = GENERATED / "f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.json"
OUT_SUMMARY = GENERATED / "f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P808, IN_F807]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F808",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p808 = load_json(IN_P808)
    f807 = load_json(IN_F807)

    if (
        p808.get("status")
        == "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_ALPHA_S_DOMAIN_INTERFACE_BLOCKED"
        and (p808.get("theorem_result") or {}).get(
            "strict_source_shannon_provider_shift_candidate_reference_lane_admitted"
        )
        is True
        and (p808.get("theorem_result") or {}).get("alpha_s_domain_interface_exported") is False
        and f807.get("status")
        == "F807_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
    ):
        status = "F808_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
    else:
        status = "F808_REQUIRES_REVIEW"

    candidate = p808.get("provider_shift_candidate_reference_lane") or {}
    shift_requirement = f807.get("exported_object") or {}

    artifact = {
        "stage": "F808",
        "packet_name": "CurrentStrictAlphaSStrictSourceShannonProviderShiftCandidateReference_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p808_provider_shift_candidate_reference_probe": rel(IN_P808),
            "f807_provider_class_shift_requirement_packet": rel(IN_F807),
        },
        "exported_object": {
            "object_id": "alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_v1",
            "goal": "Export the admitted strict-source Shannon provider-shift candidate reference lane for alpha_s while keeping alpha_s-domain interface blocked.",
            "provider_class_shift_requirement_ref": shift_requirement.get("object_id"),
            "candidate_reference_lane_id": candidate.get("candidate_id"),
            "source_upgrade_ref": candidate.get("source_upgrade_ref"),
            "candidate_grade": candidate.get("candidate_grade"),
            "alpha_s_domain_interface_status": candidate.get("alpha_s_domain_interface_status"),
            "supporting_candidate_lane_refs": candidate.get("supporting_candidate_lane_refs"),
            "scope": "strict_alpha_s_provider_shift_candidate_reference_only",
        },
        "current_honest_reading": [
            "The repo now exports one explicit strict-source Shannon provider-shift candidate reference lane for alpha_s.",
            "That object only records a candidate reference lane and its supporting strict-source evidence.",
            "It still does not export any alpha_s-domain interface and therefore does not realize provider shift.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether this strict-source Shannon provider-shift candidate reference lane exports any exact alpha_s-domain interface target without silent carrier identification.",
        "hard_limits": [
            "Does not claim that provider shift has already succeeded.",
            "Does not claim that strict-source Shannon already enters the alpha_s domain.",
            "Does not claim that alpha_s semantics are already supplied.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F808",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "alpha_s_domain_interface_status": artifact["exported_object"]["alpha_s_domain_interface_status"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
