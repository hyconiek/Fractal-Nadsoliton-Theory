#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P807 = GENERATED / "p807_current_strict_alpha_s_strict_source_shannon_reuse_as_same_domain_provider_action_source_audit_probe.json"
IN_F806 = GENERATED / "f806_current_strict_alpha_s_provider_action_continuation_boundary_packet.json"

OUT = GENERATED / "f807_current_strict_alpha_s_provider_class_shift_requirement_packet.json"
OUT_SUMMARY = GENERATED / "f807_current_strict_alpha_s_provider_class_shift_requirement_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P807, IN_F806]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F807",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p807 = load_json(IN_P807)
    f806 = load_json(IN_F806)

    if (
        p807.get("status")
        == "P807_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_REUSE_AS_SAME_DOMAIN_PROVIDER_ACTION_SOURCE_NEGATIVE_PROVIDER_SHIFT_REQUIRED"
        and (p807.get("theorem_result") or {}).get(
            "current_repo_exports_no_genuinely_new_same_domain_provider_action_source_for_alpha_s"
        )
        is True
        and (p807.get("theorem_result") or {}).get("next_honest_move_requires_provider_class_shift")
        is True
        and f806.get("status")
        == "F806_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_CONTINUATION_BOUNDARY_PACKET_NO_FALSE_PASS"
    ):
        status = "F807_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_CLASS_SHIFT_REQUIREMENT_PACKET_NO_FALSE_PASS"
    else:
        status = "F807_REQUIRES_REVIEW"

    rejected_candidate = p807.get("rejected_same_domain_source_candidate_class") or {}
    continuation_boundary = f806.get("exported_object") or {}

    artifact = {
        "stage": "F807",
        "packet_name": "CurrentStrictAlphaSProviderClassShiftRequirement_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p807_same_domain_source_audit_probe": rel(IN_P807),
            "f806_provider_action_continuation_boundary_packet": rel(IN_F806),
        },
        "exported_object": {
            "object_id": "alpha_s_provider_class_shift_requirement_v1",
            "goal": "Export the current-repo-state requirement that alpha_s continuation now proceed by provider-class shift after same-domain strict-source alpha_geo/Shannon reuse fails.",
            "continuation_boundary_ref": continuation_boundary.get("object_id"),
            "rejected_same_domain_source_candidate_id": rejected_candidate.get("candidate_id"),
            "remaining_admitted_move_class": "shift_to_a_different_provider_class_lane",
            "candidate_shift_lane_hint": "strict_source_shannon_provider_shift_lane_candidate_not_yet_alpha_s_domain_audited",
            "forbidden_move_clause": "no_same_domain_verbal_promotion_of_alpha_geo_or_shannon_candidate_language_into_alpha_s_provider_action_source",
            "scope": "strict_alpha_s_provider_class_shift_requirement_only",
        },
        "current_honest_reading": [
            "The repo now exports an explicit provider-class-shift requirement for the alpha_s lane.",
            "Same-domain reuse of strict-source alpha_geo/Shannon language has been audited and rejected as a current provider-action source.",
            "The remaining honest continuation class is a provider shift, not local rhetorical promotion.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether a strict-source Shannon lane can serve as an admissible provider-class shift candidate for alpha_s without silent domain identification.",
        "hard_limits": [
            "Does not claim that provider-class shift has already succeeded.",
            "Does not claim that strict-source Shannon already supplies alpha_s semantics.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F807",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "remaining_admitted_move_class": artifact["exported_object"]["remaining_admitted_move_class"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
