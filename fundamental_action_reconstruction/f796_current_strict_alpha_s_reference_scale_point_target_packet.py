#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P796 = GENERATED / "p796_current_strict_alpha_s_reference_scale_point_reuse_audit_probe.json"
IN_F795 = GENERATED / "f795_current_strict_alpha_s_top_boundary_semantic_role_target_packet.json"

OUT = GENERATED / "f796_current_strict_alpha_s_reference_scale_point_target_packet.json"
OUT_SUMMARY = GENERATED / "f796_current_strict_alpha_s_reference_scale_point_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P796, IN_F795]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F796",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p796 = load_json(IN_P796)
    f795 = load_json(IN_F795)

    if p796.get("status") == "P796_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_POINT_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE":
        status = "F796_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_POINT_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F796_REQUIRES_REVIEW"

    artifact = {
        "stage": "F796",
        "packet_name": "CurrentStrictAlphaSReferenceScalePointTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p796_reference_scale_reuse_audit": rel(IN_P796),
            "f795_top_boundary_semantic_role_target": rel(IN_F795),
        },
        "why_this_packet_exists": [
            "F795 freezes the larger top-boundary semantic-role target but still leaves open whether the F704 maximum can count as a reference-scale point.",
            "P796 shows that no current strict-side export already supplies that upgrade.",
        ],
        "target_object": {
            "object_id": "alpha_s_reference_scale_point_target_v1",
            "goal": "Freeze the exact reference-scale-point object still missing before the F704 maximum can count as more than a numeric extremum on the alpha_s lane.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit finite alpha_s family domain; silent domain drift is forbidden.",
                },
                {
                    "name": "supporting_strict_source_ref",
                    "required": True,
                    "hard_limit": "Must point to the strict-side source object that supports the extremum numerically without pretending that this already fixes its reference-scale role.",
                },
                {
                    "name": "numeric_extremum_ref",
                    "required": True,
                    "hard_limit": "Must identify the concrete exported extremum candidate, including the current F704 maximum, explicitly.",
                },
                {
                    "name": "reference_scale_point_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule upgrading that extremum into a strict reference-scale point without importing GeV, host identification, or external policy semantics.",
                },
                {
                    "name": "nonstrict_calibration_exclusion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off the proxy-to-GeV calibration lane from supplying the reference-scale role.",
                },
                {
                    "name": "selected_reference_point_output_schema",
                    "required": True,
                    "hard_limit": "Must output the admitted reference-scale point and its role inside the selected normalized family schema.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer whether the point 1 exists or whether the top of the F704 spectrum is numerically defined.",
            "It is the missing rule that would upgrade that extremum into a strict reference-scale point.",
            "F796 freezes that exact missing object without promoting the current family winner.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current strict-side invariant or theorem can justify treating the exported F704 maximum as a reference-scale point rather than only as the largest current proxy value.",
        "hard_limits": [
            "Does not claim that the reference-scale point already exists.",
            "Does not claim that f704_max_mode_anchor_family is already promoted.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F796",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
