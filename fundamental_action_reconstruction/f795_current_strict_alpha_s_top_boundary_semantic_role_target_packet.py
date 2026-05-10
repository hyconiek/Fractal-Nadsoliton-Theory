#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P795 = GENERATED / "p795_current_strict_alpha_s_top_boundary_semantic_role_reuse_audit_probe.json"
IN_F794 = GENERATED / "f794_current_strict_alpha_s_top_boundary_anchor_rule_target_packet.json"

OUT = GENERATED / "f795_current_strict_alpha_s_top_boundary_semantic_role_target_packet.json"
OUT_SUMMARY = GENERATED / "f795_current_strict_alpha_s_top_boundary_semantic_role_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P795, IN_F794]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F795",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p795 = load_json(IN_P795)
    f794 = load_json(IN_F794)

    if p795.get("status") == "P795_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_SEMANTIC_ROLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE":
        status = "F795_EXECUTED_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_SEMANTIC_ROLE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F795_REQUIRES_REVIEW"

    artifact = {
        "stage": "F795",
        "packet_name": "CurrentStrictAlphaSTopBoundarySemanticRoleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p795_semantic_role_reuse_audit": rel(IN_P795),
            "f794_top_boundary_anchor_target": rel(IN_F794),
        },
        "why_this_packet_exists": [
            "F794 freezes the larger top-boundary anchor target and still leaves the semantic role of the point 1 unresolved.",
            "P795 shows that no current strict-side export already supplies that role.",
        ],
        "target_object": {
            "object_id": "alpha_s_top_boundary_semantic_role_target_v1",
            "goal": "Freeze the exact semantic-role object still missing before the normalized boundary point 1 can count as anything more than a normalization artifact.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit finite alpha_s family domain; silent domain drift is forbidden.",
                },
                {
                    "name": "boundary_point_ref",
                    "required": True,
                    "hard_limit": "Must identify the explicit boundary point 1 inside the admitted normalized family schema.",
                },
                {
                    "name": "supporting_strict_source_ref",
                    "required": True,
                    "hard_limit": "Must point to the strict-side source object supporting the point numerically without pretending that this already fixes its role.",
                },
                {
                    "name": "boundary_point_semantic_role_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the strict-side rule giving semantic role to the point 1 without importing host units, GeV, or Standard Model identification.",
                },
                {
                    "name": "nonstrict_calibration_exclusion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off the proxy-to-GeV calibration lane from supplying that role.",
                },
                {
                    "name": "selected_role_output_schema",
                    "required": True,
                    "hard_limit": "Must output the admitted role of the point 1 inside the selected normalized family schema.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer the existence of the point 1 or the boundedness of the normalized grid.",
            "It is the missing strict-side semantic role of that point.",
            "F795 freezes that exact missing object without promoting the current family winner.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether the exported F704 maximum can be upgraded from a numeric extremum to a strict reference-scale point without borrowing anything from the nonstrict calibration lane.",
        "hard_limits": [
            "Does not claim that the top-boundary semantic role already exists.",
            "Does not claim that f704_max_mode_anchor_family is already promoted.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F795",
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
