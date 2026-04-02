#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P794 = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe.json"
IN_F793 = GENERATED / "f793_current_strict_alpha_s_normalization_boundary_rule_target_packet.json"

OUT = GENERATED / "f794_current_strict_alpha_s_top_boundary_anchor_rule_target_packet.json"
OUT_SUMMARY = GENERATED / "f794_current_strict_alpha_s_top_boundary_anchor_rule_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P794, IN_F793]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F794",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p794 = load_json(IN_P794)
    f793 = load_json(IN_F793)

    split = p794.get("subclause_split_audit") or {}
    if (
        p794.get("status")
        == "P794_CURRENT_STRICT_ALPHA_S_BOUNDED_GRID_CANDIDATE_SUPPORTED_TOP_BOUNDARY_ANCHOR_BLOCKED"
        and split.get("sharp_blocker_subclause") == "top_boundary_anchor"
    ):
        status = "F794_EXECUTED_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_ANCHOR_RULE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F794_REQUIRES_REVIEW"

    artifact = {
        "stage": "F794",
        "packet_name": "CurrentStrictAlphaSTopBoundaryAnchorRuleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p794_subclause_audit": rel(IN_P794),
            "f793_normalization_boundary_target": rel(IN_F793),
        },
        "why_this_packet_exists": [
            "F793 freezes the larger normalization-boundary target and includes top_boundary_anchor_rule_ref as a missing component.",
            "P794 shows that bounded-grid admissibility is no longer the sharp blocker; the sharper blocker is the missing semantic anchor for the boundary point 1.",
        ],
        "target_object": {
            "object_id": "alpha_s_top_boundary_anchor_rule_target_v1",
            "goal": "Freeze the exact top-boundary anchor object still missing before the current max-normalized family can move closer to export-grade authority.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit finite alpha_s family domain; silent domain drift is forbidden.",
                },
                {
                    "name": "top_boundary_point_ref",
                    "required": True,
                    "hard_limit": "Must identify the explicit boundary point 1 inside the admitted normalized family schema.",
                },
                {
                    "name": "boundary_point_semantic_role_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the strict-side rule assigning role to the boundary point 1 without introducing host units, GeV, or Standard Model identification.",
                },
                {
                    "name": "strict_input_chain_ref",
                    "required": True,
                    "hard_limit": "Must reuse only strict-side exported inputs and remain independent of external calibration policy objects.",
                },
                {
                    "name": "nonstrict_calibration_exclusion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off the proxy-to-GeV calibration lane from supplying the anchor meaning.",
                },
                {
                    "name": "selected_anchor_output_schema",
                    "required": True,
                    "hard_limit": "Must output the admitted boundary point and its role inside the selected normalized family schema.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer whether max-normalization yields a bounded grid.",
            "It is the missing strict-side rule that would give semantic role to the explicit top boundary point 1.",
            "F794 freezes that exact missing object without promoting the current family winner.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current strict-side export can justify semantic anchoring of the boundary point 1 without borrowing anything from the nonstrict calibration lane.",
        "hard_limits": [
            "Does not claim that the top-boundary anchor rule already exists.",
            "Does not claim that f704_max_mode_anchor_family is already promoted.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F794",
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
