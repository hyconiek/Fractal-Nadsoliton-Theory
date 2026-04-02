#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P790 = GENERATED / "p790_current_strict_alpha_s_canonical_anchor_upgrade_probe.json"
IN_F789 = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"

OUT = GENERATED / "f790_current_strict_alpha_s_canonical_anchor_selection_target_packet.json"
OUT_SUMMARY = GENERATED / "f790_current_strict_alpha_s_canonical_anchor_selection_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P790, IN_F789]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F790",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p790 = load_json(IN_P790)
    f789 = load_json(IN_F789)

    if p790.get("status") == "P790_CURRENT_STRICT_ALPHA_S_CANONICAL_ANCHOR_UPGRADE_BLOCKED_ON_CURRENT_REPO_STATE":
        status = "F790_EXECUTED_CURRENT_STRICT_ALPHA_S_CANONICAL_ANCHOR_SELECTION_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F790_REQUIRES_REVIEW"

    artifact = {
        "stage": "F790",
        "packet_name": "CurrentStrictAlphaSCanonicalAnchorSelectionTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p790_canonical_anchor_probe": rel(IN_P790),
            "f789_normalized_boundary_interface_target": rel(IN_F789),
        },
        "why_this_packet_exists": [
            "F789 freezes the missing normalized alpha_s boundary interface target.",
            "P790 narrows the strongest remaining anchor blocker to three explicit missing objects: selection principle, semantic-upgrade rule, and n_f attachment.",
        ],
        "target_object": {
            "object_id": "alpha_s_canonical_anchor_selection_target_v1",
            "goal": "Freeze the exact missing canonical-anchor selection object needed before any normalized alpha_s boundary family can be promoted into the minimal strict bridge.",
            "required_fields": [
                {
                    "name": "candidate_anchor_family_id",
                    "required": True,
                    "hard_limit": "Must identify the selected normalized family explicitly; silent choice is forbidden.",
                },
                {
                    "name": "selection_principle_ref",
                    "required": True,
                    "hard_limit": "Must point to an exported principle choosing this family over alternative normalized families.",
                },
                {
                    "name": "anchor_to_boundary_role_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule upgrading the selected anchor from internal quadratic meaning into alpha_s boundary role semantics without host matching.",
                },
                {
                    "name": "n_f_attachment_rule_ref",
                    "required": True,
                    "hard_limit": "Must attach the selected anchor to the QCD sector interface without GeV reintroduction.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny QCD closure, Standard Model identification, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The strongest dimensionless candidate family already exists, but canonical promotion remains blocked until one explicit selection object is exported.",
            "F790 does not build that object; it freezes exactly what the missing object must contain.",
            "This keeps the lane scientific and controlled: we narrow semantics before attempting any further alpha_s promotion.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether f704_max_mode_anchor_family can satisfy candidate_anchor_family_id plus a nonarbitrary selection principle without introducing host matching or GeV semantics.",
        "hard_limits": [
            "Does not claim that the canonical anchor-selection object already exists.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F790",
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
