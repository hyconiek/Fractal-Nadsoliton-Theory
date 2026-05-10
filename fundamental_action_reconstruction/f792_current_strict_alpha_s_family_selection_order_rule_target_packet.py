#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
IN_F791 = GENERATED / "f791_current_strict_alpha_s_selection_principle_target_packet.json"

OUT = GENERATED / "f792_current_strict_alpha_s_family_selection_order_rule_target_packet.json"
OUT_SUMMARY = GENERATED / "f792_current_strict_alpha_s_family_selection_order_rule_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P792, IN_F791]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F792",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p792 = load_json(IN_P792)
    f791 = load_json(IN_F791)

    if p792.get("status") == "P792_CURRENT_STRICT_ALPHA_S_PROBE_LOCAL_ORDER_RULE_UNIQUE_WINNER_NONEXPORT":
        status = "F792_EXECUTED_CURRENT_STRICT_ALPHA_S_FAMILY_SELECTION_ORDER_RULE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F792_REQUIRES_REVIEW"

    artifact = {
        "stage": "F792",
        "packet_name": "CurrentStrictAlphaSFamilySelectionOrderRuleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p792_order_rule_probe": rel(IN_P792),
            "f791_selection_principle_target": rel(IN_F791),
        },
        "why_this_packet_exists": [
            "F791 freezes the larger alpha_s selection-principle target and requires selection_objective_or_order_rule_ref.",
            "P792 shows that one explicit order-rule tuple already isolates a unique current winner, but only at probe-local and nonexport level.",
        ],
        "target_object": {
            "object_id": "alpha_s_family_selection_order_rule_target_v1",
            "goal": "Freeze the exact export-grade order-rule object needed before the current probe-local winner can fill selection_objective_or_order_rule_ref on the alpha_s lane.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit finite alpha_s family domain; silent domain changes are forbidden.",
                },
                {
                    "name": "source_authority_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule that orders basis-invariant full-spectrum sources against weaker projective pair-summary sources.",
                },
                {
                    "name": "normalization_boundary_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule governing bounded normalized grids and explicit top-boundary anchoring.",
                },
                {
                    "name": "residual_tie_break_rule_ref",
                    "required": True,
                    "hard_limit": "Must either resolve residual ties deterministically or explicitly prove that no tie remains on the admitted domain.",
                },
                {
                    "name": "nonstrict_calibration_separation_ref",
                    "required": True,
                    "hard_limit": "Must keep geometric-mean proxy-to-GeV policy language on the nonstrict side unless a new strict theorem is exported.",
                },
                {
                    "name": "selected_family_output_schema",
                    "required": True,
                    "hard_limit": "Must output the selected family identifier and the exact admitted normalization rule.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "P792 narrows the selection blocker from no-winner ambiguity to no-authority ambiguity.",
            "The current winner exists only relative to a probe-local order tuple; F792 freezes what must be exported before that can count as a real alpha_s family-selection result.",
            "This keeps the lane strict: we separate domain ranking from later semantic upgrade into boundary role or n_f attachment.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether the source-authority clause and normalization-boundary clause can be upgraded from probe-local ranking rules into exportable strict-side rules on the current P789 domain.",
        "hard_limits": [
            "Does not claim that the order rule already exists as a strict export.",
            "Does not claim that f704_max_mode_anchor_family is already promoted.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F792",
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
