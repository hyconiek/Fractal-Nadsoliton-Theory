#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P791 = GENERATED / "p791_current_strict_alpha_s_selection_principle_reuse_audit_probe.json"
IN_F790 = GENERATED / "f790_current_strict_alpha_s_canonical_anchor_selection_target_packet.json"

OUT = GENERATED / "f791_current_strict_alpha_s_selection_principle_target_packet.json"
OUT_SUMMARY = GENERATED / "f791_current_strict_alpha_s_selection_principle_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P791, IN_F790]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F791",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p791 = load_json(IN_P791)
    f790 = load_json(IN_F790)

    if p791.get("status") == "P791_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE":
        status = "F791_EXECUTED_CURRENT_STRICT_ALPHA_S_SELECTION_PRINCIPLE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F791_REQUIRES_REVIEW"

    artifact = {
        "stage": "F791",
        "packet_name": "CurrentStrictAlphaSSelectionPrincipleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p791_selection_principle_reuse_probe": rel(IN_P791),
            "f790_canonical_anchor_selection_target": rel(IN_F790),
        },
        "why_this_packet_exists": [
            "F790 freezes the larger canonical-anchor selection target and names selection_principle_ref as a required field.",
            "P791 shows that no currently exported strict selection theorem can be lawfully reused to fill that field on the alpha_s lane.",
        ],
        "target_object": {
            "object_id": "alpha_s_selection_principle_target_v1",
            "goal": "Freeze the exact alpha_s-specific strict selection-principle object required before selection_principle_ref can be filled on the canonical anchor lane.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit admissible alpha_s normalized candidate-family domain; probe-local family lists are not enough.",
                },
                {
                    "name": "selection_objective_or_order_rule_ref",
                    "required": True,
                    "hard_limit": "Must export one strict-side objective or order rule selecting among the admissible alpha_s candidate families; verbal plausibility is forbidden.",
                },
                {
                    "name": "strict_input_chain_ref",
                    "required": True,
                    "hard_limit": "Must reference only strict-side exported inputs and must not reintroduce GeV fields or host matching.",
                },
                {
                    "name": "uniqueness_or_finite_residual_rule_ref",
                    "required": True,
                    "hard_limit": "Must prove unique family selection or else name one explicit finite residual ambiguity and its deterministic fixing rule.",
                },
                {
                    "name": "selected_family_output_schema",
                    "required": True,
                    "hard_limit": "Must output the selected family identifier and the selected normalization rule explicitly.",
                },
                {
                    "name": "nontransfer_boundary_ref",
                    "required": True,
                    "hard_limit": "Must keep family selection separate from the later anchor-to-boundary semantic upgrade and from n_f attachment.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The repo has a formal strict-selection template elsewhere, but alpha_s still lacks its own admissible domain, selection rule, and uniqueness support.",
            "F791 therefore narrows the blocker from generic selection language to one exact missing alpha_s-side object.",
            "This prevents false reuse: analogy to F447/N483 is allowed as a template, not as a semantic transfer.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether the current P789 family space admits a nonarbitrary strict objective or order rule that selects one family without host matching or GeV reintroduction.",
        "hard_limits": [
            "Does not claim that the alpha_s selection principle already exists.",
            "Does not claim that selection_principle_ref is discharged.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F791",
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
