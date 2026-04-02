#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P793 = GENERATED / "p793_current_strict_alpha_s_order_rule_clause_split_audit_probe.json"
IN_F792 = GENERATED / "f792_current_strict_alpha_s_family_selection_order_rule_target_packet.json"

OUT = GENERATED / "f793_current_strict_alpha_s_normalization_boundary_rule_target_packet.json"
OUT_SUMMARY = GENERATED / "f793_current_strict_alpha_s_normalization_boundary_rule_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P793, IN_F792]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F793",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p793 = load_json(IN_P793)
    f792 = load_json(IN_F792)

    split = p793.get("clause_split_audit") or {}
    if (
        p793.get("status")
        == "P793_CURRENT_STRICT_ALPHA_S_SOURCE_AUTHORITY_CANDIDATE_SUPPORTED_NORMALIZATION_BOUNDARY_BLOCKED"
        and split.get("sharp_blocker_clause") == "normalization_boundary"
    ):
        status = "F793_EXECUTED_CURRENT_STRICT_ALPHA_S_NORMALIZATION_BOUNDARY_RULE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F793_REQUIRES_REVIEW"

    artifact = {
        "stage": "F793",
        "packet_name": "CurrentStrictAlphaSNormalizationBoundaryRuleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p793_clause_split_audit": rel(IN_P793),
            "f792_order_rule_target": rel(IN_F792),
        },
        "why_this_packet_exists": [
            "F792 freezes the larger alpha_s family-selection order-rule target and includes normalization_boundary_rule_ref as a missing component.",
            "P793 shows that normalization-boundary, not generic source-authority ranking, is now the sharp blocker on the current repo state.",
        ],
        "target_object": {
            "object_id": "alpha_s_normalization_boundary_rule_target_v1",
            "goal": "Freeze the exact normalization-boundary rule object still missing before the current probe-local family winner can move closer to export-grade authority.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit finite alpha_s family domain; silent domain drift is forbidden.",
                },
                {
                    "name": "bounded_normalized_grid_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule that admits or rejects candidate normalized grids on strict-side grounds; probe-local boundedness checks are not enough.",
                },
                {
                    "name": "top_boundary_anchor_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule assigning meaning to the explicit top boundary point 1 without reintroducing host units or GeV semantics.",
                },
                {
                    "name": "strict_input_chain_ref",
                    "required": True,
                    "hard_limit": "Must reuse only strict-side exported inputs and must not depend on external calibration policy objects.",
                },
                {
                    "name": "geometric_mean_nonstrict_separation_ref",
                    "required": True,
                    "hard_limit": "Must keep geometric-mean calibration language on the nonstrict proxy-to-GeV side unless a new strict theorem overrides that separation.",
                },
                {
                    "name": "selected_normalization_output_schema",
                    "required": True,
                    "hard_limit": "Must output the admitted normalization rule and validation-grid schema explicitly.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer generic ranking between F704 and P694.",
            "It is the missing normalization-boundary rule that would justify bounded-grid preference and top-boundary anchoring on strict-side grounds.",
            "F793 freezes that exact missing object without promoting the current family winner.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether bounded normalized-grid admissibility and top-boundary anchoring can be justified from current strict-side exports without borrowing any nonstrict calibration semantics.",
        "hard_limits": [
            "Does not claim that the normalization-boundary rule already exists.",
            "Does not claim that f704_max_mode_anchor_family is already promoted.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F793",
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
