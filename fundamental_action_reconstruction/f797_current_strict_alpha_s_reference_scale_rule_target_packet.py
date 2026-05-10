#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P797 = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe.json"
IN_F796 = GENERATED / "f796_current_strict_alpha_s_reference_scale_point_target_packet.json"

OUT = GENERATED / "f797_current_strict_alpha_s_reference_scale_rule_target_packet.json"
OUT_SUMMARY = GENERATED / "f797_current_strict_alpha_s_reference_scale_rule_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P797, IN_F796]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F797",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p797 = load_json(IN_P797)
    f796 = load_json(IN_F796)

    split = p797.get("clause_split") or {}
    if (
        p797.get("status")
        == "P797_CURRENT_STRICT_ALPHA_S_INVARIANT_EXTREMUM_CANDIDATE_SUPPORTED_REFERENCE_SCALE_RULE_BLOCKED"
        and split.get("sharp_blocker_clause") == "reference_scale_rule"
    ):
        status = "F797_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_RULE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F797_REQUIRES_REVIEW"

    artifact = {
        "stage": "F797",
        "packet_name": "CurrentStrictAlphaSReferenceScaleRuleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p797_extremum_vs_reference_scale_audit": rel(IN_P797),
            "f796_reference_scale_point_target": rel(IN_F796),
        },
        "why_this_packet_exists": [
            "F796 freezes the larger reference-scale-point target but still leaves open whether the invariant F704 maximum can be promoted beyond stable extremum status.",
            "P797 shows that extremum stability is already candidate-supported, while the actual reference-scale rule remains missing.",
        ],
        "target_object": {
            "object_id": "alpha_s_reference_scale_rule_target_v1",
            "goal": "Freeze the exact rule object still missing before the invariant F704 maximum can count as a reference-scale point.",
            "required_fields": [
                {
                    "name": "candidate_family_domain_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit finite alpha_s family domain; silent domain drift is forbidden.",
                },
                {
                    "name": "invariant_extremum_support_refs",
                    "required": True,
                    "hard_limit": "Must cite the strict-side artifacts establishing basis-invariant, gauge-stable extremum support without over-claiming semantic promotion.",
                },
                {
                    "name": "numeric_extremum_ref",
                    "required": True,
                    "hard_limit": "Must identify the concrete exported extremum candidate, including the current F704 maximum, explicitly.",
                },
                {
                    "name": "reference_scale_rule_ref",
                    "required": True,
                    "hard_limit": "Must export the rule upgrading that stable extremum into a reference-scale point without importing GeV, host identification, or external policy semantics.",
                },
                {
                    "name": "nonstrict_calibration_exclusion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off the proxy-to-GeV calibration lane from supplying the reference-scale role.",
                },
                {
                    "name": "selected_reference_scale_output_schema",
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
            "The current sharp blocker is no longer whether the F704 maximum is invariantly stable.",
            "It is the missing reference-scale rule that would promote that stable extremum into a semantic alpha_s-side point.",
            "F797 freezes that exact missing object without promoting the current family winner.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current strict-side invariant stack can be strengthened into an actual rule that promotes the invariant F704 maximum from stable extremum to reference-scale point.",
        "hard_limits": [
            "Does not claim that the reference-scale rule already exists.",
            "Does not claim that f704_max_mode_anchor_family is already promoted.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F797",
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
