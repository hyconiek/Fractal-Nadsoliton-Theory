#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P799 = GENERATED / "p799_current_strict_alpha_s_reference_scale_semantic_principle_reuse_audit_probe.json"
IN_F798 = GENERATED / "f798_current_strict_alpha_s_reference_scale_invariant_support_bundle_packet.json"

OUT = GENERATED / "f799_current_strict_alpha_s_reference_scale_semantic_principle_target_packet.json"
OUT_SUMMARY = GENERATED / "f799_current_strict_alpha_s_reference_scale_semantic_principle_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P799, IN_F798]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F799",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p799 = load_json(IN_P799)
    f798 = load_json(IN_F798)

    if (
        p799.get("status")
        == "P799_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
        and (p799.get("clause_split") or {}).get("sharp_blocker_clause") == "semantic_principle_ref"
        and f798.get("status")
        == "F798_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_INVARIANT_SUPPORT_BUNDLE_PACKET_NO_FALSE_PASS"
    ):
        status = "F799_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F799_REQUIRES_REVIEW"

    support_bundle_ref = p799.get("support_bundle_ref")

    artifact = {
        "stage": "F799",
        "packet_name": "CurrentStrictAlphaSReferenceScaleSemanticPrincipleTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p799_semantic_principle_reuse_probe": rel(IN_P799),
            "f798_support_bundle_packet": rel(IN_F798),
        },
        "why_this_packet_exists": [
            "F798 exports the strongest current same-domain invariant support bundle for the F704 maximum on the alpha_s lane.",
            "P799 shows that no current repo artifact upgrades that bundle into reference-scale semantics.",
        ],
        "target_object": {
            "object_id": "alpha_s_reference_scale_semantic_principle_target_v1",
            "goal": "Freeze the exact semantic-principle object still missing before the exported alpha_s support bundle can feed a lawful reference-scale rule.",
            "required_fields": [
                {
                    "name": "support_bundle_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit same-domain support bundle; silent domain drift is forbidden.",
                },
                {
                    "name": "same_domain_meaning_base_ref",
                    "required": True,
                    "hard_limit": "Must cite the strict same-domain meaning theorem(s) that define the current scope without silently adding new host semantics.",
                },
                {
                    "name": "semantic_principle_ref",
                    "required": True,
                    "hard_limit": "Must export the principle upgrading the support bundle into reference-scale semantics without importing GeV, host matching, or external policy objects.",
                },
                {
                    "name": "reference_scale_role_statement_ref",
                    "required": True,
                    "hard_limit": "Must state the resulting reference-scale role explicitly rather than leaving it as an inferred extremum label.",
                },
                {
                    "name": "automatic_upgrade_block_ref",
                    "required": True,
                    "hard_limit": "Must explicitly record why invariant support does not auto-upgrade into reference-scale semantics.",
                },
                {
                    "name": "nonstrict_calibration_exclusion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off the proxy-to-GeV calibration lane from supplying the semantic upgrade.",
                },
                {
                    "name": "selected_reference_scale_rule_output_schema",
                    "required": True,
                    "hard_limit": "Must output the admitted reference-scale rule and its scope on the alpha_s lane.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer same-domain invariant support for the F704 maximum.",
            "It is the missing semantic principle that would upgrade that support bundle into a reference-scale rule.",
            "F799 freezes that exact missing object without promoting the current bundle into alpha_s closure.",
        ],
        "support_bundle_ref": support_bundle_ref,
        "recommended_next_move": "Build one narrow probe testing whether any current strict-side theorem, objective, or reference-datum class can lawfully supply alpha_s_reference_scale_semantic_principle_target_v1 without domain drift.",
        "hard_limits": [
            "Does not claim that the semantic principle already exists.",
            "Does not claim that the F704 maximum is already a lawful reference-scale point.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F799",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": artifact["target_object"]["object_id"],
        "support_bundle_ref": support_bundle_ref,
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
