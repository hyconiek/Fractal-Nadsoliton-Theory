#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P800 = GENERATED / "p800_current_strict_alpha_s_reference_scale_provider_class_reuse_audit_probe.json"
IN_F799 = GENERATED / "f799_current_strict_alpha_s_reference_scale_semantic_principle_target_packet.json"

OUT = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"
OUT_SUMMARY = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P800, IN_F799]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F800",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p800 = load_json(IN_P800)
    f799 = load_json(IN_F799)

    if (
        p800.get("status")
        == "P800_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
        and (p800.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_class_ref"
        and f799.get("status")
        == "F799_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F800_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F800_REQUIRES_REVIEW"

    artifact = {
        "stage": "F800",
        "packet_name": "CurrentStrictAlphaSReferenceScaleProviderClassTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p800_provider_class_reuse_probe": rel(IN_P800),
            "f799_semantic_principle_target": rel(IN_F799),
        },
        "why_this_packet_exists": [
            "F799 freezes the missing semantic principle above the alpha_s support bundle.",
            "P800 shows that no currently exported theorem/objective/reference-datum class can lawfully supply that semantic principle on the same domain.",
        ],
        "target_object": {
            "object_id": "alpha_s_reference_scale_provider_class_target_v1",
            "goal": "Freeze the exact provider-class object still missing before the alpha_s semantic principle can be lawfully supplied without domain drift.",
            "required_fields": [
                {
                    "name": "support_bundle_ref",
                    "required": True,
                    "hard_limit": "Must point to one explicit same-domain alpha_s support bundle.",
                },
                {
                    "name": "semantic_principle_target_ref",
                    "required": True,
                    "hard_limit": "Must point to the exact unresolved alpha_s semantic-principle target object.",
                },
                {
                    "name": "same_domain_carrier_ref",
                    "required": True,
                    "hard_limit": "Must identify the carrier/domain on which the provider class acts; foreign-domain reuse is forbidden.",
                },
                {
                    "name": "provider_class_ref",
                    "required": True,
                    "hard_limit": "Must export the provider class itself rather than citing only foreign analogies.",
                },
                {
                    "name": "theorem_or_objective_grade_ref",
                    "required": True,
                    "hard_limit": "Must classify whether the provider is theorem-level, packet-level, or objective-level; probe-only classes cannot be silently promoted.",
                },
                {
                    "name": "foreign_domain_reuse_block_ref",
                    "required": True,
                    "hard_limit": "Must explicitly record why pair-plane, Z24, and per-site provider classes do not auto-transfer to the alpha_s support-bundle domain.",
                },
                {
                    "name": "probe_level_nonpromotion_ref",
                    "required": True,
                    "hard_limit": "Must explicitly fence off probe-level-only objective lanes from supplying the provider class.",
                },
                {
                    "name": "selected_semantic_principle_supply_schema",
                    "required": True,
                    "hard_limit": "Must output how the provider class supplies the alpha_s semantic principle on the same domain.",
                },
                {
                    "name": "hard_limits",
                    "required": True,
                    "hard_limit": "Must explicitly deny Standard Model identification, QCD closure, and ToE closure.",
                },
            ],
        },
        "current_honest_reading": [
            "The current sharp blocker is no longer whether the repo has any strict provider-class pattern at all.",
            "It is the missing same-domain provider class that could supply the alpha_s reference-scale semantic principle without domain drift.",
            "F800 freezes that exact missing object without promoting any foreign-domain theorem or probe into the alpha_s lane.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current same-domain F704 / OS / meaning-stack ingredient can be reorganized into alpha_s_reference_scale_provider_class_target_v1 without importing foreign-domain reference structures.",
        "hard_limits": [
            "Does not claim that the provider class already exists.",
            "Does not claim that the semantic principle already exists.",
            "Does not claim that the F704 maximum is already a lawful reference-scale point.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F800",
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
