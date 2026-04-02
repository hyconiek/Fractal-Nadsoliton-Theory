#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P798 = GENERATED / "p798_current_strict_alpha_s_reference_scale_invariant_support_bundle_admission_probe.json"
IN_F797 = GENERATED / "f797_current_strict_alpha_s_reference_scale_rule_target_packet.json"

OUT = GENERATED / "f798_current_strict_alpha_s_reference_scale_invariant_support_bundle_packet.json"
OUT_SUMMARY = GENERATED / "f798_current_strict_alpha_s_reference_scale_invariant_support_bundle_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P798, IN_F797]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F798",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p798 = load_json(IN_P798)
    f797 = load_json(IN_F797)

    if (
        p798.get("status")
        == "P798_CURRENT_STRICT_ALPHA_S_INVARIANT_SUPPORT_BUNDLE_EXPORT_ADMITTED_REFERENCE_SCALE_RULE_STILL_BLOCKED"
        and (p798.get("clause_split") or {}).get("support_bundle_clause_status") == "export_admitted"
        and (p798.get("clause_split") or {}).get("sharp_blocker_clause") == "reference_scale_rule"
        and f797.get("status")
        == "F797_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_RULE_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F798_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_INVARIANT_SUPPORT_BUNDLE_PACKET_NO_FALSE_PASS"
    else:
        status = "F798_REQUIRES_REVIEW"

    support_bundle = p798.get("support_bundle_candidate") or {}

    artifact = {
        "stage": "F798",
        "packet_name": "CurrentStrictAlphaSReferenceScaleInvariantSupportBundle_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p798_support_bundle_admission_probe": rel(IN_P798),
            "f797_reference_scale_rule_target": rel(IN_F797),
        },
        "exported_object": {
            "object_id": "alpha_s_reference_scale_invariant_support_bundle_v1",
            "goal": "Export the strongest current same-domain support bundle behind the F704 maximum on the alpha_s lane without promoting it to a reference-scale point.",
            "family_id": support_bundle.get("family_id"),
            "numeric_extremum_snapshot": support_bundle.get("numeric_extremum_snapshot"),
            "support_refs": support_bundle.get("candidate_support_refs"),
            "unresolved_reference_scale_rule_ref": "alpha_s_reference_scale_rule_target_v1",
            "scope": "strict_same_domain_invariant_support_only",
        },
        "current_honest_reading": [
            "The repo now exports one explicit same-domain support bundle for the F704 maximum on the alpha_s lane.",
            "That bundle is strong enough to support the extremum as stable, invariant, and gauge-safe inside the dimensionless OS scope.",
            "It is still not a reference-scale rule, and it does not promote the F704 maximum into a lawful alpha_s-side reference point.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether any current strict-side semantic principle can lawfully upgrade alpha_s_reference_scale_invariant_support_bundle_v1 into an actual reference-scale rule.",
        "hard_limits": [
            "Does not claim that the reference-scale rule already exists.",
            "Does not claim that the F704 maximum is already a lawful reference-scale point.",
            "Does not claim alpha_s boundary export readiness.",
            "Does not claim QCD closure.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F798",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "unresolved_reference_scale_rule_ref": artifact["exported_object"]["unresolved_reference_scale_rule_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
