#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P804 = GENERATED / "p804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_admission_probe.json"
IN_F802 = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"

OUT = GENERATED / "f804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_packet.json"
OUT_SUMMARY = GENERATED / "f804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P804, IN_F802]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F804",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p804 = load_json(IN_P804)
    f802 = load_json(IN_F802)

    if (
        p804.get("status")
        == "P804_CURRENT_STRICT_ALPHA_S_COMPOSITIONAL_ALIGNMENT_BUNDLE_EXPORT_ADMITTED_PROVIDER_ACTION_RULE_STILL_BLOCKED"
        and (p804.get("clause_split") or {}).get("alignment_bundle_clause_status") == "export_admitted"
        and (p804.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_action_rule_ref"
        and f802.get("status")
        == "F802_EXECUTED_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_RULE_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F804_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_COMPOSITIONAL_ALIGNMENT_BUNDLE_PACKET_NO_FALSE_PASS"
    else:
        status = "F804_REQUIRES_REVIEW"

    candidate = p804.get("same_domain_compositional_alignment_bundle_candidate") or {}

    artifact = {
        "stage": "F804",
        "packet_name": "CurrentStrictAlphaSSameDomainCompositionalAlignmentBundle_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p804_alignment_bundle_probe": rel(IN_P804),
            "f802_provider_action_target": rel(IN_F802),
        },
        "exported_object": {
            "object_id": "alpha_s_reference_scale_same_domain_compositional_alignment_bundle_v1",
            "goal": "Export the strongest current same-domain compositional alignment bundle beneath the missing provider action rule.",
            "relation_bundle_ref": candidate.get("relation_bundle_ref"),
            "provider_action_target_ref": candidate.get("provider_action_target_ref"),
            "shared_family_id": candidate.get("shared_family_id"),
            "shared_max_m_proxy": candidate.get("shared_max_m_proxy"),
            "shared_top_boundary_point": candidate.get("shared_top_boundary_point"),
            "alignment_refs": candidate.get("alignment_refs"),
            "scope": "strict_same_domain_compositional_alignment_only",
        },
        "current_honest_reading": [
            "The repo now exports one explicit same-domain compositional alignment bundle for the alpha_s reference-scale lane.",
            "That bundle packages one shared family winner, one shared maximum, and one shared top-boundary point on the same lane.",
            "It still does not export the provider action rule itself and therefore does not act as a provider.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether the current same-domain compositional alignment bundle can be strengthened into an actual provider action rule without importing foreign-domain reference structures or new host semantics.",
        "hard_limits": [
            "Does not claim that the provider action rule already exists.",
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
        "stage": "F804",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "provider_action_target_ref": artifact["exported_object"]["provider_action_target_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
