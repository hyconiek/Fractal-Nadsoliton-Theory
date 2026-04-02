#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P801 = GENERATED / "p801_current_strict_alpha_s_same_domain_provider_skeleton_admission_probe.json"
IN_F800 = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"

OUT = GENERATED / "f801_current_strict_alpha_s_same_domain_provider_skeleton_packet.json"
OUT_SUMMARY = GENERATED / "f801_current_strict_alpha_s_same_domain_provider_skeleton_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P801, IN_F800]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F801",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p801 = load_json(IN_P801)
    f800 = load_json(IN_F800)

    if (
        p801.get("status")
        == "P801_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_EXPORT_ADMITTED_PROVIDER_CLASS_STILL_BLOCKED"
        and (p801.get("clause_split") or {}).get("provider_skeleton_clause_status") == "export_admitted"
        and (p801.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_class_ref"
        and f800.get("status")
        == "F800_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_TARGET_PACKET_NO_FALSE_PASS"
    ):
        status = "F801_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_PACKET_NO_FALSE_PASS"
    else:
        status = "F801_REQUIRES_REVIEW"

    candidate = p801.get("same_domain_provider_skeleton_candidate") or {}

    artifact = {
        "stage": "F801",
        "packet_name": "CurrentStrictAlphaSSameDomainProviderSkeleton_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p801_provider_skeleton_probe": rel(IN_P801),
            "f800_provider_class_target": rel(IN_F800),
        },
        "exported_object": {
            "object_id": "alpha_s_reference_scale_same_domain_provider_skeleton_v1",
            "goal": "Export the strongest current same-domain provider skeleton beneath the missing alpha_s provider class.",
            "support_bundle_ref": candidate.get("support_bundle_ref"),
            "provider_class_target_ref": candidate.get("provider_class_target_ref"),
            "carrier_refs": candidate.get("carrier_refs"),
            "numeric_snapshot": candidate.get("numeric_snapshot"),
            "scope": "strict_same_domain_provider_skeleton_only",
        },
        "current_honest_reading": [
            "The repo now exports one explicit same-domain provider skeleton for the alpha_s reference-scale lane.",
            "That skeleton packages the shared OS carrier stack already present on the current repo state.",
            "It still does not export the provider class itself and therefore does not supply the missing semantic principle.",
        ],
        "recommended_next_move": "Build one narrow probe testing whether the current same-domain provider skeleton can be reorganized into an actual provider class without importing foreign-domain reference structures or new host semantics.",
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
        "stage": "F801",
        "status": status,
        "as_of": AS_OF,
        "exported_object_id": artifact["exported_object"]["object_id"],
        "provider_class_target_ref": artifact["exported_object"]["provider_class_target_ref"],
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
