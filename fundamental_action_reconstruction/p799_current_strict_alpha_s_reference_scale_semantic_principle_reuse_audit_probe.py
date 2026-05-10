#!/usr/bin/env python3
from __future__ import annotations

import json
import re
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P798 = GENERATED / "p798_current_strict_alpha_s_reference_scale_invariant_support_bundle_admission_probe.json"
IN_F798 = GENERATED / "f798_current_strict_alpha_s_reference_scale_invariant_support_bundle_packet.json"
IN_F797 = GENERATED / "f797_current_strict_alpha_s_reference_scale_rule_target_packet.json"
IN_P790 = GENERATED / "p790_current_strict_alpha_s_canonical_anchor_upgrade_probe.json"
IN_P795 = GENERATED / "p795_current_strict_alpha_s_top_boundary_semantic_role_reuse_audit_probe.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p799_current_strict_alpha_s_reference_scale_semantic_principle_reuse_audit_probe.json"
OUT_SUMMARY = GENERATED / "p799_current_strict_alpha_s_reference_scale_semantic_principle_reuse_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    lowered = (
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("‑", "-")
        .replace("–", "-")
        .replace("—", "-")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("-", " ")
        .replace("_", " ")
    )
    return " ".join(re.sub(r"[^a-z0-9]+", " ", lowered).split())


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def find_semantic_principle_artifacts() -> list[str]:
    patterns = [
        "*alpha_s*reference*scale*semantic*",
        "*alpha_s*reference*scale*principle*",
        "*alpha_s*reference*scale*rule*",
        "*reference*scale*semantic*principle*",
    ]
    exclude_tokens = [
        "P790_", "F790_", "p790_", "f790_",
        "P795_", "F795_", "p795_", "f795_",
        "P796_", "F796_", "p796_", "f796_",
        "P797_", "F797_", "p797_", "f797_",
        "P798_", "F798_", "p798_", "f798_",
        "P799_", "F799_", "p799_", "f799_",
    ]
    hits: set[str] = set()
    for pattern in patterns:
        for base in (ROOT, GENERATED):
            for path in base.glob(pattern):
                if path.suffix not in {".md", ".json", ".py"}:
                    continue
                rel_path = str(path.relative_to(REPO))
                if any(token in rel_path for token in exclude_tokens):
                    continue
                hits.add(rel_path)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P798, IN_F798, IN_F797, IN_P790, IN_P795, IN_N703, IN_P710, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P799",
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
    f798 = load_json(IN_F798)
    f797 = load_json(IN_F797)
    p790 = load_json(IN_P790)
    p795 = load_json(IN_P795)
    policy = load_json(IN_POLICY)

    n703_text = IN_N703.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    semantic_principle_artifacts = find_semantic_principle_artifacts()
    support_bundle = f798.get("exported_object") or {}
    p790_blockers = p790.get("narrowed_anchor_blockers") or []

    checks = [
        {
            "id": "f798_exports_support_bundle_but_marks_rule_unresolved",
            "pass": (
                f798.get("status")
                == "F798_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_INVARIANT_SUPPORT_BUNDLE_PACKET_NO_FALSE_PASS"
                and support_bundle.get("object_id") == "alpha_s_reference_scale_invariant_support_bundle_v1"
                and support_bundle.get("unresolved_reference_scale_rule_ref") == "alpha_s_reference_scale_rule_target_v1"
            ),
            "details": "F798 exports the same-domain invariant support bundle explicitly while still marking the reference-scale rule as unresolved.",
        },
        {
            "id": "earlier_alpha_s_audits_already_isolate_missing_semantic_layer",
            "pass": (
                p790.get("status")
                == "P790_CURRENT_STRICT_ALPHA_S_CANONICAL_ANCHOR_UPGRADE_BLOCKED_ON_CURRENT_REPO_STATE"
                and any("semantic-upgrade rule" in str(item) for item in p790_blockers)
                and p795.get("status")
                == "P795_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_SEMANTIC_ROLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
            ),
            "details": "Earlier alpha_s probes already isolate the semantic-upgrade / semantic-role layer as a real blocker rather than a missing numeric construction.",
        },
        {
            "id": "n703_same_domain_meaning_stack_still_stops_before_reference_scale_semantics",
            "pass": contains_all(
                n703_text,
                [
                    "dimensionless quadratic coefficients",
                    "not yet physical masses in gev",
                    "do not constitute a standard model match claim",
                ],
            ) and not contains_all(
                n703_text,
                [
                    "reference scale",
                    "reference-scale point",
                ],
            ),
            "details": "N703 supplies same-domain meaning discipline but still stops before any reference-scale semantics.",
        },
        {
            "id": "no_other_exported_same_domain_semantic_principle_artifact_detected",
            "pass": not semantic_principle_artifacts,
            "details": "No other exported same-domain artifact upgrades the alpha_s support bundle into a reference-scale semantic principle.",
        },
        {
            "id": "nonstrict_calibration_lane_explicitly_excluded",
            "pass": (
                policy.get("scope") == "nonstrict_physical_units_calibration_only"
                and contains_all(
                    p710_text,
                    [
                        "non-strict calibration map",
                        "no strict physical-unit map",
                    ],
                )
            ),
            "details": "The proxy-to-GeV calibration lane is explicitly nonstrict and cannot supply the missing semantic principle.",
        },
        {
            "id": "f797_still_records_reference_scale_rule_as_missing_target",
            "pass": (
                f797.get("status")
                == "F797_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and (f797.get("target_object") or {}).get("object_id") == "alpha_s_reference_scale_rule_target_v1"
            ),
            "details": "F797 still records the reference-scale rule itself as a missing target above the current support bundle.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    semantic_principle_reuse_found = bool(semantic_principle_artifacts)

    clause_split = {
        "support_bundle_status": "exported"
        if checks[0]["pass"]
        else "requires_review",
        "semantic_principle_clause_status": "blocked_nonexport"
        if checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"] and checks[5]["pass"]
        else "requires_review",
        "sharp_blocker_clause": "semantic_principle_ref"
        if checks[0]["pass"] and checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"] and checks[5]["pass"]
        else None,
    }

    if (
        p798.get("status")
        == "P798_CURRENT_STRICT_ALPHA_S_INVARIANT_SUPPORT_BUNDLE_EXPORT_ADMITTED_REFERENCE_SCALE_RULE_STILL_BLOCKED"
        and not failed_checks
        and not semantic_principle_reuse_found
    ):
        status = "P799_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
    else:
        status = "P799_REQUIRES_REVIEW"

    artifact = {
        "stage": "P799",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p798_support_bundle_probe": rel(IN_P798),
            "f798_support_bundle_packet": rel(IN_F798),
            "f797_reference_scale_rule_target": rel(IN_F797),
            "p790_canonical_anchor_upgrade_probe": rel(IN_P790),
            "p795_top_boundary_semantic_role_probe": rel(IN_P795),
            "n703_meaning_theorem": rel(IN_N703),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "semantic_principle_artifacts_detected": semantic_principle_artifacts,
        "semantic_principle_reuse_found": semantic_principle_reuse_found,
        "support_bundle_ref": support_bundle.get("object_id"),
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo now exports the strongest current same-domain invariant support bundle for the F704 maximum on the alpha_s lane.",
            "What it still does not export is any same-domain semantic principle upgrading that bundle into a lawful reference-scale rule.",
            "So the blocker is no longer invariant support; it is the missing semantic principle itself.",
        ],
        "recommended_next_packet": {
            "id": "F799_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_TARGET_PACKET",
            "goal": "Freeze the exact semantic-principle object still missing before the exported alpha_s support bundle can feed a lawful reference-scale rule.",
            "minimum_fields": [
                "support_bundle_ref",
                "same_domain_meaning_base_ref",
                "semantic_principle_ref",
                "reference_scale_role_statement_ref",
                "automatic_upgrade_block_ref",
                "nonstrict_calibration_exclusion_ref",
                "selected_reference_scale_rule_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P799",
        "status": status,
        "as_of": AS_OF,
        "semantic_principle_reuse_found": semantic_principle_reuse_found,
        "sharp_blocker_clause": clause_split["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
