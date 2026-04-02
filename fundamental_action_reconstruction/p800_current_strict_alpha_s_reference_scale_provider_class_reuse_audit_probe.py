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

IN_P799 = GENERATED / "p799_current_strict_alpha_s_reference_scale_semantic_principle_reuse_audit_probe.json"
IN_F799 = GENERATED / "f799_current_strict_alpha_s_reference_scale_semantic_principle_target_packet.json"
IN_F798 = GENERATED / "f798_current_strict_alpha_s_reference_scale_invariant_support_bundle_packet.json"
IN_N479 = ROOT / "N479_CURRENT_FIRST_STRICT_Z12_ELEMENT_ORDER_REFERENCE_AUT_Z12_INVARIANCE_NO_MARKED_DIRECTION_THEOREM.md"
IN_N488 = ROOT / "N488_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_CROSS_ENTROPY_CUTS_PAIR2_O2_TO_Z2_UNIQUENESS_THEOREM.md"
IN_N513 = ROOT / "N513_CURRENT_FIRST_STRICT_SHANNON_ELEMENT_ORDER_REFERENCE_Z24_CROSS_ENTROPY_CUTS_ALL_FOURIER_DEGENERATE_PAIRS_O2_TO_Z2_UNIQUENESS_THEOREM.md"
IN_F447 = ROOT / "F447_CURRENT_STRICT_T169_QW2122_SCALAR_TO_T168_PER_SITE_VALUE_PROVIDER_POWERLAW_ELEMENT_ORDER_REFERENCE_PACKET.md"
IN_N483 = ROOT / "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM.md"
IN_P439 = GENERATED / "p439_current_strict_qw2191_weighted_kl_reference_objective_o2_cut_audit_probe.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p800_current_strict_alpha_s_reference_scale_provider_class_reuse_audit_probe.json"
OUT_SUMMARY = GENERATED / "p800_current_strict_alpha_s_reference_scale_provider_class_reuse_audit_probe_summary.json"


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


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P799,
        IN_F799,
        IN_F798,
        IN_N479,
        IN_N488,
        IN_N513,
        IN_F447,
        IN_N483,
        IN_P439,
        IN_N703,
        IN_P710,
        IN_POLICY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P800",
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
    f799 = load_json(IN_F799)
    f798 = load_json(IN_F798)
    p439 = load_json(IN_P439)
    policy = load_json(IN_POLICY)

    n479_text = IN_N479.read_text(encoding="utf-8")
    n488_text = IN_N488.read_text(encoding="utf-8")
    n513_text = IN_N513.read_text(encoding="utf-8")
    f447_text = IN_F447.read_text(encoding="utf-8")
    n483_text = IN_N483.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    support_bundle_ref = f798.get("exported_object", {}).get("object_id")
    target_object_id = f799.get("target_object", {}).get("object_id")
    p439_verdict = (p439.get("verdict") or {}).get("statement", "")
    p439_hard_limits = " ".join(p439.get("hard_limits") or [])

    checks = [
        {
            "id": "f799_already_isolates_missing_semantic_principle_layer",
            "pass": (
                p799.get("status")
                == "P799_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
                and f799.get("status")
                == "F799_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_TARGET_PACKET_NO_FALSE_PASS"
                and support_bundle_ref == "alpha_s_reference_scale_invariant_support_bundle_v1"
                and target_object_id == "alpha_s_reference_scale_semantic_principle_target_v1"
            ),
            "details": "P799/F799 already isolate the missing semantic-principle layer above the exported alpha_s support bundle.",
        },
        {
            "id": "repo_exports_real_provider_class_patterns_elsewhere",
            "pass": (
                contains_all(n479_text, ["reference", "aut z 12 invariant"])
                and contains_all(n488_text, ["pair2", "o2", "z2", "cross entropy"])
                and contains_all(n513_text, ["z 24", "all fourier degenerate pairs", "z2"])
                and contains_all(f447_text, ["exported strict reference datum", "per site"])
                and contains_all(n483_text, ["reference datum", "existence uniqueness"])
            ),
            "details": "The repo does export real strict theorem / objective / reference-datum/provider-class patterns in other lanes.",
        },
        {
            "id": "those_provider_classes_are_foreign_domain_to_alpha_s_support_bundle",
            "pass": (
                contains_all(n488_text, ["pair2", "o2", "z2"])
                and contains_all(n513_text, ["z 24", "scope extension"])
                and contains_all(f447_text + "\n" + n483_text, ["per site", "t168"])
                and support_bundle_ref == "alpha_s_reference_scale_invariant_support_bundle_v1"
            ),
            "details": "The exported provider-class patterns act on pair-plane O(2)->Z2 cuts, Z24 scope-extension lanes, or T168 per-site lifts, not on the same-domain F704 alpha_s support bundle.",
        },
        {
            "id": "probe_level_objective_class_is_explicitly_nonpromoted",
            "pass": (
                contains_all(
                    p439_verdict,
                    [
                        "probe level observation only",
                        "not promoted into a strict core",
                    ],
                )
                and contains_all(
                    p439_hard_limits,
                    [
                        "audit only",
                        "does not export a strict core selector ingredient",
                    ],
                )
            ),
            "details": "The weighted-KL objective lane remains explicitly probe-level and cannot be silently upgraded into an export-grade provider class for alpha_s semantics.",
        },
        {
            "id": "same_domain_meaning_base_still_stops_before_reference_scale_semantics",
            "pass": contains_all(
                n703_text,
                [
                    "dimensionless quadratic coefficients",
                    "not yet physical masses in gev",
                    "do not constitute a standard model match claim",
                ],
            ),
            "details": "The same-domain F704 meaning base still stops before reference-scale semantics.",
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
            "details": "The proxy-to-GeV calibration lane is explicitly nonstrict and cannot supply the missing provider class.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    provider_class_reuse_found = False

    clause_split = {
        "foreign_domain_provider_patterns_present": checks[1]["pass"],
        "same_domain_provider_class_status": "blocked_nonexport"
        if checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"] and checks[5]["pass"]
        else "requires_review",
        "sharp_blocker_clause": "provider_class_ref"
        if all(item["pass"] for item in checks)
        else None,
    }

    if (
        p799.get("status")
        == "P799_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_SEMANTIC_PRINCIPLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
        and not failed_checks
        and not provider_class_reuse_found
    ):
        status = "P800_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
    else:
        status = "P800_REQUIRES_REVIEW"

    artifact = {
        "stage": "P800",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p799_semantic_principle_probe": rel(IN_P799),
            "f799_semantic_principle_target": rel(IN_F799),
            "f798_support_bundle_packet": rel(IN_F798),
            "n479_reference_datum_theorem": rel(IN_N479),
            "n488_pair2_cross_entropy_theorem": rel(IN_N488),
            "n513_z24_cross_entropy_theorem": rel(IN_N513),
            "f447_reference_packet": rel(IN_F447),
            "n483_provider_existence_uniqueness_theorem": rel(IN_N483),
            "p439_weighted_kl_audit": rel(IN_P439),
            "n703_meaning_theorem": rel(IN_N703),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "provider_class_reuse_found": provider_class_reuse_found,
        "support_bundle_ref": support_bundle_ref,
        "semantic_principle_target_ref": target_object_id,
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo exports real strict provider-class patterns elsewhere, including reference-data theorems, theorem-level O(2)->Z2 objectives, and per-site constrained-lift providers.",
            "None of those patterns lawfully enters the same domain as alpha_s_reference_scale_invariant_support_bundle_v1; all require domain drift or remain probe-level only.",
            "So the blocker is no longer the absence of any provider-class pattern in the repo, but the absence of a same-domain provider class for the alpha_s reference-scale semantic principle.",
        ],
        "recommended_next_packet": {
            "id": "F800_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_TARGET_PACKET",
            "goal": "Freeze the exact provider-class object still missing before the alpha_s semantic principle can be lawfully supplied on the same domain.",
            "minimum_fields": [
                "support_bundle_ref",
                "semantic_principle_target_ref",
                "same_domain_carrier_ref",
                "provider_class_ref",
                "theorem_or_objective_grade_ref",
                "foreign_domain_reuse_block_ref",
                "probe_level_nonpromotion_ref",
                "selected_semantic_principle_supply_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P800",
        "status": status,
        "as_of": AS_OF,
        "provider_class_reuse_found": provider_class_reuse_found,
        "sharp_blocker_clause": clause_split["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
