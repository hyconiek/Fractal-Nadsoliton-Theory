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

IN_P800 = GENERATED / "p800_current_strict_alpha_s_reference_scale_provider_class_reuse_audit_probe.json"
IN_F800 = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"
IN_F798 = GENERATED / "f798_current_strict_alpha_s_reference_scale_invariant_support_bundle_packet.json"
IN_F704_MD = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_F704_OBJ = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N706 = ROOT / "N706_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P696 = ROOT / "P696_CURRENT_STRICT_PHYSICAL_COMPUTABILITY_SELECTOR_ALIGNED_CHANNEL_SPECTRUM_PROXY_FROM_PROJECTIVE_SELECTOR_CLOSURE_PROBE.md"
IN_P709 = ROOT / "P709_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDIT_PROBE.md"
IN_P709_SUMMARY = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p801_current_strict_alpha_s_same_domain_provider_skeleton_admission_probe.json"
OUT_SUMMARY = GENERATED / "p801_current_strict_alpha_s_same_domain_provider_skeleton_admission_probe_summary.json"


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


def contains_any(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return any(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_P800,
        IN_F800,
        IN_F798,
        IN_F704_MD,
        IN_F704_OBJ,
        IN_N705,
        IN_N706,
        IN_N703,
        IN_P696,
        IN_P709,
        IN_P709_SUMMARY,
        IN_P710,
        IN_POLICY,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P801",
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
    f800 = load_json(IN_F800)
    f798 = load_json(IN_F798)
    f704 = load_json(IN_F704_OBJ)
    p709_summary = load_json(IN_P709_SUMMARY)
    policy = load_json(IN_POLICY)

    f704_text = IN_F704_MD.read_text(encoding="utf-8")
    n705_text = IN_N705.read_text(encoding="utf-8")
    n706_text = IN_N706.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    p696_text = IN_P696.read_text(encoding="utf-8")
    p709_text = IN_P709.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    same_domain_text = "\n".join([f704_text, n705_text, n706_text, n703_text, p696_text, p709_text])
    support_bundle_ref = (f798.get("exported_object") or {}).get("object_id")
    provider_class_target_ref = (f800.get("target_object") or {}).get("object_id")

    f704_def = f704.get("definition") or {}
    f704_outputs = f704.get("outputs") or {}
    m_proxy_ascending = f704_outputs.get("m_proxy_ascending") or []
    max_m_proxy = max(m_proxy_ascending) if m_proxy_ascending else None

    checks = [
        {
            "id": "f800_already_isolates_missing_provider_class_layer",
            "pass": (
                p800.get("status")
                == "P800_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
                and f800.get("status")
                == "F800_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_TARGET_PACKET_NO_FALSE_PASS"
                and support_bundle_ref == "alpha_s_reference_scale_invariant_support_bundle_v1"
                and provider_class_target_ref == "alpha_s_reference_scale_provider_class_target_v1"
            ),
            "details": "P800/F800 already isolate the missing provider-class layer above the same-domain alpha_s support bundle.",
        },
        {
            "id": "same_domain_carrier_stack_is_coherent",
            "pass": (
                f704.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
                and f704_def.get("basis_invariant") is True
                and contains_all(n705_text, ["basis invariant mass observable object", "tracked component"])
                and contains_all(n706_text, ["f704 basis invariant mass observable object", "gauge convention layer"])
                and contains_all(n703_text, ["psi sector hessian", "dimensionless quadratic coefficients"])
                and contains_all(p696_text, ["selector aligned", "basis level operational proxy"])
                and p709_summary.get("status") == "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED"
            ),
            "details": "F704/N705/N706/N703/P696/P709 already form one coherent same-domain carrier stack on the Release-7 OS lane.",
        },
        {
            "id": "same_domain_stack_supplies_admissible_provider_skeleton_roles",
            "pass": (
                max_m_proxy is not None
                and contains_all(f704_text, ["mass like", "basis invariant"])
                and contains_all(p696_text, ["physical computability", "strict scope"])
                and contains_all(n706_text, ["residual sign can be frozen", "without changing the exported os observable values"])
                and contains_all(n703_text, ["meaning definition", "not yet physical masses in gev"])
            ),
            "details": "The same-domain stack already supplies observable carrier, computability, gauge safety, and meaning discipline roles below any future provider class.",
        },
        {
            "id": "same_domain_stack_still_stops_before_provider_class_rule",
            "pass": not contains_any(
                same_domain_text,
                [
                    "provider class",
                    "reference scale rule",
                    "semantic principle",
                    "reference datum",
                    "selected semantic principle supply schema",
                ],
            ),
            "details": "The same-domain stack still stops before any exported provider-class rule or semantic-principle supply object.",
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
            "details": "The proxy-to-GeV calibration lane remains explicitly excluded from the same-domain provider-skeleton question.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    clause_split = {
        "provider_skeleton_clause_status": "export_admitted"
        if checks[1]["pass"] and checks[2]["pass"] and checks[4]["pass"]
        else "requires_review",
        "provider_class_clause_status": "blocked_nonexport"
        if checks[0]["pass"] and checks[3]["pass"]
        else "requires_review",
        "sharp_blocker_clause": "provider_class_ref"
        if all(item["pass"] for item in checks)
        else None,
    }

    if all(item["pass"] for item in checks):
        status = "P801_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_EXPORT_ADMITTED_PROVIDER_CLASS_STILL_BLOCKED"
    else:
        status = "P801_REQUIRES_REVIEW"

    artifact = {
        "stage": "P801",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p800_provider_class_probe": rel(IN_P800),
            "f800_provider_class_target": rel(IN_F800),
            "f798_support_bundle_packet": rel(IN_F798),
            "f704_packet": rel(IN_F704_MD),
            "f704_object": rel(IN_F704_OBJ),
            "n705_os_support_theorem": rel(IN_N705),
            "n706_sign_gauge_theorem": rel(IN_N706),
            "n703_meaning_theorem": rel(IN_N703),
            "p696_computability_probe": rel(IN_P696),
            "p709_sign_gauge_audit": rel(IN_P709),
            "p709_sign_gauge_summary": rel(IN_P709_SUMMARY),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "same_domain_provider_skeleton_candidate": {
            "support_bundle_ref": support_bundle_ref,
            "provider_class_target_ref": provider_class_target_ref,
            "carrier_refs": [
                rel(IN_F704_MD),
                rel(IN_N705),
                rel(IN_N706),
                rel(IN_N703),
                rel(IN_P696),
                rel(IN_P709),
            ],
            "numeric_snapshot": {
                "max_m_proxy": max_m_proxy,
                "max_m2_proxy": f704_outputs.get("max_m2_proxy"),
                "mode_count": f704_outputs.get("n"),
            },
        },
        "clause_split": clause_split,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current same-domain alpha_s lane already exports a coherent carrier stack below the missing provider class.",
            "That stack supplies observable carrier, computability, gauge safety, and meaning discipline on one shared Release-7 OS domain.",
            "What remains missing is the provider class itself; the current same-domain stack still does not export the rule that would supply the semantic principle.",
        ],
        "recommended_next_packet": {
            "id": "F801_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_PACKET",
            "goal": "Export the admitted same-domain provider skeleton explicitly while keeping the provider class unresolved.",
            "export_object_id": "alpha_s_reference_scale_same_domain_provider_skeleton_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P801",
        "status": status,
        "as_of": AS_OF,
        "provider_skeleton_clause_status": clause_split["provider_skeleton_clause_status"],
        "provider_class_clause_status": clause_split["provider_class_clause_status"],
        "sharp_blocker_clause": clause_split["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
