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

IN_P797 = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe.json"
IN_F797 = GENERATED / "f797_current_strict_alpha_s_reference_scale_rule_target_packet.json"
IN_F704_MD = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_F704_OBJ = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N706 = ROOT / "N706_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md"
IN_P709 = ROOT / "P709_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDIT_PROBE.md"
IN_P709_SUMMARY = GENERATED / "p709_current_strict_release_7_os_residual_sign_gauge_irrelevance_audit_probe_summary.json"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P696 = ROOT / "P696_CURRENT_STRICT_PHYSICAL_COMPUTABILITY_SELECTOR_ALIGNED_CHANNEL_SPECTRUM_PROXY_FROM_PROJECTIVE_SELECTOR_CLOSURE_PROBE.md"

OUT = GENERATED / "p798_current_strict_alpha_s_reference_scale_invariant_support_bundle_admission_probe.json"
OUT_SUMMARY = GENERATED / "p798_current_strict_alpha_s_reference_scale_invariant_support_bundle_admission_probe_summary.json"


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
        IN_P797,
        IN_F797,
        IN_F704_MD,
        IN_F704_OBJ,
        IN_N705,
        IN_N706,
        IN_P709,
        IN_P709_SUMMARY,
        IN_N703,
        IN_P696,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P798",
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
    f797 = load_json(IN_F797)
    f704 = load_json(IN_F704_OBJ)
    p709_summary = load_json(IN_P709_SUMMARY)

    f704_text = IN_F704_MD.read_text(encoding="utf-8")
    n705_text = IN_N705.read_text(encoding="utf-8")
    n706_text = IN_N706.read_text(encoding="utf-8")
    p709_text = IN_P709.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    p696_text = IN_P696.read_text(encoding="utf-8")
    same_domain_text = "\n".join([f704_text, n705_text, n706_text, p709_text, n703_text, p696_text])

    outputs = f704.get("outputs") or {}
    m_proxy = outputs.get("m_proxy_ascending") or []
    m2_proxy = outputs.get("eigenvalues_m2_proxy_ascending") or []
    max_m = max(m_proxy) if m_proxy else None
    max_m2 = max(m2_proxy) if m2_proxy else None
    unique_max = bool(m_proxy) and sum(1 for v in m_proxy if abs(float(v) - float(max_m)) <= 1e-12) == 1

    checks = [
        {
            "id": "p797_already_supports_stable_invariant_extremum_candidate",
            "pass": (
                p797.get("status")
                == "P797_CURRENT_STRICT_ALPHA_S_INVARIANT_EXTREMUM_CANDIDATE_SUPPORTED_REFERENCE_SCALE_RULE_BLOCKED"
                and (p797.get("clause_split") or {}).get("invariant_extremum_clause_status")
                == "candidate_supported_not_yet_exported"
                and (p797.get("clause_split") or {}).get("sharp_blocker_clause") == "reference_scale_rule"
            ),
            "details": "P797 already supports the F704 maximum as a stable invariant extremum candidate while keeping the reference-scale rule blocked.",
        },
        {
            "id": "f704_exports_same_domain_unique_basis_invariant_extremum",
            "pass": (
                f704.get("stage") == "F704"
                and f704.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
                and (f704.get("definition") or {}).get("basis_invariant") is True
                and unique_max
                and max_m is not None
                and max_m2 is not None
            ),
            "details": "F704 exports a same-domain basis-invariant spectrum with a unique current maximum.",
        },
        {
            "id": "n705_n706_p709_support_tracking_and_gauge_irrelevance",
            "pass": (
                contains_all(
                    n705_text,
                    [
                        "basis-invariant mass observable object",
                        "tracked component",
                    ],
                )
                and contains_all(
                    n706_text,
                    [
                        "f704 basis-invariant mass observable object",
                        "invariant under any such sign flips",
                    ],
                )
                and p709_summary.get("status") == "PASS_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_AUDITED"
                and p709_summary.get("baseline_ok") is True
                and p709_summary.get("sign_ok") is True
            ),
            "details": "N705/N706/P709 keep the F704 extremum inside the admitted tracked and gauge-stable OS stack.",
        },
        {
            "id": "n703_p696_keep_meaning_boundary_dimensionless_only",
            "pass": (
                contains_all(
                    n703_text,
                    [
                        "dimensionless quadratic coefficients",
                        "not yet physical masses in gev",
                    ],
                )
                and contains_all(
                    p696_text,
                    [
                        "basis-level operational proxy",
                        "not a claim that h psi diagonalizes",
                    ],
                )
            ),
            "details": "N703/P696 keep the lane inside dimensionless operational meaning and fence off silent physical-unit promotion.",
        },
        {
            "id": "same_domain_docs_still_stop_before_reference_scale_semantics",
            "pass": not contains_any(
                same_domain_text,
                [
                    "reference scale",
                    "reference-scale point",
                    "reference point",
                    "reference datum",
                    "alpha_s boundary anchor",
                ],
            ),
            "details": "The same-domain F704 OS lane still stops before any exported reference-scale semantics.",
        },
        {
            "id": "f797_keeps_reference_scale_rule_as_missing_object",
            "pass": (
                f797.get("status")
                == "F797_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_RULE_TARGET_PACKET_NO_FALSE_PASS"
                and (f797.get("target_object") or {}).get("object_id") == "alpha_s_reference_scale_rule_target_v1"
            ),
            "details": "F797 still records the reference-scale rule itself as the missing object above this support layer.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    support_bundle_clause_status = (
        "export_admitted"
        if checks[0]["pass"] and checks[1]["pass"] and checks[2]["pass"] and checks[3]["pass"]
        else "requires_review"
    )
    reference_scale_rule_clause_status = (
        "blocked_nonexport"
        if checks[4]["pass"] and checks[5]["pass"]
        else "requires_review"
    )

    if all(item["pass"] for item in checks):
        status = "P798_CURRENT_STRICT_ALPHA_S_INVARIANT_SUPPORT_BUNDLE_EXPORT_ADMITTED_REFERENCE_SCALE_RULE_STILL_BLOCKED"
    else:
        status = "P798_REQUIRES_REVIEW"

    artifact = {
        "stage": "P798",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p797_extremum_vs_reference_scale_audit": rel(IN_P797),
            "f797_reference_scale_rule_target": rel(IN_F797),
            "f704_packet": rel(IN_F704_MD),
            "f704_object": rel(IN_F704_OBJ),
            "n705_os_support_theorem": rel(IN_N705),
            "n706_sign_gauge_irrelevance_theorem": rel(IN_N706),
            "p709_sign_gauge_audit": rel(IN_P709),
            "p709_sign_gauge_audit_summary": rel(IN_P709_SUMMARY),
            "n703_meaning_theorem": rel(IN_N703),
            "p696_selector_aligned_proxy_probe": rel(IN_P696),
        },
        "support_bundle_candidate": {
            "family_id": "f704_max_mode_anchor_family",
            "numeric_extremum_snapshot": {
                "max_m_proxy": max_m,
                "max_m2_proxy": max_m2,
                "unique_max": unique_max,
                "mode_count": len(m_proxy),
            },
            "candidate_support_refs": [
                rel(IN_F704_MD),
                rel(IN_N705),
                rel(IN_N706),
                rel(IN_P709),
                rel(IN_N703),
                rel(IN_P696),
                rel(IN_P797),
            ],
        },
        "clause_split": {
            "support_bundle_clause_status": support_bundle_clause_status,
            "reference_scale_rule_clause_status": reference_scale_rule_clause_status,
            "sharp_blocker_clause": "reference_scale_rule"
            if support_bundle_clause_status == "export_admitted"
            and reference_scale_rule_clause_status == "blocked_nonexport"
            else None,
        },
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current same-domain F704 OS lane already supports one explicit invariant support bundle for the F704 maximum on the alpha_s route.",
            "That bundle includes basis-invariant spectrum support, tracked OS support, residual-sign gauge irrelevance, and a dimensionless meaning boundary.",
            "What remains missing is only the semantic rule that upgrades this supported extremum into a reference-scale point.",
        ],
        "recommended_next_packet": {
            "id": "F798_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_INVARIANT_SUPPORT_BUNDLE_PACKET",
            "goal": "Export the admitted invariant support bundle explicitly while leaving the reference-scale rule unresolved.",
            "export_object_id": "alpha_s_reference_scale_invariant_support_bundle_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P798",
        "status": status,
        "as_of": AS_OF,
        "support_bundle_clause_status": support_bundle_clause_status,
        "reference_scale_rule_clause_status": reference_scale_rule_clause_status,
        "sharp_blocker_clause": artifact["clause_split"]["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
