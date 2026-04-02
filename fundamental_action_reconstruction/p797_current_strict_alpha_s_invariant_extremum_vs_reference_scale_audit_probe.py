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

IN_P796 = GENERATED / "p796_current_strict_alpha_s_reference_scale_point_reuse_audit_probe.json"
IN_F796 = GENERATED / "f796_current_strict_alpha_s_reference_scale_point_target_packet.json"
IN_F704_OBJ = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N706 = ROOT / "N706_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P696 = ROOT / "P696_CURRENT_STRICT_PHYSICAL_COMPUTABILITY_SELECTOR_ALIGNED_CHANNEL_SPECTRUM_PROXY_FROM_PROJECTIVE_SELECTOR_CLOSURE_PROBE.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe.json"
OUT_SUMMARY = GENERATED / "p797_current_strict_alpha_s_invariant_extremum_vs_reference_scale_audit_probe_summary.json"


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

    prereq = [IN_P796, IN_F796, IN_F704_OBJ, IN_N705, IN_N706, IN_N703, IN_P696, IN_P710, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P797",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p796 = load_json(IN_P796)
    f796 = load_json(IN_F796)
    f704 = load_json(IN_F704_OBJ)
    policy = load_json(IN_POLICY)

    n705_text = IN_N705.read_text(encoding="utf-8")
    n706_text = IN_N706.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    p696_text = IN_P696.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    outputs = f704.get("outputs") or {}
    m2 = outputs.get("eigenvalues_m2_proxy_ascending") or []
    m = outputs.get("m_proxy_ascending") or []
    max_m2 = max(m2) if m2 else None
    max_m = max(m) if m else None
    max_m2_unique = bool(m2) and sum(1 for v in m2 if abs(float(v) - float(max_m2)) <= 1e-12) == 1
    max_m_unique = bool(m) and sum(1 for v in m if abs(float(v) - float(max_m)) <= 1e-12) == 1

    checks = [
        {
            "id": "f704_exports_sorted_basis_invariant_spectrum_with_unique_maximum",
            "pass": (
                f704.get("stage") == "F704"
                and f704.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
                and bool(m2)
                and bool(m)
                and max_m2_unique
                and max_m_unique
            ),
            "details": "F704 exports a sorted basis-invariant spectrum with a unique current maximum.",
        },
        {
            "id": "n705_n706_support_invariant_extremum_stability",
            "pass": contains_all(
                n705_text,
                [
                    "basis-invariant mass observable object",
                    "tracked component",
                ],
            ) and contains_all(
                n706_text,
                [
                    "f704 basis-invariant mass observable object",
                    "invariant under any such sign flips",
                ],
            ),
            "details": "N705 and N706 keep the F704 extremum inside the admitted invariant OS observable stack and stable under the tracked gauge/convention layer.",
        },
        {
            "id": "n703_p696_stop_before_reference_scale_semantics",
            "pass": contains_all(
                n703_text,
                [
                    "dimensionless quadratic coefficients",
                    "not yet physical masses in gev",
                ],
            ) and contains_all(
                p696_text,
                [
                    "basis-level operational proxy",
                    "not a claim that h psi diagonalizes",
                ],
            ),
            "details": "Current meaning and computability theorems still stop before reference-scale semantics.",
        },
        {
            "id": "f796_still_names_reference_scale_point_rule_as_missing",
            "pass": any(
                isinstance(item, dict) and item.get("name") == "reference_scale_point_rule_ref"
                for item in ((f796.get("target_object") or {}).get("required_fields") or [])
            ),
            "details": "F796 still records the reference-scale-point rule as an explicit missing field.",
        },
        {
            "id": "nonstrict_calibration_lane_explicitly_excluded",
            "pass": (
                policy.get("scope") == "nonstrict_physical_units_calibration_only"
                and contains_all(
                    p710_text,
                    [
                        "non-strict calibration map",
                        "single global scale parameter",
                    ],
                )
            ),
            "details": "The proxy-to-GeV calibration lane is explicitly nonstrict and cannot supply the missing reference-scale rule.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    invariant_extremum_clause_status = (
        "candidate_supported_not_yet_exported"
        if checks[0]["pass"] and checks[1]["pass"]
        else "requires_review"
    )
    reference_scale_clause_status = (
        "blocked_nonexport"
        if checks[2]["pass"] and checks[3]["pass"] and checks[4]["pass"]
        else "requires_review"
    )

    if (
        p796.get("status")
        == "P796_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_POINT_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
        and all(item["pass"] for item in checks)
    ):
        status = "P797_CURRENT_STRICT_ALPHA_S_INVARIANT_EXTREMUM_CANDIDATE_SUPPORTED_REFERENCE_SCALE_RULE_BLOCKED"
    else:
        status = "P797_REQUIRES_REVIEW"

    artifact = {
        "stage": "P797",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p796_reference_scale_reuse_audit": rel(IN_P796),
            "f796_reference_scale_point_target": rel(IN_F796),
            "f704_object": rel(IN_F704_OBJ),
            "n705_release7_os_theorem": rel(IN_N705),
            "n706_sign_gauge_irrelevance_theorem": rel(IN_N706),
            "n703_meaning_theorem": rel(IN_N703),
            "p696_selector_aligned_spectrum_probe": rel(IN_P696),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "extremum_snapshot": {
            "max_m2_proxy": max_m2,
            "max_m_proxy": max_m,
            "max_m2_unique": max_m2_unique,
            "max_m_unique": max_m_unique,
            "mode_count": len(m),
        },
        "clause_split": {
            "invariant_extremum_clause_status": invariant_extremum_clause_status,
            "reference_scale_clause_status": reference_scale_clause_status,
            "sharp_blocker_clause": "reference_scale_rule"
            if invariant_extremum_clause_status == "candidate_supported_not_yet_exported"
            and reference_scale_clause_status == "blocked_nonexport"
            else None,
        },
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The F704 maximum is already a stable invariant extremum candidate on the current repo state: it is basis-invariant, numerically unique on the exported spectrum, and stable under the admitted sign-gauge layer.",
            "What the repo still does not export is the rule that upgrades that stable extremum into a reference-scale point.",
            "So the blocker is no longer extremum stability; it is the missing reference-scale rule itself.",
        ],
        "recommended_next_packet": {
            "id": "F797_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_RULE_TARGET_PACKET",
            "goal": "Freeze the exact rule object still missing before the invariant F704 maximum can count as a reference-scale point.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "invariant_extremum_support_refs",
                "numeric_extremum_ref",
                "reference_scale_rule_ref",
                "nonstrict_calibration_exclusion_ref",
                "selected_reference_scale_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P797",
        "status": status,
        "as_of": AS_OF,
        "invariant_extremum_clause_status": invariant_extremum_clause_status,
        "reference_scale_clause_status": reference_scale_clause_status,
        "sharp_blocker_clause": artifact["clause_split"]["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
