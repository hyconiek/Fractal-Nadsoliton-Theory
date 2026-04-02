#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P795 = GENERATED / "p795_current_strict_alpha_s_top_boundary_semantic_role_reuse_audit_probe.json"
IN_F795 = GENERATED / "f795_current_strict_alpha_s_top_boundary_semantic_role_target_packet.json"
IN_F704 = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_F704_OBJ = GENERATED / "mass_observable_diagonal_local_strict_derived_v1.json"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_F447 = ROOT / "F447_CURRENT_STRICT_T169_QW2122_SCALAR_TO_T168_PER_SITE_VALUE_PROVIDER_POWERLAW_ELEMENT_ORDER_REFERENCE_PACKET.md"
IN_N483 = ROOT / "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p796_current_strict_alpha_s_reference_scale_point_reuse_audit_probe.json"
OUT_SUMMARY = GENERATED / "p796_current_strict_alpha_s_reference_scale_point_reuse_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return " ".join(
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
        .split()
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def find_reference_scale_artifacts() -> list[str]:
    patterns = [
        "*reference*scale*point*",
        "*alpha_s*reference*point*",
        "*f704*reference*",
    ]
    hits: set[str] = set()
    for pattern in patterns:
        for base in (ROOT, GENERATED):
            for path in base.glob(pattern):
                if path.suffix not in {".md", ".json", ".py"}:
                    continue
                rel_path = str(path.relative_to(REPO))
                if any(token in rel_path for token in ["P796_", "F796_", "p796_", "f796_"]):
                    continue
                hits.add(rel_path)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P795, IN_F795, IN_F704, IN_F704_OBJ, IN_N705, IN_N703, IN_F447, IN_N483, IN_P710, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P796",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p795 = load_json(IN_P795)
    f795 = load_json(IN_F795)
    f704_obj = load_json(IN_F704_OBJ)
    policy = load_json(IN_POLICY)

    f704_text = IN_F704.read_text(encoding="utf-8")
    n705_text = IN_N705.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    f447_text = IN_F447.read_text(encoding="utf-8")
    n483_text = IN_N483.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    reference_scale_artifacts = find_reference_scale_artifacts()
    m_proxy = ((f704_obj.get("outputs") or {}).get("m_proxy_ascending")) or []
    max_value = max(m_proxy) if m_proxy else None

    checks = [
        {
            "id": "repo_exports_real_reference_datum_pattern_elsewhere",
            "pass": contains_all(
                f447_text,
                [
                    "exported strict reference datum",
                    "reference datum",
                ],
            ) and contains_all(
                n483_text,
                [
                    "reference datum",
                    "existence + uniqueness",
                ],
            ),
            "details": "The repo does export real strict reference-datum objects in other lanes, so the question is lawful and not empty.",
        },
        {
            "id": "existing_reference_datum_pattern_is_domain_specific_not_f704_scale_point",
            "pass": contains_all(
                f447_text + "\n" + n483_text,
                [
                    "z 12",
                    "ord",
                    "per site",
                ],
            ) and contains_all(
                f704_text,
                [
                    "basis-invariant",
                    "mass-like",
                    "eigenvalues",
                ],
            ),
            "details": "The existing strict reference-datum pattern lives on the Z12 / per-site constrained-lift lane, not on the F704 spectral-maximum lane.",
        },
        {
            "id": "f704_n705_n703_stop_at_whole_spectrum_operational_meaning",
            "pass": max_value is not None
            and contains_all(
                f704_text,
                [
                    "basis-invariant",
                    "dimensionless",
                ],
            )
            and contains_all(
                n705_text,
                [
                    "tracked component",
                    "basis-invariant mass observable object",
                ],
            )
            and contains_all(
                n703_text,
                [
                    "dimensionless quadratic coefficients",
                    "not yet physical masses in gev",
                ],
            ),
            "details": "Current F704/N705/N703 semantics stop at whole-spectrum operational meaning and do not upgrade the maximum into a reference-scale point.",
        },
        {
            "id": "no_current_reference_scale_point_artifact_detected",
            "pass": not reference_scale_artifacts,
            "details": "No current exported repo artifact identifies the F704 maximum as a strict alpha_s reference-scale point.",
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
            "details": "The proxy-to-GeV calibration lane is explicitly nonstrict and cannot supply the missing reference-scale semantics.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    reusable_reference_scale_point_found = bool(reference_scale_artifacts)

    if (
        p795.get("status")
        == "P795_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_SEMANTIC_ROLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
        and not failed_checks
        and not reusable_reference_scale_point_found
    ):
        status = "P796_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_POINT_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
    else:
        status = "P796_REQUIRES_REVIEW"

    artifact = {
        "stage": "P796",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p795_semantic_role_reuse_audit": rel(IN_P795),
            "f795_semantic_role_target": rel(IN_F795),
            "f704_packet": rel(IN_F704),
            "f704_object": rel(IN_F704_OBJ),
            "n705_release7_os_theorem": rel(IN_N705),
            "n703_meaning_theorem": rel(IN_N703),
            "f447_reference_datum_packet": rel(IN_F447),
            "n483_reference_datum_theorem": rel(IN_N483),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "reference_scale_artifacts_detected": reference_scale_artifacts,
        "reusable_reference_scale_point_found": reusable_reference_scale_point_found,
        "numeric_extremum_snapshot": {
            "f704_max_m_proxy": max_value,
            "f704_mode_count": len(m_proxy),
        },
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo has real strict reference-datum machinery elsewhere, but none of it acts on the F704 maximum or promotes it to a reference-scale point.",
            "Current F704/N705/N703 semantics stop at basis-invariant whole-spectrum operational meaning plus internal dimensionless scope.",
            "So the sharp blocker is now the missing rule that would upgrade the F704 maximum from numeric extremum to strict reference-scale point.",
        ],
        "recommended_next_packet": {
            "id": "F796_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_POINT_TARGET_PACKET",
            "goal": "Freeze the exact reference-scale-point object still missing before the F704 maximum can count as more than a numeric extremum on the alpha_s lane.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "supporting_strict_source_ref",
                "numeric_extremum_ref",
                "reference_scale_point_rule_ref",
                "nonstrict_calibration_exclusion_ref",
                "selected_reference_point_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P796",
        "status": status,
        "as_of": AS_OF,
        "reusable_reference_scale_point_found": reusable_reference_scale_point_found,
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
