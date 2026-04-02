#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P792 = GENERATED / "p792_current_strict_alpha_s_family_selection_order_rule_probe.json"
IN_F792 = GENERATED / "f792_current_strict_alpha_s_family_selection_order_rule_target_packet.json"
IN_F704 = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N706 = ROOT / "N706_CURRENT_STRICT_RELEASE_7_OS_RESIDUAL_SIGN_GAUGE_IRRELEVANCE_PACKAGE_THEOREM.md"
IN_P694 = ROOT / "P694_CURRENT_STRICT_PHYSICAL_COMPUTABILITY_MASS_SPECTRUM_PROXY_FROM_PROJECTIVE_SELECTOR_CLOSURE_PROBE.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p793_current_strict_alpha_s_order_rule_clause_split_audit_probe.json"
OUT_SUMMARY = GENERATED / "p793_current_strict_alpha_s_order_rule_clause_split_audit_probe_summary.json"


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


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P792, IN_F792, IN_F704, IN_N705, IN_N706, IN_P694, IN_P710, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P793",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p792 = load_json(IN_P792)
    f792 = load_json(IN_F792)
    policy = load_json(IN_POLICY)

    f704_text = IN_F704.read_text(encoding="utf-8")
    n705_text = IN_N705.read_text(encoding="utf-8")
    n706_text = IN_N706.read_text(encoding="utf-8")
    p694_text = IN_P694.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    required_fields = [
        item.get("name")
        for item in (((f792.get("target_object") or {}).get("required_fields")) or [])
        if isinstance(item, dict)
    ]

    source_authority_checks = [
        {
            "id": "f704_basis_invariant_full_spectrum_exported",
            "pass": contains_all(
                f704_text,
                [
                    "basis-invariant",
                    "mass observable",
                    "dimensionless",
                ],
            ),
            "details": "F704 exports a basis-invariant full-spectrum mass observable object in strict scope.",
        },
        {
            "id": "n705_tracks_f704_as_explicit_release7_os_component",
            "pass": contains_all(
                n705_text,
                [
                    "includes the strict basis-invariant mass observable object",
                    "f704",
                ],
            ),
            "details": "N705 explicitly upgrades F704 into the tracked Release-7 OS support bundle.",
        },
        {
            "id": "n706_distinguishes_f704_from_projective_pair_proxy_layer",
            "pass": contains_all(
                n706_text,
                [
                    "p694",
                    "selector rays",
                    "f704",
                    "basis-invariant spectrum",
                ],
            ) and contains_all(p694_text, ["probe-level", "projective"]),
            "details": "N706 plus P694 support a strict observable hierarchy in which F704 is basis-invariant full-spectrum, while P694 remains a projective pair-summary proxy witness.",
        },
    ]

    normalization_boundary_checks = [
        {
            "id": "p792_rule_remains_probe_local_only",
            "pass": (
                ((p792.get("probe_local_order_rule_tuple") or {}).get("export_status"))
                == "probe_local_only"
            ),
            "details": "The current normalization-boundary preference exists only inside the probe-local ranking tuple from P792.",
        },
        {
            "id": "f792_still_names_normalization_boundary_rule_as_missing",
            "pass": "normalization_boundary_rule_ref" in required_fields,
            "details": "F792 still records normalization_boundary_rule_ref as an explicit missing field.",
        },
        {
            "id": "geometric_mean_language_is_explicitly_nonstrict_only",
            "pass": (
                policy.get("scope") == "nonstrict_physical_units_calibration_only"
                and "geometric_mean" in str(policy.get("policy_id"))
                and contains_all(
                    p710_text,
                    [
                        "non-strict calibration map",
                        "single global scale parameter",
                    ],
                )
            ),
            "details": "The only explicit geometric-mean role currently exported is the nonstrict proxy-to-GeV calibration policy/interface.",
        },
    ]

    failed_source = [item["id"] for item in source_authority_checks if not item["pass"]]
    failed_boundary = [item["id"] for item in normalization_boundary_checks if not item["pass"]]

    source_authority_clause_status = (
        "candidate_supported_not_yet_exported"
        if not failed_source
        else "requires_review"
    )
    normalization_boundary_clause_status = (
        "blocked_nonexport"
        if not failed_boundary
        else "requires_review"
    )

    if not failed_source and not failed_boundary:
        status = "P793_CURRENT_STRICT_ALPHA_S_SOURCE_AUTHORITY_CANDIDATE_SUPPORTED_NORMALIZATION_BOUNDARY_BLOCKED"
    else:
        status = "P793_REQUIRES_REVIEW"

    artifact = {
        "stage": "P793",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p792_order_rule_probe": rel(IN_P792),
            "f792_order_rule_target": rel(IN_F792),
            "f704_mass_observable_packet": rel(IN_F704),
            "n705_release7_os_theorem": rel(IN_N705),
            "n706_sign_gauge_irrelevance_theorem": rel(IN_N706),
            "p694_projective_proxy_probe": rel(IN_P694),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "clause_split_audit": {
            "source_authority_clause_status": source_authority_clause_status,
            "normalization_boundary_clause_status": normalization_boundary_clause_status,
            "sharp_blocker_clause": "normalization_boundary"
            if source_authority_clause_status == "candidate_supported_not_yet_exported"
            and normalization_boundary_clause_status == "blocked_nonexport"
            else None,
        },
        "source_authority_checks": source_authority_checks,
        "normalization_boundary_checks": normalization_boundary_checks,
        "failed_source_checks": failed_source,
        "failed_normalization_boundary_checks": failed_boundary,
        "current_honest_reading": [
            "The source-authority side now has strong strict-side support: F704 is a basis-invariant full-spectrum object tracked in Release-7 OS support, while P694 remains a projective pair-summary proxy witness.",
            "The normalization-boundary side is weaker: the winning bounded-grid preference still lives only inside the probe-local P792 ranking and has no exported strict rule.",
            "Geometric-mean language is explicitly fenced onto the nonstrict proxy-to-GeV calibration lane, so it cannot serve as strict normalization authority on the current repo state.",
        ],
        "recommended_next_packet": {
            "id": "F793_CURRENT_STRICT_ALPHA_S_NORMALIZATION_BOUNDARY_RULE_TARGET_PACKET",
            "goal": "Freeze the exact normalization-boundary rule object still missing before the current probe-local winner can move closer to export-grade authority.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "bounded_normalized_grid_rule_ref",
                "top_boundary_anchor_rule_ref",
                "strict_input_chain_ref",
                "geometric_mean_nonstrict_separation_ref",
                "selected_normalization_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P793",
        "status": status,
        "as_of": AS_OF,
        "source_authority_clause_status": source_authority_clause_status,
        "normalization_boundary_clause_status": normalization_boundary_clause_status,
        "sharp_blocker_clause": artifact["clause_split_audit"]["sharp_blocker_clause"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
