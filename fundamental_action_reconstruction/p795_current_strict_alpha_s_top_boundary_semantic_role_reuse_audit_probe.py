#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P794 = GENERATED / "p794_current_strict_alpha_s_normalization_boundary_subclause_audit_probe.json"
IN_F794 = GENERATED / "f794_current_strict_alpha_s_top_boundary_anchor_rule_target_packet.json"
IN_F704_PACKET = ROOT / "F704_CURRENT_STRICT_INVARIANT_MASS_OBSERVABLE_FROM_DIAGONAL_LOCAL_PSI_HESSIAN_EIGENSYSTEM_EXPORT_PACKET.md"
IN_F704_SUMMARY = GENERATED / "f704_current_strict_invariant_mass_observable_from_diagonal_local_psi_hessian_eigensystem_export_packet_summary.json"
IN_N705 = ROOT / "N705_CURRENT_FIRST_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_WITH_INVARIANT_MASS_OBSERVABLE_THEOREM.md"
IN_N703 = ROOT / "N703_CURRENT_STRICT_QUADRATIC_MASS_PROXY_MEANING_DEFINITION_THEOREM.md"
IN_P710 = ROOT / "P710_CURRENT_NONSTRICT_PROXY_TO_GEV_CALIBRATION_MAP_FROM_F704_EIGENSPECTRUM_PROBE.md"
IN_POLICY = ROOT / "external_data" / "proxy_to_gev_calibration_policy_v1.json"

OUT = GENERATED / "p795_current_strict_alpha_s_top_boundary_semantic_role_reuse_audit_probe.json"
OUT_SUMMARY = GENERATED / "p795_current_strict_alpha_s_top_boundary_semantic_role_reuse_audit_probe_summary.json"


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


def find_semantic_role_artifacts() -> list[str]:
    patterns = [
        "*top*boundary*semantic*role*",
        "*boundary*point*semantic*role*",
        "*alpha_s*top*boundary*anchor*",
    ]
    hits: set[str] = set()
    for pattern in patterns:
        for base in (ROOT, GENERATED):
            for path in base.glob(pattern):
                if path.suffix not in {".md", ".json", ".py"}:
                    continue
                rel_path = str(path.relative_to(REPO))
                if any(token in rel_path for token in ["P794_", "F794_", "p794_", "f794_", "P795_", "F795_", "p795_", "f795_"]):
                    continue
                hits.add(rel_path)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P794, IN_F794, IN_F704_PACKET, IN_F704_SUMMARY, IN_N705, IN_N703, IN_P710, IN_POLICY]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P795",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    p794 = load_json(IN_P794)
    f794 = load_json(IN_F794)
    f704_summary = load_json(IN_F704_SUMMARY)
    policy = load_json(IN_POLICY)

    f704_packet_text = IN_F704_PACKET.read_text(encoding="utf-8")
    n705_text = IN_N705.read_text(encoding="utf-8")
    n703_text = IN_N703.read_text(encoding="utf-8")
    p710_text = IN_P710.read_text(encoding="utf-8")

    semantic_role_artifacts = find_semantic_role_artifacts()

    checks = [
        {
            "id": "f704_tracks_max_numerically_but_exports_whole_spectrum_object",
            "pass": (
                f704_summary.get("stage") == "F704"
                and f704_summary.get("status") == "PASS_EXPORTED_STRICT_INVARIANT_MASS_OBSERVABLE_OBJECT"
                and f704_summary.get("max_m2_proxy") is not None
                and contains_all(
                    f704_packet_text,
                    [
                        "observer-limit observable",
                        "basis-invariant",
                        "mass-like",
                    ],
                )
            ),
            "details": "F704 tracks extrema numerically, but the exported strict object is the whole basis-invariant spectrum rather than a privileged top-boundary role object.",
        },
        {
            "id": "n705_tracks_f704_as_os_component_not_as_boundary_role",
            "pass": contains_all(
                n705_text,
                [
                    "includes the strict basis-invariant mass observable object",
                    "tracked component",
                ],
            ) and not contains_all(n705_text, ["top boundary", "point 1"]),
            "details": "N705 upgrades F704 into the Release-7 OS bundle, but not into any explicit top-boundary semantic anchor.",
        },
        {
            "id": "n703_blocks_host_or_particle_overread_of_boundary_point",
            "pass": contains_all(
                n703_text,
                [
                    "dimensionless quadratic coefficients",
                    "not yet physical masses in gev",
                    "do not constitute a standard model match claim",
                ],
            ),
            "details": "N703 keeps the whole proxy layer in internal dimensionless meaning scope, preventing silent upgrade of the normalized point 1 into host semantics.",
        },
        {
            "id": "no_other_exported_semantic_role_artifact_detected",
            "pass": not semantic_role_artifacts,
            "details": "No other exported repo artifact assigns a distinct alpha_s-side semantic role to the top-boundary point 1.",
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
            "details": "The proxy-to-GeV calibration lane is explicitly nonstrict and cannot supply the missing semantic role.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    reusable_semantic_role_found = bool(semantic_role_artifacts)

    if (
        p794.get("status")
        == "P794_CURRENT_STRICT_ALPHA_S_BOUNDED_GRID_CANDIDATE_SUPPORTED_TOP_BOUNDARY_ANCHOR_BLOCKED"
        and not failed_checks
        and not reusable_semantic_role_found
    ):
        status = "P795_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_SEMANTIC_ROLE_REUSE_NEGATIVE_ON_CURRENT_REPO_STATE"
    else:
        status = "P795_REQUIRES_REVIEW"

    artifact = {
        "stage": "P795",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "p794_subclause_audit": rel(IN_P794),
            "f794_top_boundary_anchor_target": rel(IN_F794),
            "f704_packet": rel(IN_F704_PACKET),
            "f704_summary": rel(IN_F704_SUMMARY),
            "n705_release7_os_theorem": rel(IN_N705),
            "n703_meaning_theorem": rel(IN_N703),
            "p710_nonstrict_proxy_to_gev_probe": rel(IN_P710),
            "proxy_to_gev_policy": rel(IN_POLICY),
        },
        "semantic_role_artifacts_detected": semantic_role_artifacts,
        "reusable_semantic_role_found": reusable_semantic_role_found,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The current repo tracks the top of the F704 spectrum numerically, but only inside a whole-spectrum strict object.",
            "No current strict-side export upgrades the normalized point 1 into a semantic alpha_s boundary role.",
            "So the sharp blocker is now the missing semantic-role rule itself, not the existence of the point or the arithmetic boundedness of the grid.",
        ],
        "recommended_next_packet": {
            "id": "F795_CURRENT_STRICT_ALPHA_S_TOP_BOUNDARY_SEMANTIC_ROLE_TARGET_PACKET",
            "goal": "Freeze the exact semantic-role object still missing before the normalized boundary point 1 can count as anything more than a normalization artifact.",
            "minimum_fields": [
                "candidate_family_domain_ref",
                "boundary_point_ref",
                "supporting_strict_source_ref",
                "boundary_point_semantic_role_rule_ref",
                "nonstrict_calibration_exclusion_ref",
                "selected_role_output_schema",
                "hard_limits",
            ],
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P795",
        "status": status,
        "as_of": AS_OF,
        "reusable_semantic_role_found": reusable_semantic_role_found,
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
