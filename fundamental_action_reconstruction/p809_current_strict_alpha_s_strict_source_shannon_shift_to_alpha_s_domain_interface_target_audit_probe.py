#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F789 = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"
IN_F800 = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"
IN_F808 = GENERATED / "f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.json"
IN_F805 = GENERATED / "f805_current_strict_alpha_s_acting_input_bundle_packet.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_P748 = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
IN_P749 = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"
IN_P808 = GENERATED / "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.json"

OUT = GENERATED / "p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe.json"
OUT_SUMMARY = GENERATED / "p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe_summary.json"


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


def find_repo_hits(candidate_id: str, target_needles: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.json",
        "generated/f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet_summary.json",
        "generated/p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.json",
        "generated/p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe_summary.json",
        "F808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET.md",
        "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMISSION_PROBE.md",
        "f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.py",
        "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.py",
        "generated/p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe.json",
        "generated/p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe_summary.json",
        "P809_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_AUDIT_PROBE.md",
        "p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe.py",
        "generated/f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet.json",
        "generated/f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet_summary.json",
        "F809_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_PACKET.md",
        "f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet.py",
    ]
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in {".md", ".json", ".py"}:
            continue
        relpath = rel(path)
        # Exclude the current reference-lane exports and the present P809/F809 files.
        if any(relpath.endswith(excluded) for excluded in excluded_suffixes):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        hay = normalize(text)
        if normalize(candidate_id) in hay and all(normalize(needle) in hay for needle in target_needles):
            hits.append(relpath)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F789,
        IN_F800,
        IN_F808,
        IN_F805,
        IN_P788,
        IN_P748,
        IN_P749,
        IN_P808,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P809",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f789 = load_json(IN_F789)
    f800 = load_json(IN_F800)
    f808 = load_json(IN_F808)
    f805 = load_json(IN_F805)
    p788 = load_json(IN_P788)
    p748 = load_json(IN_P748)
    p749 = load_json(IN_P749)
    p808 = load_json(IN_P808)

    target_interface = f789.get("target_interface") or {}
    provider_target = f800.get("target_object") or {}
    candidate_lane = f808.get("exported_object") or {}
    acting_input = f805.get("exported_object") or {}
    p808_theorem = p808.get("theorem_result") or {}

    exact_interface_target_scan_hits = find_repo_hits(
        "strict_source_shannon_provider_shift_candidate_reference_lane_v1",
        [
            "alpha_s",
            "interface target",
        ],
    )

    checks = [
        {
            "id": "f808_exports_only_reference_candidate_lane_with_domain_block",
            "pass": (
                f808.get("status")
                == "F808_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET_NO_FALSE_PASS"
                and candidate_lane.get("object_id")
                == "alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_v1"
                and candidate_lane.get("candidate_grade") == "reference_context_candidate_only"
                and candidate_lane.get("alpha_s_domain_interface_status") == "blocked_nonexport"
                and p808_theorem.get("alpha_s_domain_interface_exported") is False
            ),
            "details": "F808 exports the Shannon shift lane only at reference-candidate grade and keeps alpha_s-domain interface blocked.",
        },
        {
            "id": "f789_exports_only_generic_alpha_s_boundary_interface_target",
            "pass": (
                f789.get("status")
                == "F789_EXECUTED_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and target_interface.get("object_id") == "alpha_s_normalized_boundary_interface_target_v1"
                and "normalized alpha_s boundary interface object"
                in (p788.get("recommended_next_packet") or {}).get("goal", "")
            ),
            "details": "F789 names only the generic normalized alpha_s boundary interface target, not the shift-side domain-entry interface from the Shannon candidate lane.",
        },
        {
            "id": "strict_source_shannon_route_remains_unbridged_and_nonentering_on_its_own_lanes",
            "pass": (
                p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
                is False
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry")
                is False
            ),
            "details": "The strict-source Shannon route remains unbridged and nonentering even on its own declared future carriers.",
        },
        {
            "id": "alpha_s_provider_target_still_requires_exact_carrier_safe_interface",
            "pass": (
                provider_target.get("object_id") == "alpha_s_reference_scale_provider_class_target_v1"
                and acting_input.get("object_id") == "alpha_s_reference_scale_acting_input_bundle_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "same_domain_carrier_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "foreign_domain_reuse_block_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
            ),
            "details": "The alpha_s provider-class target still demands exact carrier-safe same-domain admission and explicitly fences off foreign-domain reuse.",
        },
        {
            "id": "no_exact_shift_to_alpha_s_domain_interface_target_exported",
            "pass": exact_interface_target_scan_hits == [],
            "details": "No exported object currently names the exact Shannon shift-candidate -> alpha_s domain interface target.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]

    theorem_result = {
        "generic_alpha_s_boundary_interface_target_exists": checks[1]["pass"],
        "exact_strict_source_shannon_shift_to_alpha_s_domain_interface_target_exported": False
        if all(item["pass"] for item in checks)
        else None,
        "next_honest_move_requires_freezing_exact_shift_to_domain_interface_target": all(
            item["pass"] for item in checks
        ),
        "no_false_pass": True,
    }

    if all(item["pass"] for item in checks):
        status = "P809_CURRENT_STRICT_ALPHA_S_NO_EXACT_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
    else:
        status = "P809_REQUIRES_REVIEW"

    artifact = {
        "stage": "P809",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f789_generic_alpha_s_interface_target": rel(IN_F789),
            "f800_provider_class_target": rel(IN_F800),
            "f805_acting_input_bundle": rel(IN_F805),
            "f808_provider_shift_candidate_reference_packet": rel(IN_F808),
            "p788_alpha_s_replacement_route_probe": rel(IN_P788),
            "p748_shannon_pair12_bridge_boundary": rel(IN_P748),
            "p749_shannon_pair_indexed_entry_boundary": rel(IN_P749),
            "p808_shift_candidate_reference_probe": rel(IN_P808),
        },
        "exact_missing_interface_target_candidate": {
            "candidate_id": "alpha_s_strict_source_shannon_shift_to_domain_interface_target_missing_v1",
            "generic_alpha_s_interface_target_ref": target_interface.get("object_id"),
            "provider_shift_candidate_reference_lane_ref": candidate_lane.get("object_id"),
            "current_alpha_s_acting_input_bundle_ref": acting_input.get("object_id"),
            "repo_scan_hits_for_exact_interface_target": exact_interface_target_scan_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo does export generic alpha_s interface-target language via F789.",
            "But it still does not export the exact interface target that would carry the strict-source Shannon provider-shift candidate lane into the alpha_s domain.",
            "So the next honest move is to freeze that exact missing target object, not to pretend F789 already solves the shift-lane entry problem.",
        ],
        "recommended_next_packet": {
            "id": "F809_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_PACKET",
            "goal": "Freeze the exact missing target object for the Shannon provider-shift candidate lane to alpha_s domain interface, without claiming interface realization.",
            "export_object_id": "alpha_s_strict_source_shannon_shift_to_domain_interface_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P809",
        "status": status,
        "as_of": AS_OF,
        "generic_alpha_s_boundary_interface_target_exists": theorem_result[
            "generic_alpha_s_boundary_interface_target_exists"
        ],
        "exact_strict_source_shannon_shift_to_alpha_s_domain_interface_target_exported": theorem_result[
            "exact_strict_source_shannon_shift_to_alpha_s_domain_interface_target_exported"
        ],
        "next_honest_move_requires_freezing_exact_shift_to_domain_interface_target": theorem_result[
            "next_honest_move_requires_freezing_exact_shift_to_domain_interface_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
