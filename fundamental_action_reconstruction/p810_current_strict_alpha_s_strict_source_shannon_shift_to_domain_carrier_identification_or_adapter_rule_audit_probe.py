#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F809 = GENERATED / "f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet.json"
IN_F808 = GENERATED / "f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.json"
IN_F800 = GENERATED / "f800_current_strict_alpha_s_reference_scale_provider_class_target_packet.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_P748 = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
IN_P749 = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"

OUT = GENERATED / "p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe.json"
OUT_SUMMARY = GENERATED / "p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe_summary.json"


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


def check_pass(artifact: dict[str, Any], check_id: str) -> bool:
    for item in artifact.get("checks", []):
        if item.get("id") == check_id:
            return bool(item.get("pass"))
    return False


def find_exact_rule_hits(interface_target_id: str, rule_terms: list[str]) -> list[str]:
    hits: list[str] = []
    excluded_suffixes = [
        "generated/f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.json",
        "generated/f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet_summary.json",
        "generated/p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.json",
        "generated/p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe_summary.json",
        "generated/f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet.json",
        "generated/f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet_summary.json",
        "generated/p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe.json",
        "generated/p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe_summary.json",
        "P808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMISSION_PROBE.md",
        "F808_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET.md",
        "p808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_admission_probe.py",
        "f808_current_strict_alpha_s_strict_source_shannon_provider_shift_candidate_reference_packet.py",
        "P809_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_AUDIT_PROBE.md",
        "F809_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_PACKET.md",
        "p809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_audit_probe.py",
        "f809_current_strict_alpha_s_strict_source_shannon_shift_to_alpha_s_domain_interface_target_packet.py",
        "generated/p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe.json",
        "generated/p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe_summary.json",
        "generated/f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.json",
        "generated/f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet_summary.json",
        "P810_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_AUDIT_PROBE.md",
        "F810_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET.md",
        "p810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_audit_probe.py",
        "f810_current_strict_alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_packet.py",
    ]
    norm_interface = normalize(interface_target_id)
    norm_terms = [normalize(term) for term in rule_terms]
    for path in ROOT.rglob("*"):
        if not path.is_file():
            continue
        if path.suffix not in {".md", ".json", ".py"}:
            continue
        relpath = rel(path)
        if any(relpath.endswith(excluded) for excluded in excluded_suffixes):
            continue
        try:
            text = path.read_text(encoding="utf-8")
        except UnicodeDecodeError:
            continue
        hay = normalize(text)
        if norm_interface in hay and any(term in hay for term in norm_terms):
            hits.append(relpath)
    return sorted(hits)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [
        IN_F809,
        IN_F808,
        IN_F800,
        IN_P788,
        IN_P748,
        IN_P749,
    ]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P810",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f809 = load_json(IN_F809)
    f808 = load_json(IN_F808)
    f800 = load_json(IN_F800)
    p788 = load_json(IN_P788)
    p748 = load_json(IN_P748)
    p749 = load_json(IN_P749)

    target_object = f809.get("target_object") or {}
    target_refs = f809.get("target_refs") or {}
    candidate_lane = f808.get("exported_object") or {}
    provider_target = f800.get("target_object") or {}

    exact_rule_scan_hits = find_exact_rule_hits(
        "alpha_s_strict_source_shannon_shift_to_domain_interface_target_v1",
        [
            "carrier identification",
            "adapter rule",
            "transport rule",
            "adapter",
        ],
    )

    checks = [
        {
            "id": "f809_freezes_exact_shift_to_domain_interface_target_and_requires_adapter_rule",
            "pass": (
                f809.get("status")
                == "F809_EXECUTED_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_ALPHA_S_DOMAIN_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
                and target_object.get("object_id")
                == "alpha_s_strict_source_shannon_shift_to_domain_interface_target_v1"
                and any(
                    isinstance(item, dict)
                    and item.get("name") == "carrier_identification_or_adapter_rule_ref"
                    for item in (target_object.get("required_fields") or [])
                )
            ),
            "details": "F809 already freezes the exact shift-to-domain interface target and names the still-missing carrier-identification-or-adapter rule slot.",
        },
        {
            "id": "p788_still_reports_no_exported_alpha_s_adapter",
            "pass": (
                p788.get("status")
                == "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"
                and check_pass(p788, "no_exported_dimensionless_alpha_s_adapter_detected")
                and (p788.get("adapter_token_hits") or {}) == {}
            ),
            "details": "P788 still audits the canonical alpha_s lane as lacking any exported adapter from strict normalized observables into the alpha_s boundary schema.",
        },
        {
            "id": "strict_source_shannon_route_remains_unbridged_and_nonentering_on_current_exports",
            "pass": (
                p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
                is False
                and p749.get("current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry")
                is False
            ),
            "details": "The strict-source Shannon route remains unbridged and nonentering even on its own current exported lanes.",
        },
        {
            "id": "alpha_s_provider_lane_still_demands_same_domain_carrier_discipline",
            "pass": (
                f800.get("status")
                == "F800_EXECUTED_CURRENT_STRICT_ALPHA_S_REFERENCE_SCALE_PROVIDER_CLASS_TARGET_PACKET_NO_FALSE_PASS"
                and candidate_lane.get("object_id")
                == "alpha_s_strict_source_shannon_provider_shift_candidate_reference_lane_v1"
                and any(
                    isinstance(item, dict) and item.get("name") == "same_domain_carrier_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
                and any(
                    isinstance(item, dict) and item.get("name") == "foreign_domain_reuse_block_ref"
                    for item in (provider_target.get("required_fields") or [])
                )
            ),
            "details": "The alpha_s provider lane still requires exact same-domain carrier discipline and explicitly blocks foreign-domain reuse.",
        },
        {
            "id": "no_exported_rule_names_exact_shannon_shift_to_domain_carrier_adapter",
            "pass": exact_rule_scan_hits == [],
            "details": "No current export names one exact carrier-identification or adapter rule for the frozen strict-source Shannon -> alpha_s domain interface target.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "exact_shift_to_domain_interface_target_frozen": checks[0]["pass"],
        "exported_alpha_s_adapter_rule_present": not checks[1]["pass"],
        "strict_source_shannon_side_bridge_or_entry_present": not checks[2]["pass"],
        "exact_carrier_identification_or_adapter_rule_exported_for_f809_target": False if all_pass else None,
        "next_honest_move_requires_freezing_exact_carrier_identification_or_adapter_rule_target": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P810_CURRENT_STRICT_ALPHA_S_NO_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_FOR_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_INTERFACE_TARGET_EXPORTED_TARGET_FREEZE_REQUIRED"
        if all_pass
        else "P810_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P810",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f809_shift_to_domain_interface_target_packet": rel(IN_F809),
            "f808_provider_shift_candidate_reference_packet": rel(IN_F808),
            "f800_provider_class_target_packet": rel(IN_F800),
            "p788_alpha_s_adapter_route_probe": rel(IN_P788),
            "p748_strict_source_shannon_bridge_boundary": rel(IN_P748),
            "p749_strict_source_shannon_entry_boundary": rel(IN_P749),
        },
        "exact_missing_rule_target_candidate": {
            "candidate_id": "alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_missing_v1",
            "shift_to_domain_interface_target_ref": target_object.get("object_id"),
            "provider_shift_candidate_reference_lane_ref": target_refs.get("provider_shift_candidate_reference_lane_ref"),
            "target_alpha_s_acting_input_bundle_ref": target_refs.get("target_alpha_s_acting_input_bundle_ref"),
            "repo_scan_hits_for_exact_rule": exact_rule_scan_hits,
        },
        "theorem_result": theorem_result,
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The repo now freezes the exact strict-source Shannon -> alpha_s domain interface target via F809.",
            "But the canonical alpha_s lane still exports no adapter rule, and the strict-source Shannon lane still exports no bridge or entry even on its own side.",
            "So the next honest move is to freeze the exact missing carrier-identification or adapter-rule target, not to pretend that lawful domain admission already exists.",
        ],
        "recommended_next_packet": {
            "id": "F810_CURRENT_STRICT_ALPHA_S_STRICT_SOURCE_SHANNON_SHIFT_TO_DOMAIN_CARRIER_IDENTIFICATION_OR_ADAPTER_RULE_TARGET_PACKET",
            "goal": "Freeze the exact missing carrier-identification or adapter-rule target required to instantiate the frozen F809 shift-to-domain interface target without silent domain identification.",
            "export_object_id": "alpha_s_strict_source_shannon_shift_to_domain_carrier_identification_or_adapter_rule_target_v1",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P810",
        "status": status,
        "as_of": AS_OF,
        "exact_shift_to_domain_interface_target_frozen": theorem_result[
            "exact_shift_to_domain_interface_target_frozen"
        ],
        "exact_carrier_identification_or_adapter_rule_exported_for_f809_target": theorem_result[
            "exact_carrier_identification_or_adapter_rule_exported_for_f809_target"
        ],
        "next_honest_move_requires_freezing_exact_carrier_identification_or_adapter_rule_target": theorem_result[
            "next_honest_move_requires_freezing_exact_carrier_identification_or_adapter_rule_target"
        ],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
