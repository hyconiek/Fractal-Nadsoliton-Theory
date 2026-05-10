#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F787 = GENERATED / "f787_current_minimal_strict_sm_gr_bridge_export_refinement_packet.json"
IN_F788 = GENERATED / "f788_current_minimal_strict_sm_gr_bridge_authority_refinement_packet.json"
IN_P788 = GENERATED / "p788_current_alpha_s_dimensionless_or_normalized_replacement_route_probe.json"
IN_WORKING_NOTE = REPO / "WORKING_NOTE_LEGACY_KEEP_CUT_AND_MINIMAL_STRICT_SM_GR_BRIDGE.md"

OUT = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet.json"
OUT_SUMMARY = GENERATED / "f789_current_strict_alpha_s_normalized_boundary_interface_target_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def normalize(text: str) -> str:
    return (
        text.lower()
        .replace("“", '"')
        .replace("”", '"')
        .replace("’", "'")
        .replace("->", " ")
        .replace("→", " ")
        .replace("/", " ")
        .replace("-", " ")
        .replace("_", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F787, IN_F788, IN_P788, IN_WORKING_NOTE]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F789",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        OUT_SUMMARY.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT_SUMMARY)
        return

    f787 = load_json(IN_F787)
    f788 = load_json(IN_F788)
    p788 = load_json(IN_P788)
    working_note = IN_WORKING_NOTE.read_text(encoding="utf-8")

    alpha_s_refined = (f787.get("refined_targets") or {}).get("alpha_s_boundary_mu0_alpha0") or {}
    authority_ok = (
        ((f788.get("authority_decision") or {}).get("authoritative_current_minimal_bridge_packet") == "F787")
        and (((f788.get("overlap_audit") or {}).get("alpha_s_export_conflict")) is True)
    )
    alpha_s_demoted = (
        alpha_s_refined.get("status_label") == "open"
        and alpha_s_refined.get("bridge_ready") is False
    )
    unit_rule_present = contains_all(
        working_note,
        [
            "dimensionless or explicitly normalized observables",
            "proxy -> GeV",
        ],
    )

    target_fields = [
        {
            "name": "mu0_tilde",
            "type": "dimensionless_or_explicitly_normalized_boundary_scale",
            "required": True,
            "hard_limit": "Must replace mu0_gev inside the canonical strict alpha_s boundary interface.",
        },
        {
            "name": "alpha_s_mu0",
            "type": "dimensionless_coupling_value",
            "required": True,
            "hard_limit": "Must remain dimensionless; no unit-calibration semantics follow from this field alone.",
        },
        {
            "name": "n_f_active_at_mu0",
            "type": "integer_sector_selector",
            "required": True,
            "hard_limit": "May classify the running sector but does not by itself close QCD.",
        },
        {
            "name": "normalized_validation_points",
            "type": "list_of_dimensionless_or_explicitly_normalized_scales",
            "required": True,
            "hard_limit": "Must replace mu_gev validation samples in the canonical interface.",
        },
        {
            "name": "normalization_rule_ref",
            "type": "explicit_reference_to_exported_normalization_rule_or_adapter",
            "required": True,
            "hard_limit": "No hidden GeV translation is allowed inside this rule.",
        },
        {
            "name": "strict_input_chain",
            "type": "ordered_list_of_exported_strict_source_objects",
            "required": True,
            "hard_limit": "Must not silently promote QW-2063 mass-chain objects into strict first-principles if the repo still marks them below that level.",
        },
        {
            "name": "hard_limits",
            "type": "list_of_nonpromotion_boundaries",
            "required": True,
            "hard_limit": "Must explicitly deny QCD closure, SM identification, and ToE closure.",
        },
    ]

    target_rules = [
        "No field named mu0_gev or any equivalent GeV boundary scale may appear in the canonical interface object.",
        "No GeV-only validation grid may appear in the canonical interface object.",
        "A real adapter from strict mass proxies or explicitly normalized strict observables into the alpha_s boundary schema must be exported explicitly.",
        "Any downstream running consumer must accept the normalized boundary directly or be fenced as an external translation layer.",
        "The interface must preserve the current no-false-pass boundary from F787/F788.",
    ]

    if (
        authority_ok
        and alpha_s_demoted
        and unit_rule_present
        and p788.get("status") == "P788_CURRENT_ALPHA_S_DIMENSIONLESS_OR_NORMALIZED_REPLACEMENT_ROUTE_BLOCKED_ON_CURRENT_REPO_STATE"
    ):
        status = "F789_EXECUTED_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_TARGET_PACKET_NO_FALSE_PASS"
    else:
        status = "F789_REQUIRES_REVIEW"

    artifact = {
        "stage": "F789",
        "packet_name": "CurrentStrictAlphaSNormalizedBoundaryInterfaceTarget_v1",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f787_export_refinement": rel(IN_F787),
            "f788_authority_refinement": rel(IN_F788),
            "p788_replacement_route_probe": rel(IN_P788),
            "working_note": rel(IN_WORKING_NOTE),
        },
        "renumbering_note": {
            "p788_suggested_packet_id": "F788_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_TARGET_PACKET",
            "actual_packet_id": "F789_CURRENT_STRICT_ALPHA_S_NORMALIZED_BOUNDARY_INTERFACE_TARGET_PACKET",
            "reason": "Slot F788 is already occupied by the minimal bridge authority refinement packet on current repo state.",
        },
        "why_this_packet_exists": [
            "F787 demotes alpha_s_boundary_mu0_alpha0 out of the minimal export set.",
            "F788 fixes F787 as the authoritative current minimal bridge packet.",
            "P788 shows that no dimensionless or explicitly normalized replacement route is currently exported.",
        ],
        "blocked_by_current_repo_state": p788.get("missing_interfaces") or [],
        "target_interface": {
            "object_id": "alpha_s_normalized_boundary_interface_target_v1",
            "goal": "Specify the exact canonical interface that must exist before alpha_s can re-enter the minimal strict bridge without GeV-level boundary semantics.",
            "required_fields": target_fields,
            "target_rules": target_rules,
        },
        "current_honest_reading": [
            "The repo is no longer blocked on vague alpha_s language; it is now blocked on one explicit normalized boundary interface target.",
            "This packet does not claim the interface already exists; it only freezes the exact object that must be built next.",
            "Using F789 is methodologically safer than reopening kernel parameter fitting, because it removes one missing interface rather than compensating that gap by retuning.",
        ],
        "recommended_next_move": "Build one candidate interface probe that instantiates mu0_tilde plus normalized_validation_points from exported strict objects, without introducing GeV-level fields into the canonical alpha_s boundary object.",
        "hard_limits": [
            "Does not claim that a normalized alpha_s boundary interface is already exported.",
            "Does not claim QCD closure.",
            "Does not claim Standard Model identification.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F789",
        "status": status,
        "as_of": AS_OF,
        "actual_packet_id": artifact["renumbering_note"]["actual_packet_id"],
        "blocked_interface_count": len(artifact["blocked_by_current_repo_state"]),
        "recommended_next_move": artifact["recommended_next_move"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
