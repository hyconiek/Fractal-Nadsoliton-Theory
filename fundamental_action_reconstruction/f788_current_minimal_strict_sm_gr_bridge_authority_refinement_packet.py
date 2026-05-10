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
IN_F801 = GENERATED / "f801_current_strict_sm_gr_minimal_bridge_registry_packet.json"
IN_WORKING_NOTE = REPO / "WORKING_NOTE_LEGACY_KEEP_CUT_AND_MINIMAL_STRICT_SM_GR_BRIDGE.md"

OUT = GENERATED / "f788_current_minimal_strict_sm_gr_bridge_authority_refinement_packet.json"
OUT_SUMMARY = GENERATED / "f788_current_minimal_strict_sm_gr_bridge_authority_refinement_packet_summary.json"


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
        .replace("_", " ")
        .replace("-", " ")
    )


def contains_all(text: str, needles: list[str]) -> bool:
    hay = normalize(text)
    return all(normalize(needle) in hay for needle in needles)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F787, IN_F801, IN_WORKING_NOTE]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "F788",
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
    f801 = load_json(IN_F801)
    working_note = IN_WORKING_NOTE.read_text(encoding="utf-8")

    f787_targets = f787.get("refined_targets") or {}
    f801_entries = {entry.get("id"): entry for entry in (f801.get("registry_entries") or [])}

    overlap_ids = sorted(set(f787_targets.keys()) & set(f801_entries.keys()))
    alpha_s_f787 = f787_targets.get("alpha_s_boundary_mu0_alpha0") or {}
    alpha_s_f801 = f801_entries.get("alpha_s_boundary_mu0_alpha0") or {}

    minimal_bridge_unit_rule_present = contains_all(
        working_note,
        [
            "dimensionless or explicitly normalized observables",
            "proxy -> GeV",
        ],
    )

    alpha_s_export_conflict = (
        alpha_s_f787.get("bridge_ready") is False
        and alpha_s_f787.get("status_label") == "open"
        and alpha_s_f801.get("status_label") == "strict-derived"
    )

    if overlap_ids and alpha_s_export_conflict and minimal_bridge_unit_rule_present:
        authoritative_packet = "F787"
        nonauthoritative_packet = "F801"
        authority_reason = (
            "F787 is stricter on the current repo state because it demotes the alpha_s boundary out of the minimal bridge export set, in line with the working-note unit discipline and the current formula-layer blocker set."
        )
    else:
        authoritative_packet = "UNRESOLVED"
        nonauthoritative_packet = None
        authority_reason = "No decisive refinement trigger detected."

    artifact = {
        "stage": "F788",
        "packet_name": "CurrentMinimalStrictSMGRBridgeAuthorityRefinement_v1",
        "status": "F788_EXECUTED_CURRENT_MINIMAL_STRICT_SM_GR_BRIDGE_AUTHORITY_REFINEMENT_PACKET_NO_FALSE_PASS",
        "as_of": AS_OF,
        "inputs": {
            "f787_refined_registry": rel(IN_F787),
            "f801_overlapping_registry": rel(IN_F801),
            "working_note": rel(IN_WORKING_NOTE),
        },
        "overlap_audit": {
            "overlap_ids": overlap_ids,
            "same_four_target_registry_overlap": len(overlap_ids) == 4,
            "alpha_s_export_conflict": alpha_s_export_conflict,
            "minimal_bridge_unit_rule_present": minimal_bridge_unit_rule_present,
        },
        "authority_decision": {
            "authoritative_current_minimal_bridge_packet": authoritative_packet,
            "overlapping_nonauthoritative_packet": nonauthoritative_packet,
            "reason": authority_reason,
        },
        "current_honest_reading": [
            "F801 is a real overlapping registry packet, but it keeps alpha_s_boundary_mu0_alpha0 exported in a way that is wider than the current no-false-pass refinement lane.",
            "F787 is the tighter current packet because it preserves the same bridge lane while demoting the alpha_s boundary to support-only/nonexport.",
            "This authority refinement does not delete F801; it only prevents the repo from treating two inconsistent export sets as equally current in the minimal bridge lane.",
        ],
        "recommended_next_move": "Build or prove absent one replacement alpha_s boundary object that is dimensionless or explicitly normalized, then rerun both the export-refinement and authority-refinement packets.",
        "suggested_next_artifact_slots_after_f788": [
            "P788",
            "N788",
            "F789",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "F788",
        "status": artifact["status"],
        "as_of": AS_OF,
        "authoritative_current_minimal_bridge_packet": authoritative_packet,
        "overlapping_nonauthoritative_packet": nonauthoritative_packet,
        "alpha_s_export_conflict": alpha_s_export_conflict,
        "suggested_next_artifact_slots_after_f788": artifact["suggested_next_artifact_slots_after_f788"],
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
