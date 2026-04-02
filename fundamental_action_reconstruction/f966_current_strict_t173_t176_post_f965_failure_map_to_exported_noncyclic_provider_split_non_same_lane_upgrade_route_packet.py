#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-28"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1091 = GENERATED / "p1091_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_decision_audit_probe_summary.json"
OUT_JSON = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet.json"
OUT_SUMMARY = GENERATED / "f966_current_strict_t173_t176_post_f965_failure_map_to_exported_noncyclic_provider_split_non_same_lane_upgrade_route_packet_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    if not IN_P1091.exists():
        artifact = {
            "stage": "F966",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_P1091.relative_to(REPO))],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1091 = load_json(IN_P1091)
    admitted = p1091.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_DECISION_AUDITED"
    status = (
        "F966_EXECUTED_CURRENT_STRICT_T173_T176_POST_F965_FAILURE_MAP_TO_EXPORTED_NONCYCLIC_PROVIDER_SPLIT_NON_SAME_LANE_UPGRADE_ROUTE_PACKET_NO_FALSE_PASS"
        if admitted
        else "F966_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_ROUTE_STATE"
    )

    artifact = {
        "stage": "F966",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1091_route_decision_audit_probe_summary": str(IN_P1091.relative_to(REPO)),
        },
        "exported_route_packet": {
            "object_id": "PostF965FailureMapConstrainedToExportedNoncyclicProviderSplitNonSameLaneUpgradeRoute_v1",
            "preferred_search_family": "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1",
            "allowed_next_move_contract": "search_one_genuinely_new_non_same_lane_upgrade_route_within_exported_noncyclic_provider_split_family_or_one_genuinely_new_inversion_sensitive_source_side_provider_class",
            "active_missing_bridge": p1091.get("active_missing_bridge"),
            "typed_meaning": (
                "one strict route packet freezing that the post-F965 honest search space stays on the exported noncyclic provider-split family, "
                "with all already stopped same-lane reentries disallowed, so continuation must be a genuinely new non-same-lane upgrade route or a genuinely new inversion-sensitive source-side provider class"
            ),
        },
        "admission_properties": {
            "process_failure_map_frozen": True,
            "pair12_entry_same_lane_reentry_disallowed_as_primary_move": bool(p1091.get("pair12_entry_same_lane_reentry_disallowed_as_primary_move")),
            "pair_side_sharper_same_lane_reentry_disallowed_as_primary_move": bool(p1091.get("pair_side_sharper_same_lane_reentry_disallowed_as_primary_move")),
            "feeder_sharper_same_lane_reentry_disallowed_as_primary_move": bool(p1091.get("feeder_sharper_same_lane_reentry_disallowed_as_primary_move")),
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "hard_limits": [
            "Does not export a lawful supplier.",
            "Does not export a strict physical orientation datum.",
            "Does not discharge T183.",
            "Does not discharge T176.",
            "Does not discharge QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "route_object_id": artifact["exported_route_packet"]["object_id"],
        "preferred_search_family": artifact["exported_route_packet"]["preferred_search_family"],
        "allowed_next_move_contract": artifact["exported_route_packet"]["allowed_next_move_contract"],
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
