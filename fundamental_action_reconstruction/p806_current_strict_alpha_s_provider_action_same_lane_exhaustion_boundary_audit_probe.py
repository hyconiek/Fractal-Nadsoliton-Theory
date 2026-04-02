#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-19"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P802 = GENERATED / "p802_current_strict_alpha_s_provider_class_reorganization_audit_probe.json"
IN_P803 = GENERATED / "p803_current_strict_alpha_s_same_domain_relation_bundle_admission_probe.json"
IN_P804 = GENERATED / "p804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_admission_probe.json"
IN_P805 = GENERATED / "p805_current_strict_alpha_s_acting_input_bundle_admission_probe.json"
IN_F801 = GENERATED / "f801_current_strict_alpha_s_same_domain_provider_skeleton_packet.json"
IN_F802 = GENERATED / "f802_current_strict_alpha_s_provider_action_rule_target_packet.json"
IN_F803 = GENERATED / "f803_current_strict_alpha_s_same_domain_relation_bundle_packet.json"
IN_F804 = GENERATED / "f804_current_strict_alpha_s_same_domain_compositional_alignment_bundle_packet.json"
IN_F805 = GENERATED / "f805_current_strict_alpha_s_acting_input_bundle_packet.json"
IN_S2 = ROOT / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md"

OUT_JSON = GENERATED / "p806_current_strict_alpha_s_provider_action_same_lane_exhaustion_boundary_audit_probe.json"
OUT_SUMMARY = GENERATED / "p806_current_strict_alpha_s_provider_action_same_lane_exhaustion_boundary_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_text(path: Path) -> str:
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_P802,
        IN_P803,
        IN_P804,
        IN_P805,
        IN_F801,
        IN_F802,
        IN_F803,
        IN_F804,
        IN_F805,
        IN_S2,
    ]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P806",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p802 = load_json(IN_P802)
    p803 = load_json(IN_P803)
    p804 = load_json(IN_P804)
    p805 = load_json(IN_P805)
    f801 = load_json(IN_F801)
    f802 = load_json(IN_F802)
    f803 = load_json(IN_F803)
    f804 = load_json(IN_F804)
    f805 = load_json(IN_F805)
    s2_text = load_text(IN_S2)

    checks: list[dict[str, Any]] = []
    blocking_checks: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking_checks.append(check_id)

    passive_same_lane_stack_already_exported = (
        f801.get("status")
        == "F801_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_PROVIDER_SKELETON_PACKET_NO_FALSE_PASS"
        and (f801.get("exported_object") or {}).get("object_id")
        == "alpha_s_reference_scale_same_domain_provider_skeleton_v1"
        and f803.get("status")
        == "F803_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_RELATION_BUNDLE_PACKET_NO_FALSE_PASS"
        and (f803.get("exported_object") or {}).get("object_id")
        == "alpha_s_reference_scale_same_domain_relation_bundle_v1"
        and f804.get("status")
        == "F804_EXECUTED_CURRENT_STRICT_ALPHA_S_SAME_DOMAIN_COMPOSITIONAL_ALIGNMENT_BUNDLE_PACKET_NO_FALSE_PASS"
        and (f804.get("exported_object") or {}).get("object_id")
        == "alpha_s_reference_scale_same_domain_compositional_alignment_bundle_v1"
        and f805.get("status")
        == "F805_EXECUTED_CURRENT_STRICT_ALPHA_S_ACTING_INPUT_BUNDLE_PACKET_NO_FALSE_PASS"
        and (f805.get("exported_object") or {}).get("object_id")
        == "alpha_s_reference_scale_acting_input_bundle_v1"
    )

    blocker_stays_unchanged_across_same_lane_split = (
        (p802.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_action_rule_ref"
        and (p803.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_action_rule_ref"
        and (p804.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_action_rule_ref"
        and (p805.get("clause_split") or {}).get("sharp_blocker_clause") == "provider_action_rule_ref"
        and (f802.get("target_object") or {}).get("object_id")
        == "alpha_s_reference_scale_provider_action_rule_target_v1"
    )

    current_alpha_s_provider_action_same_lane_split_is_exhaustive = (
        passive_same_lane_stack_already_exported
        and blocker_stays_unchanged_across_same_lane_split
        and (p803.get("clause_split") or {}).get("relation_bundle_clause_status") == "export_admitted"
        and (p804.get("clause_split") or {}).get("alignment_bundle_clause_status") == "export_admitted"
        and (p805.get("clause_split") or {}).get("acting_input_bundle_clause_status") == "export_admitted"
    )

    current_alpha_s_provider_action_same_lane_remains_nonexport = (
        (p802.get("clause_split") or {}).get("active_provider_action_rule_clause_status")
        == "blocked_nonexport"
        and (p803.get("clause_split") or {}).get("provider_action_rule_clause_status")
        == "blocked_nonexport"
        and (p804.get("clause_split") or {}).get("provider_action_rule_clause_status")
        == "blocked_nonexport"
        and (p805.get("clause_split") or {}).get("provider_action_rule_clause_status")
        == "blocked_nonexport"
    )

    current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole = not (
        current_alpha_s_provider_action_same_lane_split_is_exhaustive
        and current_alpha_s_provider_action_same_lane_remains_nonexport
    )

    s2_noncyclic_continuation_discipline_applies = (
        "strict-core ToE closure using only strict-side sources" in s2_text
        and "new provider class and noncyclic anchor, not a repetition of L5/L12." in s2_text
    )

    same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move = (
        current_alpha_s_provider_action_same_lane_split_is_exhaustive
        and current_alpha_s_provider_action_same_lane_remains_nonexport
        and not current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole
        and s2_noncyclic_continuation_discipline_applies
    )

    next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift = (
        same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move
    )

    add_check(
        "passive_same_lane_stack_already_exported",
        passive_same_lane_stack_already_exported,
        True,
        "F801/F803/F804/F805 already export the passive same-domain alpha_s chain beneath the missing provider action rule.",
    )
    add_check(
        "blocker_stays_unchanged_across_same_lane_split",
        blocker_stays_unchanged_across_same_lane_split,
        True,
        "P802/P803/P804/P805 all keep the same sharp blocker: provider_action_rule_ref.",
    )
    add_check(
        "current_alpha_s_provider_action_same_lane_split_is_exhaustive",
        current_alpha_s_provider_action_same_lane_split_is_exhaustive,
        True,
        "The current same-domain alpha_s lane has already exported the passive layers below provider_action_rule_ref as far as the current honest chain goes.",
    )
    add_check(
        "current_alpha_s_provider_action_same_lane_remains_nonexport",
        current_alpha_s_provider_action_same_lane_remains_nonexport,
        True,
        "Even after exporting the acting input bundle, the missing provider action rule remains nonexport on the same lane.",
    )
    add_check(
        "current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole",
        current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole,
        False,
        "No residual passive same-lane loophole remains below the current provider_action_rule_ref blocker.",
    )
    add_check(
        "s2_noncyclic_continuation_discipline_applies",
        s2_noncyclic_continuation_discipline_applies,
        True,
        "S2 already freezes the continuation discipline: repetition under the same blocker-cut must give way to a genuinely new provider class or noncyclic anchor.",
    )
    add_check(
        "same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move",
        same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move,
        True,
        "Therefore further same-level passive continuation inside the current alpha_s provider-action lane is no longer an admitted primary move.",
    )
    add_check(
        "next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift",
        next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift,
        True,
        "So the next honest move must either export a genuinely new same-domain provider-action source or shift provider class rather than keep splitting the same passive lane.",
    )

    status = (
        "P806_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_SAME_LANE_EXHAUSTION_BOUNDARY_AUDITED"
        if not blocking_checks
        else "P806_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P806",
        "status": status,
        "as_of": AS_OF,
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "theorem_result": {
            "boundary_name": "CurrentStrictAlphaSProviderActionSameLaneExhaustionBoundary_strict_v1",
            "boundary_exported_on_current_repo_state": status
            == "P806_CURRENT_STRICT_ALPHA_S_PROVIDER_ACTION_SAME_LANE_EXHAUSTION_BOUNDARY_AUDITED",
            "passive_same_lane_stack_already_exported": passive_same_lane_stack_already_exported,
            "blocker_stays_unchanged_across_same_lane_split": blocker_stays_unchanged_across_same_lane_split,
            "current_alpha_s_provider_action_same_lane_split_is_exhaustive": current_alpha_s_provider_action_same_lane_split_is_exhaustive,
            "current_alpha_s_provider_action_same_lane_remains_nonexport": current_alpha_s_provider_action_same_lane_remains_nonexport,
            "current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole": current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole,
            "s2_noncyclic_continuation_discipline_applies": s2_noncyclic_continuation_discipline_applies,
            "same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move": same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move,
            "next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift": next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift,
            "no_false_pass": True,
        },
        "inputs": {
            "p802_provider_class_reorganization_probe": rel(IN_P802),
            "p803_relation_bundle_probe": rel(IN_P803),
            "p804_alignment_bundle_probe": rel(IN_P804),
            "p805_acting_input_bundle_probe": rel(IN_P805),
            "f801_provider_skeleton_packet": rel(IN_F801),
            "f802_provider_action_rule_target_packet": rel(IN_F802),
            "f803_relation_bundle_packet": rel(IN_F803),
            "f804_alignment_bundle_packet": rel(IN_F804),
            "f805_acting_input_bundle_packet": rel(IN_F805),
            "s2_priority_packet": rel(IN_S2),
        },
        "checks": checks,
        "blocking_checks": blocking_checks,
        "current_honest_reading": [
            "The current alpha_s same-domain passive lane is now exported all the way down to one explicit acting input bundle.",
            "The sharp blocker has not changed: provider_action_rule_ref remains missing on the same lane.",
            "Under current noncyclic discipline, further same-lane passive decomposition is no longer the admitted primary move.",
        ],
        "recommended_next_move_classes": [
            "export_one_genuinely_new_same_domain_provider_action_source",
            "shift_to_a_different_provider_class_lane",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P806",
        "status": status,
        "as_of": AS_OF,
        "boundary_name": artifact["theorem_result"]["boundary_name"],
        "boundary_exported_on_current_repo_state": artifact["theorem_result"][
            "boundary_exported_on_current_repo_state"
        ],
        "current_alpha_s_provider_action_same_lane_split_is_exhaustive": artifact[
            "theorem_result"
        ]["current_alpha_s_provider_action_same_lane_split_is_exhaustive"],
        "current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole": artifact[
            "theorem_result"
        ]["current_alpha_s_provider_action_same_lane_has_residual_unsplit_loophole"],
        "same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move": artifact[
            "theorem_result"
        ]["same_level_alpha_s_provider_action_lane_continuation_no_longer_admitted_primary_move"],
        "next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift": artifact[
            "theorem_result"
        ]["next_honest_move_requires_genuinely_new_provider_action_source_or_provider_shift"],
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
