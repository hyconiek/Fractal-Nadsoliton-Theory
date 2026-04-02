#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P745 = GENERATED / "p745_current_strict_t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_F160 = GENERATED / "f160_first_actual_nonstrict_declared_scope_selector_closure_packet_summary.json"
IN_N270 = GENERATED / "n270_current_first_nonstrict_declared_scope_selector_closure_theorem_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_F160_MD = ROOT / "F160_FIRST_ACTUAL_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_PACKET.md"

OUT_JSON = GENERATED / "p746_current_strict_t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p746_current_strict_t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_nonexport_audit_probe_summary.json"

CLOSURE_TARGET_NAME = "axiom_augmented_declared_scope_selector_closure_target_v1"
PAIR12_TYPED_CARRIER_NAME = "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
CURRENT_THEOREM_FILE = "N742_CURRENT_STRICT_T200_AXIOM_AUGMENTED_DECLARED_SCOPE_SELECTOR_CLOSURE_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_STRICT_CORE_UPGRADE_BRIDGE_NONEXPORT_BOUNDARY_THEOREM.md"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def scan_positive_bridge_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "f*.py", "n*.py", "t*.py")
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path.name == CURRENT_THEOREM_FILE or path in seen:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if CLOSURE_TARGET_NAME in text and PAIR12_TYPED_CARRIER_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P745, IN_F160, IN_N270, IN_F301, IN_F160_MD]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P746",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p731 = load_json(IN_P731)
    p745 = load_json(IN_P745)
    f160 = load_json(IN_F160)
    n270 = load_json(IN_N270)
    f301 = load_json(IN_F301)
    f160_md_text = IN_F160_MD.read_text(encoding="utf-8")

    pair_index_set = f301.get("pair_index_set") or []
    f301_notes = f301.get("notes") or []
    f301_current_absence = f301.get("current_absence") or []
    n270_result = n270.get("theorem_result") or {}
    positive_bridge_candidates = scan_positive_bridge_candidates()

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

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
            blocking.append(check_id)

    current_actual_nonstrict_declared_scope_selector_closure_exported = (
        f160.get("closure_witness_name") == "C_sel_declared_scope_nonstrict_actual_witness_v1"
        and f160.get("closure_target_name") == CLOSURE_TARGET_NAME
        and bool(f160.get("nonstrict_declared_scope_selector_closure_exported"))
        and bool(n270_result.get("nonstrict_declared_scope_selector_closure_exported"))
        and bool(n270_result.get("discharged"))
    )

    current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier = (
        current_actual_nonstrict_declared_scope_selector_closure_exported
        and "tau_src_candidate_v1" in f160_md_text
        and CLOSURE_TARGET_NAME in f160_md_text
        and f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
    )

    current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged = (
        current_actual_nonstrict_declared_scope_selector_closure_exported
        and f160.get("accepted_scope") == "axiom_augmented_only"
        and n270_result.get("accepted_scope") == "axiom_augmented_only"
        and not bool(f160.get("strict_core_changed"))
        and not bool(n270_result.get("strict_core_changed"))
        and not bool(f160.get("strict_core_selector_closure_claimed"))
        and not bool(n270_result.get("actual_strict_core_selector_closure"))
        and not bool(f160.get("global_selector_closure_claimed"))
        and not bool(n270_result.get("actual_global_selector_closure"))
        and not bool(f160.get("global_qw2191_discharge_claimed"))
        and not bool(n270_result.get("actual_global_qw2191_discharge"))
        and not bool(f160.get("toe_closure_claimed"))
        and not bool(n270_result.get("toe_closure_claimed"))
    )

    surviving_pair12_residual_datum_carrier_remains_selector_neutral = (
        f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and contains_token(f301_notes, "selector-neutral")
        and contains_token(f301_current_absence, "no strict-core selector closure")
        and contains_token(f301_current_absence, "no QW-2191 discharge")
    )

    current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge = (
        current_actual_nonstrict_declared_scope_selector_closure_exported
        and current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and len(positive_bridge_candidates) > 0
        and not current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged
    )

    current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade = (
        current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged
        and current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and not current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge
        and len(positive_bridge_candidates) == 0
    )

    p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge
    )

    add_check(
        "p731_pair12_witness_split_already_opposite_and_unpromoted",
        {
            "split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "pair1_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_sign": p731.get("pair2_w_break_branch_score_sign"),
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
        },
        {
            "split_separated": True,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
            "t185_exported": False,
        },
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed source-side promotion bridge remains unexported.",
    )
    add_check(
        "p745_declared_scope_theorem_target_already_unbridged_to_pair12_typed_carrier",
        {
            "target_exported": bool(p745.get("current_declared_scope_source_topology_selector_theorem_target_exported")),
            "same_tau_src_packet": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "remains_unbridged": bool(
                p745.get("current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier")
            ),
            "t199_exported": bool(p745.get("t199_target_exported_on_current_repo_state")),
        },
        {
            "target_exported": True,
            "same_tau_src_packet": True,
            "remains_unbridged": True,
            "t199_exported": False,
        },
        "P745 already proves that the declared-scope theorem target itself remains unbridged to the surviving pair1/pair2 typed residual-datum carrier.",
    )
    add_check(
        "current_actual_nonstrict_declared_scope_selector_closure_exported",
        {
            "closure_witness_name": f160.get("closure_witness_name"),
            "closure_target_name": f160.get("closure_target_name"),
            "f160_exported": bool(f160.get("nonstrict_declared_scope_selector_closure_exported")),
            "n270_discharged": bool(n270_result.get("discharged")),
            "n270_exported": bool(n270_result.get("nonstrict_declared_scope_selector_closure_exported")),
        },
        {
            "closure_witness_name": "C_sel_declared_scope_nonstrict_actual_witness_v1",
            "closure_target_name": CLOSURE_TARGET_NAME,
            "f160_exported": True,
            "n270_discharged": True,
            "n270_exported": True,
        },
        "F160/N270 already export one real axiom-augmented declared-scope selector-closure witness below strict core.",
    )
    add_check(
        "current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier",
        {
            "f160_mentions_tau_src_candidate_v1": "tau_src_candidate_v1" in f160_md_text,
            "f160_mentions_closure_target": CLOSURE_TARGET_NAME in f160_md_text,
            "carrier_source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
        },
        {
            "f160_mentions_tau_src_candidate_v1": True,
            "f160_mentions_closure_target": True,
            "carrier_source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
        },
        "That axiom-augmented declared-scope closure also lives on the same tau_src_candidate_v1 packet as the surviving F301 pair1/pair2 carrier.",
    )
    add_check(
        "current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged",
        {
            "f160_accepted_scope": f160.get("accepted_scope"),
            "n270_accepted_scope": n270_result.get("accepted_scope"),
            "f160_strict_core_changed": bool(f160.get("strict_core_changed")),
            "n270_strict_core_changed": bool(n270_result.get("strict_core_changed")),
            "f160_strict_core_selector_closure_claimed": bool(f160.get("strict_core_selector_closure_claimed")),
            "n270_actual_strict_core_selector_closure": bool(n270_result.get("actual_strict_core_selector_closure")),
            "f160_global_qw2191_discharge_claimed": bool(f160.get("global_qw2191_discharge_claimed")),
            "n270_actual_global_qw2191_discharge": bool(n270_result.get("actual_global_qw2191_discharge")),
            "f160_toe_closure_claimed": bool(f160.get("toe_closure_claimed")),
            "n270_toe_closure_claimed": bool(n270_result.get("toe_closure_claimed")),
        },
        {
            "f160_accepted_scope": "axiom_augmented_only",
            "n270_accepted_scope": "axiom_augmented_only",
            "f160_strict_core_changed": False,
            "n270_strict_core_changed": False,
            "f160_strict_core_selector_closure_claimed": False,
            "n270_actual_strict_core_selector_closure": False,
            "f160_global_qw2191_discharge_claimed": False,
            "n270_actual_global_qw2191_discharge": False,
            "f160_toe_closure_claimed": False,
            "n270_toe_closure_claimed": False,
        },
        "That closure remains explicitly axiom-augmented only and leaves strict core unchanged; it does not claim strict-core selector closure, global QW-2191 discharge, or ToE closure.",
    )
    add_check(
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral",
        {
            "source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
            "selector_neutral_note_present": contains_token(f301_notes, "selector-neutral"),
            "strict_core_selector_closure_absent": contains_token(
                f301_current_absence, "no strict-core selector closure"
            ),
            "qw2191_discharge_absent": contains_token(f301_current_absence, "no QW-2191 discharge"),
        },
        {
            "source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
            "selector_neutral_note_present": True,
            "strict_core_selector_closure_absent": True,
            "qw2191_discharge_absent": True,
        },
        "The surviving F301 residual-datum carrier remains selector-neutral and still sits below strict-core selector closure or QW-2191 discharge.",
    )
    add_check(
        "positive_packet_theorem_or_spec_bridge_candidates_from_axiom_augmented_declared_scope_selector_closure_to_pair12_typed_carrier",
        positive_bridge_candidates,
        [],
        "No current positive packet/theorem/spec export explicitly bridges the axiom-augmented declared-scope selector-closure target to the surviving pair1/pair2 typed residual-datum carrier.",
    )
    add_check(
        "current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade",
        current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade,
        True,
        "So the current axiom-augmented declared-scope selector closure remains non-strict and does not upgrade to one typed pair1/pair2 residual-datum strict-core bridge.",
    )
    add_check(
        "current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge",
        current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge,
        False,
        "No current export turns the axiom-augmented declared-scope selector closure into one pair1/pair2 typed-carrier strict-core upgrade bridge.",
    )
    add_check(
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure",
        p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not upgrade into one strict-core branch distinction through the current axiom-augmented declared-scope selector closure.",
    )
    add_check(
        "t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_exported",
        False,
        False,
        "So the axiom-augmented declared-scope selector closure to residual-datum pair1/pair2 typed-carrier strict-core upgrade bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_AXIOM_AUGMENTED_DECLARED_SCOPE_SELECTOR_CLOSURE_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_STRICT_CORE_UPGRADE_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P746_REQUIRES_REVIEW_CHANGED_NONSTRICT_DECLARED_SCOPE_SELECTOR_CLOSURE_TO_PAIR12_STRICT_CORE_UPGRADE_STATE"
    )

    artifact = {
        "stage": "P746",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t200_axiom_augmented_declared_scope_selector_closure_to_residual_datum_pair12_typed_carrier_strict_core_upgrade_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P745": str(IN_P745.relative_to(REPO)),
            "F160": str(IN_F160.relative_to(REPO)),
            "N270": str(IN_N270.relative_to(REPO)),
            "F160_markdown": str(IN_F160_MD.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t200_target_name": "AxiomAugmentedDeclaredScopeSelectorClosureToResidualDatumPair12TypedCarrierStrictCoreUpgradeBridge_global_C_v1_strict_v1",
        "t200_target_exported_on_current_repo_state": False,
        "current_actual_nonstrict_declared_scope_selector_closure_exported": (
            current_actual_nonstrict_declared_scope_selector_closure_exported
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier": (
            current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged": (
            current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_axiom_augmented_declared_scope_selector_closure_to_pair12_typed_carrier": (
            positive_bridge_candidates
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade": (
            current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade
        ),
        "current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge": (
            current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge
        ),
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure": (
            p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure
        ),
        "audit_conclusion": {
            "current_repo_already_exports_nonstrict_declared_scope_selector_closure": (
                current_actual_nonstrict_declared_scope_selector_closure_exported
            ),
            "current_repo_already_exports_t200_target": False,
            "next_honest_move": (
                "leave_the_current_axiom_augmented_declared_scope_selector_closure_lane_and_attack_a_genuinely_new_strict_chart_sensitive_source_side_provider_class_or_explicitly_freeze_the_current_nonstrict_lane_as_nonupgrading"
            ),
        },
        "hard_limits": [
            "No T200 discharge claim.",
            "No claim that the current axiom-augmented declared-scope selector closure changes strict core.",
            "No claim that the current non-strict closure selects one raw pair1/pair2 branch in strict core.",
            "No strict-core selector closure claim.",
            "No global QW-2191 discharge claim.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P746",
        "status": status,
        "as_of": AS_OF,
        "t200_target_name": "AxiomAugmentedDeclaredScopeSelectorClosureToResidualDatumPair12TypedCarrierStrictCoreUpgradeBridge_global_C_v1_strict_v1",
        "t200_target_exported_on_current_repo_state": False,
        "current_actual_nonstrict_declared_scope_selector_closure_exported": (
            current_actual_nonstrict_declared_scope_selector_closure_exported
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier": (
            current_actual_nonstrict_declared_scope_selector_closure_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged": (
            current_actual_nonstrict_declared_scope_selector_closure_remains_axiom_augmented_only_and_strict_core_unchanged
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_axiom_augmented_declared_scope_selector_closure_to_pair12_typed_carrier": (
            positive_bridge_candidates
        ),
        "current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade": (
            current_axiom_augmented_declared_scope_selector_closure_remains_nonstrict_not_pair12_typed_strict_core_upgrade
        ),
        "current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge": (
            current_axiom_augmented_declared_scope_selector_closure_has_exported_pair12_typed_carrier_strict_core_upgrade_bridge
        ),
        "p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure": (
            p731_pair12_witness_split_upgrades_to_strict_core_via_current_axiom_augmented_declared_scope_selector_closure
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
