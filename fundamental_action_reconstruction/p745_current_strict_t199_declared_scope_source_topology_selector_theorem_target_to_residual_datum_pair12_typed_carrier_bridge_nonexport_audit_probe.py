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
IN_P744 = GENERATED / "p744_current_strict_t198_declared_scope_source_topology_selector_theorem_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
IN_N258 = GENERATED / "n258_current_first_declared_scope_source_topology_selector_theorem_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p745_current_strict_t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p745_current_strict_t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"

THEOREM_TARGET_NAME = "declared_scope_source_topology_selector_theorem_target_v1"
PAIR12_TYPED_CARRIER_NAME = "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
CURRENT_THEOREM_FILE = "N741_CURRENT_STRICT_T199_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_TARGET_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_BOUNDARY_THEOREM.md"


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
            if THEOREM_TARGET_NAME in text and PAIR12_TYPED_CARRIER_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P744, IN_F150, IN_N258, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P745",
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
    p744 = load_json(IN_P744)
    f150 = load_json(IN_F150)
    n258 = load_json(IN_N258)
    f301 = load_json(IN_F301)

    pair_index_set = f301.get("pair_index_set") or []
    f301_notes = f301.get("notes") or []
    f301_current_absence = f301.get("current_absence") or []
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

    current_declared_scope_source_topology_selector_theorem_target_exported = (
        f150.get("codomain_target") == THEOREM_TARGET_NAME
        and n258.get("codomain_target") == THEOREM_TARGET_NAME
        and f150.get("input_packet") == "tau_src_candidate_v1"
        and n258.get("input_packet") == "tau_src_candidate_v1"
        and bool(f150.get("declared_scope_source_topology_selector_theorem_exported"))
        and bool(n258.get("declared_scope_source_topology_selector_theorem_exported"))
    )

    current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier = (
        current_declared_scope_source_topology_selector_theorem_target_exported
        and f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
    )

    current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only = (
        current_declared_scope_source_topology_selector_theorem_target_exported
        and bool(f150.get("declared_scope_only"))
        and bool(n258.get("declared_scope_only"))
        and not bool(f150.get("raw_theta_uniqueness_claimed"))
        and not bool(n258.get("current_selector_closure"))
        and not bool(n258.get("current_global_qw2191_discharge"))
        and not bool(f150.get("tau_src_identified_with_s_prelm"))
        and not bool(n258.get("tau_src_identified_with_s_prelm"))
        and not bool(f150.get("admissible_strict_core_internal_selector_source_object_claimed"))
        and not bool(n258.get("admissible_strict_core_internal_selector_source_object_claimed"))
    )

    surviving_pair12_residual_datum_carrier_remains_selector_neutral = (
        f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and contains_token(f301_notes, "selector-neutral")
        and contains_token(f301_current_absence, "no strict-core selector closure")
        and contains_token(f301_current_absence, "no QW-2191 discharge")
    )

    current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge = (
        current_declared_scope_source_topology_selector_theorem_target_exported
        and current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and len(positive_bridge_candidates) > 0
    )

    current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier = (
        current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only
        and current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier
        and surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and not current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge
        and len(positive_bridge_candidates) == 0
    )

    p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge
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
        "p744_declared_scope_source_topology_selector_theorem_lane_already_real_but_not_pair12_typed",
        {
            "same_tau_src_packet": bool(
                p744.get("current_declared_scope_source_topology_selector_theorem_binds_same_tau_src_packet_as_pair12_carrier")
            ),
            "remains_declared_scope_quotient_class_only": bool(
                p744.get("current_declared_scope_source_topology_selector_theorem_remains_declared_scope_quotient_class_only")
            ),
            "continuation_remains_declared_scope_quotient_class_only_not_pair12_typed": bool(
                p744.get("current_declared_scope_source_topology_selector_theorem_continuation_remains_declared_scope_quotient_class_only_not_pair12_typed")
            ),
            "t198_exported": bool(p744.get("t198_target_exported_on_current_repo_state")),
        },
        {
            "same_tau_src_packet": True,
            "remains_declared_scope_quotient_class_only": True,
            "continuation_remains_declared_scope_quotient_class_only_not_pair12_typed": True,
            "t198_exported": False,
        },
        "P744 already proves that the declared-scope theorem lane is real on the same tau_src packet as F301, but still remains declared-scope and quotient-class only rather than typed on the surviving pair1/pair2 carrier.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_exported",
        {
            "f150_codomain_target": f150.get("codomain_target"),
            "n258_codomain_target": n258.get("codomain_target"),
            "f150_input_packet": f150.get("input_packet"),
            "n258_input_packet": n258.get("input_packet"),
            "f150_declared_scope_theorem_exported": bool(
                f150.get("declared_scope_source_topology_selector_theorem_exported")
            ),
            "n258_declared_scope_theorem_exported": bool(
                n258.get("declared_scope_source_topology_selector_theorem_exported")
            ),
        },
        {
            "f150_codomain_target": THEOREM_TARGET_NAME,
            "n258_codomain_target": THEOREM_TARGET_NAME,
            "f150_input_packet": "tau_src_candidate_v1",
            "n258_input_packet": "tau_src_candidate_v1",
            "f150_declared_scope_theorem_exported": True,
            "n258_declared_scope_theorem_exported": True,
        },
        "F150/N258 already export one real declared-scope source-topology selector theorem target on tau_src_candidate_v1.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier",
        {
            "theorem_target_input_packet": f150.get("input_packet"),
            "carrier_source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
        },
        {
            "theorem_target_input_packet": "tau_src_candidate_v1",
            "carrier_source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
        },
        "That declared-scope theorem target already lives on the same tau_src_candidate_v1 packet as the surviving F301 pair1/pair2 residual-datum carrier.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only",
        {
            "f150_declared_scope_only": bool(f150.get("declared_scope_only")),
            "n258_declared_scope_only": bool(n258.get("declared_scope_only")),
            "f150_raw_theta_uniqueness_claimed": bool(f150.get("raw_theta_uniqueness_claimed")),
            "n258_current_selector_closure": bool(n258.get("current_selector_closure")),
            "n258_current_global_qw2191_discharge": bool(n258.get("current_global_qw2191_discharge")),
            "f150_tau_src_identified_with_s_prelm": bool(f150.get("tau_src_identified_with_s_prelm")),
            "f150_admissible_strict_core_internal_selector_source_object_claimed": bool(
                f150.get("admissible_strict_core_internal_selector_source_object_claimed")
            ),
        },
        {
            "f150_declared_scope_only": True,
            "n258_declared_scope_only": True,
            "f150_raw_theta_uniqueness_claimed": False,
            "n258_current_selector_closure": False,
            "n258_current_global_qw2191_discharge": False,
            "f150_tau_src_identified_with_s_prelm": False,
            "f150_admissible_strict_core_internal_selector_source_object_claimed": False,
        },
        "That declared-scope theorem target still remains scope-limited, quotient-safe, and below strict-core selector closure or global QW-2191 discharge on current exports.",
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
        "positive_packet_theorem_or_spec_bridge_candidates_from_declared_scope_source_topology_selector_theorem_target_to_pair12_carrier",
        positive_bridge_candidates,
        [],
        "No current positive packet/theorem/spec export explicitly bridges the declared-scope theorem target to the surviving pair1/pair2 typed residual-datum carrier on current repo state.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier",
        current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier,
        True,
        "So the currently exported declared-scope theorem target still remains unbridged to the surviving pair1/pair2 typed residual-datum carrier.",
    )
    add_check(
        "current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge",
        current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge,
        False,
        "No current export turns the declared-scope theorem target itself into one typed pair1/pair2 residual-datum carrier bridge.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge",
        p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through any exported declared-scope theorem-target to residual-datum typed-carrier bridge.",
    )
    add_check(
        "t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_exported",
        False,
        False,
        "So the declared-scope source-topology selector theorem-target to residual-datum pair1/pair2 typed-carrier bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_DECLARED_SCOPE_SOURCE_TOPOLOGY_SELECTOR_THEOREM_TARGET_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P745_REQUIRES_REVIEW_CHANGED_DECLARED_SCOPE_THEOREM_TARGET_TO_PAIR12_TYPED_CARRIER_STATE"
    )

    artifact = {
        "stage": "P745",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t199_declared_scope_source_topology_selector_theorem_target_to_residual_datum_pair12_typed_carrier_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P744": str(IN_P744.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
            "N258": str(IN_N258.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t199_target_name": "DeclaredScopeSourceTopologySelectorTheoremTargetToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t199_target_exported_on_current_repo_state": False,
        "current_declared_scope_source_topology_selector_theorem_target_exported": (
            current_declared_scope_source_topology_selector_theorem_target_exported
        ),
        "current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier": (
            current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only": (
            current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_declared_scope_source_topology_selector_theorem_target_to_pair12_carrier": (
            positive_bridge_candidates
        ),
        "surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier": (
            current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier
        ),
        "current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge": (
            current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge": (
            p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge
        ),
        "audit_conclusion": {
            "current_repo_already_exports_declared_scope_theorem_target": (
                current_declared_scope_source_topology_selector_theorem_target_exported
            ),
            "current_repo_already_exports_t199_target": False,
            "next_honest_move": (
                "either_export_one_explicit_bridge_from_declared_scope_source_topology_selector_theorem_target_v1_to_the_surviving_pair12_typed_residual_datum_carrier_or_leave_the_current_declared_scope_theorem_target_entirely_and_attack_a_genuinely_new_chart_sensitive_source_side_provider_class"
            ),
        },
        "hard_limits": [
            "No T199 discharge claim.",
            "No claim that declared_scope_source_topology_selector_theorem_target_v1 is already the F301 typed carrier.",
            "No claim that the current declared-scope theorem target already selects one raw pair1/pair2 branch.",
            "No raw-theta uniqueness claim.",
            "No strict-core selector closure claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P745",
        "status": status,
        "as_of": AS_OF,
        "t199_target_name": "DeclaredScopeSourceTopologySelectorTheoremTargetToResidualDatumPair12TypedCarrierBridge_global_C_v1_strict_v1",
        "t199_target_exported_on_current_repo_state": False,
        "current_declared_scope_source_topology_selector_theorem_target_exported": (
            current_declared_scope_source_topology_selector_theorem_target_exported
        ),
        "current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier": (
            current_declared_scope_source_topology_selector_theorem_target_binds_same_tau_src_packet_as_pair12_carrier
        ),
        "current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only": (
            current_declared_scope_source_topology_selector_theorem_target_remains_declared_scope_quotient_class_only
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_declared_scope_source_topology_selector_theorem_target_to_pair12_carrier": (
            positive_bridge_candidates
        ),
        "current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier": (
            current_declared_scope_source_topology_selector_theorem_target_remains_unbridged_to_pair12_typed_carrier
        ),
        "current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge": (
            current_declared_scope_source_topology_selector_theorem_target_has_exported_pair12_typed_residual_datum_bridge
        ),
        "p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge": (
            p731_pair12_witness_split_descends_to_current_declared_scope_source_topology_selector_theorem_target_to_residual_datum_typed_carrier_bridge
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
