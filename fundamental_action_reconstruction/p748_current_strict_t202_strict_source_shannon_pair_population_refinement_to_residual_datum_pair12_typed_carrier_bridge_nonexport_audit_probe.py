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

IN_ALPHA = GENERATED / "alpha_geo_strict_derived_v1.json"
IN_SIGMA = GENERATED / "sigma_int_strict_derived_v1.json"
IN_P728 = GENERATED / "p728_current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_audit_probe_summary.json"
IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_F319_MD = ROOT / "F319_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_NONEQUALITY_FEEDER_LAW_REFINEMENT_CANDIDATE_PACKET.md"
IN_F320_MD = ROOT / "F320_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_THETA_EXPORT_REFINEMENT_CANDIDATE_PACKET.md"
IN_F321_MD = ROOT / "F321_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_REFINEMENT_CANDIDATE_PACKET.md"
IN_T26_MD = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"

STRICT_SOURCE_SHANNON_PAIR_POPULATION_OBJECT = (
    "BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_population_refinement_candidate_v1"
)
PAIR12_TYPED_CARRIER_NAME = "Omicron_residual_datum_bridge_export_map_object_support_carrier_v1"
PAIR12_ANCHOR_TARGET_NAME = "pair_indexed_population_anchor_target_v1"
CURRENT_THEOREM_FILE = (
    "N744_CURRENT_STRICT_T202_STRICT_SOURCE_SHANNON_PAIR_POPULATION_REFINEMENT_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_BOUNDARY_THEOREM.md"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


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
            if STRICT_SOURCE_SHANNON_PAIR_POPULATION_OBJECT in text and PAIR12_TYPED_CARRIER_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [
        IN_ALPHA,
        IN_SIGMA,
        IN_P728,
        IN_P729,
        IN_P731,
        IN_F301,
        IN_F319_MD,
        IN_F320_MD,
        IN_F321_MD,
        IN_T26_MD,
    ]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P748",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    alpha = load_json(IN_ALPHA)
    sigma = load_json(IN_SIGMA)
    p728 = load_json(IN_P728)
    p729 = load_json(IN_P729)
    p731 = load_json(IN_P731)
    f301 = load_json(IN_F301)
    f319_text = IN_F319_MD.read_text(encoding="utf-8")
    f320_text = IN_F320_MD.read_text(encoding="utf-8")
    f321_text = IN_F321_MD.read_text(encoding="utf-8")
    t26_text = IN_T26_MD.read_text(encoding="utf-8")
    positive_bridge_candidates = scan_positive_bridge_candidates()

    pair_index_set = f301.get("pair_index_set") or []
    f301_notes = f301.get("notes") or []
    f301_current_absence = f301.get("current_absence") or []

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

    current_strict_source_shannon_source_upgrades_exported = (
        alpha.get("object") == "alpha_geo_strict_derived_v1"
        and alpha.get("value") == "4 ln(2)"
        and sigma.get("object") == "sigma_int_strict_derived_v1"
        and sigma.get("value") == -1
    )

    current_t26_pair12_noncyclic_anchor_target_frozen = (
        "strict_source_orientation_seed_target_v1" in t26_text
        and PAIR12_ANCHOR_TARGET_NAME in t26_text
        and "strict_internal_orientation_provider_target_v1" in t26_text
        and "[pair1, pair2]" in t26_text
    )

    current_strict_source_shannon_refinement_stack_exported = (
        "Shannon4ln2_strict_nad12_sigma_residual_nonequality_feeder_law_refinement_candidate_v1" in f319_text
        and "ThetaPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_export_refinement_candidate_v1"
        in f320_text
        and STRICT_SOURCE_SHANNON_PAIR_POPULATION_OBJECT in f321_text
    )

    current_strict_source_shannon_pair_population_refinement_candidate_exported = (
        current_strict_source_shannon_refinement_stack_exported
        and "actual packaged strict-source Shannon-weighted pair-population refinement candidate" in f321_text
        and "strict_weight: alpha_geo_strict_derived_v1" in f321_text
        and "pair_indexed_slot_status = present_via_R1" in f321_text
    )

    current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed = (
        current_strict_source_shannon_pair_population_refinement_candidate_exported
        and "still below actual pair population" in f321_text
        and "still below actual theta export" in f321_text
        and "actual pair population" in f321_text
        and "actual theta" in f321_text
        and "tau_src_candidate_v1" not in f321_text
        and "pair1" not in f321_text
        and "pair2" not in f321_text
        and PAIR12_TYPED_CARRIER_NAME not in f321_text
    )

    current_surviving_pair12_residual_datum_carrier_exported = (
        f301.get("object") == PAIR12_TYPED_CARRIER_NAME
        and f301.get("source_domain") == "tau_src_candidate_v1"
        and sorted(pair_index_set) == ["pair1", "pair2"]
        and bool(p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane"))
        and sorted(p728.get("residual_datum_source_side_supported_positive_charts") or []) == ["pair1", "pair2"]
    )

    current_surviving_pair12_residual_datum_carrier_remains_selector_neutral = (
        current_surviving_pair12_residual_datum_carrier_exported
        and any("selector-neutral" in str(entry) for entry in f301_notes)
        and any("no strict-core selector closure" in str(entry) for entry in f301_current_absence)
        and any("no QW-2191 discharge" in str(entry) for entry in f301_current_absence)
    )

    current_pair12_split_remains_unresolved_on_current_frontier = (
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions"))
        and p729.get("pair1_orbit_branch_kind") == "delta_k_positive_index_branch"
        and p729.get("pair2_orbit_branch_kind") == "delta_minus_k_negative_index_branch"
        and bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and p731.get("pair1_w_break_branch_score_sign") == "negative"
        and p731.get("pair2_w_break_branch_score_sign") == "positive"
    )

    current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge = (
        current_strict_source_shannon_source_upgrades_exported
        and current_t26_pair12_noncyclic_anchor_target_frozen
        and current_strict_source_shannon_pair_population_refinement_candidate_exported
        and current_surviving_pair12_residual_datum_carrier_exported
        and len(positive_bridge_candidates) > 0
        and not current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed
    )

    current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier = (
        current_strict_source_shannon_source_upgrades_exported
        and current_t26_pair12_noncyclic_anchor_target_frozen
        and current_strict_source_shannon_pair_population_refinement_candidate_exported
        and current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed
        and current_surviving_pair12_residual_datum_carrier_remains_selector_neutral
        and current_pair12_split_remains_unresolved_on_current_frontier
        and len(positive_bridge_candidates) == 0
        and not current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge
    )

    p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge
    )

    add_check(
        "current_strict_source_shannon_source_upgrades_exported",
        {
            "alpha_object": alpha.get("object"),
            "alpha_value": alpha.get("value"),
            "sigma_object": sigma.get("object"),
            "sigma_value": sigma.get("value"),
        },
        {
            "alpha_object": "alpha_geo_strict_derived_v1",
            "alpha_value": "4 ln(2)",
            "sigma_object": "sigma_int_strict_derived_v1",
            "sigma_value": -1,
        },
        "The current repo already exports the strict-source Shannon upgrade data alpha_geo_strict_derived_v1 = 4 ln(2) and sigma_int_strict_derived_v1 = -1.",
    )
    add_check(
        "t26_pair12_noncyclic_anchor_target_frozen",
        {
            "strict_source_orientation_seed_target_present": "strict_source_orientation_seed_target_v1" in t26_text,
            "pair_indexed_population_anchor_target_present": PAIR12_ANCHOR_TARGET_NAME in t26_text,
            "strict_internal_orientation_provider_target_present": (
                "strict_internal_orientation_provider_target_v1" in t26_text
            ),
            "minimal_designated_pair_family_pair12_present": "[pair1, pair2]" in t26_text,
        },
        {
            "strict_source_orientation_seed_target_present": True,
            "pair_indexed_population_anchor_target_present": True,
            "strict_internal_orientation_provider_target_present": True,
            "minimal_designated_pair_family_pair12_present": True,
        },
        "T26 already freezes the genuinely new noncyclic anchor target class, with minimal designated pair family [pair1, pair2], as a future-only strict-side route.",
    )
    add_check(
        "current_strict_source_shannon_refinement_stack_exported",
        {
            "feeder_refinement_packet_present": (
                "Shannon4ln2_strict_nad12_sigma_residual_nonequality_feeder_law_refinement_candidate_v1"
                in f319_text
            ),
            "theta_refinement_packet_present": (
                "ThetaPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_export_refinement_candidate_v1"
                in f320_text
            ),
            "pair_population_refinement_packet_present": STRICT_SOURCE_SHANNON_PAIR_POPULATION_OBJECT in f321_text,
        },
        {
            "feeder_refinement_packet_present": True,
            "theta_refinement_packet_present": True,
            "pair_population_refinement_packet_present": True,
        },
        "F319/F320/F321 already export the strict-source Shannon refinement stack up through pair-population refinement-candidate language.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_candidate_exported",
        {
            "pair_population_object_present": STRICT_SOURCE_SHANNON_PAIR_POPULATION_OBJECT in f321_text,
            "candidate_packet_phrase_present": (
                "actual packaged strict-source Shannon-weighted pair-population refinement candidate" in f321_text
            ),
            "strict_weight_field_present": "strict_weight: alpha_geo_strict_derived_v1" in f321_text,
            "pair_indexed_slot_status_present": "pair_indexed_slot_status = present_via_R1" in f321_text,
        },
        {
            "pair_population_object_present": True,
            "candidate_packet_phrase_present": True,
            "strict_weight_field_present": True,
            "pair_indexed_slot_status_present": True,
        },
        "F321 already packages one real strict-source Shannon pair-population refinement candidate with strict weight and pair-indexed slot language.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed",
        {
            "still_below_actual_pair_population": "still below actual pair population" in f321_text,
            "still_below_actual_theta_export": "still below actual theta export" in f321_text,
            "tau_src_candidate_v1_mentioned": "tau_src_candidate_v1" in f321_text,
            "pair1_or_pair2_mentioned": ("pair1" in f321_text or "pair2" in f321_text),
            "pair12_typed_carrier_mentioned": PAIR12_TYPED_CARRIER_NAME in f321_text,
        },
        {
            "still_below_actual_pair_population": True,
            "still_below_actual_theta_export": True,
            "tau_src_candidate_v1_mentioned": False,
            "pair1_or_pair2_mentioned": False,
            "pair12_typed_carrier_mentioned": False,
        },
        "That strict-source Shannon pair-population route still remains candidate-only, below actual pair population/theta export, and is not yet typed on tau_src_candidate_v1 or the surviving pair1/pair2 carrier.",
    )
    add_check(
        "current_surviving_pair12_residual_datum_carrier_exported",
        {
            "carrier_object": f301.get("object"),
            "source_domain": f301.get("source_domain"),
            "pair_index_set": pair_index_set,
            "supported_positive_charts": p728.get("residual_datum_source_side_supported_positive_charts"),
            "corridor_reduced_to_pair12": bool(
                p728.get("current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane")
            ),
        },
        {
            "carrier_object": PAIR12_TYPED_CARRIER_NAME,
            "source_domain": "tau_src_candidate_v1",
            "pair_index_set": ["pair1", "pair2"],
            "supported_positive_charts": ["pair1", "pair2"],
            "corridor_reduced_to_pair12": True,
        },
        "The surviving source-side frontier is already reduced to the selector-neutral pair1/pair2 residual-datum carrier on tau_src_candidate_v1 (P728/F301).",
    )
    add_check(
        "current_pair12_split_remains_unresolved_on_current_frontier",
        {
            "orbit_split_localized": bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
            "pair1_orbit_branch_kind": p729.get("pair1_orbit_branch_kind"),
            "pair2_orbit_branch_kind": p729.get("pair2_orbit_branch_kind"),
            "w_break_split_separated": bool(
                p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")
            ),
            "pair1_w_break_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_w_break_sign": p731.get("pair2_w_break_branch_score_sign"),
        },
        {
            "orbit_split_localized": True,
            "pair1_orbit_branch_kind": "delta_k_positive_index_branch",
            "pair2_orbit_branch_kind": "delta_minus_k_negative_index_branch",
            "w_break_split_separated": True,
            "pair1_w_break_sign": "negative",
            "pair2_w_break_sign": "positive",
        },
        "P729/P731 already freeze the surviving pair1/pair2 frontier as an unresolved opposite orbit-direction split with antisymmetric witness scores.",
    )
    add_check(
        "positive_packet_theorem_or_spec_bridge_candidates_from_strict_source_shannon_pair_population_refinement_to_pair12_typed_carrier",
        positive_bridge_candidates,
        [],
        "No current packet, theorem, or spec simultaneously exports the strict-source Shannon pair-population refinement object together with the surviving pair1/pair2 typed carrier as a positive bridge lane.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier",
        current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier,
        True,
        "So the current strict-source Shannon pair-population refinement route still remains unbridged to the surviving pair1/pair2 residual-datum carrier on current exports.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge",
        current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge,
        False,
        "No current export carries the strict-source Shannon pair-population refinement route into the surviving pair1/pair2 typed carrier as an explicit source-side bridge.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge",
        p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through any exported strict-source Shannon pair-population refinement to pair1/pair2 typed-carrier bridge.",
    )
    add_check(
        "t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_exported",
        False,
        False,
        "So the strict-source Shannon pair-population refinement to residual-datum pair1/pair2 typed-carrier bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_STRICT_SOURCE_SHANNON_PAIR_POPULATION_REFINEMENT_TO_RESIDUAL_DATUM_PAIR12_TYPED_CARRIER_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P748_REQUIRES_REVIEW_CHANGED_STRICT_SOURCE_SHANNON_PAIR12_BRIDGE_STATE"
    )

    artifact = {
        "stage": "P748",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_boundary_only",
        "inputs": {
            "alpha_geo_strict_derived": str(IN_ALPHA.relative_to(REPO)),
            "sigma_int_strict_derived": str(IN_SIGMA.relative_to(REPO)),
            "P728": str(IN_P728.relative_to(REPO)),
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "F301_carrier": str(IN_F301.relative_to(REPO)),
            "F319": str(IN_F319_MD.relative_to(REPO)),
            "F320": str(IN_F320_MD.relative_to(REPO)),
            "F321": str(IN_F321_MD.relative_to(REPO)),
            "T26": str(IN_T26_MD.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t202_target_name": "StrictSourceShannonPairPopulationRefinementToResidualDatumPair12TypedCarrierBridge_strict_v1",
        "t202_target_exported_on_current_repo_state": False,
        "current_strict_source_shannon_source_upgrades_exported": (
            current_strict_source_shannon_source_upgrades_exported
        ),
        "current_t26_pair12_noncyclic_anchor_target_frozen": current_t26_pair12_noncyclic_anchor_target_frozen,
        "current_strict_source_shannon_refinement_stack_exported": (
            current_strict_source_shannon_refinement_stack_exported
        ),
        "current_strict_source_shannon_pair_population_refinement_candidate_exported": (
            current_strict_source_shannon_pair_population_refinement_candidate_exported
        ),
        "current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed": (
            current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed
        ),
        "current_surviving_pair12_residual_datum_carrier_exported": (
            current_surviving_pair12_residual_datum_carrier_exported
        ),
        "current_surviving_pair12_residual_datum_carrier_remains_selector_neutral": (
            current_surviving_pair12_residual_datum_carrier_remains_selector_neutral
        ),
        "current_pair12_split_remains_unresolved_on_current_frontier": (
            current_pair12_split_remains_unresolved_on_current_frontier
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_strict_source_shannon_pair_population_refinement_to_pair12_typed_carrier": (
            positive_bridge_candidates
        ),
        "current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier": (
            current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier
        ),
        "current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge": (
            current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge
        ),
        "p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge": (
            p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge
        ),
        "frontier": {
            "current_repo_already_exports_strict_source_shannon_pair_population_refinement_candidate": (
                current_strict_source_shannon_pair_population_refinement_candidate_exported
            ),
            "current_repo_already_exports_surviving_pair12_typed_carrier": (
                current_surviving_pair12_residual_datum_carrier_exported
            ),
            "current_repo_already_exports_strict_source_shannon_pair_population_refinement_to_pair12_typed_carrier_bridge": (
                current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge
            ),
            "next_honest_move": (
                "attempt_a_genuinely_new_strict_source_bridge_from_the_shannon_weighted_pair_population_refinement_route_into_the_surviving_pair1_pair2_source_side_carrier_without_silently_identifying_candidate_pair_population_syntax_with_tau_src_or_with_the_current_residual_datum_carrier"
            ),
        },
        "hard_limits": [
            "No claim that the strict-source Shannon pair-population refinement candidate is already an actual pair population.",
            "No claim that the strict-source Shannon route is already typed on tau_src_candidate_v1.",
            "No claim that pair-indexed slot syntax already identifies one surviving pair1/pair2 branch.",
            "No claim that the current surviving residual-datum carrier is already identified with the strict-source Shannon pair-population route.",
            "No kernel-alone/global QW-2191 discharge.",
            "No strict-core selector closure.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t202_target_name": artifact["t202_target_name"],
        "t202_target_exported_on_current_repo_state": artifact["t202_target_exported_on_current_repo_state"],
        "current_strict_source_shannon_source_upgrades_exported": (
            artifact["current_strict_source_shannon_source_upgrades_exported"]
        ),
        "current_t26_pair12_noncyclic_anchor_target_frozen": artifact[
            "current_t26_pair12_noncyclic_anchor_target_frozen"
        ],
        "current_strict_source_shannon_pair_population_refinement_candidate_exported": (
            artifact["current_strict_source_shannon_pair_population_refinement_candidate_exported"]
        ),
        "current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed": (
            artifact["current_strict_source_shannon_pair_population_refinement_candidate_remains_candidate_only_not_pair12_typed"]
        ),
        "current_surviving_pair12_residual_datum_carrier_exported": (
            artifact["current_surviving_pair12_residual_datum_carrier_exported"]
        ),
        "current_pair12_split_remains_unresolved_on_current_frontier": (
            artifact["current_pair12_split_remains_unresolved_on_current_frontier"]
        ),
        "current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier": (
            artifact["current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier"]
        ),
        "current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge": (
            artifact["current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge"]
        ),
        "p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge": (
            artifact["p731_pair12_witness_split_descends_to_current_strict_source_shannon_pair_population_bridge"]
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
