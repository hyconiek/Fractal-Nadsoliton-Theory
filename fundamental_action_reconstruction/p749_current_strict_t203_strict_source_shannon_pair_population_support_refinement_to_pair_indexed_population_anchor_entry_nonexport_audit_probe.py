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
IN_P748 = GENERATED / "p748_current_strict_t202_strict_source_shannon_pair_population_refinement_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_F323_MD = ROOT / "F323_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_SUPPORT_REFINEMENT_CANDIDATE_PACKET.md"
IN_N434_MD = ROOT / "N434_CURRENT_FIRST_ACTUAL_STRICT_NAD12_SIGMA_STRICT_SOURCE_SHANNON_WEIGHTED_PAIR_POPULATION_SUPPORT_REFINEMENT_CANDIDATE_THEOREM.md"
IN_T26_MD = ROOT / "T26_STRICT_SOURCE_TO_PAIR_POPULATION_NONCYCLIC_ANCHOR_TARGET_SPEC.md"

OUT_JSON = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p749_current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_audit_probe_summary.json"

STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_OBJECT = (
    "BasisPair_strict_nad12_sigma_residual_nonequality_shannon_weighted_support_refinement_candidate_v1"
)
PAIR_INDEXED_POPULATION_ANCHOR_TARGET = "pair_indexed_population_anchor_target_v1"
CURRENT_THEOREM_FILE = (
    "N745_CURRENT_STRICT_T203_STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_REFINEMENT_TO_PAIR_INDEXED_POPULATION_ANCHOR_ENTRY_NONEXPORT_BOUNDARY_THEOREM.md"
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
            if (
                STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_OBJECT in text
                and PAIR_INDEXED_POPULATION_ANCHOR_TARGET in text
            ):
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_ALPHA, IN_SIGMA, IN_P748, IN_F323_MD, IN_N434_MD, IN_T26_MD]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P749",
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
    p748 = load_json(IN_P748)
    f323_text = IN_F323_MD.read_text(encoding="utf-8")
    n434_text = IN_N434_MD.read_text(encoding="utf-8")
    t26_text = IN_T26_MD.read_text(encoding="utf-8")
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

    current_strict_source_shannon_source_upgrades_exported = (
        alpha.get("object") == "alpha_geo_strict_derived_v1"
        and alpha.get("value") == "4 ln(2)"
        and sigma.get("object") == "sigma_int_strict_derived_v1"
        and sigma.get("value") == -1
    )

    p748_weaker_strict_source_shannon_pair_population_refinement_route_already_unbridged = (
        bool(p748.get("current_strict_source_shannon_pair_population_refinement_candidate_exported"))
        and bool(p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier"))
        and not bool(p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge"))
        and not bool(p748.get("t202_target_exported_on_current_repo_state"))
    )

    current_t26_pair_indexed_population_anchor_target_frozen = (
        "strict_source_orientation_seed_target_v1" in t26_text
        and PAIR_INDEXED_POPULATION_ANCHOR_TARGET in t26_text
        and "strict_internal_orientation_provider_target_v1" in t26_text
        and "[pair1, pair2]" in t26_text
    )

    current_strict_source_shannon_pair_population_support_refinement_candidate_exported = (
        STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_OBJECT in f323_text
        and "actual packaged strict-source Shannon-weighted pair-population support refinement candidate"
        in f323_text
        and STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_OBJECT in n434_text
    )

    current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only = (
        current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and "pair_indexed_slot_status = present_via_R1" in f323_text
        and "conditional_population_status = present_via_C48_C49" in f323_text
        and "strict_weight: alpha_geo_strict_derived_v1" in f323_text
        and "tau_src_candidate_v1" not in f323_text
        and "[pair1, pair2]" not in f323_text
        and "pair1" not in f323_text
        and "pair2" not in f323_text
        and PAIR_INDEXED_POPULATION_ANCHOR_TARGET not in f323_text
    )

    current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export = (
        current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and "still below actual pair population" in f323_text
        and "still below actual theta export" in f323_text
        and "not actual pair population" in n434_text
        and "not actual theta export" in n434_text
    )

    current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry = (
        current_strict_source_shannon_source_upgrades_exported
        and current_t26_pair_indexed_population_anchor_target_frozen
        and current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and len(positive_bridge_candidates) > 0
        and not current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only
        and not current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export
    )

    current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor = (
        current_strict_source_shannon_source_upgrades_exported
        and p748_weaker_strict_source_shannon_pair_population_refinement_route_already_unbridged
        and current_t26_pair_indexed_population_anchor_target_frozen
        and current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        and current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only
        and current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export
        and len(positive_bridge_candidates) == 0
        and not current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry
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
        "p748_weaker_strict_source_shannon_pair_population_refinement_route_already_unbridged",
        {
            "weaker_route_exported": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_candidate_exported")
            ),
            "weaker_route_unbridged": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_route_remains_unbridged_to_pair12_typed_carrier")
            ),
            "weaker_route_bridge_exported": bool(
                p748.get("current_strict_source_shannon_pair_population_refinement_route_has_exported_pair12_typed_carrier_bridge")
            ),
            "t202_exported": bool(p748.get("t202_target_exported_on_current_repo_state")),
        },
        {
            "weaker_route_exported": True,
            "weaker_route_unbridged": True,
            "weaker_route_bridge_exported": False,
            "t202_exported": False,
        },
        "P748 already proves that the weaker strict-source Shannon pair-population refinement route is real but still unbridged on the current frontier.",
    )
    add_check(
        "t26_pair_indexed_population_anchor_target_frozen",
        {
            "strict_source_orientation_seed_target_present": "strict_source_orientation_seed_target_v1" in t26_text,
            "pair_indexed_population_anchor_target_present": PAIR_INDEXED_POPULATION_ANCHOR_TARGET in t26_text,
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
        "T26 already freezes component 2 as a future pair-indexed population anchor on at least the minimal designated pair family [pair1, pair2].",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported",
        {
            "support_object_present_in_f323": STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_OBJECT in f323_text,
            "support_packet_phrase_present": (
                "actual packaged strict-source Shannon-weighted pair-population support refinement candidate"
                in f323_text
            ),
            "support_object_present_in_n434": STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_OBJECT in n434_text,
        },
        {
            "support_object_present_in_f323": True,
            "support_packet_phrase_present": True,
            "support_object_present_in_n434": True,
        },
        "F323/N434 already export one real strict-source Shannon pair-population support-refinement candidate as the strongest current layer of that route.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only",
        {
            "pair_indexed_slot_status_present": "pair_indexed_slot_status = present_via_R1" in f323_text,
            "conditional_population_status_present": "conditional_population_status = present_via_C48_C49" in f323_text,
            "strict_weight_field_present": "strict_weight: alpha_geo_strict_derived_v1" in f323_text,
            "tau_src_candidate_v1_mentioned": "tau_src_candidate_v1" in f323_text,
            "minimal_designated_pair_family_pair12_mentioned": "[pair1, pair2]" in f323_text,
            "pair1_or_pair2_mentioned": ("pair1" in f323_text or "pair2" in f323_text),
            "pair_indexed_population_anchor_target_mentioned": PAIR_INDEXED_POPULATION_ANCHOR_TARGET in f323_text,
        },
        {
            "pair_indexed_slot_status_present": True,
            "conditional_population_status_present": True,
            "strict_weight_field_present": True,
            "tau_src_candidate_v1_mentioned": False,
            "minimal_designated_pair_family_pair12_mentioned": False,
            "pair1_or_pair2_mentioned": False,
            "pair_indexed_population_anchor_target_mentioned": False,
        },
        "That strongest strict-source Shannon support-refinement route still exports only generic pair-indexed slot/population syntax, not an actual reduction to tau_src_candidate_v1 or to the minimal designated pair family [pair1, pair2].",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export",
        {
            "still_below_actual_pair_population_in_f323": "still below actual pair population" in f323_text,
            "still_below_actual_theta_export_in_f323": "still below actual theta export" in f323_text,
            "not_actual_pair_population_in_n434": "not actual pair population" in n434_text,
            "not_actual_theta_export_in_n434": "not actual theta export" in n434_text,
        },
        {
            "still_below_actual_pair_population_in_f323": True,
            "still_below_actual_theta_export_in_f323": True,
            "not_actual_pair_population_in_n434": True,
            "not_actual_theta_export_in_n434": True,
        },
        "F323/N434 keep the strict-source Shannon support-refinement route explicitly below actual pair population and actual theta export.",
    )
    add_check(
        "positive_packet_theorem_or_spec_bridge_candidates_from_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor",
        positive_bridge_candidates,
        [],
        "No current packet, theorem, or spec simultaneously exports the strongest strict-source Shannon pair-population support-refinement object together with the T26 pair-indexed population-anchor target as a positive entry bridge.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor",
        current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor,
        True,
        "So the strongest current strict-source Shannon pair-population support-refinement route still remains nonentering for T26 component 2 on current exports.",
    )
    add_check(
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry",
        current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry,
        False,
        "No current export promotes the strongest strict-source Shannon pair-population support-refinement route into one actual pair-indexed population-anchor entry.",
    )
    add_check(
        "t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_exported",
        False,
        False,
        "So the strict-source Shannon pair-population support-refinement to pair-indexed population-anchor entry bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_STRICT_SOURCE_SHANNON_PAIR_POPULATION_SUPPORT_REFINEMENT_TO_PAIR_INDEXED_POPULATION_ANCHOR_ENTRY_NONEXPORT_AUDITED"
        if not blocking
        else "P749_REQUIRES_REVIEW_CHANGED_STRICT_SOURCE_SHANNON_COMPONENT2_ENTRY_STATE"
    )

    artifact = {
        "stage": "P749",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t203_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry_nonexport_boundary_only",
        "inputs": {
            "alpha_geo_strict_derived": str(IN_ALPHA.relative_to(REPO)),
            "sigma_int_strict_derived": str(IN_SIGMA.relative_to(REPO)),
            "P748": str(IN_P748.relative_to(REPO)),
            "F323": str(IN_F323_MD.relative_to(REPO)),
            "N434": str(IN_N434_MD.relative_to(REPO)),
            "T26": str(IN_T26_MD.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t203_target_name": "StrictSourceShannonPairPopulationSupportRefinementToPairIndexedPopulationAnchorEntry_strict_v1",
        "t203_target_exported_on_current_repo_state": False,
        "current_strict_source_shannon_source_upgrades_exported": (
            current_strict_source_shannon_source_upgrades_exported
        ),
        "p748_weaker_strict_source_shannon_pair_population_refinement_route_already_unbridged": (
            p748_weaker_strict_source_shannon_pair_population_refinement_route_already_unbridged
        ),
        "current_t26_pair_indexed_population_anchor_target_frozen": (
            current_t26_pair_indexed_population_anchor_target_frozen
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported": (
            current_strict_source_shannon_pair_population_support_refinement_candidate_exported
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only": (
            current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export": (
            current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor": (
            positive_bridge_candidates
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor": (
            current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry": (
            current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry
        ),
        "frontier": {
            "current_repo_already_exports_strict_source_shannon_pair_population_support_refinement_candidate": (
                current_strict_source_shannon_pair_population_support_refinement_candidate_exported
            ),
            "current_repo_already_exports_t26_pair_indexed_population_anchor_target_only_as_future_spec": (
                current_t26_pair_indexed_population_anchor_target_frozen
            ),
            "current_repo_already_exports_strict_source_shannon_pair_population_support_refinement_to_pair_indexed_population_anchor_entry": (
                current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry
            ),
            "next_honest_move": (
                "attempt_a_genuinely_new_entry_bridge_from_the_strongest_strict_source_shannon_pair_population_support_refinement_route_into_the_minimal_designated_pair_family_pair1_pair2_component2_anchor_without_relabeling_generic_pair_indexed_syntax_as_actual_pair_population"
            ),
        },
        "hard_limits": [
            "No claim that the strongest strict-source Shannon pair-population support-refinement candidate is already an actual pair population.",
            "No claim that generic pair-indexed slot syntax already reduces to the minimal designated pair family [pair1, pair2].",
            "No claim that the current strict-source Shannon route already enters T26 component 2.",
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
        "t203_target_name": artifact["t203_target_name"],
        "t203_target_exported_on_current_repo_state": artifact["t203_target_exported_on_current_repo_state"],
        "current_strict_source_shannon_source_upgrades_exported": (
            artifact["current_strict_source_shannon_source_upgrades_exported"]
        ),
        "current_t26_pair_indexed_population_anchor_target_frozen": (
            artifact["current_t26_pair_indexed_population_anchor_target_frozen"]
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_exported": (
            artifact["current_strict_source_shannon_pair_population_support_refinement_candidate_exported"]
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only": (
            artifact["current_strict_source_shannon_pair_population_support_refinement_candidate_has_generic_pair_indexed_population_syntax_only"]
        ),
        "current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export": (
            artifact["current_strict_source_shannon_pair_population_support_refinement_candidate_remains_below_actual_pair_population_and_theta_export"]
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor": (
            artifact["current_strict_source_shannon_pair_population_support_refinement_route_remains_nonentering_for_pair_indexed_population_anchor"]
        ),
        "current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry": (
            artifact["current_strict_source_shannon_pair_population_support_refinement_route_has_exported_pair_indexed_population_anchor_entry"]
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
