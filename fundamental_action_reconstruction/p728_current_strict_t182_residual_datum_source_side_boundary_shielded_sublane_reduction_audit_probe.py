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

IN_P724 = GENERATED / "p724_current_strict_t178_positive_source_polarity_atlas_entry_corridor_reduction_audit_probe_summary.json"
IN_P727 = GENERATED / "p727_current_strict_t181_positive_corridor_excluded_negative_boundary_adjacency_chart_selection_bridge_nonexport_audit_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_F148 = GENERATED / "f148_first_actual_source_topology_basis_independent_promotion_witness_packet_summary.json"
IN_F150 = GENERATED / "f150_first_actual_declared_scope_source_topology_selector_theorem_packet_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"

OUT_JSON = GENERATED / "p728_current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_audit_probe.json"
OUT_SUMMARY = GENERATED / "p728_current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P724, IN_P727, IN_F147, IN_F148, IN_F150, IN_F301]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P728",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p724 = load_json(IN_P724)
    p727 = load_json(IN_P727)
    f147 = load_json(IN_F147)
    f148 = load_json(IN_F148)
    f150 = load_json(IN_F150)
    f301 = load_json(IN_F301)

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

    positive_corridor = list(p724.get("atlas_entry_roots_compatible_with_current_positive_source_polarity") or [])
    positive_boundary_adjacent_charts = list(p727.get("positive_boundary_adjacent_charts") or [])
    positive_boundary_shielded_charts = list(p727.get("positive_boundary_shielded_charts") or [])
    pair_index_set = list(f301.get("pair_index_set") or [])
    derived_parameters = f301.get("derived_parameters") or {}
    nad12_depth = f301.get("nad12_depth") or {}
    pair1_weights = list(nad12_depth.get("pair1_weights") or [])
    pair2_weights = list(nad12_depth.get("pair2_weights") or [])
    notes = list(f301.get("notes") or [])
    current_absence = list(f301.get("current_absence") or [])

    residual_datum_support_reduces_positive_corridor_to_boundary_shielded_sublane = (
        pair_index_set == positive_boundary_shielded_charts
    )
    residual_datum_support_intersection_with_boundary_adjacent = [
        pair for pair in pair_index_set if pair in positive_boundary_adjacent_charts
    ]
    residual_datum_pair12_support_is_symmetric = (
        float(derived_parameters.get("a") or 0.0) > 0.0
        and float(derived_parameters.get("b") or 0.0) > 0.0
        and abs(float(derived_parameters.get("a") or 0.0) - float(derived_parameters.get("b") or 0.0)) < 1.0e-12
        and pair1_weights == pair2_weights
    )
    tau_src_identified = {
        "F147": bool(f147.get("tau_src_identified_with_s_prelm")),
        "F148": bool(f148.get("tau_src_identified_with_s_prelm")),
        "F150": bool(f150.get("tau_src_identified_with_s_prelm")),
    }
    selector_neutral_boundary_intact = {
        "selector_neutral_note_present": any("selector-neutral" in note for note in notes),
        "no_selector_closure_claim_present": any("no strict-core selector closure" in item for item in current_absence),
        "no_qw2191_discharge_claim_present": any("no QW-2191 discharge" in item for item in current_absence),
    }

    add_check(
        "positive_corridor_still_present_after_p724",
        positive_corridor,
        ["pair1", "pair2", "pair3", "pair5"],
        "P724 still fixes the current positive atlas-entry corridor to four charts: pair1,pair2,pair3,pair5.",
    )
    add_check(
        "p727_boundary_split_already_localized",
        {
            "positive_boundary_adjacent_charts": positive_boundary_adjacent_charts,
            "positive_boundary_shielded_charts": positive_boundary_shielded_charts,
        },
        {
            "positive_boundary_adjacent_charts": ["pair3", "pair5"],
            "positive_boundary_shielded_charts": ["pair1", "pair2"],
        },
        "P727 already localizes the surviving positive corridor as boundary-adjacent charts pair3,pair5 versus boundary-shielded charts pair1,pair2.",
    )
    add_check(
        "f301_source_domain_is_tau_src_candidate",
        f301.get("source_domain"),
        "tau_src_candidate_v1",
        "The exported residual-datum support carrier F301 remains a genuine source-side carrier on tau_src_candidate_v1.",
    )
    add_check(
        "f301_pair_index_set_matches_boundary_shielded_sublane",
        pair_index_set,
        ["pair1", "pair2"],
        "The already-exported F301 residual-datum carrier names exactly pair1 and pair2.",
    )
    add_check(
        "f301_support_reduces_positive_corridor_to_boundary_shielded_sublane",
        residual_datum_support_reduces_positive_corridor_to_boundary_shielded_sublane,
        True,
        "That F301 pair-index support exactly matches the positive boundary-shielded sublane localized by P727, giving the first non-geometric source-side reduction below that split.",
    )
    add_check(
        "f301_support_avoids_positive_boundary_adjacent_charts",
        residual_datum_support_intersection_with_boundary_adjacent,
        [],
        "The same F301 support excludes the positive boundary-adjacent charts pair3 and pair5 on current exports.",
    )
    add_check(
        "f301_pair12_support_is_symmetric",
        residual_datum_pair12_support_is_symmetric,
        True,
        "The current F301 carrier is still exactly symmetric between pair1 and pair2: induced parameters satisfy a=b>0 and the exported pair1/pair2 weight profiles match.",
    )
    add_check(
        "f301_selector_neutral_boundary_intact",
        selector_neutral_boundary_intact,
        {
            "selector_neutral_note_present": True,
            "no_selector_closure_claim_present": True,
            "no_qw2191_discharge_claim_present": True,
        },
        "F301 remains selector-neutral and explicitly below strict-core selector closure or QW-2191 discharge.",
    )
    add_check(
        "tau_src_still_not_identified_with_current_selector_carrier",
        tau_src_identified,
        {
            "F147": False,
            "F148": False,
            "F150": False,
        },
        "Even after importing the F301 residual-datum carrier, tau_src is still not identified with the current selector carrier on the source-topology lane.",
    )
    add_check(
        "unique_chart_selected_within_boundary_shielded_sublane",
        len(pair_index_set) == 1,
        False,
        "Current residual-datum source-side support still does not select one unique chart inside the surviving boundary-shielded sublane.",
    )
    add_check(
        "t182_pair12_chart_selection_bridge_exported",
        False,
        False,
        "No current export provides the finer residual-datum source-side rule selecting between pair1 and pair2 inside the boundary-shielded sublane.",
    )

    status = (
        "PARTIAL_RESIDUAL_DATUM_SOURCE_SIDE_REDUCES_POSITIVE_CORRIDOR_TO_BOUNDARY_SHIELDED_SUBLANE_ONLY"
        if not blocking
        else "P728_REQUIRES_REVIEW_CHANGED_RESIDUAL_DATUM_BOUNDARY_SHIELDED_REDUCTION_STATE"
    )

    artifact = {
        "stage": "P728",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t182_residual_datum_source_side_boundary_shielded_sublane_reduction_only",
        "inputs": {
            "P724": str(IN_P724.relative_to(REPO)),
            "P727": str(IN_P727.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "F148": str(IN_F148.relative_to(REPO)),
            "F150": str(IN_F150.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t182_target_name": "ResidualDatumSourceSideBoundaryShieldedPair12ChartSelectionBridge_global_C_v1_strict_v1",
        "t182_target_exported_on_current_repo_state": False,
        "positive_corridor_roots": positive_corridor,
        "positive_boundary_adjacent_charts": positive_boundary_adjacent_charts,
        "positive_boundary_shielded_charts": positive_boundary_shielded_charts,
        "residual_datum_source_side_supported_positive_charts": pair_index_set,
        "current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane": (
            residual_datum_support_reduces_positive_corridor_to_boundary_shielded_sublane
        ),
        "residual_datum_source_side_support_excludes_boundary_adjacent_positive_charts": (
            residual_datum_support_intersection_with_boundary_adjacent == []
        ),
        "residual_datum_source_side_pair12_support_is_symmetric": residual_datum_pair12_support_is_symmetric,
        "unique_chart_selected_within_boundary_shielded_sublane": False,
        "audit_conclusion": {
            "current_repo_already_exports_first_non_geometric_source_side_positive_corridor_reduction": True,
            "current_source_side_reduction_lands_on_boundary_shielded_sublane": (
                residual_datum_support_reduces_positive_corridor_to_boundary_shielded_sublane
            ),
            "current_repo_already_exports_t182_target": False,
            "next_honest_move": (
                "export_or_attack_a_finer_residual_datum_source_side_rule_selecting_between_pair1_and_pair2"
            ),
        },
        "hard_limits": [
            "No T182 discharge claim.",
            "No claim that the current source-topology packet family itself discharges T181.",
            "No unique chart-seed selection claim.",
            "No identification of tau_src with the current selector carrier.",
            "No strict physical orientation datum claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P728",
        "status": status,
        "as_of": AS_OF,
        "t182_target_name": "ResidualDatumSourceSideBoundaryShieldedPair12ChartSelectionBridge_global_C_v1_strict_v1",
        "t182_target_exported_on_current_repo_state": False,
        "residual_datum_source_side_supported_positive_charts": pair_index_set,
        "current_residual_datum_source_side_support_reduces_positive_corridor_to_boundary_shielded_sublane": (
            residual_datum_support_reduces_positive_corridor_to_boundary_shielded_sublane
        ),
        "unique_chart_selected_within_boundary_shielded_sublane": False,
        "next_honest_move": "finer_residual_datum_source_side_rule_selecting_between_pair1_and_pair2",
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
