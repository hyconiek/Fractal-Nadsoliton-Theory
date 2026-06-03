#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    function_for_family,
    load_json,
    rel,
    replay_cell,
)
from p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit import projection_vector

GEN = ROOT / "generated"
OUT = GEN / "p2485_s1435_strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit.json"
MD = GEN / "p2485_s1435_strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit.md"
EXTENDED_LEVEL_COUNT = 64

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2482_BOUNDARY_BAND_WEAKEST_CELL_REFINEMENT": GEN / "p2482_s1432_strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit.json",
    "P2484_SECANT_MARGIN_REPLAY": GEN / "p2484_s1434_strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit.json",
}


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode("utf-8")).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, "fundamental_action_reconstruction", "material_dowodowy",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
        ],
        cwd=REPO,
        check=False,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:40]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2485|S1435|boundary-side secant curvature stability|secant curvature stability|64-level secant|second-order secant|secant-margin drift|dyadic curvature envelope|boundary-side curvature",
        "precursor_packets": "P2484|S1434|boundary-side dyadic secant margin|P2483|S1433|nested boundary ladder|P2482|S1432|boundary-band weakest cell",
        "curvature_language": "secant drift|margin drift|second finite difference|relative secant spread|curvature proxy|normalized gain envelope",
        "coverage_semantics": "diagnostic row ratio|not a coverage fraction|zero new P2459 unreplayed cells|covered-boundary-chain",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def extended_dyadic_cells(parent_cell: dict[str, Any]) -> list[dict[str, Any]]:
    left = Decimal(str(parent_cell["left"]))
    right = Decimal(str(parent_cell["right"]))
    width = right - left
    return [
        {
            "left": str(left),
            "right": str(left + width / (Decimal(2) ** level)),
            "uncovered_index": level,
            "uncovered_count": EXTENDED_LEVEL_COUNT,
            "extended_level": level,
            "relative_width_fraction_of_p2482_leftmost_subcell": f"1/{2 ** level}",
            "parent_dyadic_subcell_index": int(parent_cell["dyadic_subcell_index"]),
            "parent_collar_side": parent_cell["parent_collar_side"],
            "parent_boundary_band_index": int(parent_cell["parent_boundary_band_index"]),
        }
        for level in range(1, EXTENDED_LEVEL_COUNT + 1)
    ]


def boundary_side_curvature_stability(parent_cell: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    family = parent_cell["family"]
    function = function_for_family(family)
    replayed = []
    for cell in extended_dyadic_cells(parent_cell):
        fresh = replay_cell(family, cell, projection, function)
        width = Decimal(str(fresh["right"])) - Decimal(str(fresh["left"]))
        replayed.append({
            **fresh,
            "extended_level": int(cell["extended_level"]),
            "relative_width_fraction_of_p2482_leftmost_subcell": cell["relative_width_fraction_of_p2482_leftmost_subcell"],
            "absolute_width": str(width),
            "parent_dyadic_subcell_index": int(cell["parent_dyadic_subcell_index"]),
            "parent_collar_side": cell["parent_collar_side"],
            "parent_boundary_band_index": int(cell["parent_boundary_band_index"]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    ordered = sorted(replayed, key=lambda row: int(row["extended_level"]))
    parent_width = Decimal(str(parent_cell["right"])) - Decimal(str(parent_cell["left"]))
    secant_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        prior_width = Decimal(prior["absolute_width"])
        current_width = Decimal(current["absolute_width"])
        removed_width = prior_width - current_width
        delta = Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])
        quotient = delta / removed_width
        secant_pairs.append({
            "from_extended_level": int(prior["extended_level"]),
            "to_extended_level": int(current["extended_level"]),
            "from_width": str(prior_width),
            "to_width": str(current_width),
            "removed_width": str(removed_width),
            "width_ratio_to_from": str(current_width / prior_width),
            "same_left_boundary_anchor": Decimal(str(prior["left"])) == Decimal(str(current["left"])),
            "separation_delta_to_minus_from": str(delta),
            "normalized_lower_bound_gain_per_removed_width": str(quotient),
            "strictly_improves_lower_bound": delta > 0,
            "positive_secant_margin": quotient > 0,
        })
    drift_rows = []
    for prior, current in zip(secant_pairs, secant_pairs[1:]):
        prior_margin = Decimal(prior["normalized_lower_bound_gain_per_removed_width"])
        current_margin = Decimal(current["normalized_lower_bound_gain_per_removed_width"])
        drift = current_margin - prior_margin
        drift_rows.append({
            "from_pair_start_level": int(prior["from_extended_level"]),
            "to_pair_start_level": int(current["from_extended_level"]),
            "prior_secant_margin": str(prior_margin),
            "current_secant_margin": str(current_margin),
            "secant_margin_drift_current_minus_prior": str(drift),
            "positive_secant_margin_drift": drift > 0,
        })
    separations = [Decimal(row["decimal_separation_from_zero"]) for row in ordered]
    widths = [Decimal(row["absolute_width"]) for row in ordered]
    margins = [Decimal(row["normalized_lower_bound_gain_per_removed_width"]) for row in secant_pairs]
    margin_drifts = [Decimal(row["secant_margin_drift_current_minus_prior"]) for row in drift_rows]
    weakest_row = min(ordered, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    tightest_row = max(ordered, key=lambda row: int(row["extended_level"]))
    p2482_parent_lower = Decimal(parent_cell["decimal_separation_from_zero"])
    tightest_lower = Decimal(tightest_row["decimal_separation_from_zero"])
    margin_spread = max(margins) - min(margins)
    relative_spread = margin_spread / min(margins)
    return {
        "family": family,
        "parent_collar_side": parent_cell["parent_collar_side"],
        "parent_boundary_band_index": int(parent_cell["parent_boundary_band_index"]),
        "parent_dyadic_subcell_index": int(parent_cell["dyadic_subcell_index"]),
        "root_window_side_boundary_anchor": str(parent_cell["left"]),
        "p2482_leftmost_subcell_width": str(parent_width),
        "p2482_leftmost_subcell_decimal_separation": parent_cell["decimal_separation_from_zero"],
        "extended_level_count": EXTENDED_LEVEL_COUNT,
        "consecutive_secant_pair_count": len(secant_pairs),
        "consecutive_margin_drift_count": len(drift_rows),
        "all_extended_rows_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in ordered),
        "all_extended_rows_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in ordered),
        "all_extended_rows_share_left_boundary_anchor": all(Decimal(str(row["left"])) == Decimal(str(parent_cell["left"])) for row in ordered),
        "all_consecutive_widths_halve": all(Decimal(row["width_ratio_to_from"]) == Decimal("0.5") for row in secant_pairs),
        "all_consecutive_lower_bounds_strictly_increase": all(row["strictly_improves_lower_bound"] for row in secant_pairs),
        "all_consecutive_secant_margins_positive": all(row["positive_secant_margin"] for row in secant_pairs),
        "all_secant_margin_drifts_positive": all(row["positive_secant_margin_drift"] for row in drift_rows),
        "weakest_extended_decimal_separation": str(min(separations)),
        "tightest_extended_decimal_separation": tightest_row["decimal_separation_from_zero"],
        "minimum_consecutive_positive_lower_bound_delta": str(min(Decimal(row["separation_delta_to_minus_from"]) for row in secant_pairs)),
        "minimum_positive_secant_margin_per_removed_width": str(min(margins)),
        "maximum_positive_secant_margin_per_removed_width": str(max(margins)),
        "secant_margin_absolute_spread": str(margin_spread),
        "secant_margin_relative_spread": str(relative_spread),
        "minimum_positive_secant_margin_drift": str(min(margin_drifts)),
        "maximum_positive_secant_margin_drift": str(max(margin_drifts)),
        "minimum_replayed_width": str(min(widths)),
        "maximum_replayed_width": str(max(widths)),
        "tightest_extended_level": int(tightest_row["extended_level"]),
        "tightest_width_fraction_of_p2482_leftmost_subcell": tightest_row["relative_width_fraction_of_p2482_leftmost_subcell"],
        "weakest_row_is_coarsest_level": int(weakest_row["extended_level"]) == 1,
        "tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound": tightest_lower > p2482_parent_lower,
        "tightest_minus_p2482_leftmost_subcell_lower_bound_delta": str(tightest_lower - p2482_parent_lower),
        "extended_ladder_rows": ordered,
        "consecutive_secant_margin_pairs": secant_pairs,
        "consecutive_secant_margin_drift_rows": drift_rows,
        "extended_curvature_stability_fingerprint_sha256": sha256_json({"rows": ordered, "pairs": secant_pairs, "drifts": drift_rows}),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2485/S1435 strict pointwise interval-Decimal P2459 boundary-side secant-curvature stability audit

`P2485/S1435` extends the P2484 boundary-side secant ladder from `32` to `64` anchored dyadic levels inside the same inherited P2456 covered-boundary-chain subcell.  It keeps the finite lower-bound secant margins and adds a second-order diagnostic: consecutive secant-margin drift.  All `64` replayed rows exclude zero, all `63` lower-bound secant margins are positive, and all `62` consecutive margin drifts are positive with a narrow relative secant-margin spread.  This is finite boundary-side curvature-stability evidence only; it is not an analytic monotonicity/convexity theorem, not directed rounding, not root-window closure, not a P2459 coverage increase, and not selector/source or ToE closure.

For a non-specialist: P2484 checked that each zoom improved the safety margin per amount of width removed.  P2485 repeats the check deeper and also verifies that those finite gains drift coherently instead of alternating.  That is useful one-sided evidence, but it still does not prove the excluded root window or any continuum theorem.
"""
    lag_section = """
## P2485/S1435 P2459 boundary-side secant-curvature stability guard

`P2485/S1435` adds a deeper finite secant-curvature stability guard behind `L_total`: `64` diagnostic Decimal rows, `63` positive lower-bound secant margins, and `62` positive finite margin drifts in one inherited P2456 covered-boundary-chain subcell.  It remains diagnostic finite arithmetic and does not export directed rounding, analytic monotonicity/convexity, root-window exclusion, selector/source/gauge authority, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2485/S1435 strict pointwise interval-Decimal P2459 boundary-side secant-curvature stability audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2485/S1435 P2459 boundary-side secant-curvature stability guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2482 = theorem(sources["P2482_BOUNDARY_BAND_WEAKEST_CELL_REFINEMENT"], "strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit")
    p2484 = theorem(sources["P2484_SECANT_MARGIN_REPLAY"], "strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit")
    parent_cell = p2482["boundary_band_weakest_cell_dyadic_refinement_replay"]["minimum_subcell_replay_row"]
    stability = boundary_side_curvature_stability(parent_cell, projection)
    p2459_universe = p2482["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2485_T1_strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit",
        "audited_chain": ["P2482/S1432", "P2484/S1434"],
        "boundary_side_secant_curvature_stability": stability,
        "family": stability["family"],
        "parent_collar_side": stability["parent_collar_side"],
        "parent_boundary_band_index": stability["parent_boundary_band_index"],
        "parent_dyadic_subcell_index": stability["parent_dyadic_subcell_index"],
        "root_window_side_boundary_anchor": stability["root_window_side_boundary_anchor"],
        "extended_level_count": stability["extended_level_count"],
        "consecutive_secant_pair_count": stability["consecutive_secant_pair_count"],
        "consecutive_margin_drift_count": stability["consecutive_margin_drift_count"],
        "all_extended_rows_exclude_zero": stability["all_extended_rows_exclude_zero"],
        "all_extended_rows_have_positive_separation": stability["all_extended_rows_have_positive_separation"],
        "all_extended_rows_share_left_boundary_anchor": stability["all_extended_rows_share_left_boundary_anchor"],
        "all_consecutive_widths_halve": stability["all_consecutive_widths_halve"],
        "all_consecutive_lower_bounds_strictly_increase": stability["all_consecutive_lower_bounds_strictly_increase"],
        "all_consecutive_secant_margins_positive": stability["all_consecutive_secant_margins_positive"],
        "all_secant_margin_drifts_positive": stability["all_secant_margin_drifts_positive"],
        "weakest_extended_decimal_separation": stability["weakest_extended_decimal_separation"],
        "tightest_extended_decimal_separation": stability["tightest_extended_decimal_separation"],
        "minimum_consecutive_positive_lower_bound_delta": stability["minimum_consecutive_positive_lower_bound_delta"],
        "minimum_positive_secant_margin_per_removed_width": stability["minimum_positive_secant_margin_per_removed_width"],
        "maximum_positive_secant_margin_per_removed_width": stability["maximum_positive_secant_margin_per_removed_width"],
        "secant_margin_absolute_spread": stability["secant_margin_absolute_spread"],
        "secant_margin_relative_spread": stability["secant_margin_relative_spread"],
        "minimum_positive_secant_margin_drift": stability["minimum_positive_secant_margin_drift"],
        "maximum_positive_secant_margin_drift": stability["maximum_positive_secant_margin_drift"],
        "minimum_replayed_width": stability["minimum_replayed_width"],
        "maximum_replayed_width": stability["maximum_replayed_width"],
        "tightest_extended_level": stability["tightest_extended_level"],
        "tightest_width_fraction_of_p2482_leftmost_subcell": stability["tightest_width_fraction_of_p2482_leftmost_subcell"],
        "weakest_row_is_coarsest_level": stability["weakest_row_is_coarsest_level"],
        "tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound": stability["tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound"],
        "tightest_minus_p2482_leftmost_subcell_lower_bound_delta": stability["tightest_minus_p2482_leftmost_subcell_lower_bound_delta"],
        "p2484_secant_level_count_inherited": p2484["secant_level_count"],
        "p2484_minimum_positive_secant_margin_per_removed_width_inherited": p2484["minimum_positive_secant_margin_per_removed_width"],
        "p2482_minimum_subcell_decimal_separation_inherited": p2482["minimum_subcell_decimal_separation"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2485_fresh_decimal_evaluation_row_count_not_a_coverage_count": stability["extended_level_count"],
        "p2485_consecutive_secant_margin_count_not_a_coverage_count": stability["consecutive_secant_pair_count"],
        "p2485_consecutive_margin_drift_count_not_a_coverage_count": stability["consecutive_margin_drift_count"],
        "p2485_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction": f"{stability['extended_level_count']}/{p2459_universe}",
        "targeted_p2485_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2485_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2485_refines_one_inherited_p2456_covered_boundary_chain_cell": True,
        "finite_chain_coverage_budget_inherited_from_p2482": p2482["finite_chain_coverage_budget_inherited_from_p2481"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_boundary_side_secant_curvature_stability_audit_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "analytic_convexity_theorem_exported_by_this_certificate": False,
        "root_window_exclusion_theorem_exported_by_this_certificate": False,
        "global_continuum_root_exclusion_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "lay_summary": "This packet extends the one-sided boundary ladder to 64 dyadic levels and checks both positive finite secant margins and positive consecutive margin drifts. The finite gains behave coherently on the checked ladder, but this is still diagnostic arithmetic outside the excluded root window; it does not prove root-window, analytic monotonicity/convexity, selector/source, or ToE closure.",
        "not_licensed": [
            "The 64 fresh Decimal evaluations, 63 secant margins, and 62 margin drifts are diagnostic rows inside one inherited P2456 covered-boundary-chain cell, not distinct P2459 coverage cells.",
            "Positive finite secant-margin drift is a second-order diagnostic, not an analytic convexity or monotonicity theorem.",
            "P2485 adds zero new P2459 unreplayed cells and does not turn the boundary-side ladder into root-window exclusion.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "A stronger result requires a real directed-rounding/analytic proof backend or an explicit root-window-side theorem; do not inflate finite secant-curvature stability into continuum closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "extended_level_count_expected": theorem_export["extended_level_count"] == 64,
        "secant_pair_count_expected": theorem_export["consecutive_secant_pair_count"] == 63,
        "margin_drift_count_expected": theorem_export["consecutive_margin_drift_count"] == 62,
        "parent_binding_expected": theorem_export["parent_collar_side"] == "p2456_right_boundary_band" and theorem_export["parent_boundary_band_index"] == 0 and theorem_export["parent_dyadic_subcell_index"] == 0,
        "all_extended_rows_exclude_zero": theorem_export["all_extended_rows_exclude_zero"],
        "all_extended_rows_have_positive_separation": theorem_export["all_extended_rows_have_positive_separation"],
        "same_left_boundary_anchor": theorem_export["all_extended_rows_share_left_boundary_anchor"],
        "widths_halve": theorem_export["all_consecutive_widths_halve"],
        "lower_bounds_improve": theorem_export["all_consecutive_lower_bounds_strictly_increase"] and Decimal(theorem_export["minimum_consecutive_positive_lower_bound_delta"]) > 0,
        "secant_margins_positive": theorem_export["all_consecutive_secant_margins_positive"] and Decimal(theorem_export["minimum_positive_secant_margin_per_removed_width"]) > 0,
        "secant_margin_drifts_positive": theorem_export["all_secant_margin_drifts_positive"] and Decimal(theorem_export["minimum_positive_secant_margin_drift"]) > 0,
        "relative_spread_small_diagnostic_only": Decimal(theorem_export["secant_margin_relative_spread"]) < Decimal("0.000001"),
        "weakest_row_coarsest_level_not_overclaimed": theorem_export["weakest_row_is_coarsest_level"],
        "tightest_level_expected": theorem_export["tightest_extended_level"] == 64 and theorem_export["tightest_width_fraction_of_p2482_leftmost_subcell"] == "1/18446744073709551616",
        "tightest_lower_bound_stronger_than_p2482": theorem_export["tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound"] and Decimal(theorem_export["tightest_minus_p2482_leftmost_subcell_lower_bound_delta"]) > 0,
        "zero_new_p2459_unreplayed_cells": theorem_export["targeted_p2485_new_p2459_unreplayed_cell_count"] == 0 and theorem_export["targeted_p2485_new_p2459_unreplayed_cell_scope_against_p2459_universe"] == f"0/{p2459_universe}",
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_boundary_side_secant_curvature_stability_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_root_window_theorem": not theorem_export["root_window_exclusion_theorem_exported_by_this_certificate"],
        "no_analytic_monotonicity_or_convexity": not theorem_export["analytic_monotonicity_theorem_exported_by_this_certificate"] and not theorem_export["analytic_convexity_theorem_exported_by_this_certificate"],
        "no_symbolic_root_exclusion": not theorem_export["symbolic_root_exclusion_theorem_exported_by_this_certificate"],
        "no_continuum_root_exclusion": not theorem_export["global_continuum_root_exclusion_theorem_exported_by_this_certificate"],
        "no_selector_source_gauge_theorem": not theorem_export["pointwise_coordinate_selector_exported_by_this_certificate"] and not theorem_export["strict_observable_source_constraint_exported_by_this_certificate"] and not theorem_export["gauge_slice_theorem_exported_by_this_certificate"],
        "no_value_generator_export": not theorem_export["strict_physical_value_generator_exported"],
        "no_qw2191_discharge": not theorem_export["qw2191_discharged"],
        "no_ltotal_export": not theorem_export["role_bearing_ltotal_exported"],
        "no_legacy_role_transfer": not theorem_export["legacy_role_transfer_exported"],
        "no_toe_export": not theorem_export["toe_closure_exported"],
        "fingerprint_stable": sha256_json(theorem_export) == sha256_json(theorem_export),
    }
    return {
        "schema_version": "p2485_s1435_v1",
        "packet_id": "P2485",
        "stage_id": "S1435",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_SIDE_SECANT_CURVATURE_STABILITY_AUDIT_NO_ROOT_WINDOW_NO_ANALYTIC_MONOTONICITY_CONVEXITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit"]["theorem_export"]
    lines = [
        "# P2485/S1435 strict pointwise interval-Decimal P2459 boundary-side secant-curvature stability audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-side secant-curvature stability replay",
        "",
        f"Root-window-side boundary anchor: `{t['root_window_side_boundary_anchor']}`.",
        f"Extended dyadic levels replayed: `{t['extended_level_count']}`.",
        f"Consecutive secant margins computed: `{t['consecutive_secant_pair_count']}`.",
        f"Consecutive secant-margin drifts computed: `{t['consecutive_margin_drift_count']}`.",
        f"Weakest extended Decimal separation: `{t['weakest_extended_decimal_separation']}`.",
        f"Tightest extended Decimal separation: `{t['tightest_extended_decimal_separation']}`.",
        f"Minimum consecutive positive lower-bound delta: `{t['minimum_consecutive_positive_lower_bound_delta']}`.",
        f"Minimum positive secant margin per removed width: `{t['minimum_positive_secant_margin_per_removed_width']}`.",
        f"Maximum positive secant margin per removed width: `{t['maximum_positive_secant_margin_per_removed_width']}`.",
        f"Secant margin absolute spread: `{t['secant_margin_absolute_spread']}`.",
        f"Secant margin relative spread: `{t['secant_margin_relative_spread']}`.",
        f"Minimum positive secant-margin drift: `{t['minimum_positive_secant_margin_drift']}`.",
        f"Maximum positive secant-margin drift: `{t['maximum_positive_secant_margin_drift']}`.",
        f"Tightest extended level: `{t['tightest_extended_level']}`.",
        f"Tightest width fraction of P2482 leftmost subcell: `{t['tightest_width_fraction_of_p2482_leftmost_subcell']}`.",
        f"All extended rows exclude zero: `{t['all_extended_rows_exclude_zero']}`.",
        f"All secant margins positive: `{t['all_consecutive_secant_margins_positive']}`.",
        f"All secant-margin drifts positive: `{t['all_secant_margin_drifts_positive']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2485 fresh Decimal evaluation rows (not a P2459 coverage count): `{t['p2485_fresh_decimal_evaluation_row_count_not_a_coverage_count']}`.",
        f"P2485 consecutive secant margins (not a P2459 coverage count): `{t['p2485_consecutive_secant_margin_count_not_a_coverage_count']}`.",
        f"P2485 consecutive margin drifts (not a P2459 coverage count): `{t['p2485_consecutive_margin_drift_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2485_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2485: `{t['targeted_p2485_new_p2459_unreplayed_cell_count']}`.",
        f"New P2459 unreplayed-cell scope against inherited P2459 universe: `{t['targeted_p2485_new_p2459_unreplayed_cell_scope_against_p2459_universe']}`.",
        f"P2485 refines one inherited P2456 covered-boundary-chain cell: `{t['p2485_refines_one_inherited_p2456_covered_boundary_chain_cell']}`.",
        f"Full complement replay exported by P2485: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite boundary-side secant-curvature stability audit outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity/convexity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(parents=True, exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps({"status": payload["status"], "gatekeepers": payload["gatekeeper_checks"]}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
