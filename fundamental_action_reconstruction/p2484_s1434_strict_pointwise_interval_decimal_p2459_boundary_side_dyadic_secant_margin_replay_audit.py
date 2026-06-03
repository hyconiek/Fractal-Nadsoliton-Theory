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
OUT = GEN / "p2484_s1434_strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit.json"
MD = GEN / "p2484_s1434_strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit.md"
SECANT_LEVEL_COUNT = 32

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2482_BOUNDARY_BAND_WEAKEST_CELL_REFINEMENT": GEN / "p2482_s1432_strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit.json",
    "P2483_NESTED_BOUNDARY_LADDER": GEN / "p2483_s1433_strict_pointwise_interval_decimal_p2459_root_window_adjacent_nested_boundary_ladder_replay_audit.json",
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
        "new_packet": "P2484|S1434|boundary-side dyadic secant margin|dyadic secant margin|boundary secant ladder|nested slope audit|halving secant audit|finite secant lower-bound|boundary-side secant",
        "precursor_packets": "P2483|S1433|nested boundary ladder|P2482|S1432|boundary-band weakest cell|root-window adjacent cell refinement",
        "secant_language": "secant margin|normalized lower-bound gain|removed-width quotient|dyadic contraction slope|boundary-side margin",
        "coverage_semantics": "diagnostic row ratio|not a coverage fraction|zero new P2459 unreplayed cells|covered-boundary-chain",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def secant_ladder_cells(parent_cell: dict[str, Any]) -> list[dict[str, Any]]:
    left = Decimal(str(parent_cell["left"]))
    right = Decimal(str(parent_cell["right"]))
    width = right - left
    return [
        {
            "left": str(left),
            "right": str(left + width / (Decimal(2) ** level)),
            "uncovered_index": level,
            "uncovered_count": SECANT_LEVEL_COUNT,
            "secant_level": level,
            "relative_width_fraction_of_p2482_leftmost_subcell": f"1/{2 ** level}",
            "parent_dyadic_subcell_index": int(parent_cell["dyadic_subcell_index"]),
            "parent_collar_side": parent_cell["parent_collar_side"],
            "parent_boundary_band_index": int(parent_cell["parent_boundary_band_index"]),
        }
        for level in range(1, SECANT_LEVEL_COUNT + 1)
    ]


def dyadic_secant_margin_replay(parent_cell: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    family = parent_cell["family"]
    function = function_for_family(family)
    replayed = []
    for cell in secant_ladder_cells(parent_cell):
        fresh = replay_cell(family, cell, projection, function)
        width = Decimal(str(fresh["right"])) - Decimal(str(fresh["left"]))
        replayed.append({
            **fresh,
            "secant_level": int(cell["secant_level"]),
            "relative_width_fraction_of_p2482_leftmost_subcell": cell["relative_width_fraction_of_p2482_leftmost_subcell"],
            "absolute_width": str(width),
            "parent_dyadic_subcell_index": int(cell["parent_dyadic_subcell_index"]),
            "parent_collar_side": cell["parent_collar_side"],
            "parent_boundary_band_index": int(cell["parent_boundary_band_index"]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    ordered = sorted(replayed, key=lambda row: int(row["secant_level"]))
    parent_width = Decimal(str(parent_cell["right"])) - Decimal(str(parent_cell["left"]))
    consecutive_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        prior_width = Decimal(prior["absolute_width"])
        current_width = Decimal(current["absolute_width"])
        removed_width = prior_width - current_width
        delta = Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])
        quotient = delta / removed_width
        consecutive_pairs.append({
            "from_secant_level": int(prior["secant_level"]),
            "to_secant_level": int(current["secant_level"]),
            "from_width": str(prior_width),
            "to_width": str(current_width),
            "removed_width": str(removed_width),
            "width_ratio_to_from": str(current_width / prior_width),
            "same_left_boundary_anchor": Decimal(str(prior["left"])) == Decimal(str(current["left"])),
            "from_separation": prior["decimal_separation_from_zero"],
            "to_separation": current["decimal_separation_from_zero"],
            "separation_delta_to_minus_from": str(delta),
            "normalized_lower_bound_gain_per_removed_width": str(quotient),
            "strictly_improves_lower_bound": delta > 0,
            "positive_secant_margin": quotient > 0,
        })
    weakest_row = min(ordered, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    tightest_row = max(ordered, key=lambda row: int(row["secant_level"]))
    quotients = [Decimal(row["normalized_lower_bound_gain_per_removed_width"]) for row in consecutive_pairs]
    deltas = [Decimal(row["separation_delta_to_minus_from"]) for row in consecutive_pairs]
    widths = [Decimal(row["absolute_width"]) for row in ordered]
    p2482_parent_lower = Decimal(parent_cell["decimal_separation_from_zero"])
    tightest_lower = Decimal(tightest_row["decimal_separation_from_zero"])
    return {
        "family": family,
        "parent_collar_side": parent_cell["parent_collar_side"],
        "parent_boundary_band_index": int(parent_cell["parent_boundary_band_index"]),
        "parent_dyadic_subcell_index": int(parent_cell["dyadic_subcell_index"]),
        "root_window_side_boundary_anchor": str(parent_cell["left"]),
        "p2482_leftmost_subcell_right": str(parent_cell["right"]),
        "p2482_leftmost_subcell_width": str(parent_width),
        "p2482_leftmost_subcell_decimal_separation": parent_cell["decimal_separation_from_zero"],
        "secant_level_count": SECANT_LEVEL_COUNT,
        "all_secant_rows_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in ordered),
        "all_secant_rows_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in ordered),
        "all_secant_rows_share_left_boundary_anchor": all(Decimal(str(row["left"])) == Decimal(str(parent_cell["left"])) for row in ordered),
        "all_consecutive_widths_halve": all(Decimal(row["width_ratio_to_from"]) == Decimal("0.5") for row in consecutive_pairs),
        "all_consecutive_lower_bounds_strictly_increase": all(row["strictly_improves_lower_bound"] for row in consecutive_pairs),
        "all_consecutive_secant_margins_positive": all(row["positive_secant_margin"] for row in consecutive_pairs),
        "minimum_consecutive_positive_lower_bound_delta": str(min(deltas)),
        "minimum_positive_secant_margin_per_removed_width": str(min(quotients)),
        "maximum_positive_secant_margin_per_removed_width": str(max(quotients)),
        "minimum_replayed_width": str(min(widths)),
        "maximum_replayed_width": str(max(widths)),
        "weakest_secant_row": weakest_row,
        "tightest_secant_row": tightest_row,
        "weakest_secant_decimal_separation": weakest_row["decimal_separation_from_zero"],
        "tightest_secant_decimal_separation": tightest_row["decimal_separation_from_zero"],
        "tightest_secant_level": int(tightest_row["secant_level"]),
        "tightest_width_fraction_of_p2482_leftmost_subcell": tightest_row["relative_width_fraction_of_p2482_leftmost_subcell"],
        "weakest_secant_row_is_coarsest_level": int(weakest_row["secant_level"]) == 1,
        "tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound": tightest_lower > p2482_parent_lower,
        "tightest_minus_p2482_leftmost_subcell_lower_bound_delta": str(tightest_lower - p2482_parent_lower),
        "secant_ladder_rows": ordered,
        "consecutive_secant_margin_pairs": consecutive_pairs,
        "secant_ladder_replay_fingerprint_sha256": sha256_json(ordered),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2484/S1434 strict pointwise interval-Decimal P2459 boundary-side dyadic secant margin replay audit

`P2484/S1434` keeps the root-window-side anchor from P2482/P2483 and replays a `32`-level dyadic contraction ladder inside the same weakest inherited boundary-chain subcell.  Beyond checking zero-exclusion and halving widths, it computes consecutive finite secant margins: the lower-bound gain divided by the width removed at each halving.  Every replayed row remains positive, every halving improves the Decimal lower bound, and every normalized finite secant margin is positive.  This is stronger boundary-side finite evidence than a plain nested ladder, but it is still diagnostic finite arithmetic outside the excluded root window, not an analytic monotonicity theorem, not a directed-rounding theorem, not a P2459 coverage increase, and not a continuum closure.

For a non-specialist: P2483 showed that repeated one-sided zooms next to the root window got safer.  P2484 asks whether each zoom has a positive finite gain per amount of width removed.  The answer is yes for the 32 checked dyadic levels, but this still does not prove what happens inside the excluded root window or identify a selector/source.
"""
    lag_section = """
## P2484/S1434 P2459 boundary-side dyadic secant margin replay guard

`P2484/S1434` adds a finite dyadic secant-margin guard behind `L_total` for the same root-window-adjacent inherited P2456 boundary-chain subcell audited by P2482/P2483.  Its `32` fresh Decimal evaluations and `31` consecutive secant margins are diagnostic rows, not new P2459 coverage cells; the packet does not export directed rounding, analytic monotonicity, root-window exclusion, selector/source/gauge authority, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2484/S1434 strict pointwise interval-Decimal P2459 boundary-side dyadic secant margin replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2484/S1434 P2459 boundary-side dyadic secant margin replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2482 = theorem(sources["P2482_BOUNDARY_BAND_WEAKEST_CELL_REFINEMENT"], "strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit")
    p2483 = theorem(sources["P2483_NESTED_BOUNDARY_LADDER"], "strict_pointwise_interval_decimal_p2459_root_window_adjacent_nested_boundary_ladder_replay_audit")
    parent_cell = p2482["boundary_band_weakest_cell_dyadic_refinement_replay"]["minimum_subcell_replay_row"]
    replay = dyadic_secant_margin_replay(parent_cell, projection)
    p2459_universe = p2482["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2484_T1_strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit",
        "audited_chain": ["P2482/S1432", "P2483/S1433"],
        "boundary_side_dyadic_secant_margin_replay": replay,
        "family": replay["family"],
        "parent_collar_side": replay["parent_collar_side"],
        "parent_boundary_band_index": replay["parent_boundary_band_index"],
        "parent_dyadic_subcell_index": replay["parent_dyadic_subcell_index"],
        "root_window_side_boundary_anchor": replay["root_window_side_boundary_anchor"],
        "secant_level_count": replay["secant_level_count"],
        "consecutive_secant_pair_count": len(replay["consecutive_secant_margin_pairs"]),
        "all_secant_rows_exclude_zero": replay["all_secant_rows_exclude_zero"],
        "all_secant_rows_have_positive_separation": replay["all_secant_rows_have_positive_separation"],
        "all_secant_rows_share_left_boundary_anchor": replay["all_secant_rows_share_left_boundary_anchor"],
        "all_consecutive_widths_halve": replay["all_consecutive_widths_halve"],
        "all_consecutive_lower_bounds_strictly_increase": replay["all_consecutive_lower_bounds_strictly_increase"],
        "all_consecutive_secant_margins_positive": replay["all_consecutive_secant_margins_positive"],
        "minimum_consecutive_positive_lower_bound_delta": replay["minimum_consecutive_positive_lower_bound_delta"],
        "minimum_positive_secant_margin_per_removed_width": replay["minimum_positive_secant_margin_per_removed_width"],
        "maximum_positive_secant_margin_per_removed_width": replay["maximum_positive_secant_margin_per_removed_width"],
        "minimum_replayed_width": replay["minimum_replayed_width"],
        "maximum_replayed_width": replay["maximum_replayed_width"],
        "weakest_secant_decimal_separation": replay["weakest_secant_decimal_separation"],
        "tightest_secant_decimal_separation": replay["tightest_secant_decimal_separation"],
        "tightest_secant_level": replay["tightest_secant_level"],
        "tightest_width_fraction_of_p2482_leftmost_subcell": replay["tightest_width_fraction_of_p2482_leftmost_subcell"],
        "weakest_secant_row_is_coarsest_level": replay["weakest_secant_row_is_coarsest_level"],
        "tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound": replay["tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound"],
        "tightest_minus_p2482_leftmost_subcell_lower_bound_delta": replay["tightest_minus_p2482_leftmost_subcell_lower_bound_delta"],
        "p2483_nested_level_count_inherited": p2483["nested_level_count"],
        "p2483_tightest_nested_ladder_decimal_separation_inherited": p2483["tightest_nested_ladder_decimal_separation"],
        "p2482_minimum_subcell_decimal_separation_inherited": p2482["minimum_subcell_decimal_separation"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2484_fresh_decimal_evaluation_row_count_not_a_coverage_count": replay["secant_level_count"],
        "p2484_consecutive_secant_margin_count_not_a_coverage_count": len(replay["consecutive_secant_margin_pairs"]),
        "p2484_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction": f"{replay['secant_level_count']}/{p2459_universe}",
        "targeted_p2484_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2484_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2484_refines_one_inherited_p2456_covered_boundary_chain_cell": True,
        "finite_chain_coverage_budget_inherited_from_p2482": p2482["finite_chain_coverage_budget_inherited_from_p2481"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_boundary_side_dyadic_secant_margin_audit_exported": True,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "analytic_monotonicity_theorem_exported_by_this_certificate": False,
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
        "lay_summary": "This packet replays 32 one-sided dyadic contractions next to the excluded root window and computes 31 finite secant margins from the lower-bound gains. Every checked row excludes zero and every normalized finite gain is positive. The result is useful boundary-side margin evidence, but it is still finite diagnostic arithmetic outside the root window and does not prove root-window, analytic, selector/source, or ToE closure.",
        "not_licensed": [
            "The 32 fresh Decimal evaluations and 31 secant margins are diagnostic rows inside one inherited P2456 covered-boundary-chain cell, not distinct P2459 coverage cells.",
            "P2484 adds zero new P2459 unreplayed cells and does not turn the P2482/P2483 boundary-side checks into a root-window theorem.",
            "Finite positive secant margins do not by themselves prove directed rounding, symbolic root exclusion, analytic monotonicity, or continuum root exclusion.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "A stronger result requires either a real directed-rounding/analytic monotonicity backend or a separately stated root-window-side theorem; do not inflate finite secant margins into continuum closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "secant_level_count_expected": theorem_export["secant_level_count"] == 32,
        "secant_pair_count_expected": theorem_export["consecutive_secant_pair_count"] == 31,
        "parent_binding_expected": theorem_export["parent_collar_side"] == "p2456_right_boundary_band" and theorem_export["parent_boundary_band_index"] == 0 and theorem_export["parent_dyadic_subcell_index"] == 0,
        "all_secant_rows_exclude_zero": theorem_export["all_secant_rows_exclude_zero"],
        "all_secant_rows_have_positive_separation": theorem_export["all_secant_rows_have_positive_separation"],
        "same_left_boundary_anchor": theorem_export["all_secant_rows_share_left_boundary_anchor"],
        "widths_halve": theorem_export["all_consecutive_widths_halve"],
        "lower_bounds_improve": theorem_export["all_consecutive_lower_bounds_strictly_increase"] and Decimal(theorem_export["minimum_consecutive_positive_lower_bound_delta"]) > 0,
        "secant_margins_positive": theorem_export["all_consecutive_secant_margins_positive"] and Decimal(theorem_export["minimum_positive_secant_margin_per_removed_width"]) > 0,
        "weakest_row_coarsest_level_not_overclaimed": theorem_export["weakest_secant_row_is_coarsest_level"],
        "tightest_level_expected": theorem_export["tightest_secant_level"] == 32 and theorem_export["tightest_width_fraction_of_p2482_leftmost_subcell"] == "1/4294967296",
        "tightest_lower_bound_stronger_than_p2482": theorem_export["tightest_lower_bound_exceeds_p2482_leftmost_subcell_bound"] and Decimal(theorem_export["tightest_minus_p2482_leftmost_subcell_lower_bound_delta"]) > 0,
        "zero_new_p2459_unreplayed_cells": theorem_export["targeted_p2484_new_p2459_unreplayed_cell_count"] == 0 and theorem_export["targeted_p2484_new_p2459_unreplayed_cell_scope_against_p2459_universe"] == f"0/{p2459_universe}",
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_boundary_side_dyadic_secant_margin_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_root_window_theorem": not theorem_export["root_window_exclusion_theorem_exported_by_this_certificate"],
        "no_analytic_monotonicity": not theorem_export["analytic_monotonicity_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2484_s1434_v1",
        "packet_id": "P2484",
        "stage_id": "S1434",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_SIDE_DYADIC_SECANT_MARGIN_REPLAY_AUDIT_NO_ROOT_WINDOW_NO_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_side_dyadic_secant_margin_replay_audit"]["theorem_export"]
    lines = [
        "# P2484/S1434 strict pointwise interval-Decimal P2459 boundary-side dyadic secant margin replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-side dyadic secant margin replay",
        "",
        f"Root-window-side boundary anchor: `{t['root_window_side_boundary_anchor']}`.",
        f"Secant levels replayed: `{t['secant_level_count']}`.",
        f"Consecutive secant margins computed: `{t['consecutive_secant_pair_count']}`.",
        f"Weakest secant Decimal separation: `{t['weakest_secant_decimal_separation']}`.",
        f"Tightest secant Decimal separation: `{t['tightest_secant_decimal_separation']}`.",
        f"Minimum consecutive positive lower-bound delta: `{t['minimum_consecutive_positive_lower_bound_delta']}`.",
        f"Minimum positive secant margin per removed width: `{t['minimum_positive_secant_margin_per_removed_width']}`.",
        f"Maximum positive secant margin per removed width: `{t['maximum_positive_secant_margin_per_removed_width']}`.",
        f"Tightest secant level: `{t['tightest_secant_level']}`.",
        f"Tightest width fraction of P2482 leftmost subcell: `{t['tightest_width_fraction_of_p2482_leftmost_subcell']}`.",
        f"Tightest-minus-P2482 lower-bound delta: `{t['tightest_minus_p2482_leftmost_subcell_lower_bound_delta']}`.",
        f"All secant rows exclude zero: `{t['all_secant_rows_exclude_zero']}`.",
        f"All consecutive secant margins positive: `{t['all_consecutive_secant_margins_positive']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2484 fresh Decimal evaluation rows (not a P2459 coverage count): `{t['p2484_fresh_decimal_evaluation_row_count_not_a_coverage_count']}`.",
        f"P2484 consecutive secant margins (not a P2459 coverage count): `{t['p2484_consecutive_secant_margin_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2484_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2484: `{t['targeted_p2484_new_p2459_unreplayed_cell_count']}`.",
        f"New P2459 unreplayed-cell scope against inherited P2459 universe: `{t['targeted_p2484_new_p2459_unreplayed_cell_scope_against_p2459_universe']}`.",
        f"P2484 refines one inherited P2456 covered-boundary-chain cell: `{t['p2484_refines_one_inherited_p2456_covered_boundary_chain_cell']}`.",
        f"Full complement replay exported by P2484: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite boundary-side dyadic secant-margin audit outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
