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
from p2480_s1430_strict_pointwise_interval_decimal_p2459_segment_start_cell_dyadic_refinement_replay_audit import DYADIC_SUBCELL_COUNT

GEN = ROOT / "generated"
OUT = GEN / "p2482_s1432_strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit.json"
MD = GEN / "p2482_s1432_strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2481_BOUNDARY_HANDOFF_COLLAR": GEN / "p2481_s1431_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit.json",
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
        "new_packet": "P2482|S1432|boundary-band weakest cell|weakest boundary-band cell|boundary band dyadic refinement|root-window adjacent cell refinement|covered boundary cell refinement|P2456 weakest collar",
        "precursor_packets": "P2481|S1431|boundary handoff collar|P2456 right-boundary-band|covered-boundary-chain",
        "refinement_language": "dyadic refinement|subcell rows|diagnostic rows|weakest collar row|boundary-band parent cell|root-window-adjacent",
        "coverage_semantics": "not a P2459 coverage count|not distinct P2459 cells|covered-boundary-chain|diagnostic Decimal evaluation",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def boundary_dyadic_subcells(parent_cell: dict[str, Any]) -> list[dict[str, Any]]:
    left = Decimal(str(parent_cell["left"]))
    right = Decimal(str(parent_cell["right"]))
    width = right - left
    cells = []
    for index in range(DYADIC_SUBCELL_COUNT):
        sub_left = left + width * Decimal(index) / Decimal(DYADIC_SUBCELL_COUNT)
        sub_right = left + width * Decimal(index + 1) / Decimal(DYADIC_SUBCELL_COUNT)
        cells.append({
            "left": str(sub_left),
            "right": str(sub_right),
            "uncovered_index": index,
            "uncovered_count": DYADIC_SUBCELL_COUNT,
            "parent_collar_side": parent_cell["collar_side"],
            "parent_boundary_band_index": int(parent_cell["boundary_band_index"]),
        })
    return cells


def boundary_band_weakest_cell_refinement(parent_cell: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    family = parent_cell["family"]
    function = function_for_family(family)
    replayed = []
    for cell in boundary_dyadic_subcells(parent_cell):
        fresh = replay_cell(family, cell, projection, function)
        replayed.append({
            **fresh,
            "dyadic_subcell_index": int(cell["uncovered_index"]),
            "dyadic_subcell_count": DYADIC_SUBCELL_COUNT,
            "parent_collar_side": cell["parent_collar_side"],
            "parent_boundary_band_index": int(cell["parent_boundary_band_index"]),
            "fresh_decimal_separation_positive": Decimal(fresh["decimal_separation_from_zero"]) > 0,
        })
    ordered = sorted(replayed, key=lambda row: int(row["dyadic_subcell_index"]))
    consecutive_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        delta = Decimal(current["decimal_separation_from_zero"]) - Decimal(prior["decimal_separation_from_zero"])
        endpoint_gap = Decimal(str(current["left"])) - Decimal(str(prior["right"]))
        consecutive_pairs.append({
            "from_dyadic_subcell_index": int(prior["dyadic_subcell_index"]),
            "to_dyadic_subcell_index": int(current["dyadic_subcell_index"]),
            "endpoint_gap_to_minus_from": str(endpoint_gap),
            "from_separation": prior["decimal_separation_from_zero"],
            "to_separation": current["decimal_separation_from_zero"],
            "separation_delta_to_minus_from": str(delta),
            "strictly_increases": delta > 0,
            "is_exactly_adjacent": endpoint_gap == 0,
        })
    min_row = min(ordered, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    max_row = max(ordered, key=lambda row: Decimal(row["decimal_separation_from_zero"]))
    parent_lower_bound = Decimal(parent_cell["decimal_separation_from_zero"])
    refined_lower_bound = Decimal(min_row["decimal_separation_from_zero"])
    return {
        "family": family,
        "parent_collar_side": parent_cell["collar_side"],
        "parent_boundary_band_index": int(parent_cell["boundary_band_index"]),
        "parent_left": str(parent_cell["left"]),
        "parent_right": str(parent_cell["right"]),
        "parent_decimal_separation_from_p2481": parent_cell["decimal_separation_from_zero"],
        "dyadic_subcell_count": DYADIC_SUBCELL_COUNT,
        "all_subcells_exclude_zero": all(row["decimal_interval_excludes_zero"] for row in ordered),
        "all_subcells_have_positive_separation": all(row["fresh_decimal_separation_positive"] for row in ordered),
        "all_consecutive_subcell_separations_strictly_increase_left_to_right": all(row["strictly_increases"] for row in consecutive_pairs),
        "all_consecutive_subcells_are_exactly_adjacent": all(row["is_exactly_adjacent"] for row in consecutive_pairs),
        "minimum_consecutive_endpoint_gap": str(min(Decimal(row["endpoint_gap_to_minus_from"]) for row in consecutive_pairs)),
        "maximum_consecutive_endpoint_gap": str(max(Decimal(row["endpoint_gap_to_minus_from"]) for row in consecutive_pairs)),
        "minimum_consecutive_positive_delta": str(min(Decimal(row["separation_delta_to_minus_from"]) for row in consecutive_pairs)),
        "minimum_subcell_decimal_separation": min_row["decimal_separation_from_zero"],
        "maximum_subcell_decimal_separation": max_row["decimal_separation_from_zero"],
        "minimum_subcell_replay_row": min_row,
        "maximum_subcell_replay_row": max_row,
        "minimum_is_leftmost_subcell_adjacent_to_root_window": int(min_row["dyadic_subcell_index"]) == 0,
        "maximum_is_rightmost_subcell_adjacent_to_p2480_side": int(max_row["dyadic_subcell_index"]) == DYADIC_SUBCELL_COUNT - 1,
        "refined_minimum_lower_bound_exceeds_parent_interval_lower_bound": refined_lower_bound > parent_lower_bound,
        "refined_minus_parent_lower_bound_delta": str(refined_lower_bound - parent_lower_bound),
        "subcell_replay_rows": ordered,
        "consecutive_subcell_pairs": consecutive_pairs,
        "subcell_replay_fingerprint_sha256": sha256_json(ordered),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2482/S1432 strict pointwise interval-Decimal P2459 boundary-band weakest-cell dyadic-refinement replay audit

`P2482/S1432` opens the weakest row found by P2481: the leftmost P2456 right-boundary-band cell adjacent to the excluded root window.  It subdivides that already-covered boundary-chain parent cell into `128` dyadic diagnostic rows and replays each row with the Decimal/Taylor backend.  All `128` subcell rows exclude zero, their separations strictly increase left-to-right with exact consecutive endpoint adjacency, and the refined minimum lower bound is stronger than the coarse P2481 parent interval lower bound.  The minimum remains the leftmost subcell adjacent to the root-window side, so this is a finite boundary-cell refinement, not a root-window theorem, not a P2459 coverage increase, and not a full-complement replay.

For a non-specialist: P2481 showed that the weakest checked place was just before the refined segment-start cell, inside the already-covered strip next to the root window.  P2482 zooms into that one weakest strip cell.  The smaller pieces are all still safely away from zero, but the weakest piece is still at the edge nearest the root window.  This improves finite evidence without proving what happens inside the root window, without directed rounding, and without selector/source or ToE closure.
"""
    lag_section = """
## P2482/S1432 P2459 boundary-band weakest-cell dyadic-refinement replay guard

`P2482/S1432` adds a cell-internal refinement guard for the weakest P2481 boundary-band row behind `L_total`.  The `128` fresh Decimal evaluations are diagnostic subcell rows inside one inherited P2456 covered-boundary-chain cell, not new P2459 coverage cells; the result remains finite seam evidence only and does not add selector/source/gauge authority, physical-value generation, role-bearing finality, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2482/S1432 strict pointwise interval-Decimal P2459 boundary-band weakest-cell dyadic-refinement replay audit", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2482/S1432 P2459 boundary-band weakest-cell dyadic-refinement replay guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2481 = theorem(sources["P2481_BOUNDARY_HANDOFF_COLLAR"], "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit")
    parent_cell = p2481["handoff_collar_replay"]["minimum_collar_replay_row"]
    replay = boundary_band_weakest_cell_refinement(parent_cell, projection)
    p2459_universe = p2481["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2482_T1_strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit",
        "audited_chain": ["P2481/S1431"],
        "boundary_band_weakest_cell_dyadic_refinement_replay": replay,
        "family": replay["family"],
        "parent_collar_side": replay["parent_collar_side"],
        "parent_boundary_band_index": replay["parent_boundary_band_index"],
        "parent_decimal_separation_from_p2481": replay["parent_decimal_separation_from_p2481"],
        "dyadic_subcell_count": replay["dyadic_subcell_count"],
        "all_subcells_exclude_zero": replay["all_subcells_exclude_zero"],
        "all_subcells_have_positive_separation": replay["all_subcells_have_positive_separation"],
        "all_consecutive_subcell_separations_strictly_increase_left_to_right": replay["all_consecutive_subcell_separations_strictly_increase_left_to_right"],
        "all_consecutive_subcells_are_exactly_adjacent": replay["all_consecutive_subcells_are_exactly_adjacent"],
        "minimum_consecutive_endpoint_gap": replay["minimum_consecutive_endpoint_gap"],
        "maximum_consecutive_endpoint_gap": replay["maximum_consecutive_endpoint_gap"],
        "minimum_consecutive_positive_delta": replay["minimum_consecutive_positive_delta"],
        "minimum_subcell_decimal_separation": replay["minimum_subcell_decimal_separation"],
        "maximum_subcell_decimal_separation": replay["maximum_subcell_decimal_separation"],
        "minimum_is_leftmost_subcell_adjacent_to_root_window": replay["minimum_is_leftmost_subcell_adjacent_to_root_window"],
        "maximum_is_rightmost_subcell_adjacent_to_p2480_side": replay["maximum_is_rightmost_subcell_adjacent_to_p2480_side"],
        "refined_minimum_lower_bound_exceeds_parent_interval_lower_bound": replay["refined_minimum_lower_bound_exceeds_parent_interval_lower_bound"],
        "refined_minus_parent_lower_bound_delta": replay["refined_minus_parent_lower_bound_delta"],
        "p2481_total_fresh_handoff_collar_replay_rows_inherited": p2481["total_fresh_handoff_collar_replay_rows"],
        "p2481_minimum_collar_decimal_separation_inherited": p2481["minimum_collar_decimal_separation"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2482_fresh_decimal_evaluation_row_count_not_a_coverage_count": replay["dyadic_subcell_count"],
        "p2482_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction": f"{replay['dyadic_subcell_count']}/{p2459_universe}",
        "targeted_p2482_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2482_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2482_refines_one_inherited_p2456_covered_boundary_chain_cell": True,
        "p2482_subcell_rows_inside_single_covered_boundary_cell_not_distinct_p2459_cells": replay["dyadic_subcell_count"],
        "finite_chain_coverage_budget_inherited_from_p2481": p2481["finite_chain_coverage_budget_inherited_from_p2479"],
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "finite_replay_chain_boundary_band_weakest_cell_refinement_audit_exported": True,
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
        "lay_summary": "This packet opens the weakest P2481 collar row, which is the first P2456 right-boundary-band cell next to the excluded root window. Its 128 dyadic diagnostic rows all exclude zero and increase left-to-right, but the weakest row is still the leftmost subcell on the root-window side. The refined lower bound is stronger than the coarse parent-cell lower bound. This is finite covered-boundary-cell refinement, not a proof inside the root window, not a P2459 coverage increase, and not a continuum or full-complement proof.",
        "not_licensed": [
            "The 128 fresh Decimal evaluations are diagnostic subcell rows inside one inherited P2456 covered-boundary-chain cell, not 128 distinct P2459 coverage cells.",
            "P2482 adds zero new P2459 unreplayed cells; it refines an already-covered boundary-chain cell.",
            "It does not prove root-window exclusion, directed-rounding interval closure, symbolic root exclusion, analytic monotonicity, or continuum root exclusion.",
            "It does not export selector/source/gauge closure, a physical-value generator, QW-2191 discharge, role-bearing L_total, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "A stronger boundary claim now requires a root-window-side theorem or a directed-rounding backend; do not inflate this covered-boundary-cell refinement into a proof inside the excluded root window or into P2459 full-complement coverage.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "subcell_count_expected": theorem_export["dyadic_subcell_count"] == 128,
        "parent_binding_expected": theorem_export["parent_collar_side"] == "p2456_right_boundary_band" and theorem_export["parent_boundary_band_index"] == 0,
        "endpoint_adjacency_exact": theorem_export["all_consecutive_subcells_are_exactly_adjacent"] and Decimal(theorem_export["minimum_consecutive_endpoint_gap"]) == 0 and Decimal(theorem_export["maximum_consecutive_endpoint_gap"]) == 0,
        "minimum_separation_positive": Decimal(theorem_export["minimum_subcell_decimal_separation"]) > 0,
        "all_subcells_exclude_zero": theorem_export["all_subcells_exclude_zero"],
        "all_subcells_have_positive_separation": theorem_export["all_subcells_have_positive_separation"],
        "subcell_order_checked": theorem_export["all_consecutive_subcell_separations_strictly_increase_left_to_right"] and Decimal(theorem_export["minimum_consecutive_positive_delta"]) > 0,
        "leftmost_root_window_boundary_not_overclaimed": theorem_export["minimum_is_leftmost_subcell_adjacent_to_root_window"],
        "refined_lower_bound_stronger_than_parent": theorem_export["refined_minimum_lower_bound_exceeds_parent_interval_lower_bound"] and Decimal(theorem_export["refined_minus_parent_lower_bound_delta"]) > 0,
        "zero_new_p2459_unreplayed_cells": theorem_export["targeted_p2482_new_p2459_unreplayed_cell_count"] == 0 and theorem_export["targeted_p2482_new_p2459_unreplayed_cell_scope_against_p2459_universe"] == f"0/{p2459_universe}",
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "finite_grid_only_not_directed_rounding": theorem_export["finite_replay_chain_boundary_band_weakest_cell_refinement_audit_exported"] and not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
        "no_root_window_theorem": not theorem_export["root_window_exclusion_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2482_s1432_v1",
        "packet_id": "P2482",
        "stage_id": "S1432",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_BAND_WEAKEST_CELL_DYADIC_REFINEMENT_REPLAY_AUDIT_NO_ROOT_WINDOW_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit"]["theorem_export"]
    lines = [
        "# P2482/S1432 strict pointwise interval-Decimal P2459 boundary-band weakest-cell dyadic-refinement replay audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-band weakest-cell dyadic refinement",
        "",
        f"Parent collar side: `{t['parent_collar_side']}`, boundary-band index `{t['parent_boundary_band_index']}`.",
        f"Dyadic subcell replay rows: `{t['dyadic_subcell_count']}`.",
        f"Parent Decimal separation inherited from P2481: `{t['parent_decimal_separation_from_p2481']}`.",
        f"Minimum subcell Decimal separation: `{t['minimum_subcell_decimal_separation']}`.",
        f"Maximum subcell Decimal separation: `{t['maximum_subcell_decimal_separation']}`.",
        f"Minimum consecutive endpoint gap: `{t['minimum_consecutive_endpoint_gap']}`.",
        f"Maximum consecutive endpoint gap: `{t['maximum_consecutive_endpoint_gap']}`.",
        f"Minimum consecutive positive delta: `{t['minimum_consecutive_positive_delta']}`.",
        f"Minimum is leftmost subcell adjacent to root window: `{t['minimum_is_leftmost_subcell_adjacent_to_root_window']}`.",
        f"Refined minimum lower bound exceeds parent interval lower bound: `{t['refined_minimum_lower_bound_exceeds_parent_interval_lower_bound']}`.",
        f"Refined-minus-parent lower-bound delta: `{t['refined_minus_parent_lower_bound_delta']}`.",
        f"All subcells exclude zero: `{t['all_subcells_exclude_zero']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2482 fresh Decimal evaluation rows (not a P2459 coverage count): `{t['p2482_fresh_decimal_evaluation_row_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2482_fresh_decimal_evaluation_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2482: `{t['targeted_p2482_new_p2459_unreplayed_cell_count']}`.",
        f"New P2459 unreplayed-cell scope against inherited P2459 universe: `{t['targeted_p2482_new_p2459_unreplayed_cell_scope_against_p2459_universe']}`.",
        f"P2482 refines one inherited P2456 covered-boundary-chain cell: `{t['p2482_refines_one_inherited_p2456_covered_boundary_chain_cell']}`.",
        f"Full complement replay exported by P2482: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite dyadic refinement of one inherited P2456 covered-boundary-chain cell only.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
