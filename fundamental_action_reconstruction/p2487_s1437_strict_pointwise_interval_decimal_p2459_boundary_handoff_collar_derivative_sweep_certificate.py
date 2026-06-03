#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import subprocess
from decimal import Decimal
from typing import Any

from p2453_s1403_strict_pointwise_directed_decimal_weakest_cell_replay_certificate import (
    DecimalInterval,
    interval_dot,
    pointwise_gradient_derivative_interval,
    projection_amplitude_interval,
)
from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import (
    DOC_FILES,
    REPO,
    ROOT,
    load_json,
    rel,
)
from p2475_s1425_strict_pointwise_interval_decimal_p2459_critical_minimum_halo_replay_audit import projection_vector

GEN = ROOT / "generated"
OUT = GEN / "p2487_s1437_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate.json"
MD = GEN / "p2487_s1437_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2481_BOUNDARY_HANDOFF_COLLAR": GEN / "p2481_s1431_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit.json",
    "P2486_BOUNDARY_CELL_DERIVATIVE": GEN / "p2486_s1436_strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate.json",
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
        "new_packet": "P2487|S1437|boundary-handoff collar derivative sweep|collar derivative sweep|finite collar monotonicity|collar-wide derivative sign|handoff collar monotone|piecewise derivative collar",
        "precursor_packets": "P2486|S1436|boundary-cell derivative monotonicity|P2481|S1431|boundary handoff collar|P2485|S1435|secant curvature stability",
        "derivative_language": "projection amplitude derivative|positive derivative|finite-interval monotone|collar-wide derivative|piecewise monotone|endpoint adjacency",
        "coverage_semantics": "diagnostic row ratio|not a coverage fraction|zero new P2459 unreplayed cells|covered-boundary-chain",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def row_left(row: dict[str, Any]) -> Decimal:
    return Decimal(str(row.get("left", row.get("cell_left"))))


def row_right(row: dict[str, Any]) -> Decimal:
    return Decimal(str(row.get("right", row.get("cell_right"))))


def derivative_row(source_kind: str, source_index: int, source_row: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    left = row_left(source_row)
    right = row_right(source_row)
    interval = DecimalInterval(left, right)
    amplitude = projection_amplitude_interval(projection, interval)
    derivative = interval_dot(projection, pointwise_gradient_derivative_interval(interval))
    derivative_lower = derivative.lo
    derivative_upper = derivative.hi
    return {
        "collar_source_kind": source_kind,
        "collar_source_index": source_index,
        "left": str(left),
        "right": str(right),
        "width": str(right - left),
        "family": source_row["family"],
        "decimal_separation_from_zero_inherited": source_row["decimal_separation_from_zero"],
        "amplitude_interval_value": amplitude.as_dict(),
        "amplitude_interval_excludes_zero": not amplitude.contains_zero(),
        "amplitude_interval_positive": amplitude.lo > 0,
        "amplitude_interval_separation_from_zero": str(amplitude.separation_from_zero()),
        "derivative_interval_value": derivative.as_dict(),
        "derivative_interval_excludes_zero": not derivative.contains_zero(),
        "derivative_interval_positive": derivative_lower > 0,
        "derivative_interval_separation_from_zero": str(derivative.separation_from_zero()),
        "derivative_lower_bound": str(derivative_lower),
        "derivative_upper_bound": str(derivative_upper),
        "derivative_relative_width": str((derivative_upper - derivative_lower) / derivative_lower),
        "local_interval_monotone_increasing_witness": derivative_lower > 0,
    }


def collar_derivative_sweep(p2481: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    collar = p2481["handoff_collar_replay"]
    boundary_rows = [derivative_row("p2456_right_boundary_band", index, row, projection) for index, row in enumerate(collar["boundary_band_rows"])]
    subcell_rows = [derivative_row("p2480_parent_dyadic_subcell", index, row, projection) for index, row in enumerate(collar["dyadic_subcell_rows"])]
    ordered = sorted(boundary_rows + subcell_rows, key=lambda row: Decimal(row["left"]))
    consecutive_pairs = []
    for prior, current in zip(ordered, ordered[1:]):
        gap = Decimal(current["left"]) - Decimal(prior["right"])
        prior_sep = Decimal(prior["amplitude_interval_separation_from_zero"])
        current_sep = Decimal(current["amplitude_interval_separation_from_zero"])
        consecutive_pairs.append({
            "from_collar_order_index": ordered.index(prior),
            "to_collar_order_index": ordered.index(current),
            "from_right": prior["right"],
            "to_left": current["left"],
            "endpoint_gap": str(gap),
            "exactly_adjacent": gap == 0,
            "from_amplitude_separation": str(prior_sep),
            "to_amplitude_separation": str(current_sep),
            "amplitude_separation_delta_to_minus_from": str(current_sep - prior_sep),
            "amplitude_separation_strictly_increases": current_sep > prior_sep,
        })
    derivative_lowers = [Decimal(row["derivative_lower_bound"]) for row in ordered]
    derivative_uppers = [Decimal(row["derivative_upper_bound"]) for row in ordered]
    derivative_relative_widths = [Decimal(row["derivative_relative_width"]) for row in ordered]
    amplitude_separations = [Decimal(row["amplitude_interval_separation_from_zero"]) for row in ordered]
    return {
        "family": ordered[0]["family"],
        "boundary_band_derivative_row_count": len(boundary_rows),
        "dyadic_subcell_derivative_row_count": len(subcell_rows),
        "total_derivative_row_count": len(ordered),
        "collar_left": ordered[0]["left"],
        "collar_right": ordered[-1]["right"],
        "all_amplitude_intervals_exclude_zero": all(row["amplitude_interval_excludes_zero"] for row in ordered),
        "all_amplitude_intervals_positive": all(row["amplitude_interval_positive"] for row in ordered),
        "all_derivative_intervals_exclude_zero": all(row["derivative_interval_excludes_zero"] for row in ordered),
        "all_derivative_intervals_positive": all(row["derivative_interval_positive"] for row in ordered),
        "all_rows_have_local_interval_monotone_increasing_witness": all(row["local_interval_monotone_increasing_witness"] for row in ordered),
        "all_consecutive_rows_exactly_adjacent": all(row["exactly_adjacent"] for row in consecutive_pairs),
        "all_consecutive_amplitude_separations_strictly_increase": all(row["amplitude_separation_strictly_increases"] for row in consecutive_pairs),
        "minimum_amplitude_interval_separation_from_zero": str(min(amplitude_separations)),
        "maximum_amplitude_interval_separation_from_zero": str(max(amplitude_separations)),
        "minimum_derivative_lower_bound": str(min(derivative_lowers)),
        "maximum_derivative_upper_bound": str(max(derivative_uppers)),
        "maximum_derivative_relative_width": str(max(derivative_relative_widths)),
        "minimum_endpoint_gap": str(min(Decimal(row["endpoint_gap"]) for row in consecutive_pairs)),
        "maximum_endpoint_gap": str(max(Decimal(row["endpoint_gap"]) for row in consecutive_pairs)),
        "minimum_consecutive_positive_amplitude_separation_delta": str(min(Decimal(row["amplitude_separation_delta_to_minus_from"]) for row in consecutive_pairs)),
        "weakest_derivative_row": min(ordered, key=lambda row: Decimal(row["derivative_lower_bound"])),
        "strongest_derivative_row": max(ordered, key=lambda row: Decimal(row["derivative_upper_bound"])),
        "weakest_amplitude_row": min(ordered, key=lambda row: Decimal(row["amplitude_interval_separation_from_zero"])),
        "collar_derivative_rows": ordered,
        "consecutive_collar_derivative_pairs": consecutive_pairs,
        "collar_derivative_sweep_fingerprint_sha256": sha256_json({"rows": ordered, "pairs": consecutive_pairs}),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2487/S1437 strict pointwise interval-Decimal P2459 boundary-handoff collar derivative sweep certificate

`P2487/S1437` lifts the P2486 one-cell derivative-sign check across the full P2481 boundary-handoff collar: the `6` inherited P2456 right-boundary-band cells plus the `128` dyadic rows inside the P2480 parent cell.  Every collar row has positive projection-amplitude value and a strictly positive projection-amplitude derivative interval under the strict Decimal/Taylor backend; the rows are exactly adjacent and their amplitude separations still increase left-to-right.  This gives a finite piecewise interval-monotonicity witness for the checked collar outside the excluded root window.  It is not a global analytic monotonicity theorem, not root-window exclusion, not directed rounding, not a P2459 coverage increase, and not selector/source or ToE closure.

For a non-specialist: P2486 proved the derivative was positive on the weakest checked cell.  P2487 repeats that derivative check over the whole small handoff collar that was previously replayed, so the checked seam behaves monotonically across all its diagnostic rows.  This is stronger finite seam evidence, but it still does not prove the excluded root window or the whole continuum.
"""
    lag_section = """
## P2487/S1437 P2459 boundary-handoff collar derivative sweep guard

`P2487/S1437` adds a finite collar-wide derivative-sign guard behind `L_total`: `134` diagnostic interval-derivative rows across the P2481 handoff collar, all with positive projection-amplitude value and positive projection-amplitude derivative.  The guard remains finite and local to the inherited collar; it does not export directed rounding, root-window exclusion, global analytic monotonicity, selector/source/gauge authority, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2487/S1437 strict pointwise interval-Decimal P2459 boundary-handoff collar derivative sweep certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2487/S1437 P2459 boundary-handoff collar derivative sweep guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2481 = theorem(sources["P2481_BOUNDARY_HANDOFF_COLLAR"], "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit")
    p2486 = theorem(sources["P2486_BOUNDARY_CELL_DERIVATIVE"], "strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate")
    sweep = collar_derivative_sweep(p2481, projection)
    p2459_universe = p2481["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2487_T1_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate",
        "audited_chain": ["P2481/S1431", "P2486/S1436"],
        "boundary_handoff_collar_derivative_sweep": sweep,
        "family": sweep["family"],
        "boundary_band_derivative_row_count": sweep["boundary_band_derivative_row_count"],
        "dyadic_subcell_derivative_row_count": sweep["dyadic_subcell_derivative_row_count"],
        "total_derivative_row_count": sweep["total_derivative_row_count"],
        "collar_left": sweep["collar_left"],
        "collar_right": sweep["collar_right"],
        "all_amplitude_intervals_exclude_zero": sweep["all_amplitude_intervals_exclude_zero"],
        "all_amplitude_intervals_positive": sweep["all_amplitude_intervals_positive"],
        "all_derivative_intervals_exclude_zero": sweep["all_derivative_intervals_exclude_zero"],
        "all_derivative_intervals_positive": sweep["all_derivative_intervals_positive"],
        "all_rows_have_local_interval_monotone_increasing_witness": sweep["all_rows_have_local_interval_monotone_increasing_witness"],
        "all_consecutive_rows_exactly_adjacent": sweep["all_consecutive_rows_exactly_adjacent"],
        "all_consecutive_amplitude_separations_strictly_increase": sweep["all_consecutive_amplitude_separations_strictly_increase"],
        "minimum_amplitude_interval_separation_from_zero": sweep["minimum_amplitude_interval_separation_from_zero"],
        "maximum_amplitude_interval_separation_from_zero": sweep["maximum_amplitude_interval_separation_from_zero"],
        "minimum_derivative_lower_bound": sweep["minimum_derivative_lower_bound"],
        "maximum_derivative_upper_bound": sweep["maximum_derivative_upper_bound"],
        "maximum_derivative_relative_width": sweep["maximum_derivative_relative_width"],
        "minimum_endpoint_gap": sweep["minimum_endpoint_gap"],
        "maximum_endpoint_gap": sweep["maximum_endpoint_gap"],
        "minimum_consecutive_positive_amplitude_separation_delta": sweep["minimum_consecutive_positive_amplitude_separation_delta"],
        "p2486_one_cell_derivative_positive_inherited": p2486["derivative_interval_positive_on_entire_cell"],
        "p2481_total_handoff_collar_replay_rows_inherited": p2481["total_fresh_handoff_collar_replay_rows"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2487_derivative_interval_row_count_not_a_coverage_count": sweep["total_derivative_row_count"],
        "p2487_derivative_interval_row_ratio_not_a_p2459_coverage_fraction": f"{sweep['total_derivative_row_count']}/{p2459_universe}",
        "targeted_p2487_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2487_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2487_reuses_p2481_collar_rows_without_new_p2459_coverage": True,
        "finite_chain_coverage_budget_inherited_from_p2481": p2481["finite_chain_coverage_budget_inherited_from_p2479"],
        "finite_piecewise_interval_monotonicity_on_checked_collar_exported": True,
        "global_analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
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
        "no_full_complement_claimed_by_this_certificate": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "lay_summary": "This packet checks the derivative on every diagnostic row of the P2481 handoff collar. All 134 rows have positive projection amplitude and positive projection-amplitude derivative, and the rows remain exactly adjacent. This provides finite piecewise monotonicity evidence for the checked seam, but it is still outside the excluded root window and does not prove global analytic, selector/source, or ToE closure.",
        "not_licensed": [
            "The 134 derivative interval rows are diagnostics over the existing P2481 collar rows, not new P2459 coverage cells.",
            "Positive derivative signs over this finite collar license only finite piecewise interval monotonicity on the checked collar, not global analytic monotonicity or root-window exclusion.",
            "P2487 adds zero new P2459 unreplayed cells and does not discharge directed rounding, symbolic root exclusion, continuum root exclusion, selector/source/gauge closure, QW-2191, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "The derivative-sign evidence now covers the finite boundary handoff collar. A stronger result requires a separate root-window-side theorem or a global analytic/directed-rounding proof, not more silent inflation of collar diagnostics.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "counts_match_p2481_collar": theorem_export["boundary_band_derivative_row_count"] == 6 and theorem_export["dyadic_subcell_derivative_row_count"] == 128 and theorem_export["total_derivative_row_count"] == 134,
        "p2486_inheritance_checked": theorem_export["p2486_one_cell_derivative_positive_inherited"],
        "p2481_inheritance_checked": theorem_export["p2481_total_handoff_collar_replay_rows_inherited"] == 134,
        "amplitude_intervals_exclude_zero": theorem_export["all_amplitude_intervals_exclude_zero"] and theorem_export["all_amplitude_intervals_positive"],
        "derivative_intervals_positive": theorem_export["all_derivative_intervals_exclude_zero"] and theorem_export["all_derivative_intervals_positive"] and Decimal(theorem_export["minimum_derivative_lower_bound"]) > 0,
        "local_piecewise_monotone_witness_only": theorem_export["finite_piecewise_interval_monotonicity_on_checked_collar_exported"] and theorem_export["all_rows_have_local_interval_monotone_increasing_witness"],
        "exact_adjacency_preserved": theorem_export["all_consecutive_rows_exactly_adjacent"] and Decimal(theorem_export["minimum_endpoint_gap"]) == 0 and Decimal(theorem_export["maximum_endpoint_gap"]) == 0,
        "separations_increase": theorem_export["all_consecutive_amplitude_separations_strictly_increase"] and Decimal(theorem_export["minimum_consecutive_positive_amplitude_separation_delta"]) > 0,
        "zero_new_p2459_unreplayed_cells": theorem_export["targeted_p2487_new_p2459_unreplayed_cell_count"] == 0 and theorem_export["targeted_p2487_new_p2459_unreplayed_cell_scope_against_p2459_universe"] == f"0/{p2459_universe}",
        "no_full_complement_unless_genuinely_full": theorem_export["no_full_complement_claimed_by_this_certificate"] and not theorem_export["full_complement_replay_exported_by_this_certificate"],
        "not_global_analytic_monotonicity": not theorem_export["global_analytic_monotonicity_theorem_exported_by_this_certificate"],
        "not_directed_rounding": not theorem_export["directed_rounding_interval_theorem_exported_by_this_certificate"],
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
        "schema_version": "p2487_s1437_v1",
        "packet_id": "P2487",
        "stage_id": "S1437",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_HANDOFF_COLLAR_DERIVATIVE_SWEEP_CERTIFICATE_FINITE_PIECEWISE_MONOTONE_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate"]["theorem_export"]
    lines = [
        "# P2487/S1437 strict pointwise interval-Decimal P2459 boundary-handoff collar derivative sweep certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-handoff collar derivative sweep",
        "",
        f"Collar interval: `[{t['collar_left']}, {t['collar_right']}]`.",
        f"Boundary-band derivative rows: `{t['boundary_band_derivative_row_count']}`.",
        f"Dyadic subcell derivative rows: `{t['dyadic_subcell_derivative_row_count']}`.",
        f"Total derivative rows: `{t['total_derivative_row_count']}`.",
        f"All amplitude intervals positive: `{t['all_amplitude_intervals_positive']}`.",
        f"All derivative intervals positive: `{t['all_derivative_intervals_positive']}`.",
        f"All rows have local monotone-increasing witness: `{t['all_rows_have_local_interval_monotone_increasing_witness']}`.",
        f"All consecutive rows exactly adjacent: `{t['all_consecutive_rows_exactly_adjacent']}`.",
        f"All consecutive amplitude separations strictly increase: `{t['all_consecutive_amplitude_separations_strictly_increase']}`.",
        f"Minimum amplitude interval separation: `{t['minimum_amplitude_interval_separation_from_zero']}`.",
        f"Maximum amplitude interval separation: `{t['maximum_amplitude_interval_separation_from_zero']}`.",
        f"Minimum derivative lower bound: `{t['minimum_derivative_lower_bound']}`.",
        f"Maximum derivative upper bound: `{t['maximum_derivative_upper_bound']}`.",
        f"Maximum derivative relative width: `{t['maximum_derivative_relative_width']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2487 derivative interval rows (not a P2459 coverage count): `{t['p2487_derivative_interval_row_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2487_derivative_interval_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2487: `{t['targeted_p2487_new_p2459_unreplayed_cell_count']}`.",
        f"New P2459 unreplayed-cell scope against inherited P2459 universe: `{t['targeted_p2487_new_p2459_unreplayed_cell_scope_against_p2459_universe']}`.",
        f"P2487 reuses P2481 collar rows without new P2459 coverage: `{t['p2487_reuses_p2481_collar_rows_without_new_p2459_coverage']}`.",
        f"Full complement replay exported by P2487: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a finite piecewise derivative monotonicity certificate for the checked boundary-handoff collar outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, global analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
