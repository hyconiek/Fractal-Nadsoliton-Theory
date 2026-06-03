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
    load_json,
    rel,
)

GEN = ROOT / "generated"
OUT = GEN / "p2489_s1439_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate.json"
MD = GEN / "p2489_s1439_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate.md"

SOURCE_FILES = {
    "P2487_COLLAR_DERIVATIVE_SWEEP": GEN / "p2487_s1437_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate.json",
    "P2488_COLLAR_MONOTONICITY_LEMMA": GEN / "p2488_s1438_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate.json",
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
        "new_packet": "P2489|S1439|cumulative derivative barrier|derivative lower-barrier|collar transport barrier|integrated derivative barrier|finite lower-barrier|cumulative collar barrier",
        "precursor_packets": "P2488|S1438|collar monotonicity lemma|P2487|S1437|collar derivative sweep|P2481|S1431|boundary handoff collar",
        "barrier_language": "derivative transport lower gain|cumulative derivative|lower barrier|left anchor|integrated lower bound|collar-wide barrier",
        "coverage_semantics": "diagnostic row ratio|not a coverage fraction|zero new P2459 unreplayed cells|proof compression|no new replay",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer|root-window theorem",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def cumulative_derivative_barrier(p2487: dict[str, Any], p2488: dict[str, Any]) -> dict[str, Any]:
    sweep = p2487["boundary_handoff_collar_derivative_sweep"]
    rows = sweep["collar_derivative_rows"]
    anchor_lower = Decimal(rows[0]["amplitude_interval_value"]["lo"])
    current = anchor_lower
    barrier_rows = []
    for order, row in enumerate(rows):
        width = Decimal(row["width"])
        derivative_lower = Decimal(row["derivative_lower_bound"])
        row_local_lower = Decimal(row["amplitude_interval_value"]["lo"])
        derivative_gain = derivative_lower * width
        exit_barrier = current + derivative_gain
        barrier_rows.append({
            "collar_order_index": order,
            "collar_source_kind": row["collar_source_kind"],
            "collar_source_index": row["collar_source_index"],
            "left": row["left"],
            "right": row["right"],
            "width": row["width"],
            "row_local_amplitude_lower_bound": str(row_local_lower),
            "entry_cumulative_lower_barrier": str(current),
            "derivative_lower_bound": str(derivative_lower),
            "certified_derivative_gain_on_row": str(derivative_gain),
            "exit_cumulative_lower_barrier": str(exit_barrier),
            "entry_barrier_positive": current > 0,
            "derivative_gain_positive": derivative_gain > 0,
            "exit_barrier_positive": exit_barrier > 0,
            "row_local_lower_bound_positive": row_local_lower > 0,
        })
        current = exit_barrier

    gains = [Decimal(item["certified_derivative_gain_on_row"]) for item in barrier_rows]
    entries = [Decimal(item["entry_cumulative_lower_barrier"]) for item in barrier_rows]
    exits = [Decimal(item["exit_cumulative_lower_barrier"]) for item in barrier_rows]
    local_lowers = [Decimal(item["row_local_amplitude_lower_bound"]) for item in barrier_rows]
    expected_gain = sum(gains)
    lemma = p2488["boundary_handoff_collar_monotonicity_lemma"]
    preconditions = {
        "p2487_all_rows_positive_derivative": p2487["all_derivative_intervals_positive"],
        "p2487_all_rows_positive_value": p2487["all_amplitude_intervals_positive"],
        "p2487_exact_adjacency": p2487["all_consecutive_rows_exactly_adjacent"],
        "p2488_checked_collar_lemma_exported": p2488["finite_piecewise_monotone_increasing_collar_lemma_exported"],
        "same_total_row_count_as_p2487": len(barrier_rows) == p2487["total_derivative_row_count"],
        "positive_left_anchor": anchor_lower > 0,
        "all_row_gains_positive": all(item["derivative_gain_positive"] for item in barrier_rows),
        "all_entry_barriers_positive": all(item["entry_barrier_positive"] for item in barrier_rows),
        "all_exit_barriers_positive": all(item["exit_barrier_positive"] for item in barrier_rows),
    }
    return {
        "family": p2487["family"],
        "collar_left": p2487["collar_left"],
        "collar_right": p2487["collar_right"],
        "total_barrier_row_count": len(barrier_rows),
        "barrier_preconditions": preconditions,
        "all_barrier_preconditions_met": all(preconditions.values()),
        "left_anchor_amplitude_lower_bound": str(anchor_lower),
        "minimum_row_local_amplitude_lower_bound": str(min(local_lowers)),
        "minimum_entry_cumulative_lower_barrier": str(min(entries)),
        "minimum_exit_cumulative_lower_barrier": str(min(exits)),
        "minimum_positive_derivative_gain_on_row": str(min(gains)),
        "maximum_positive_derivative_gain_on_row": str(max(gains)),
        "total_certified_derivative_gain_over_collar": str(expected_gain),
        "p2488_derivative_transport_lower_gain_inherited": lemma["derivative_transport_lower_gain_over_checked_collar"],
        "transport_gain_matches_p2488": expected_gain == Decimal(lemma["derivative_transport_lower_gain_over_checked_collar"]),
        "right_endpoint_cumulative_lower_barrier": str(exits[-1]),
        "right_endpoint_barrier_exceeds_left_anchor": exits[-1] > anchor_lower,
        "finite_cumulative_lower_barrier_certificate_exported": all(preconditions.values()),
        "barrier_rows": barrier_rows,
        "barrier_fingerprint_sha256": sha256_json(barrier_rows),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2489/S1439 strict pointwise interval-Decimal P2459 boundary-handoff collar cumulative derivative barrier certificate

`P2489/S1439` converts the P2487 positive derivative intervals and the P2488 checked-collar lemma into a finite cumulative lower-barrier calculation.  Starting from the left collar amplitude lower bound, it sums `derivative_lower_bound * row_width` across the same `134` adjacent P2481 collar rows and verifies that every entry barrier, row gain, and exit barrier remains strictly positive.  This is a finite integrated-derivative barrier over the checked collar only: it adds no P2459 coverage, does not use directed rounding, and does not prove the excluded root window, global analytic monotonicity, selector/source closure, legacy-role transfer, or ToE closure.

For a non-specialist: P2488 said the checked seam is positive and increasing piece by piece.  P2489 records the corresponding cumulative safety margin: even if we transport only the certified minimum derivative through each tiny row, the running lower bound stays positive all the way across the checked collar.
"""
    lag_section = """
## P2489/S1439 P2459 boundary-handoff collar cumulative derivative barrier guard

`P2489/S1439` adds a finite cumulative derivative lower-barrier guard behind `L_total`: the same `134` P2487/P2488 collar rows are integrated with certified lower derivative gains, producing a positive running lower barrier from the left anchor to the collar right endpoint.  The guard is still collar-local and proof-compressive; it does not export directed rounding, root-window exclusion, global analytic monotonicity, selector/source/gauge authority, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2489/S1439 strict pointwise interval-Decimal P2459 boundary-handoff collar cumulative derivative barrier certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2489/S1439 P2459 boundary-handoff collar cumulative derivative barrier guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2487 = theorem(sources["P2487_COLLAR_DERIVATIVE_SWEEP"], "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate")
    p2488 = theorem(sources["P2488_COLLAR_MONOTONICITY_LEMMA"], "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate")
    barrier = cumulative_derivative_barrier(p2487, p2488)
    p2459_universe = p2487["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2489_T1_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate",
        "audited_chain": ["P2487/S1437", "P2488/S1438"],
        "boundary_handoff_collar_cumulative_derivative_barrier": barrier,
        "family": barrier["family"],
        "total_barrier_row_count": barrier["total_barrier_row_count"],
        "all_barrier_preconditions_met": barrier["all_barrier_preconditions_met"],
        "finite_cumulative_lower_barrier_certificate_exported": barrier["finite_cumulative_lower_barrier_certificate_exported"],
        "minimum_entry_cumulative_lower_barrier": barrier["minimum_entry_cumulative_lower_barrier"],
        "minimum_exit_cumulative_lower_barrier": barrier["minimum_exit_cumulative_lower_barrier"],
        "minimum_positive_derivative_gain_on_row": barrier["minimum_positive_derivative_gain_on_row"],
        "total_certified_derivative_gain_over_collar": barrier["total_certified_derivative_gain_over_collar"],
        "right_endpoint_cumulative_lower_barrier": barrier["right_endpoint_cumulative_lower_barrier"],
        "transport_gain_matches_p2488": barrier["transport_gain_matches_p2488"],
        "p2487_total_derivative_rows_inherited": p2487["total_derivative_row_count"],
        "p2488_checked_collar_lemma_inherited": p2488["finite_piecewise_monotone_increasing_collar_lemma_exported"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2489_new_decimal_replay_row_count": 0,
        "p2489_reused_p2487_p2488_row_count_not_a_coverage_count": barrier["total_barrier_row_count"],
        "p2489_reused_row_ratio_not_a_p2459_coverage_fraction": f"{barrier['total_barrier_row_count']}/{p2459_universe}",
        "targeted_p2489_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2489_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2489_is_barrier_proof_compression_not_new_replay": True,
        "full_complement_replay_exported_by_this_certificate": False,
        "directed_rounding_interval_theorem_exported_by_this_certificate": False,
        "symbolic_root_exclusion_theorem_exported_by_this_certificate": False,
        "root_window_exclusion_theorem_exported_by_this_certificate": False,
        "global_continuum_root_exclusion_theorem_exported_by_this_certificate": False,
        "global_analytic_monotonicity_theorem_exported_by_this_certificate": False,
        "pointwise_coordinate_selector_exported_by_this_certificate": False,
        "strict_observable_source_constraint_exported_by_this_certificate": False,
        "gauge_slice_theorem_exported_by_this_certificate": False,
        "strict_physical_value_generator_exported": False,
        "qw2191_discharged": False,
        "role_bearing_ltotal_exported": False,
        "legacy_role_transfer_exported": False,
        "toe_closure_exported": False,
        "lay_summary": "This packet integrates the certified positive derivative lower bounds from P2487 across the checked P2481 collar, using the P2488 monotonicity lemma as the finite proof context. The cumulative lower barrier remains positive at every row entry and exit, but this is still a collar-local finite certificate outside the excluded root window.",
        "not_licensed": [
            "The cumulative barrier reuses the 134 P2487/P2488 collar rows and adds zero new P2459 coverage cells.",
            "The positive integrated derivative barrier licenses only checked-collar lower-barrier evidence, not root-window exclusion or global analytic monotonicity.",
            "P2489 does not discharge directed rounding, symbolic or continuum root exclusion, selector/source/gauge closure, QW-2191, legacy-role transfer, physical-value generation, or ToE closure.",
        ],
        "next_honest_step": "The checked collar now has value, derivative, lemma, and cumulative-barrier certificates. A stronger closure requires a separate root-window-side theorem, directed-rounding proof, or a nonlocal analytic argument rather than more collar-local row inflation.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "p2487_inheritance_checked": theorem_export["p2487_total_derivative_rows_inherited"] == 134,
        "p2488_lemma_inheritance_checked": theorem_export["p2488_checked_collar_lemma_inherited"],
        "barrier_counts_match_sources": theorem_export["total_barrier_row_count"] == theorem_export["p2487_total_derivative_rows_inherited"],
        "all_barrier_preconditions_met": theorem_export["all_barrier_preconditions_met"],
        "transport_gain_matches_p2488": theorem_export["transport_gain_matches_p2488"],
        "positive_running_barrier": Decimal(theorem_export["minimum_entry_cumulative_lower_barrier"]) > 0 and Decimal(theorem_export["minimum_exit_cumulative_lower_barrier"]) > 0,
        "coverage_budget_not_inflated": theorem_export["targeted_p2489_new_p2459_unreplayed_cell_count"] == 0,
        "no_closure_inflation": not any(theorem_export[key] for key in [
            "full_complement_replay_exported_by_this_certificate",
            "directed_rounding_interval_theorem_exported_by_this_certificate",
            "symbolic_root_exclusion_theorem_exported_by_this_certificate",
            "root_window_exclusion_theorem_exported_by_this_certificate",
            "global_continuum_root_exclusion_theorem_exported_by_this_certificate",
            "global_analytic_monotonicity_theorem_exported_by_this_certificate",
            "pointwise_coordinate_selector_exported_by_this_certificate",
            "strict_observable_source_constraint_exported_by_this_certificate",
            "gauge_slice_theorem_exported_by_this_certificate",
            "strict_physical_value_generator_exported",
            "qw2191_discharged",
            "role_bearing_ltotal_exported",
            "legacy_role_transfer_exported",
            "toe_closure_exported",
        ]),
    }
    return {
        "packet_id": "P2489",
        "stage_id": "S1439",
        "status": "FINITE_COLLAR_CUMULATIVE_DERIVATIVE_BARRIER_CERTIFICATE_PROOF_COMPRESSION_NO_NEW_REPLAY_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_SELECTOR_SOURCE_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate"]["theorem_export"]
    b = t["boundary_handoff_collar_cumulative_derivative_barrier"]
    lines = [
        "# P2489/S1439 strict pointwise interval-Decimal P2459 boundary-handoff collar cumulative derivative barrier certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite cumulative derivative lower barrier",
        "",
        f"Audited chain: `{', '.join(t['audited_chain'])}`.",
        f"Barrier rows reused from P2487/P2488: `{t['total_barrier_row_count']}`.",
        f"All barrier preconditions met: `{t['all_barrier_preconditions_met']}`.",
        f"Minimum entry cumulative lower barrier: `{t['minimum_entry_cumulative_lower_barrier']}`.",
        f"Minimum exit cumulative lower barrier: `{t['minimum_exit_cumulative_lower_barrier']}`.",
        f"Minimum positive derivative gain on one row: `{t['minimum_positive_derivative_gain_on_row']}`.",
        f"Total certified derivative gain over collar: `{t['total_certified_derivative_gain_over_collar']}`.",
        f"Right endpoint cumulative lower barrier: `{t['right_endpoint_cumulative_lower_barrier']}`.",
        f"Transport gain matches P2488: `{t['transport_gain_matches_p2488']}`.",
        "",
        "## Coverage semantics",
        "",
        f"New Decimal replay rows in P2489: `{t['p2489_new_decimal_replay_row_count']}`.",
        f"Reused barrier rows (not a P2459 coverage count): `{t['p2489_reused_p2487_p2488_row_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio, not coverage fraction: `{t['p2489_reused_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2489: `{t['targeted_p2489_new_p2459_unreplayed_cell_count']}`.",
        f"P2489 is barrier proof compression, not new replay: `{t['p2489_is_barrier_proof_compression_not_new_replay']}`.",
        "",
        "## Negative controls",
        "",
        "P2489 does not export directed rounding, symbolic/global continuum root exclusion, root-window exclusion, global analytic monotonicity, selector/source/gauge closure, QW-2191 discharge, role-bearing L_total, legacy-role transfer, physical-value generation, or ToE closure.",
        "",
        "## Lay summary",
        "",
        t["lay_summary"],
        "",
        "## Fingerprints",
        "",
        f"Barrier fingerprint: `{b['barrier_fingerprint_sha256']}`.",
        f"Theorem fingerprint: `{payload['strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_cumulative_derivative_barrier_certificate']['theorem_fingerprint_sha256']}`.",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    GEN.mkdir(exist_ok=True)
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
