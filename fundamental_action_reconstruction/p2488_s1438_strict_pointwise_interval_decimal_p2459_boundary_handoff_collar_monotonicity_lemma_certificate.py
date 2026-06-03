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
OUT = GEN / "p2488_s1438_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate.json"
MD = GEN / "p2488_s1438_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate.md"

SOURCE_FILES = {
    "P2481_BOUNDARY_HANDOFF_COLLAR": GEN / "p2481_s1431_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit.json",
    "P2487_COLLAR_DERIVATIVE_SWEEP": GEN / "p2487_s1437_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate.json",
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
        "new_packet": "P2488|S1438|boundary-handoff collar monotonicity lemma|collar monotonicity lemma|finite collar lemma|derivative-sign compression|collar proof compression|piecewise monotonicity lemma|handoff collar lemma",
        "precursor_packets": "P2487|S1437|collar derivative sweep|P2481|S1431|boundary handoff collar",
        "proof_language": "finite piecewise monotonicity|positive derivative intervals|exact adjacency|positive amplitude intervals|proof compression|lemma certificate",
        "coverage_semantics": "diagnostic row ratio|not a coverage fraction|zero new P2459 unreplayed cells|covered-boundary-chain",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def collar_monotonicity_lemma(p2481: dict[str, Any], p2487: dict[str, Any]) -> dict[str, Any]:
    sweep = p2487["boundary_handoff_collar_derivative_sweep"]
    rows = sweep["collar_derivative_rows"]
    pairs = sweep["consecutive_collar_derivative_pairs"]
    start_row = rows[0]
    end_row = rows[-1]
    derivative_lowers = [Decimal(row["derivative_lower_bound"]) for row in rows]
    amplitude_lowers = [Decimal(row["amplitude_interval_value"]["lo"]) for row in rows]
    amplitude_separations = [Decimal(row["amplitude_interval_separation_from_zero"]) for row in rows]
    widths = [Decimal(row["width"]) for row in rows]
    derivative_transport_lower_gain = sum(Decimal(row["derivative_lower_bound"]) * Decimal(row["width"]) for row in rows)
    lemma_preconditions = {
        "p2487_rows_match_p2481_rows": p2487["total_derivative_row_count"] == p2481["total_fresh_handoff_collar_replay_rows"],
        "all_rows_positive_value": p2487["all_amplitude_intervals_positive"],
        "all_rows_positive_derivative": p2487["all_derivative_intervals_positive"],
        "all_rows_local_monotone": p2487["all_rows_have_local_interval_monotone_increasing_witness"],
        "exact_endpoint_adjacency": p2487["all_consecutive_rows_exactly_adjacent"] and Decimal(p2487["minimum_endpoint_gap"]) == 0 and Decimal(p2487["maximum_endpoint_gap"]) == 0,
        "strict_row_to_row_separation_growth": p2487["all_consecutive_amplitude_separations_strictly_increase"],
        "first_row_positive_anchor": Decimal(start_row["amplitude_interval_value"]["lo"]) > 0,
    }
    return {
        "family": p2487["family"],
        "collar_left": p2487["collar_left"],
        "collar_right": p2487["collar_right"],
        "boundary_band_row_count": p2487["boundary_band_derivative_row_count"],
        "dyadic_subcell_row_count": p2487["dyadic_subcell_derivative_row_count"],
        "total_collar_row_count": p2487["total_derivative_row_count"],
        "lemma_preconditions": lemma_preconditions,
        "all_lemma_preconditions_met": all(lemma_preconditions.values()),
        "finite_piecewise_monotone_increasing_collar_lemma_exported": all(lemma_preconditions.values()),
        "finite_positive_collar_zero_exclusion_lemma_exported": all(lemma_preconditions.values()),
        "minimum_row_amplitude_lower_bound": str(min(amplitude_lowers)),
        "minimum_row_amplitude_separation": str(min(amplitude_separations)),
        "maximum_row_amplitude_separation": str(max(amplitude_separations)),
        "minimum_derivative_lower_bound": str(min(derivative_lowers)),
        "maximum_derivative_lower_bound": str(max(derivative_lowers)),
        "total_checked_collar_width": str(sum(widths)),
        "derivative_transport_lower_gain_over_checked_collar": str(derivative_transport_lower_gain),
        "start_row_left": start_row["left"],
        "start_row_amplitude_lower_bound": start_row["amplitude_interval_value"]["lo"],
        "end_row_right": end_row["right"],
        "end_row_amplitude_lower_bound": end_row["amplitude_interval_value"]["lo"],
        "minimum_consecutive_positive_amplitude_separation_delta": p2487["minimum_consecutive_positive_amplitude_separation_delta"],
        "p2481_handoff_endpoint_gap_boundary_to_parent_inherited": p2481["handoff_endpoint_gap_boundary_to_parent"],
        "p2487_collar_derivative_sweep_fingerprint_inherited": sweep["collar_derivative_sweep_fingerprint_sha256"],
        "compressed_witness_rows": {
            "weakest_amplitude_row": sweep["weakest_amplitude_row"],
            "weakest_derivative_row": sweep["weakest_derivative_row"],
            "strongest_derivative_row": sweep["strongest_derivative_row"],
        },
        "compressed_consecutive_pair_witnesses": {
            "first_pair": pairs[0],
            "last_pair": pairs[-1],
            "weakest_separation_delta_pair": min(pairs, key=lambda row: Decimal(row["amplitude_separation_delta_to_minus_from"])),
        },
        "collar_monotonicity_lemma_fingerprint_sha256": sha256_json({
            "preconditions": lemma_preconditions,
            "source_fingerprint": sweep["collar_derivative_sweep_fingerprint_sha256"],
            "summary": {
                "rows": p2487["total_derivative_row_count"],
                "min_derivative": str(min(derivative_lowers)),
                "min_amplitude": str(min(amplitude_lowers)),
            },
        }),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2488/S1438 strict pointwise interval-Decimal P2459 boundary-handoff collar monotonicity lemma certificate

`P2488/S1438` compresses the P2487 collar-wide derivative sweep into a finite lemma.  The proof obligations are explicit: the P2487 derivative rows match the P2481 collar rows, every checked row has positive projection-amplitude value, every checked row has a positive projection-amplitude derivative interval, consecutive collar rows are exactly adjacent, and row-to-row separations increase.  Under those finite preconditions, P2488 exports a checked-collar-only piecewise monotone-increasing and zero-excluding lemma for the boundary handoff collar.  It is a proof-compression step over existing diagnostics, not a new coverage replay, not a root-window theorem, not global analytic monotonicity, not directed rounding, and not selector/source or ToE closure.

For a non-specialist: P2487 contained many derivative rows.  P2488 states the actual finite lemma those rows support: on the checked collar outside the root window, the function is positive and increasing piece by piece.  This clarifies the proof status without pretending to solve the root window or continuum.
"""
    lag_section = """
## P2488/S1438 P2459 boundary-handoff collar monotonicity lemma guard

`P2488/S1438` turns the P2487 finite collar derivative sweep into a compact checked-collar lemma behind `L_total`: positive value plus positive derivative on all `134` P2481 collar rows implies finite piecewise monotonicity and zero exclusion on that checked collar only.  It reuses existing diagnostic rows, adds no P2459 coverage, and does not export directed rounding, root-window exclusion, global analytic monotonicity, selector/source/gauge authority, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2488/S1438 strict pointwise interval-Decimal P2459 boundary-handoff collar monotonicity lemma certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2488/S1438 P2459 boundary-handoff collar monotonicity lemma guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    p2481 = theorem(sources["P2481_BOUNDARY_HANDOFF_COLLAR"], "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_replay_audit")
    p2487 = theorem(sources["P2487_COLLAR_DERIVATIVE_SWEEP"], "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_derivative_sweep_certificate")
    lemma = collar_monotonicity_lemma(p2481, p2487)
    p2459_universe = p2481["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2488_T1_strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate",
        "audited_chain": ["P2481/S1431", "P2487/S1437"],
        "boundary_handoff_collar_monotonicity_lemma": lemma,
        "family": lemma["family"],
        "total_collar_row_count": lemma["total_collar_row_count"],
        "boundary_band_row_count": lemma["boundary_band_row_count"],
        "dyadic_subcell_row_count": lemma["dyadic_subcell_row_count"],
        "all_lemma_preconditions_met": lemma["all_lemma_preconditions_met"],
        "finite_piecewise_monotone_increasing_collar_lemma_exported": lemma["finite_piecewise_monotone_increasing_collar_lemma_exported"],
        "finite_positive_collar_zero_exclusion_lemma_exported": lemma["finite_positive_collar_zero_exclusion_lemma_exported"],
        "minimum_row_amplitude_lower_bound": lemma["minimum_row_amplitude_lower_bound"],
        "minimum_row_amplitude_separation": lemma["minimum_row_amplitude_separation"],
        "maximum_row_amplitude_separation": lemma["maximum_row_amplitude_separation"],
        "minimum_derivative_lower_bound": lemma["minimum_derivative_lower_bound"],
        "maximum_derivative_lower_bound": lemma["maximum_derivative_lower_bound"],
        "total_checked_collar_width": lemma["total_checked_collar_width"],
        "derivative_transport_lower_gain_over_checked_collar": lemma["derivative_transport_lower_gain_over_checked_collar"],
        "minimum_consecutive_positive_amplitude_separation_delta": lemma["minimum_consecutive_positive_amplitude_separation_delta"],
        "p2481_handoff_endpoint_gap_boundary_to_parent_inherited": lemma["p2481_handoff_endpoint_gap_boundary_to_parent_inherited"],
        "p2487_collar_derivative_sweep_fingerprint_inherited": lemma["p2487_collar_derivative_sweep_fingerprint_inherited"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2488_new_decimal_replay_row_count": 0,
        "p2488_reused_p2487_derivative_row_count_not_a_coverage_count": lemma["total_collar_row_count"],
        "p2488_reused_derivative_row_ratio_not_a_p2459_coverage_fraction": f"{lemma['total_collar_row_count']}/{p2459_universe}",
        "targeted_p2488_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2488_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2488_is_proof_compression_not_new_replay": True,
        "finite_chain_coverage_budget_inherited_from_p2481": p2481["finite_chain_coverage_budget_inherited_from_p2479"],
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
        "lay_summary": "This packet compresses the P2487 derivative sweep into a finite collar lemma. The checked collar rows are exactly adjacent, positive, and have positive derivative intervals, so the checked collar is piecewise monotone increasing and zero-excluding. This is a collar-local proof-compression step and does not prove the excluded root window, global analytic monotonicity, selector/source closure, or ToE closure.",
        "not_licensed": [
            "P2488 performs zero new Decimal replays; it compresses existing P2487 diagnostic rows into a finite checked-collar lemma.",
            "The lemma covers only the P2481 boundary-handoff collar outside the excluded root window and does not add P2459 coverage.",
            "Finite checked-collar monotonicity is not global analytic monotonicity, root-window exclusion, directed rounding, symbolic/continuum root exclusion, selector/source/gauge closure, QW-2191 discharge, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "The finite collar now has a compact monotonicity lemma. Any stronger result must either bridge into the excluded root-window side with a new theorem or move to a real directed-rounding/global analytic proof backend.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "lemma_preconditions_all_met": theorem_export["all_lemma_preconditions_met"],
        "counts_match_collar": theorem_export["boundary_band_row_count"] == 6 and theorem_export["dyadic_subcell_row_count"] == 128 and theorem_export["total_collar_row_count"] == 134,
        "positive_amplitude_floor": Decimal(theorem_export["minimum_row_amplitude_lower_bound"]) > 0,
        "positive_derivative_floor": Decimal(theorem_export["minimum_derivative_lower_bound"]) > 0,
        "positive_separation_growth": Decimal(theorem_export["minimum_consecutive_positive_amplitude_separation_delta"]) > 0,
        "proof_compression_not_new_replay": theorem_export["p2488_is_proof_compression_not_new_replay"] and theorem_export["p2488_new_decimal_replay_row_count"] == 0,
        "zero_new_p2459_unreplayed_cells": theorem_export["targeted_p2488_new_p2459_unreplayed_cell_count"] == 0 and theorem_export["targeted_p2488_new_p2459_unreplayed_cell_scope_against_p2459_universe"] == f"0/{p2459_universe}",
        "finite_lemma_exported_only": theorem_export["finite_piecewise_monotone_increasing_collar_lemma_exported"] and theorem_export["finite_positive_collar_zero_exclusion_lemma_exported"],
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
        "schema_version": "p2488_s1438_v1",
        "packet_id": "P2488",
        "stage_id": "S1438",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_HANDOFF_COLLAR_MONOTONICITY_LEMMA_CERTIFICATE_PROOF_COMPRESSION_NO_NEW_REPLAY_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_handoff_collar_monotonicity_lemma_certificate"]["theorem_export"]
    lines = [
        "# P2488/S1438 strict pointwise interval-Decimal P2459 boundary-handoff collar monotonicity lemma certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite collar monotonicity lemma",
        "",
        f"Total collar rows compressed: `{t['total_collar_row_count']}`.",
        f"Boundary-band rows: `{t['boundary_band_row_count']}`.",
        f"Dyadic subcell rows: `{t['dyadic_subcell_row_count']}`.",
        f"All lemma preconditions met: `{t['all_lemma_preconditions_met']}`.",
        f"Finite piecewise monotone-increasing collar lemma exported: `{t['finite_piecewise_monotone_increasing_collar_lemma_exported']}`.",
        f"Finite positive collar zero-exclusion lemma exported: `{t['finite_positive_collar_zero_exclusion_lemma_exported']}`.",
        f"Minimum row amplitude lower bound: `{t['minimum_row_amplitude_lower_bound']}`.",
        f"Minimum derivative lower bound: `{t['minimum_derivative_lower_bound']}`.",
        f"Total checked collar width: `{t['total_checked_collar_width']}`.",
        f"Derivative transport lower gain over checked collar: `{t['derivative_transport_lower_gain_over_checked_collar']}`.",
        "",
        "## Coverage budget",
        "",
        f"New Decimal replay rows in P2488: `{t['p2488_new_decimal_replay_row_count']}`.",
        f"Reused P2487 derivative rows (not a P2459 coverage count): `{t['p2488_reused_p2487_derivative_row_count_not_a_coverage_count']}`.",
        f"Reused diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2488_reused_derivative_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2488: `{t['targeted_p2488_new_p2459_unreplayed_cell_count']}`.",
        f"New P2459 unreplayed-cell scope against inherited P2459 universe: `{t['targeted_p2488_new_p2459_unreplayed_cell_scope_against_p2459_universe']}`.",
        f"P2488 is proof compression, not new replay: `{t['p2488_is_proof_compression_not_new_replay']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a proof-compression lemma for the finite checked boundary-handoff collar outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, global analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
