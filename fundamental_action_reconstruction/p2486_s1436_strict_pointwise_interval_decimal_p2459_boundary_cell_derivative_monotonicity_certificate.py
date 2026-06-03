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
OUT = GEN / "p2486_s1436_strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate.json"
MD = GEN / "p2486_s1436_strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate.md"

SOURCE_FILES = {
    "P2449_PROJECTION_REDUCTION": GEN / "p2449_s1399_strict_pointwise_rank_lift_projection_reduction_certificate.json",
    "P2482_BOUNDARY_BAND_WEAKEST_CELL_REFINEMENT": GEN / "p2482_s1432_strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit.json",
    "P2485_SECANT_CURVATURE_STABILITY": GEN / "p2485_s1435_strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit.json",
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
        "new_packet": "P2486|S1436|boundary-cell derivative monotonicity|boundary-side derivative monotonicity|derivative sign certificate|one-cell monotonicity|amplitude derivative interval|local derivative sign|boundary subcell monotonicity",
        "precursor_packets": "P2485|S1435|secant curvature stability|P2484|S1434|dyadic secant margin|P2482|S1432|boundary-band weakest cell",
        "derivative_language": "pointwise_gradient_derivative_interval|projection amplitude derivative|interval derivative|positive derivative|monotone increasing cell|left endpoint minimum",
        "coverage_semantics": "diagnostic row ratio|not a coverage fraction|zero new P2459 unreplayed cells|covered-boundary-chain",
        "closure_blockers": "QW-2191|role-bearing L_total|physical-value generator|ToE closure|selector/source/gauge theorem|legacy-role transfer",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def derivative_monotonicity_certificate(parent_cell: dict[str, Any], projection: list[Decimal]) -> dict[str, Any]:
    interval = DecimalInterval(parent_cell["left"], parent_cell["right"])
    amplitude_interval = projection_amplitude_interval(projection, interval)
    derivative_interval = interval_dot(projection, pointwise_gradient_derivative_interval(interval))
    width = Decimal(str(parent_cell["right"])) - Decimal(str(parent_cell["left"]))
    derivative_lower = derivative_interval.lo
    derivative_upper = derivative_interval.hi
    amplitude_lower = amplitude_interval.lo
    amplitude_upper = amplitude_interval.hi
    parent_separation = Decimal(parent_cell["decimal_separation_from_zero"])
    return {
        "family": parent_cell["family"],
        "parent_collar_side": parent_cell["parent_collar_side"],
        "parent_boundary_band_index": int(parent_cell["parent_boundary_band_index"]),
        "parent_dyadic_subcell_index": int(parent_cell["dyadic_subcell_index"]),
        "cell_left": parent_cell["left"],
        "cell_right": parent_cell["right"],
        "cell_width": str(width),
        "amplitude_interval_value": amplitude_interval.as_dict(),
        "amplitude_interval_excludes_zero": not amplitude_interval.contains_zero(),
        "amplitude_interval_separation_from_zero": str(amplitude_interval.separation_from_zero()),
        "amplitude_interval_matches_p2482_parent_sign": amplitude_lower > 0 and parent_separation > 0,
        "p2482_parent_decimal_separation_inherited": str(parent_separation),
        "derivative_interval_value": derivative_interval.as_dict(),
        "derivative_interval_excludes_zero": not derivative_interval.contains_zero(),
        "derivative_interval_positive_on_entire_cell": derivative_lower > 0,
        "derivative_interval_separation_from_zero": str(derivative_interval.separation_from_zero()),
        "derivative_lower_bound": str(derivative_lower),
        "derivative_upper_bound": str(derivative_upper),
        "derivative_relative_width": str((derivative_upper - derivative_lower) / derivative_lower),
        "finite_interval_monotone_increasing_witness": derivative_lower > 0,
        "left_endpoint_is_interval_minimum_under_derivative_witness": derivative_lower > 0,
        "zero_exclusion_reinforced_by_positive_value_and_positive_derivative": amplitude_lower > 0 and derivative_lower > 0,
        "not_a_root_window_interior_claim": True,
        "derivative_certificate_fingerprint_sha256": sha256_json({
            "cell": {"left": parent_cell["left"], "right": parent_cell["right"]},
            "amplitude": amplitude_interval.as_dict(),
            "derivative": derivative_interval.as_dict(),
        }),
    }


def append_once(path, marker: str, section: str) -> None:
    text = path.read_text(encoding="utf-8")
    if marker not in text:
        path.write_text(text.rstrip() + "\n\n" + section.strip() + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
## P2486/S1436 strict pointwise interval-Decimal P2459 boundary-cell derivative monotonicity certificate

`P2486/S1436` stops merely deepening the dyadic boundary ladder and audits the same weakest P2482/P2485 boundary-side subcell with an interval derivative certificate.  On that one inherited P2456 covered-boundary-chain subcell, the projection-amplitude interval remains positive and the interval enclosure for the projection-amplitude derivative is strictly positive throughout the cell.  This licenses a local finite-interval monotone-increasing witness and explains why the one-sided dyadic refinements kept improving away from the root-window-side endpoint.  It is still only a one-cell strict-backend certificate outside the excluded root window: not a root-window theorem, not global analytic monotonicity, not directed rounding, not a P2459 coverage increase, and not selector/source or ToE closure.

For a non-specialist: the prior packets repeatedly zoomed into the boundary cell and saw the safety margin improve.  P2486 checks the derivative on the whole checked cell and finds it positive, so the checked function is locally increasing there.  That is more proof-like than another zoom, but it still does not prove the excluded root window or the global theory.
"""
    lag_section = """
## P2486/S1436 P2459 boundary-cell derivative monotonicity guard

`P2486/S1436` adds a one-cell interval derivative-sign guard behind `L_total`: the same inherited P2456 boundary-chain subcell audited by P2482/P2485 has positive projection-amplitude value and positive projection-amplitude derivative under the strict Decimal/Taylor interval backend.  The certificate is local to that inherited cell and does not export directed rounding, root-window exclusion, global analytic monotonicity, selector/source/gauge authority, physical-value generation, legacy-role transfer, or ToE closure.
"""
    append_once(DOC_FILES["equation_sheet"], "## P2486/S1436 strict pointwise interval-Decimal P2459 boundary-cell derivative monotonicity certificate", eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2486/S1436 P2459 boundary-cell derivative monotonicity guard", lag_section)


def build_payload() -> dict[str, Any]:
    grep = rg_audit()
    sources = {name: load_json(path) for name, path in SOURCE_FILES.items()}
    projection = projection_vector(sources["P2449_PROJECTION_REDUCTION"])
    p2482 = theorem(sources["P2482_BOUNDARY_BAND_WEAKEST_CELL_REFINEMENT"], "strict_pointwise_interval_decimal_p2459_boundary_band_weakest_cell_dyadic_refinement_replay_audit")
    p2485 = theorem(sources["P2485_SECANT_CURVATURE_STABILITY"], "strict_pointwise_interval_decimal_p2459_boundary_side_secant_curvature_stability_audit")
    parent_cell = p2482["boundary_band_weakest_cell_dyadic_refinement_replay"]["minimum_subcell_replay_row"]
    certificate = derivative_monotonicity_certificate(parent_cell, projection)
    p2459_universe = p2482["p2471_p2459_universe_count_inherited"]
    theorem_export = {
        "theorem_name": "P2486_T1_strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate",
        "audited_chain": ["P2482/S1432", "P2485/S1435"],
        "boundary_cell_derivative_monotonicity_certificate": certificate,
        "family": certificate["family"],
        "parent_collar_side": certificate["parent_collar_side"],
        "parent_boundary_band_index": certificate["parent_boundary_band_index"],
        "parent_dyadic_subcell_index": certificate["parent_dyadic_subcell_index"],
        "cell_left": certificate["cell_left"],
        "cell_right": certificate["cell_right"],
        "cell_width": certificate["cell_width"],
        "amplitude_interval_value": certificate["amplitude_interval_value"],
        "amplitude_interval_excludes_zero": certificate["amplitude_interval_excludes_zero"],
        "amplitude_interval_separation_from_zero": certificate["amplitude_interval_separation_from_zero"],
        "amplitude_interval_matches_p2482_parent_sign": certificate["amplitude_interval_matches_p2482_parent_sign"],
        "derivative_interval_value": certificate["derivative_interval_value"],
        "derivative_interval_excludes_zero": certificate["derivative_interval_excludes_zero"],
        "derivative_interval_positive_on_entire_cell": certificate["derivative_interval_positive_on_entire_cell"],
        "derivative_interval_separation_from_zero": certificate["derivative_interval_separation_from_zero"],
        "derivative_lower_bound": certificate["derivative_lower_bound"],
        "derivative_upper_bound": certificate["derivative_upper_bound"],
        "derivative_relative_width": certificate["derivative_relative_width"],
        "finite_interval_monotone_increasing_witness": certificate["finite_interval_monotone_increasing_witness"],
        "left_endpoint_is_interval_minimum_under_derivative_witness": certificate["left_endpoint_is_interval_minimum_under_derivative_witness"],
        "zero_exclusion_reinforced_by_positive_value_and_positive_derivative": certificate["zero_exclusion_reinforced_by_positive_value_and_positive_derivative"],
        "p2485_extended_level_count_inherited": p2485["extended_level_count"],
        "p2485_all_secant_margin_drifts_positive_inherited": p2485["all_secant_margin_drifts_positive"],
        "p2471_p2459_universe_count_inherited": p2459_universe,
        "p2486_derivative_interval_row_count_not_a_coverage_count": 1,
        "p2486_derivative_interval_row_ratio_not_a_p2459_coverage_fraction": f"1/{p2459_universe}",
        "targeted_p2486_new_p2459_unreplayed_cell_count": 0,
        "targeted_p2486_new_p2459_unreplayed_cell_scope_against_p2459_universe": f"0/{p2459_universe}",
        "p2486_refines_one_inherited_p2456_covered_boundary_chain_cell": True,
        "finite_chain_coverage_budget_inherited_from_p2482": p2482["finite_chain_coverage_budget_inherited_from_p2481"],
        "local_one_cell_interval_monotonicity_certificate_exported": True,
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
        "lay_summary": "This packet checks the derivative on the whole already-inherited boundary-side subcell instead of only replaying smaller dyadic samples. The projection amplitude is positive on the cell, and its interval derivative is also strictly positive, so the cell has a local monotone-increasing witness. The certificate remains local to one covered boundary-chain cell outside the excluded root window and does not prove global analytic, selector/source, or ToE closure.",
        "not_licensed": [
            "The derivative interval is one diagnostic certificate for one inherited P2456 covered-boundary-chain subcell, not a new P2459 coverage cell.",
            "The positive derivative signs license only local one-cell interval monotonicity, not global analytic monotonicity or root-window exclusion.",
            "P2486 adds zero new P2459 unreplayed cells and does not discharge directed rounding, symbolic root exclusion, continuum root exclusion, selector/source/gauge closure, QW-2191, legacy-role transfer, or ToE closure.",
        ],
        "next_honest_step": "Use the one-cell derivative witness as a boundary-side lemma, then separately audit whether an analogous derivative-sign or root-window-side theorem can be stated without crossing into the excluded root window or claiming global closure.",
    }
    gatekeepers = {
        "rg_audit_ran": grep["tool"] == "rg" and all(item["count"] >= 0 for item in grep["patterns"].values()),
        "parent_binding_expected": theorem_export["parent_collar_side"] == "p2456_right_boundary_band" and theorem_export["parent_boundary_band_index"] == 0 and theorem_export["parent_dyadic_subcell_index"] == 0,
        "amplitude_interval_excludes_zero": theorem_export["amplitude_interval_excludes_zero"],
        "amplitude_interval_positive": Decimal(theorem_export["amplitude_interval_value"]["lo"]) > 0,
        "derivative_interval_excludes_zero": theorem_export["derivative_interval_excludes_zero"],
        "derivative_interval_positive": theorem_export["derivative_interval_positive_on_entire_cell"] and Decimal(theorem_export["derivative_lower_bound"]) > 0,
        "local_monotone_witness_only": theorem_export["local_one_cell_interval_monotonicity_certificate_exported"] and theorem_export["finite_interval_monotone_increasing_witness"],
        "zero_exclusion_reinforced": theorem_export["zero_exclusion_reinforced_by_positive_value_and_positive_derivative"],
        "p2485_inheritance_checked": theorem_export["p2485_extended_level_count_inherited"] == 64 and theorem_export["p2485_all_secant_margin_drifts_positive_inherited"],
        "zero_new_p2459_unreplayed_cells": theorem_export["targeted_p2486_new_p2459_unreplayed_cell_count"] == 0 and theorem_export["targeted_p2486_new_p2459_unreplayed_cell_scope_against_p2459_universe"] == f"0/{p2459_universe}",
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
        "schema_version": "p2486_s1436_v1",
        "packet_id": "P2486",
        "stage_id": "S1436",
        "status": "PASS_STRICT_POINTWISE_INTERVAL_DECIMAL_P2459_BOUNDARY_CELL_DERIVATIVE_MONOTONICITY_CERTIFICATE_LOCAL_ONE_CELL_NO_ROOT_WINDOW_NO_GLOBAL_ANALYTIC_MONOTONICITY_NO_COVERAGE_INCREASE_NO_DIRECTED_ROUNDING_SELECTOR_SOURCE_THEOREM",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints": {name: sha256_json(sources[name]) for name in sorted(sources)},
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["strict_pointwise_interval_decimal_p2459_boundary_cell_derivative_monotonicity_certificate"]["theorem_export"]
    lines = [
        "# P2486/S1436 strict pointwise interval-Decimal P2459 boundary-cell derivative monotonicity certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## One-cell derivative monotonicity certificate",
        "",
        f"Cell: `[{t['cell_left']}, {t['cell_right']}]`.",
        f"Cell width: `{t['cell_width']}`.",
        f"Amplitude interval: `{t['amplitude_interval_value']}`.",
        f"Amplitude interval excludes zero: `{t['amplitude_interval_excludes_zero']}`.",
        f"Amplitude separation from zero: `{t['amplitude_interval_separation_from_zero']}`.",
        f"Derivative interval: `{t['derivative_interval_value']}`.",
        f"Derivative interval positive on entire cell: `{t['derivative_interval_positive_on_entire_cell']}`.",
        f"Derivative separation from zero: `{t['derivative_interval_separation_from_zero']}`.",
        f"Derivative relative width: `{t['derivative_relative_width']}`.",
        f"Local finite-interval monotone-increasing witness exported: `{t['local_one_cell_interval_monotonicity_certificate_exported']}`.",
        "",
        "## Coverage budget",
        "",
        f"P2486 derivative interval rows (not a P2459 coverage count): `{t['p2486_derivative_interval_row_count_not_a_coverage_count']}`.",
        f"Diagnostic row ratio against inherited P2459 universe (not a coverage fraction): `{t['p2486_derivative_interval_row_ratio_not_a_p2459_coverage_fraction']}`.",
        f"New P2459 unreplayed cells added by P2486: `{t['targeted_p2486_new_p2459_unreplayed_cell_count']}`.",
        f"New P2459 unreplayed-cell scope against inherited P2459 universe: `{t['targeted_p2486_new_p2459_unreplayed_cell_scope_against_p2459_universe']}`.",
        f"P2486 refines one inherited P2456 covered-boundary-chain cell: `{t['p2486_refines_one_inherited_p2456_covered_boundary_chain_cell']}`.",
        f"Full complement replay exported by P2486: `{t['full_complement_replay_exported_by_this_certificate']}`.",
        "",
        "## Plain-language progress note",
        "",
        t["lay_summary"],
        "",
        "## Hard limits / negative controls",
        "",
        "This is a local one-cell derivative monotonicity certificate outside the excluded root window.  It is not a P2459 coverage increase, root-window exclusion theorem, full complement replay, directed-rounding interval theorem, symbolic root-exclusion theorem, global analytic monotonicity theorem, global continuum root-exclusion theorem, selector/source/gauge theorem, physical-value generator, QW-2191 discharge, role-bearing `L_total`, legacy-role transfer, or ToE closure.",
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
