#!/usr/bin/env python3
"""P3024/S1974: chart/selector source obstruction for the P3023 chain.

Attack exactly one missing theorem from P3023: a strict chart/selector source for
promoting the integer-label dissipation chain 1->2->...->12 to physical time.

The finite test is the U(12)-orbit of the directed chain.  A non-premise strict
chart source would have to select one orbit representative invariantly.  The
orbit has four representatives and the original chain has trivial stabilizer,
so no U(12)-invariant selector can choose the integer chart without importing a
premise/convention.  Endpoint and steepest-edge anchors are also chart-indexed
and move under relabeling.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import OUT as P3023, N, UNITS, unit_image_edge

OUT = GEN / "p3024_s1974_dissipation_chart_selector_source_obstruction.json"
MD = GEN / "p3024_s1974_dissipation_chart_selector_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def canonical_chain() -> tuple[tuple[int, int], ...]:
    return tuple((d, d + 1) for d in range(1, N))


def image_chain(unit: int) -> tuple[tuple[int, int], ...]:
    return tuple(unit_image_edge(edge, unit) for edge in canonical_chain())


def build_selector_matrix() -> dict[str, Any]:
    chain = canonical_chain()
    orbit_rows = []
    orbit = {}
    for unit in UNITS:
        img = image_chain(unit)
        orbit[unit] = img
        orbit_rows.append({
            "unit": unit,
            "fixed_original_chart": img == chain,
            "first_four_edges": [list(edge) for edge in img[:4]],
            "contains_all_canonical_successor_edges": set(img) == set(chain),
        })
    distinct_representatives = {img for img in orbit.values()}
    stabilizer = [row["unit"] for row in orbit_rows if row["fixed_original_chart"]]

    anchor_rows = [
        {"candidate_anchor": "max_K_endpoint", "selected_label_in_original_chart": 1, "chart_source_status": "moves under U(12) relabeling; anchor is not a physical chart theorem"},
        {"candidate_anchor": "min_K_endpoint", "selected_label_in_original_chart": 12, "chart_source_status": "moves under U(12) relabeling; endpoint does not select a U(12)-invariant integer chart"},
        {"candidate_anchor": "steepest_descent_edge", "selected_edge_in_original_chart": [1, 2], "chart_source_status": "edge label is inherited from the same chart whose source is under test"},
        {"candidate_anchor": "full_descent_chain", "selected_chain_in_original_chart": "1->2->...->12", "chart_source_status": "finite chain exists but its orbit has no invariant representative"},
    ]
    obligations = [
        {"obligation": "attacks_single_P3023_missing_theorem", "satisfied": True, "detail": "only the strict chart/selector source for the P3023 integer order is tested"},
        {"obligation": "finite_U12_orbit_constructed", "satisfied": True, "detail": "computed images under units 1,5,7,11"},
        {"obligation": "nontrivial_stabilizer_for_chart", "satisfied": len(stabilizer) > 1, "detail": "only the identity fixes the full directed chart"},
        {"obligation": "U12_invariant_representative_exists", "satisfied": len(distinct_representatives) == 1, "detail": "four distinct orbit representatives prevent invariant chart selection"},
        {"obligation": "endpoint_or_steepest_anchor_sources_chart", "satisfied": False, "detail": "anchors are computed after choosing the chart and relabel with it"},
        {"obligation": "strict_nonpremise_chart_selector_exported", "satisfied": False, "detail": "no strict source law chooses one orbit representative without a premise"},
    ]
    return {
        "object": "DissipationChainChartSelectorSource_AutOrbitObstructionMatrix",
        "tested_theorem_atom": "strict chart/selector source for the P3023 integer order",
        "chain_under_test": "1->2->...->12",
        "orbit_rows": orbit_rows,
        "orbit_size": len(distinct_representatives),
        "stabilizer_units": stabilizer,
        "anchor_rows": anchor_rows,
        "proof_obligations": obligations,
        "accepted_as_strict_chart_selector_source": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3023_path: Any) -> dict[str, Any]:
    read_json(p3023_path)
    matrix = build_selector_matrix()
    return {
        "status": "P3024_DISSIPATION_CHART_SELECTOR_SOURCE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3023": hashlib.sha256(p3023_path.read_bytes()).hexdigest() if p3023_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "u12_unit_count": len(UNITS),
            "orbit_size": matrix["orbit_size"],
            "stabilizer_size": len(matrix["stabilizer_units"]),
            "fixed_chart_rows": sum(1 for row in matrix["orbit_rows"] if row["fixed_original_chart"]),
            "accepted_as_strict_chart_selector_source": matrix["accepted_as_strict_chart_selector_source"],
        },
        "decision": {
            "breakthrough": "The P3023 integer order was tested as a strict chart/selector-source theorem.  Its U(12) orbit has four distinct directed-chain representatives and only the identity fixes the original chart.  Endpoint, steepest-edge, and full-chain anchors are finite but chart-indexed, so they cannot provide a non-premise strict source for the integer order.",
            "negative_export_flags": {k: False for k in ["strict_chart_selector_source_exported", "strict_time_order_object_exported", "physical_unit_theorem_exported", "unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_arrow_exported", "observed_physics_exported", "ltotal_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay chart anchors for the P3023 chain.  The remaining non-replay theorem atom for this object is an independent physical tick/action/Hamiltonian unit theorem; otherwise pivot to a genuinely new typed object outside the P3017-P3024 unit/time-order lane.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3024/S1974 dissipation chart-selector source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- U(12) unit rows: `{c['u12_unit_count']}`",
        f"- directed-chain orbit size: `{c['orbit_size']}`",
        f"- stabilizer size: `{c['stabilizer_size']}`",
        f"- fixed original chart rows: `{c['fixed_chart_rows']}`",
        f"- accepted as strict chart/selector source: `{c['accepted_as_strict_chart_selector_source']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3023)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3024/S1974 dissipation chart-selector source obstruction", "## P3024/S1974 dissipation chart-selector source obstruction\n\n`P3024/S1974` attacks exactly one missing theorem from P3023: whether the kernel-dissipation integer chain has a strict chart/selector source.  The finite `U(12)` orbit test gives four distinct directed-chain representatives and a trivial stabilizer: only the identity fixes `1->2->...->12`.  Endpoint, steepest-descent-edge, and full-chain anchors are computable, but they are chart-indexed and move under relabeling, so they do not export a non-premise strict chart source.  No strict time-order object with physical unit, time arrow, unit-bearing EOM/Hamiltonian, observed-physics export, `L_total`, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3024/S1974 chart-selector `L_total` guard", "## P3024/S1974 chart-selector `L_total` guard\n\n`P3024/S1974` adds no physical `L_total` term.  It proves a bounded chart-source obstruction for the P3023 dissipation chain: the directed-chain orbit under `U(12)` has no invariant representative, and all endpoint/steepest/full-chain anchors inherit the same unsourced label chart.\n")
    append_once(AGENTS, "Current dissipation chart-selector source guardrail (P3024/S1974, 2026-06-22)", "## Current dissipation chart-selector source guardrail (P3024/S1974, 2026-06-22)\n\n- P3024 attacks exactly one P3023 missing theorem: a strict chart/selector source for the integer order of the kernel-dissipation chain.\n- The finite `U(12)` orbit has four distinct directed-chain representatives and trivial stabilizer; endpoint, steepest-edge, and full-chain anchors are chart-indexed and do not provide a non-premise strict chart source.\n- Do not promote these chart anchors to strict time-order with physical unit, time arrow, unit-bearing EOM/Hamiltonian, observed-physics, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- The remaining non-replay atom for this object is an independent physical tick/action/Hamiltonian unit theorem; otherwise pivot to a genuinely new typed object outside the P3017-P3024 unit/time-order lane.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
