#!/usr/bin/env python3
"""P3033/S1983: Z12 block coarse-graining classical-limit obstruction.

Attack exactly one P3028 foundation atom outside the closed P3029-P3032 matter
spectral lane: a classical coarse-graining limit.  We construct a finite
block-averaging/RG candidate on sampled K_strict_gate over Z12 and audit whether
it supplies a chart-independent continuum/classical limit.  It does not: aligned
block averages are computable, but the scale tower is finite, translation-chart
sensitive, and lacks physical length/unit and continuum convergence data.
"""
from __future__ import annotations

import hashlib, itertools, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice import OUT as P3028
from p3032_s1982_matter_spectral_lane_reconciliation_no_new_live_frontier import OUT as P3032

OUT = GEN / "p3033_s1983_z12_block_coarse_graining_classical_limit_obstruction.json"
MD = GEN / "p3033_s1983_z12_block_coarse_graining_classical_limit_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

BLOCK_SIZES = [1, 2, 3, 4, 6, 12]


def kernel_vector(shift: int = 0) -> list[float]:
    values = [k_strict(label) for label in range(1, N + 1)]
    return [values[(j - shift) % N] for j in range(N)]


def block_average(values: list[float], block_size: int) -> list[float]:
    if len(values) % block_size != 0:
        raise ValueError("block size must divide vector length")
    return [sum(values[i:i + block_size]) / block_size for i in range(0, len(values), block_size)]


def l1_delta(a: list[float], b: list[float]) -> float:
    return round(sum(abs(x - y) for x, y in zip(a, b)), 12)


def composition_rows(values: list[float]) -> list[dict[str, Any]]:
    rows = []
    for first, second in itertools.product(BLOCK_SIZES, BLOCK_SIZES):
        if first == 1 or second == 1:
            continue
        product = first * second
        if N % product != 0:
            continue
        direct = block_average(values, product)
        staged = block_average(block_average(values, first), second)
        rows.append({
            "first_block": first,
            "second_block": second,
            "product_block": product,
            "composition_matches_direct": all(abs(x - y) <= 1e-12 for x, y in zip(direct, staged)),
            "l1_delta": l1_delta(direct, staged),
        })
    return rows


def translation_sensitivity_rows() -> list[dict[str, Any]]:
    base = kernel_vector(0)
    rows = []
    for block_size in BLOCK_SIZES:
        base_avg = block_average(base, block_size)
        deltas = []
        for shift in range(1, N):
            shifted_avg = block_average(kernel_vector(shift), block_size)
            deltas.append(l1_delta(base_avg, shifted_avg))
        rows.append({
            "block_size": block_size,
            "coarse_cells": len(base_avg),
            "max_translation_l1_delta": round(max(deltas) if deltas else 0.0, 12),
            "translation_invariant": all(delta <= 1e-12 for delta in deltas),
        })
    return rows


def build_matrix() -> dict[str, Any]:
    values = kernel_vector()
    coarse_rows = []
    for block_size in BLOCK_SIZES:
        avg = block_average(values, block_size)
        coarse_rows.append({
            "block_size": block_size,
            "coarse_cells": len(avg),
            "mean": round(sum(avg) / len(avg), 12),
            "variation": round(max(avg) - min(avg), 12) if len(avg) > 1 else 0.0,
        })
    comp = composition_rows(values)
    translation_rows = translation_sensitivity_rows()
    obligations = [
        {"obligation": "attacks_single_P3028_foundation_atom", "satisfied": True, "detail": "only classical coarse-graining limit is tested"},
        {"obligation": "explicit_coarse_graining_operator", "satisfied": True, "detail": "aligned block averaging on Z12 divisors is constructed"},
        {"obligation": "finite_rg_composition_law", "satisfied": all(row["composition_matches_direct"] for row in comp), "detail": "aligned block averages compose exactly where product blocks divide 12"},
        {"obligation": "chart_translation_independent_limit", "satisfied": False, "detail": "nontrivial block averages depend on the chosen origin/chart"},
        {"obligation": "infinite_refinement_or_continuum_parameter", "satisfied": False, "detail": "Z12 supplies only six divisor block scales, not a continuum/refinement limit"},
        {"obligation": "physical_length_unit_for_scale", "satisfied": False, "detail": "no external physical length/unit theorem labels the block scale"},
        {"obligation": "observer_independent_classical_readout", "satisfied": False, "detail": "coarse cells remain pre-physical internal averages, not classical observables"},
    ]
    return {
        "object": "Z12BlockCoarseGraining_ClassicalLimitObstructionMatrix",
        "tested_foundation_atom": "classical_coarse_graining_limit",
        "coarse_graining_operator": "aligned block average C_b: R^12 -> R^(12/b) for b in {1,2,3,4,6,12}",
        "coarse_rows": coarse_rows,
        "composition_rows": comp,
        "translation_sensitivity_rows": translation_rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "block_scales": len(BLOCK_SIZES),
            "composition_rows": len(comp),
            "composition_rows_passed": sum(1 for row in comp if row["composition_matches_direct"]),
            "translation_sensitive_nontrivial_scales": sum(1 for row in translation_rows if row["block_size"] not in (1, 12) and not row["translation_invariant"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "classical_coarse_graining_limit_exported": all(row["satisfied"] for row in obligations),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3028)
    read_json(P3032)
    matrix = build_matrix()
    return {
        "status": "P3033_Z12_BLOCK_COARSE_GRAINING_CLASSICAL_LIMIT_OBSTRUCTION_NO_CLASSICAL_EXPORT",
        "input_hashes": {
            "P3028": hashlib.sha256(P3028.read_bytes()).hexdigest() if P3028.exists() else None,
            "P3032": hashlib.sha256(P3032.read_bytes()).hexdigest() if P3032.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "A concrete Z12 block-averaging coarse-graining operator exists and has exact finite composition rows, but it does not export a classical coarse-graining limit: nontrivial block averages are chart/translation sensitive, the scale tower is finite, no physical length/unit labels the scale, and the readout remains an internal average rather than an observer-independent classical observable.",
            "negative_export_flags": {k: False for k in ["classical_coarse_graining_limit_exported", "classical_transition_exported", "observed_spacetime_exported", "matter_sector_exported", "observer_readout_exported", "unit_bearing_action_eom_source_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay aligned Z12 block averages as a classical limit.  The next proof-grade move must supply a genuinely new chart-independent refinement/continuum object with a physical length unit, or pivot to another single P3028 foundation atom such as strict selector/branch source or unit-bearing action/EOM/Hamiltonian.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3033/S1983 Z12 block coarse-graining classical-limit obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- block scales: `{c['block_scales']}`",
        f"- composition rows: `{c['composition_rows']}`",
        f"- composition rows passed: `{c['composition_rows_passed']}`",
        f"- translation-sensitive nontrivial scales: `{c['translation_sensitive_nontrivial_scales']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- classical coarse-graining limit exported: `{c['classical_coarse_graining_limit_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3033/S1983 Z12 block coarse-graining classical-limit obstruction", "## P3033/S1983 Z12 block coarse-graining classical-limit obstruction\n\n`P3033/S1983` pivots from the closed P3029-P3032 matter spectral lane to exactly one P3028 foundation atom: classical coarse-graining limit.  It constructs aligned block-averaging operators on sampled `K_strict_gate` over `Z12` for block sizes `{1,2,3,4,6,12}`.  The finite positive is real: aligned composition rows pass where product blocks divide `12`.  The bounded obstruction is that nontrivial block averages are chart/translation sensitive, the scale tower is finite rather than a continuum/refinement limit, no physical length/unit theorem labels scale, and no observer-independent classical readout is exported.  No classical transition, observed physics, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3033/S1983 coarse-graining classical-limit `L_total` guard", "## P3033/S1983 coarse-graining classical-limit `L_total` guard\n\n`P3033/S1983` adds no physical `L_total` term.  The Z12 block-averaging operator is a finite internal coarse-graining scaffold, not a chart-independent continuum limit with a physical length unit or unit-bearing action/EOM insertion.\n")
    append_once(AGENTS, "Current Z12 block coarse-graining classical-limit guardrail (P3033/S1983, 2026-06-23)", "## Current Z12 block coarse-graining classical-limit guardrail (P3033/S1983, 2026-06-23)\n\n- P3033 attacks exactly one P3028 foundation atom outside the closed P3029-P3032 matter spectral lane: classical coarse-graining limit.\n- Aligned Z12 block averages are computable and compose on divisor scales, but nontrivial coarse rows are chart/translation sensitive, finite-scale only, and lack a physical length/unit theorem.\n- Do not promote block averages, finite divisor-scale RG scaffolds, or internal coarse cells to classical transition, observed physics, unit-bearing action/EOM, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- A next move must supply a genuinely new chart-independent refinement/continuum object with physical units, or pivot to another single P3028 foundation atom.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
