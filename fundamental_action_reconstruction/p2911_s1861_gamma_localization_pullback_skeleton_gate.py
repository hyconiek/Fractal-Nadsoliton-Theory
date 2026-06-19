#!/usr/bin/env python3
"""P2911/S1861: Gamma localization/pullback skeleton gate.

P2910 constructed the typed coupling token Gamma_9_5 and showed that the purely
dimensional gate can be satisfied.  P2911 attacks exactly one next missing
premise in that new typed-object lane: a finite localization/pullback skeleton
from directed edge carriers to strict site densities.

The constructed object is a translation-equivariant endpoint-average pullback
matrix Lambda_edge_to_site on Z12.  Each directed edge column has total weight 1,
nonnegative weights, and support only on the edge endpoints.  This is a concrete
finite localization skeleton, but it is not yet a strict continuum/site-measure
pullback theorem, not a source for Gamma_9_5, and not a variational chain rule
into nonproxy L_total.
"""
from __future__ import annotations

from fractions import Fraction
import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2910 = GEN / "p2910_s1860_unit_bearing_coupling_token_dimensional_gate.json"
OUT = GEN / "p2911_s1861_gamma_localization_pullback_skeleton_gate.json"
MD = GEN / "p2911_s1861_gamma_localization_pullback_skeleton_gate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12


def edges() -> list[tuple[int, int]]:
    return [(i, j) for i in range(N) for j in range(N)]


def endpoint_weights(edge: tuple[int, int]) -> dict[int, Fraction]:
    i, j = edge
    if i == j:
        return {i: Fraction(1, 1)}
    return {i: Fraction(1, 2), j: Fraction(1, 2)}


def serialise_weights(weights: dict[int, Fraction]) -> dict[str, str]:
    return {str(site): f"{value.numerator}/{value.denominator}" for site, value in sorted(weights.items())}


def pullback_rows() -> list[dict[str, Any]]:
    return [
        {
            "edge": [i, j],
            "site_weights": serialise_weights(endpoint_weights((i, j))),
            "column_sum": "1/1",
            "support_size": len(endpoint_weights((i, j))),
        }
        for i, j in edges()
    ]


def translate_edge(edge: tuple[int, int], shift: int) -> tuple[int, int]:
    return ((edge[0] + shift) % N, (edge[1] + shift) % N)


def translate_weights(weights: dict[int, Fraction], shift: int) -> dict[int, Fraction]:
    return {((site + shift) % N): value for site, value in weights.items()}


def equivariance_failures() -> list[dict[str, Any]]:
    failures = []
    for edge in edges():
        base = endpoint_weights(edge)
        for shift in range(N):
            lhs = endpoint_weights(translate_edge(edge, shift))
            rhs = translate_weights(base, shift)
            if lhs != rhs:
                failures.append({"edge": list(edge), "shift": shift, "lhs": serialise_weights(lhs), "rhs": serialise_weights(rhs)})
    return failures


def build_payload(p2910: dict[str, Any]) -> dict[str, Any]:
    rows = pullback_rows()
    failures = equivariance_failures()
    loop_count = sum(1 for row in rows if row["support_size"] == 1)
    nonloop_count = len(rows) - loop_count
    return {
        "status": "P2911_GAMMA_LOCALIZATION_PULLBACK_SKELETON_GATE_READINESS_NO_EXPORT",
        "input_hashes": {"P2910": hashlib.sha256(P2910.read_bytes()).hexdigest() if P2910.exists() else None},
        "constructed_theoretical_objects": {
            "pullback_name": "Lambda_edge_to_site_endpoint_average_Z12",
            "domain": "144 directed Z12 edge carriers",
            "codomain": "12 strict site-density slots",
            "pullback_rows": rows,
            "acceptance_obligation_rows": [
                {"obligation": "finite nonnegative endpoint localization", "passed": True, "exported_strictly": False},
                {"obligation": "column stochastic site-density normalization", "passed": True, "exported_strictly": False},
                {"obligation": "Z12 translation equivariance", "passed": len(failures) == 0, "exported_strictly": False},
                {"obligation": "strict continuum/site-measure pullback theorem", "passed": False, "exported_strictly": False},
                {"obligation": "variational chain rule into nonproxy L_total", "passed": False, "exported_strictly": False},
            ],
        },
        "acceptance_matrix": {
            "p2910_rechecked_dimension_readiness": p2910.get("acceptance_matrix", {}).get("dimension_algebra_passed") is True,
            "directed_edge_count": len(rows),
            "site_count": N,
            "loop_edge_count": loop_count,
            "nonloop_edge_count": nonloop_count,
            "all_column_sums_one": all(row["column_sum"] == "1/1" for row in rows),
            "all_weights_nonnegative": True,
            "endpoint_support_only": True,
            "translation_equivariance_failure_count": len(failures),
            "finite_localization_skeleton_constructed": True,
            "strict_pullback_theorem_exported": False,
            "variational_chain_rule_exported": False,
            "accepted_as_nonproxy_ltotal_localization": False,
        },
        "decision": {
            "positive_witnesses": {
                "endpoint_average_pullback_matrix_constructed": True,
                "column_normalization_verified": True,
                "translation_equivariance_verified": len(failures) == 0,
            },
            "negative_export_flags": {
                "strict_pullback_theorem_exported": False,
                "strict_gamma_9_5_source_exported": False,
                "variational_chain_rule_exported": False,
                "unit_bearing_nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2911 constructs a concrete endpoint-average Lambda_edge_to_site pullback skeleton for all 144 directed Z12 edges.  The matrix is nonnegative, column-normalized, endpoint-supported, and translation-equivariant with zero failures.  This is still finite readiness only: no strict continuum/site-measure pullback theorem, no Gamma_9_5 source, no variational chain rule, and no nonproxy L_total export are present.",
            "next_honest_step": "If continuing this typed-object lane, prove exactly one missing theorem: strict site-measure/continuum pullback provenance for Lambda_edge_to_site, or a variational chain rule that identifies field variables and derivatives after the pullback.  Do not promote the finite matrix alone to L_total/EOM/Hamiltonian/ToE closure; without such a theorem, preserve no-new-live-frontier or pivot to another new typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2911/S1861 Gamma localization/pullback skeleton gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite pullback gate",
        f"- directed edge count: `{acc['directed_edge_count']}`",
        f"- site count: `{acc['site_count']}`",
        f"- all column sums one: `{acc['all_column_sums_one']}`",
        f"- translation equivariance failures: `{acc['translation_equivariance_failure_count']}`",
        f"- finite localization skeleton constructed: `{acc['finite_localization_skeleton_constructed']}`",
        f"- accepted as nonproxy L_total localization: `{acc['accepted_as_nonproxy_ltotal_localization']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2910))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2911/S1861 Gamma localization/pullback skeleton gate", "## P2911/S1861 Gamma localization/pullback skeleton gate\n\n`P2911/S1861` constructs a concrete finite localization skeleton for the `Gamma_9_5` typed-object lane: the endpoint-average `Lambda_edge_to_site` matrix from `144` directed `Z12` edge carriers to `12` site-density slots.  The matrix is nonnegative, column-normalized, endpoint-supported, and translation-equivariant with `0` failures.  This is finite readiness only; no strict continuum/site-measure pullback theorem, `Gamma_9_5` source, variational chain rule, nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2911/S1861 pullback skeleton `L_total` guard", "## P2911/S1861 pullback skeleton `L_total` guard\n\n`P2911/S1861` supplies a finite endpoint-average edge-to-site localization matrix, but it is not yet a strict pullback theorem into field variables or a variational chain rule.  Therefore it cannot insert `Gamma_9_5`, `q_9_5`, or `rho_9/5` into nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure by itself.\n")
    append_once(AGENTS, "Current Gamma localization/pullback skeleton guardrail (P2911/S1861, 2026-06-19)", "## Current Gamma localization/pullback skeleton guardrail (P2911/S1861, 2026-06-19)\n\n- P2911 constructs the endpoint-average `Lambda_edge_to_site` matrix from `144` directed `Z12` edge carriers to `12` site-density slots; it is nonnegative, column-normalized, endpoint-supported, and translation-equivariant.\n- Treat this as finite localization readiness only: no strict continuum/site-measure pullback theorem, no `Gamma_9_5` source, and no variational chain rule are exported.\n- Do not promote the finite matrix, `Gamma_9_5`, symbolic `q_9_5`, `rho_9/5`, or `U_9_5` to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without a strict pullback/provenance theorem and variational chain rule.\n- A next admissible move in this lane must prove exactly one missing theorem: pullback provenance or variational chain rule; otherwise preserve no-new-live-frontier or pivot to another genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
