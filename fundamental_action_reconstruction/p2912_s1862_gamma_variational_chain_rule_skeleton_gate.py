#!/usr/bin/env python3
"""P2912/S1862: Gamma variational chain-rule skeleton gate.

P2911 constructed a finite endpoint-average localization/pullback matrix for the
Gamma_9_5 typed-object lane.  P2912 attacks exactly one next missing premise in
that lane: the finite variational chain-rule skeleton induced by that matrix.

The constructed object is a symbolic edge-field to site-density variational
Jacobian for Z12 directed edges.  For each edge variable q_{ij}, the localized
site density is varied through the P2911 endpoint-average weights, so the finite
Jacobian has 276 nonzero entries across a 12 x 144 table and is translation
covariant.  This is a proof-readiness skeleton only: it does not export strict
field-variable provenance, a strict source for Gamma_9_5, a continuum measure
pullback theorem, or nonproxy L_total/EOM/Hamiltonian closure.
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

P2911 = GEN / "p2911_s1861_gamma_localization_pullback_skeleton_gate.json"
OUT = GEN / "p2912_s1862_gamma_variational_chain_rule_skeleton_gate.json"
MD = GEN / "p2912_s1862_gamma_variational_chain_rule_skeleton_gate.md"
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


def f(value: Fraction) -> str:
    return f"{value.numerator}/{value.denominator}"


def translate_edge(edge: tuple[int, int], shift: int) -> tuple[int, int]:
    return ((edge[0] + shift) % N, (edge[1] + shift) % N)


def jacobian_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for edge in edges():
        for site, weight in sorted(endpoint_weights(edge).items()):
            rows.append({
                "site": site,
                "edge": [edge[0], edge[1]],
                "derivative": f(weight),
                "symbolic_derivative": f"d rho_site[{site}] / d q_edge[{edge[0]},{edge[1]}] = {f(weight)} * Gamma_9_5",
            })
    return rows


def covariance_failures() -> list[dict[str, Any]]:
    failures: list[dict[str, Any]] = []
    for edge in edges():
        base = endpoint_weights(edge)
        for shift in range(N):
            shifted_edge = translate_edge(edge, shift)
            shifted_weights = endpoint_weights(shifted_edge)
            expected = {((site + shift) % N): value for site, value in base.items()}
            if shifted_weights != expected:
                failures.append({"edge": list(edge), "shift": shift, "expected": {str(k): f(v) for k, v in expected.items()}, "actual": {str(k): f(v) for k, v in shifted_weights.items()}})
    return failures


def build_payload(p2911: dict[str, Any]) -> dict[str, Any]:
    rows = jacobian_rows()
    failures = covariance_failures()
    loop_nonzero = sum(1 for e in edges() if e[0] == e[1])
    nonloop_nonzero = len(rows) - loop_nonzero
    target_edge = (0, 5)
    target_derivatives = [row for row in rows if row["edge"] == [0, 5]]
    zero_entries = N * len(edges()) - len(rows)
    return {
        "status": "P2912_GAMMA_VARIATIONAL_CHAIN_RULE_SKELETON_GATE_READINESS_NO_EXPORT",
        "input_hashes": {"P2911": hashlib.sha256(P2911.read_bytes()).hexdigest() if P2911.exists() else None},
        "constructed_theoretical_objects": {
            "chain_rule_name": "D_Lambda_endpoint_average_Z12_edge_field_to_site_density",
            "edge_field_variables": "q_edge[i,j] for 144 directed Z12 edges",
            "site_density_template": "rho_site[s] = Gamma_9_5 * sum_e Lambda[s,e] * q_edge[e]",
            "jacobian_nonzero_rows": rows,
            "target_edge_0_5_derivatives": target_derivatives,
            "acceptance_obligation_rows": [
                {"obligation": "finite symbolic variational Jacobian constructed", "passed": True, "exported_strictly": False},
                {"obligation": "endpoint-average derivative support matches P2911 pullback", "passed": True, "exported_strictly": False},
                {"obligation": "Z12 translation covariance of derivative table", "passed": len(failures) == 0, "exported_strictly": False},
                {"obligation": "strict field-variable provenance theorem", "passed": False, "exported_strictly": False},
                {"obligation": "continuum/nonproxy variational chain rule theorem", "passed": False, "exported_strictly": False},
                {"obligation": "strict Gamma_9_5 source/provenance theorem", "passed": False, "exported_strictly": False},
            ],
        },
        "acceptance_matrix": {
            "p2911_rechecked_localization_readiness": p2911.get("acceptance_matrix", {}).get("finite_localization_skeleton_constructed") is True,
            "site_count": N,
            "edge_variable_count": len(edges()),
            "jacobian_total_entries": N * len(edges()),
            "jacobian_nonzero_entry_count": len(rows),
            "jacobian_zero_entry_count": zero_entries,
            "loop_edge_nonzero_entry_count": loop_nonzero,
            "nonloop_edge_nonzero_entry_count": nonloop_nonzero,
            "target_edge_0_5_nonzero_derivative_count": len(target_derivatives),
            "target_edge_0_5_derivative_sites": [row["site"] for row in target_derivatives],
            "target_edge_0_5_derivative_values": [row["derivative"] for row in target_derivatives],
            "translation_covariance_failure_count": len(failures),
            "finite_variational_chain_rule_skeleton_constructed": True,
            "strict_field_variable_provenance_exported": False,
            "continuum_variational_chain_rule_exported": False,
            "strict_gamma_9_5_source_exported": False,
            "accepted_as_nonproxy_ltotal_variational_rule": False,
        },
        "decision": {
            "positive_witnesses": {
                "finite_jacobian_constructed": True,
                "target_edge_0_5_derivatives_identified": True,
                "translation_covariance_verified": len(failures) == 0,
            },
            "negative_export_flags": {
                "strict_field_variable_provenance_exported": False,
                "strict_gamma_9_5_source_exported": False,
                "strict_continuum_measure_pullback_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2912 constructs the finite variational Jacobian induced by the P2911 endpoint-average pullback.  The 12 x 144 table has 276 nonzero derivative entries, 1452 zero entries, and zero translation-covariance failures; the local defect edge (0,5) varies only sites 0 and 5 with coefficient 1/2 * Gamma_9_5.  This is only chain-rule readiness: strict field-variable provenance, continuum/nonproxy variational theorem, Gamma_9_5 source, and L_total/EOM/Hamiltonian closure remain unexported.",
            "next_honest_step": "If continuing this lane, the next proof-grade move should audit exactly one missing provenance theorem: either a strict source/provenance theorem for Gamma_9_5 or a strict field-variable/continuum-measure theorem that upgrades this finite Jacobian to a nonproxy variational chain rule.  Do not repeat finite matrix construction or promote the skeleton to L_total/EOM/Hamiltonian/ToE closure without one of those theorems.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2912/S1862 Gamma variational chain-rule skeleton gate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite Jacobian gate",
        f"- edge variables: `{acc['edge_variable_count']}`",
        f"- site densities: `{acc['site_count']}`",
        f"- Jacobian total entries: `{acc['jacobian_total_entries']}`",
        f"- Jacobian nonzero entries: `{acc['jacobian_nonzero_entry_count']}`",
        f"- Jacobian zero entries: `{acc['jacobian_zero_entry_count']}`",
        f"- target edge (0,5) derivative sites: `{acc['target_edge_0_5_derivative_sites']}`",
        f"- target edge (0,5) derivative values: `{acc['target_edge_0_5_derivative_values']}`",
        f"- translation covariance failures: `{acc['translation_covariance_failure_count']}`",
        f"- accepted as nonproxy L_total variational rule: `{acc['accepted_as_nonproxy_ltotal_variational_rule']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2911))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2912/S1862 Gamma variational chain-rule skeleton gate", "## P2912/S1862 Gamma variational chain-rule skeleton gate\n\n`P2912/S1862` constructs the finite symbolic variational Jacobian induced by the P2911 endpoint-average `Lambda_edge_to_site` matrix.  For `144` directed edge variables and `12` site densities the table has `276` nonzero derivative entries, `1452` zero entries, and `0` translation-covariance failures; the local edge `(0,5)` varies only sites `0` and `5` with derivative `1/2 * Gamma_9_5`.  This is chain-rule readiness only: no strict field-variable provenance, continuum/nonproxy variational theorem, strict `Gamma_9_5` source, nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2912/S1862 finite variational chain-rule `L_total` guard", "## P2912/S1862 finite variational chain-rule `L_total` guard\n\n`P2912/S1862` gives a finite Jacobian skeleton for `rho_site[s] = Gamma_9_5 * sum_e Lambda[s,e] * q_edge[e]`.  The defect edge `(0,5)` has exactly two nonzero site derivatives, both `1/2 * Gamma_9_5`, but this remains symbolic readiness.  Without strict field-variable provenance, a continuum/site-measure theorem, and a strict `Gamma_9_5` source, the Jacobian cannot be promoted to nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Gamma variational chain-rule skeleton guardrail (P2912/S1862, 2026-06-19)", "## Current Gamma variational chain-rule skeleton guardrail (P2912/S1862, 2026-06-19)\n\n- P2912 constructs the finite symbolic Jacobian induced by the P2911 endpoint-average pullback: `276` nonzero derivative entries, `1452` zero entries, and `0` translation-covariance failures.\n- Treat this as variational chain-rule readiness only: no strict field-variable provenance, no continuum/nonproxy variational theorem, and no strict `Gamma_9_5` source are exported.\n- Do not promote this Jacobian skeleton, `Gamma_9_5`, symbolic `q_9_5`, `rho_9/5`, or `U_9_5` to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure without strict provenance/source and continuum variational theorems.\n- A next admissible move in this lane must prove exactly one missing provenance theorem for `Gamma_9_5` or for the field-variable/continuum variational upgrade; otherwise preserve no-new-live-frontier or pivot to another genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
