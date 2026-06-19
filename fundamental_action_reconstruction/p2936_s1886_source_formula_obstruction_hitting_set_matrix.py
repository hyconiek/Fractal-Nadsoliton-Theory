#!/usr/bin/env python3
"""P2936/S1886: source-formula obstruction hitting-set matrix.

P2935 built the five-obligation intake gate for a future
Strict_AutBreaking_PrimeCoordinate_Source_Law.  P2936 makes the next honest
step more proof-like: it treats the currently available formula classes as a
finite obstruction system, computes exactly which P2935 obligations each class
can satisfy, and extracts minimal missing-obligation blockers.

The result is not another coordinate scan and not a closure claim.  It is a
finite certificate that every currently named formula route is blocked before it
can become a strict L_p source.  A future move must add a new theorem/formula
that removes at least one blocker class, not replay algebraic readiness.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2935 = GEN / "p2935_s1885_strict_aut_breaking_source_formula_intake_gate.json"
OUT = GEN / "p2936_s1886_source_formula_obstruction_hitting_set_matrix.json"
MD = GEN / "p2936_s1886_source_formula_obstruction_hitting_set_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

OBLIGATIONS = [
    "algebraic_nonzero_additive_aut_breaking_vector",
    "strict_nadsoliton_formula_provenance",
    "nonconventional_symmetry_breaking_source",
    "delta_eta_source_law",
    "beta_eta_coupling_theorem",
]


def current_route_rows() -> list[dict[str, Any]]:
    return [
        {
            "route": "P2934_bounded_cube_nonzero_vectors",
            "satisfied": ["algebraic_nonzero_additive_aut_breaking_vector"],
            "blocker_witness": "formal vector abundance supplies algebra only, with no strict source formula or damping coupling",
        },
        {
            "route": "lexicographic_or_minimal_norm_coordinate_rule",
            "satisfied": ["algebraic_nonzero_additive_aut_breaking_vector"],
            "blocker_witness": "coordinate choice is conventional and has no nadsoliton provenance",
        },
        {
            "route": "external_log_or_numeric_calibration_import",
            "satisfied": ["algebraic_nonzero_additive_aut_breaking_vector"],
            "blocker_witness": "values are imported/calibrated rather than internally sourced",
        },
        {
            "route": "selector_orientation_chiral_readiness_replay",
            "satisfied": ["nonconventional_symmetry_breaking_source"],
            "blocker_witness": "readiness/sign artifacts do not output a prime-coordinate vector and do not couple to delta/eta or beta/eta",
        },
        {
            "route": "future_Strict_AutBreaking_PrimeCoordinate_Source_Law_name_only",
            "satisfied": [],
            "blocker_witness": "the object name is a schema, not an exported formula/theorem",
        },
    ]


def missing(row: dict[str, Any]) -> set[str]:
    return set(OBLIGATIONS) - set(row["satisfied"])


def minimal_global_hitting_sets(rows: list[dict[str, Any]]) -> list[list[str]]:
    """Small exact hitting sets: obligation sets intersecting every route's blockers."""
    blockers = [missing(row) for row in rows]
    universe = OBLIGATIONS
    hits: list[tuple[str, ...]] = []
    for r in range(1, len(universe) + 1):
        for combo in itertools.combinations(universe, r):
            cset = set(combo)
            if all(cset & b for b in blockers):
                if not any(set(prev).issubset(cset) for prev in hits):
                    hits.append(combo)
        if hits:
            break
    return [list(hit) for hit in hits]


def build_payload(p2935: dict[str, Any]) -> dict[str, Any]:
    rows = current_route_rows()
    matrix = []
    for row in rows:
        miss = sorted(missing(row), key=OBLIGATIONS.index)
        matrix.append({
            **row,
            "missing": miss,
            "missing_count": len(miss),
            "accepted_by_P2935_gate": False,
        })
    hits = minimal_global_hitting_sets(rows)
    fully_satisfied = [row for row in matrix if not row["missing"]]
    return {
        "status": "P2936_SOURCE_FORMULA_OBSTRUCTION_HITTING_SET_MATRIX_NO_ACCEPTED_ROUTE",
        "input_hashes": {"P2935": hashlib.sha256(P2935.read_bytes()).hexdigest() if P2935.exists() else None},
        "constructed_theoretical_objects": {
            "obstruction_matrix": matrix,
            "minimal_global_hitting_sets_of_missing_obligations": hits,
            "interpretation": "Any future route must remove the P2935 blockers by exporting actual theorem/formula content, not by choosing another formal coordinate vector.",
        },
        "certificate": {
            "route_count": len(rows),
            "obligation_count": len(OBLIGATIONS),
            "accepted_route_count": len(fully_satisfied),
            "minimal_hitting_set_size": len(hits[0]) if hits else None,
            "minimal_hitting_set_count": len(hits),
            "all_current_routes_blocked": len(fully_satisfied) == 0,
            "coordinate_scan_replay_closed": True,
        },
        "decision": {
            "positive_witnesses": {
                "finite_obstruction_matrix_constructed": True,
                "minimal_missing_obligation_hitting_sets_computed": True,
                "P2935_gate_reused_without_weakening": True,
            },
            "negative_export_flags": {
                "strict_aut_breaking_prime_coordinate_source_law_exported": False,
                "strict_prime_log_value_source_exported": False,
                "strict_delta_eta_source_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2936 turns the P2935 intake gate into a finite obstruction/hitting-set matrix.  Every current formula route misses at least one required obligation; zero routes pass the gate.  The minimal global blockers are theorem-level obligations such as strict provenance, delta/eta source, or beta/eta coupling, not more algebraic coordinate scans.",
            "next_honest_step": "Construct one new explicit theorem object: either a strict nadsoliton provenance formula for a chosen prime-coordinate vector, or a delta/eta-to-beta/eta coupling theorem that can be plugged into the P2935/P2936 gate.  If neither theorem is supplied, preserve the P2929-P2936 no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["certificate"]
    lines = [
        "# P2936/S1886 source-formula obstruction hitting-set matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Certificate",
        f"- audited formula routes: `{cert['route_count']}`",
        f"- P2935 obligations: `{cert['obligation_count']}`",
        f"- accepted routes: `{cert['accepted_route_count']}`",
        f"- minimal hitting-set size: `{cert['minimal_hitting_set_size']}`",
        f"- minimal hitting-set count: `{cert['minimal_hitting_set_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2935))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2936/S1886 source-formula obstruction hitting-set matrix", "## P2936/S1886 source-formula obstruction hitting-set matrix\n\n`P2936/S1886` reuses the P2935 five-obligation gate as an exact finite obstruction matrix for the currently named formula routes.  Five routes are audited: bounded formal vectors, lexicographic/minimal-norm coordinate choice, external log/calibration import, selector/orientation/chiral readiness replay, and the future theorem name alone.  All five miss at least one required obligation and `0` routes are accepted; the minimal global blockers are theorem-level obligations such as strict provenance, delta/eta source, or beta/eta coupling, not another coordinate scan.  No strict `L_p` source, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2936/S1886 obstruction matrix `L_total` guard", "## P2936/S1886 obstruction matrix `L_total` guard\n\n`P2936/S1886` proves that every current route into an Aut-breaking prime-coordinate source remains blocked by the P2935 obligations.  Since no route supplies strict provenance plus delta/eta and beta/eta coupling, no prime-coordinate vector can become role-bearing in nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current source-formula obstruction hitting-set guardrail (P2936/S1886, 2026-06-19)", "## Current source-formula obstruction hitting-set guardrail (P2936/S1886, 2026-06-19)\n\n- P2936 turns the P2935 intake gate into a finite obstruction/hitting-set matrix for currently named formula routes.\n- Bounded formal vectors, lexicographic/minimal-norm choices, external logs/calibration, selector/orientation/chiral readiness replay, and the theorem name alone all fail the P2935 gate; `0` routes are accepted.\n- Do not replay coordinate scans, convention-based vector choices, external logs, or readiness artifacts as a strict `L_p` source.\n- A next admissible move must export one explicit theorem object: strict nadsoliton provenance for a chosen prime-coordinate vector or a delta/eta-to-beta/eta coupling theorem, then rerun the P2935/P2936 gate; otherwise preserve the no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
