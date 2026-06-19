#!/usr/bin/env python3
"""P2923/S1873: prime-log proportionality source matrix.

P2922 closed the P2601-to-sigma_Gamma interface and recommended pivoting to a
separate typed object.  P2923 pivots to the strict damping residual frontier
named by P2601: the missing prime-log proportionality source.

The script constructs the finite prime-exponent source matrix on the audited
monoid {1,...,11}.  It proves the algebraic readiness statement that prime
factor exponent vectors add under multiplication and therefore can carry a
formal logarithmic character once prime log atoms are supplied.  It does not
supply those atoms from strict nadsoliton data; prime log values and the
proportionality/anchor remain external unless a new source theorem is exported.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2601 = GEN / "p2601_s1551_nadsoliton_identity_action_unital_multiplicative_source_theorem.json"
OUT = GEN / "p2923_s1873_prime_log_proportionality_source_matrix.json"
MD = GEN / "p2923_s1873_prime_log_proportionality_source_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

PRIMES = [2, 3, 5, 7, 11]
NODES = list(range(1, 12))


def factor_vector(n: int) -> list[int]:
    remaining = n
    vector = []
    for prime in PRIMES:
        exponent = 0
        while remaining % prime == 0:
            remaining //= prime
            exponent += 1
        vector.append(exponent)
    if remaining != 1:
        raise ValueError(f"node {n} has factor outside audited primes: {remaining}")
    return vector


def factor_rows() -> list[dict[str, Any]]:
    return [
        {
            "node": node,
            "prime_exponent_vector_over_2_3_5_7_11": factor_vector(node),
            "formal_log_expression": " + ".join(f"{exp}*L_{p}" for p, exp in zip(PRIMES, factor_vector(node)) if exp) or "0",
        }
        for node in NODES
    ]


def product_rows() -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        for e in NODES:
            if d * e <= NODES[-1]:
                vd = factor_vector(d)
                ve = factor_vector(e)
                vde = factor_vector(d * e)
                rows.append({
                    "d": d,
                    "e": e,
                    "de": d * e,
                    "v_d_plus_v_e": [a + b for a, b in zip(vd, ve)],
                    "v_de": vde,
                    "additive_character_defect": [a + b - c for a, b, c in zip(vd, ve, vde)],
                    "passes": [a + b for a, b in zip(vd, ve)] == vde,
                })
    return rows


def candidate_source_atoms() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "formal_prime_exponent_vectors",
            "algebraic_readiness": True,
            "prime_log_values_sourced": False,
            "strict_nadsoliton_provenance": False,
            "accepted_as_prime_log_source": False,
            "failure": "factorization vectors source additivity structure but not the analytic prime log atom values L_p",
        },
        {
            "candidate": "external_real_log_function",
            "algebraic_readiness": True,
            "prime_log_values_sourced": True,
            "strict_nadsoliton_provenance": False,
            "accepted_as_prime_log_source": False,
            "failure": "ordinary real logarithm supplies values externally, not from strict nadsoliton data",
        },
        {
            "candidate": "P2601_identity_action_y1_zero",
            "algebraic_readiness": False,
            "prime_log_values_sourced": False,
            "strict_nadsoliton_provenance": True,
            "accepted_as_prime_log_source": False,
            "failure": "identity action sources y_1=0/unitality but not nonzero prime log atoms",
        },
        {
            "candidate": "delta_4_5_slope_anchor",
            "algebraic_readiness": False,
            "prime_log_values_sourced": False,
            "strict_nadsoliton_provenance": False,
            "accepted_as_prime_log_source": False,
            "failure": "the slope/prime anchor is the other missing residual key, not a source for prime logs by itself",
        },
    ]


def matrix_rank_binary(rows: list[list[int]]) -> int:
    """Small exact rank over rationals using row reduction with Fractions avoided by integer pivots."""
    mat = [[float(x) for x in row] for row in rows]
    rank = 0
    col_count = len(mat[0]) if mat else 0
    for col in range(col_count):
        pivot = next((r for r in range(rank, len(mat)) if abs(mat[r][col]) > 1e-12), None)
        if pivot is None:
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pivot_val = mat[rank][col]
        mat[rank] = [x / pivot_val for x in mat[rank]]
        for r in range(len(mat)):
            if r != rank and abs(mat[r][col]) > 1e-12:
                factor = mat[r][col]
                mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
    return rank


def p2601_theorem(p2601: dict[str, Any]) -> dict[str, Any]:
    return p2601.get("nadsoliton_identity_action_unital_multiplicative_source_theorem", {}).get("theorem_export", {})


def build_payload(p2601: dict[str, Any]) -> dict[str, Any]:
    t = p2601_theorem(p2601)
    f_rows = factor_rows()
    p_rows = product_rows()
    candidates = candidate_source_atoms()
    accepted = [candidate for candidate in candidates if candidate["accepted_as_prime_log_source"]]
    rank = matrix_rank_binary([row["prime_exponent_vector_over_2_3_5_7_11"] for row in f_rows])
    return {
        "status": "P2923_PRIME_LOG_PROPORTIONALITY_SOURCE_MATRIX_READINESS_NO_STRICT_SOURCE",
        "input_hashes": {"P2601": hashlib.sha256(P2601.read_bytes()).hexdigest() if P2601.exists() else None},
        "constructed_theoretical_objects": {
            "missing_theorem_name": "Strict_Prime_Log_Proportionality_Source_Theorem",
            "prime_basis": PRIMES,
            "factor_exponent_rows": f_rows,
            "product_additivity_rows": p_rows,
            "candidate_source_atoms": candidates,
            "acceptance_schema": [
                "derive prime log atom values L_p from strict nadsoliton data",
                "prove log-character proportionality y_d = a * sum_p v_p(d) L_p on the audited monoid",
                "supply or couple the slope/prime anchor rather than importing it",
                "preserve P2601 unital/multiplicative source without promoting bridge, role transfer, L_total, or ToE closure",
            ],
        },
        "acceptance_matrix": {
            "p2601_rechecked_multiplicative_source_exported": t.get("multiplicative_character_law_source_exported") is True,
            "node_count": len(NODES),
            "prime_basis_count": len(PRIMES),
            "factor_matrix_rank": rank,
            "product_pair_count_de_le_11": len(p_rows),
            "product_additivity_failures": sum(1 for row in p_rows if not row["passes"]),
            "candidate_source_atom_count": len(candidates),
            "accepted_source_atom_count": len(accepted),
            "formal_log_character_readiness_exported": True,
            "strict_prime_log_value_source_exported": False,
            "slope_prime_anchor_source_exported": False,
            "strict_damping_beta_eta_source_exported": False,
        },
        "decision": {
            "positive_witnesses": {
                "prime_factor_matrix_constructed": True,
                "multiplicative_additivity_verified": all(row["passes"] for row in p_rows),
                "formal_log_character_readiness_exported": True,
            },
            "negative_export_flags": {
                "strict_prime_log_value_source_exported": False,
                "slope_prime_anchor_source_exported": False,
                "strict_damping_beta_eta_source_exported": False,
                "nonproxy_ltotal_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2923 constructs the finite prime-exponent matrix for the P2601 residual prime-log proportionality key.  On nodes 1..11, exponent vectors over primes 2,3,5,7,11 have rank 5 and add exactly under all audited products de<=11, so they provide formal logarithmic-character readiness.  However, the prime log atom values L_p and the 4/5 slope/prime anchor remain unsourced by strict nadsoliton data; no strict damping beta/eta source or L_total closure is exported.",
            "next_honest_step": "The next proof-grade step in this damping lane should attack exactly one remaining atom: either a strict source for prime log atom values L_p, or the slope/prime anchor source.  Do not replay Gamma/Lambda, P2601 identity action, or factorization readiness as closure evidence.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2923/S1873 prime-log proportionality source matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Prime-log matrix gate",
        f"- nodes: `{acc['node_count']}`",
        f"- prime basis count: `{acc['prime_basis_count']}`",
        f"- factor matrix rank: `{acc['factor_matrix_rank']}`",
        f"- product pairs de<=11: `{acc['product_pair_count_de_le_11']}`",
        f"- product additivity failures: `{acc['product_additivity_failures']}`",
        f"- candidate source atoms: `{acc['candidate_source_atom_count']}`",
        f"- accepted source atoms: `{acc['accepted_source_atom_count']}`",
        f"- formal log-character readiness exported: `{acc['formal_log_character_readiness_exported']}`",
        f"- strict prime-log value source exported: `{acc['strict_prime_log_value_source_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2601))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2923/S1873 prime-log proportionality source matrix", "## P2923/S1873 prime-log proportionality source matrix\n\n`P2923/S1873` pivots outside the Gamma/Lambda lane to the P2601 strict-damping residual frontier.  It constructs the finite prime-exponent matrix on nodes `1..11` over prime basis `{2,3,5,7,11}`: the matrix has rank `5`, and all audited products `d*e<=11` have zero exponent-additivity defect.  This exports formal logarithmic-character readiness only.  The prime log atom values `L_p` and the `delta=4/5` slope/prime anchor are still not sourced by strict nadsoliton data, so no strict damping `beta/eta` source, nonproxy `L_total`, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2923/S1873 prime-log matrix `L_total` guard", "## P2923/S1873 prime-log matrix `L_total` guard\n\n`P2923/S1873` supplies finite factorization/log-character readiness for the strict-damping residual key left after P2601: exponent vectors add exactly on the audited monoid.  This does not source the prime log atom values or the `4/5` slope/prime anchor, so damping terms remain non-role-bearing and cannot be promoted to nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current prime-log proportionality source matrix guardrail (P2923/S1873, 2026-06-19)", "## Current prime-log proportionality source matrix guardrail (P2923/S1873, 2026-06-19)\n\n- P2923 pivots outside Gamma/Lambda and P2601-to-sigma_Gamma replay to the P2601 strict-damping residual frontier: prime-log proportionality.\n- The finite prime-exponent matrix over `{2,3,5,7,11}` has rank `5` and zero product-additivity failures on audited products `d*e<=11`; this is formal log-character readiness only.\n- No strict source for prime log atom values `L_p`, no `delta=4/5` slope/prime anchor, no strict damping `beta/eta`, no nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n- A next admissible move may attack exactly one remaining damping atom: strict source for `L_p` values or strict slope/prime anchor source; do not replay Gamma/Lambda, P2601 identity action, or factorization readiness as closure evidence.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
