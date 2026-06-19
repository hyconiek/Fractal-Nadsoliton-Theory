#!/usr/bin/env python3
"""P2937/S1887: Z12 multiplication-kernel prime-coordinate source candidate.

P2936 closed replay of merely naming formula classes.  This file therefore tests
one explicit new theorem-object candidate: source the prime-coordinate vector
from the intrinsic finite multiplication endomorphisms m_p:x -> p*x on Z/12Z.

The candidate assigns each audited prime p the kernel-excess

    K_p = |ker(m_p)| - 1 = gcd(p, 12) - 1.

This is a genuine finite Z12 object rather than a coordinate convention.  The
result is still not accepted as a strict prime-log source: it only detects the
non-unit primes 2 and 3, gives zero values for 5, 7, and 11, and does not export
nadsoliton provenance, delta/eta source, or beta/eta coupling.  Its value is a
partial obstruction witness: Z12 multiplication alone cannot source all five
nonzero prime coordinates required by the P2935 gate.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2936 = GEN / "p2936_s1886_source_formula_obstruction_hitting_set_matrix.json"
OUT = GEN / "p2937_s1887_z12_multiplication_kernel_prime_coordinate_source_candidate.json"
MD = GEN / "p2937_s1887_z12_multiplication_kernel_prime_coordinate_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NODES = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]


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
        raise ValueError(f"node {n} has unaudited factor {remaining}")
    return vector


def kernel_size(multiplier: int) -> int:
    return sum(1 for x in range(MODULUS) if (multiplier * x) % MODULUS == 0)


def image_size(multiplier: int) -> int:
    return len({(multiplier * x) % MODULUS for x in range(MODULUS)})


def prime_kernel_rows() -> list[dict[str, Any]]:
    rows = []
    for p in PRIMES:
        k_size = kernel_size(p)
        i_size = image_size(p)
        rows.append({
            "prime": p,
            "gcd_p_12": math.gcd(p, MODULUS),
            "kernel_size": k_size,
            "image_size": i_size,
            "kernel_times_image_equals_12": k_size * i_size == MODULUS,
            "kernel_excess_coordinate": k_size - 1,
            "unit_mod_12": math.gcd(p, MODULUS) == 1,
        })
    return rows


def node_value_rows(prime_values: dict[int, int]) -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        vec = factor_vector(d)
        value = sum(exp * prime_values[p] for exp, p in zip(vec, PRIMES))
        rows.append({
            "node": d,
            "prime_exponent_vector": vec,
            "kernel_excess_value": value,
        })
    return rows


def product_rows(node_values: dict[int, int]) -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        for e in NODES:
            if d * e <= NODES[-1]:
                defect = node_values[d * e] - node_values[d] - node_values[e]
                rows.append({
                    "d": d,
                    "e": e,
                    "de": d * e,
                    "additive_defect": defect,
                    "passes_additivity": defect == 0,
                })
    return rows


def acceptance_rows(prime_rows: list[dict[str, Any]], product_audit: list[dict[str, Any]]) -> list[dict[str, Any]]:
    nonzero_count = sum(1 for row in prime_rows if row["kernel_excess_coordinate"] != 0)
    return [
        {
            "criterion": "finite_Z12_multiplication_object_constructed",
            "satisfied": True,
            "evidence": "kernel and image of x -> p*x mod 12 computed exactly for all audited primes",
        },
        {
            "criterion": "product_additive_carrier_on_nodes_1_to_11",
            "satisfied": all(row["passes_additivity"] for row in product_audit),
            "evidence": "node values are extended linearly through prime exponent vectors",
        },
        {
            "criterion": "all_five_prime_coordinates_nonzero",
            "satisfied": nonzero_count == len(PRIMES),
            "evidence": f"only {nonzero_count} of {len(PRIMES)} prime coordinates are nonzero; unit primes 5,7,11 have zero kernel excess",
        },
        {
            "criterion": "strict_nadsoliton_formula_provenance",
            "satisfied": False,
            "evidence": "Z12 multiplication endomorphism data are finite arithmetic data, not an exported strict nadsoliton source theorem",
        },
        {
            "criterion": "delta_eta_source_law",
            "satisfied": False,
            "evidence": "candidate does not source delta=4/5 or eta=9/5",
        },
        {
            "criterion": "beta_eta_coupling_theorem",
            "satisfied": False,
            "evidence": "candidate does not couple its partial coordinates to the strict damping tail",
        },
    ]


def build_payload(p2936: dict[str, Any]) -> dict[str, Any]:
    p_rows = prime_kernel_rows()
    p_values = {row["prime"]: row["kernel_excess_coordinate"] for row in p_rows}
    n_rows = node_value_rows(p_values)
    n_values = {row["node"]: row["kernel_excess_value"] for row in n_rows}
    prod_rows = product_rows(n_values)
    criteria = acceptance_rows(p_rows, prod_rows)
    accepted = all(row["satisfied"] for row in criteria)
    zero_prime_rows = [row for row in p_rows if row["kernel_excess_coordinate"] == 0]
    return {
        "status": "P2937_Z12_MULTIPLICATION_KERNEL_PRIME_COORDINATE_SOURCE_CANDIDATE_REJECTED_PARTIAL",
        "input_hashes": {"P2936": hashlib.sha256(P2936.read_bytes()).hexdigest() if P2936.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_name": "Z12_Multiplication_Kernel_Excess_PrimeCoordinate_Source_Candidate",
            "definition": "K_p = |ker(x -> p*x mod 12)| - 1 = gcd(p,12)-1 for p in {2,3,5,7,11}",
            "prime_kernel_rows": p_rows,
            "node_value_rows": n_rows,
            "product_additivity_rows": prod_rows,
            "acceptance_rows": criteria,
        },
        "candidate_certificate": {
            "prime_coordinate_vector_order_2_3_5_7_11": [p_values[p] for p in PRIMES],
            "nonzero_prime_coordinate_count": sum(1 for p in PRIMES if p_values[p] != 0),
            "zero_prime_coordinates": [row["prime"] for row in zero_prime_rows],
            "product_pair_count_de_le_11": len(prod_rows),
            "product_additivity_defect_count": sum(1 for row in prod_rows if not row["passes_additivity"]),
            "accepted_strict_prime_log_source": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "new_finite_Z12_multiplication_object_constructed": True,
                "intrinsic_kernel_excess_coordinates_computed": True,
                "product_additivity_verified": all(row["passes_additivity"] for row in prod_rows),
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
            "reason": "P2937 constructs an explicit finite Z12 multiplication-kernel candidate instead of replaying coordinate scans.  The candidate is product-additive and intrinsic to multiplication on Z/12Z, but it only sources nonzero coordinates for 2 and 3; unit primes 5, 7, and 11 receive zero.  It also lacks strict nadsoliton provenance and delta/eta plus beta/eta coupling, so it is rejected as a strict L_p source.",
            "next_honest_step": "Do not continue multiplication-kernel variants unless a new strict theorem couples unit-prime information to nonzero 5/7/11 coordinates and to delta/beta damping.  The next admissible move is a genuinely new unit-sensitive strict source object, or preservation of the P2929-P2937 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["candidate_certificate"]
    lines = [
        "# P2937/S1887 Z12 multiplication-kernel prime-coordinate source candidate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate certificate",
        f"- vector [L2,L3,L5,L7,L11]: `{cert['prime_coordinate_vector_order_2_3_5_7_11']}`",
        f"- nonzero prime coordinates: `{cert['nonzero_prime_coordinate_count']}`",
        f"- zero prime coordinates: `{cert['zero_prime_coordinates']}`",
        f"- product pairs d*e<=11: `{cert['product_pair_count_de_le_11']}`",
        f"- product additivity defects: `{cert['product_additivity_defect_count']}`",
        f"- accepted strict prime-log source: `{cert['accepted_strict_prime_log_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2936))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2937/S1887 Z12 multiplication-kernel prime-coordinate source candidate", "## P2937/S1887 Z12 multiplication-kernel prime-coordinate source candidate\n\n`P2937/S1887` tests one explicit theorem-object candidate after P2936: `K_p = |ker(x -> p*x mod 12)| - 1 = gcd(p,12)-1` for `p in {2,3,5,7,11}`.  This finite Z12 multiplication object gives the coordinate vector `[1,2,0,0,0]` in prime order `2,3,5,7,11` and has `0` product-additivity defects on all `29` audited products.  It is rejected as a strict `L_p` source because unit primes `5,7,11` receive zero and no strict nadsoliton provenance, delta/eta source, or beta/eta coupling theorem is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2937/S1887 multiplication-kernel candidate `L_total` guard", "## P2937/S1887 multiplication-kernel candidate `L_total` guard\n\n`P2937/S1887` constructs a concrete finite Z12 multiplication-kernel coordinate candidate, but the vector `[1,2,0,0,0]` is partial and unsourced on the strict damping side.  Because it does not supply nonzero unit-prime coordinates or delta/eta and beta/eta coupling, it cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 multiplication-kernel prime-coordinate guardrail (P2937/S1887, 2026-06-19)", "## Current Z12 multiplication-kernel prime-coordinate guardrail (P2937/S1887, 2026-06-19)\n\n- P2937 constructs a genuine finite candidate from Z12 multiplication kernels: `K_p = |ker(x -> p*x mod 12)| - 1 = gcd(p,12)-1`.\n- The candidate is product-additive on all `29` audited products and yields `[1,2,0,0,0]` for primes `2,3,5,7,11`, so it detects only nonunit primes and leaves unit primes `5,7,11` at zero.\n- Do not promote this partial multiplication-kernel candidate to a strict `L_p` source, beta/eta damping source, nonproxy `L_total`, EOM/Hamiltonian, bridge closure, role transfer, or ToE.\n- A next admissible move needs a genuinely new unit-sensitive strict source object with provenance and delta/beta coupling, or preservation of the P2929-P2937 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
