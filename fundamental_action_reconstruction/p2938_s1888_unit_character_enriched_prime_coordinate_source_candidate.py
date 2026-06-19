#!/usr/bin/env python3
"""P2938/S1888: unit-character-enriched prime-coordinate source candidate.

P2937 found a concrete Z12 multiplication-kernel source candidate, but it left
unit primes 5, 7, and 11 at zero.  P2938 tests the honest follow-up requested by
that failure: add one explicit unit-sensitive finite object instead of replaying
kernel-size variants.

The new object is the full {+1,-1}-character table of the unit group U(12).  For
unit primes, it adds the nontrivial-character negativity count

    C_p = #{nontrivial characters chi of U(12) with chi(p) = -1}.

The combined candidate is

    V_p = (|ker(x -> p*x mod 12)| - 1) + C_p.

This yields [1,2,2,2,2] for primes [2,3,5,7,11] and has zero product-additivity
defects when extended through prime exponent vectors.  It is still only a finite
carrier candidate: current artifacts do not export a strict nadsoliton theorem
selecting this character aggregate as a physical value source, nor a delta/eta
source or beta/eta coupling theorem.
"""
from __future__ import annotations

import hashlib
import itertools
import json
import math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2937 = GEN / "p2937_s1887_z12_multiplication_kernel_prime_coordinate_source_candidate.json"
OUT = GEN / "p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate.json"
MD = GEN / "p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NODES = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
UNITS = [1, 5, 7, 11]


def mul_mod(a: int, b: int) -> int:
    return (a * b) % MODULUS


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


def kernel_excess(p: int) -> int:
    return math.gcd(p, MODULUS) - 1


def unit_characters() -> list[dict[int, int]]:
    chars = []
    for values in itertools.product([-1, 1], repeat=len(UNITS)):
        chi = dict(zip(UNITS, values))
        if chi[1] != 1:
            continue
        if all(chi[mul_mod(a, b)] == chi[a] * chi[b] for a in UNITS for b in UNITS):
            chars.append(chi)
    return sorted(chars, key=lambda row: tuple(row[u] for u in UNITS))


def character_rows(chars: list[dict[int, int]]) -> list[dict[str, Any]]:
    rows = []
    for index, chi in enumerate(chars):
        rows.append({
            "character_index": index,
            "values_on_units_1_5_7_11": [chi[u] for u in UNITS],
            "trivial_character": all(chi[u] == 1 for u in UNITS),
            "multiplicative_on_U12": all(chi[mul_mod(a, b)] == chi[a] * chi[b] for a in UNITS for b in UNITS),
        })
    return rows


def prime_coordinate_rows(chars: list[dict[int, int]]) -> list[dict[str, Any]]:
    nontrivial = [chi for chi in chars if not all(chi[u] == 1 for u in UNITS)]
    rows = []
    for p in PRIMES:
        unit = p in UNITS
        negativity = sum(1 for chi in nontrivial if unit and chi[p] == -1)
        k = kernel_excess(p)
        rows.append({
            "prime": p,
            "kernel_excess": k,
            "unit_character_negativity_count": negativity,
            "combined_coordinate": k + negativity,
            "unit_mod_12": unit,
        })
    return rows


def node_value_rows(prime_values: dict[int, int]) -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        vec = factor_vector(d)
        rows.append({
            "node": d,
            "prime_exponent_vector": vec,
            "combined_value": sum(exp * prime_values[p] for exp, p in zip(vec, PRIMES)),
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


def acceptance_rows(prime_rows: list[dict[str, Any]], prod_rows: list[dict[str, Any]], chars: list[dict[int, int]]) -> list[dict[str, Any]]:
    return [
        {
            "criterion": "full_U12_character_table_constructed",
            "satisfied": len(chars) == 4,
            "evidence": "all homomorphisms U(12)->{+1,-1} are enumerated",
        },
        {
            "criterion": "unit_primes_receive_nonzero_coordinates",
            "satisfied": all(row["combined_coordinate"] != 0 for row in prime_rows if row["unit_mod_12"]),
            "evidence": "nontrivial-character negativity count repairs the P2937 zero rows for 5,7,11",
        },
        {
            "criterion": "all_five_prime_coordinates_nonzero",
            "satisfied": all(row["combined_coordinate"] != 0 for row in prime_rows),
            "evidence": "combined vector is nonzero on every audited prime",
        },
        {
            "criterion": "product_additive_carrier_on_nodes_1_to_11",
            "satisfied": all(row["passes_additivity"] for row in prod_rows),
            "evidence": "node values are extended through prime exponent vectors",
        },
        {
            "criterion": "strict_nadsoliton_formula_provenance",
            "satisfied": False,
            "evidence": "finite U12 character aggregation is not yet exported as a strict nadsoliton source theorem",
        },
        {
            "criterion": "delta_eta_source_law",
            "satisfied": False,
            "evidence": "candidate does not derive delta=4/5 or eta=9/5",
        },
        {
            "criterion": "beta_eta_coupling_theorem",
            "satisfied": False,
            "evidence": "candidate does not couple its vector to the strict damping/compression tail",
        },
    ]


def build_payload(p2937: dict[str, Any]) -> dict[str, Any]:
    chars = unit_characters()
    char_rows = character_rows(chars)
    p_rows = prime_coordinate_rows(chars)
    p_values = {row["prime"]: row["combined_coordinate"] for row in p_rows}
    n_rows = node_value_rows(p_values)
    n_values = {row["node"]: row["combined_value"] for row in n_rows}
    prod_rows = product_rows(n_values)
    criteria = acceptance_rows(p_rows, prod_rows, chars)
    accepted = all(row["satisfied"] for row in criteria)
    return {
        "status": "P2938_UNIT_CHARACTER_ENRICHED_PRIME_COORDINATE_SOURCE_CANDIDATE_CONDITIONAL_REJECTED",
        "input_hashes": {"P2937": hashlib.sha256(P2937.read_bytes()).hexdigest() if P2937.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_name": "UnitCharacter_Enriched_Z12_PrimeCoordinate_Source_Candidate",
            "definition": "V_p=(gcd(p,12)-1)+#{nontrivial chi:U(12)->±1 with chi(p)=-1}",
            "unit_character_rows": char_rows,
            "prime_coordinate_rows": p_rows,
            "node_value_rows": n_rows,
            "product_additivity_rows": prod_rows,
            "acceptance_rows": criteria,
        },
        "candidate_certificate": {
            "unit_group": UNITS,
            "character_count": len(chars),
            "nontrivial_character_count": sum(1 for chi in chars if not all(chi[u] == 1 for u in UNITS)),
            "prime_coordinate_vector_order_2_3_5_7_11": [p_values[p] for p in PRIMES],
            "all_prime_coordinates_nonzero": all(p_values[p] != 0 for p in PRIMES),
            "product_pair_count_de_le_11": len(prod_rows),
            "product_additivity_defect_count": sum(1 for row in prod_rows if not row["passes_additivity"]),
            "accepted_strict_prime_log_source": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "unit_sensitive_finite_object_constructed": True,
                "P2937_unit_prime_zero_defect_repaired_formally": True,
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
            "reason": "P2938 supplies a genuinely unit-sensitive finite carrier by aggregating the full U(12) character table.  It repairs the P2937 zero unit-prime rows and gives [1,2,2,2,2] with zero product-additivity defects.  However, it remains a conditional finite carrier because current artifacts do not export the required strict nadsoliton provenance theorem or the delta/eta and beta/eta coupling theorems.",
            "next_honest_step": "The next admissible move is not another finite character aggregate.  It must prove a strict provenance theorem that identifies this U12 character aggregate as nadsoliton-sourced and couples it to delta/eta and beta/eta, or else preserve the P2929-P2938 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["candidate_certificate"]
    lines = [
        "# P2938/S1888 unit-character-enriched prime-coordinate source candidate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate certificate",
        f"- U(12): `{cert['unit_group']}`",
        f"- characters: `{cert['character_count']}`",
        f"- nontrivial characters: `{cert['nontrivial_character_count']}`",
        f"- vector [L2,L3,L5,L7,L11]: `{cert['prime_coordinate_vector_order_2_3_5_7_11']}`",
        f"- all prime coordinates nonzero: `{cert['all_prime_coordinates_nonzero']}`",
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
    payload = build_payload(read_json(P2937))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2938/S1888 unit-character-enriched prime-coordinate source candidate", "## P2938/S1888 unit-character-enriched prime-coordinate source candidate\n\n`P2938/S1888` repairs the P2937 unit-prime zero defect at the finite-carrier level by adding the full `U(12)` character table.  The candidate `V_p=(gcd(p,12)-1)+#{nontrivial chi:U(12)->±1 with chi(p)=-1}` gives `[1,2,2,2,2]` for primes `2,3,5,7,11` and has `0` product-additivity defects on all `29` audited products.  This is still conditional carrier readiness only: no strict nadsoliton provenance theorem, delta/eta source law, or beta/eta coupling theorem is exported, so no strict `L_p` source, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2938/S1888 unit-character candidate `L_total` guard", "## P2938/S1888 unit-character candidate `L_total` guard\n\n`P2938/S1888` constructs a unit-sensitive finite carrier with vector `[1,2,2,2,2]`, but it remains unsourced on the strict nadsoliton and damping-coupling sides.  Until a provenance theorem and delta/eta plus beta/eta coupling are exported, this carrier cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current unit-character-enriched prime-coordinate guardrail (P2938/S1888, 2026-06-19)", "## Current unit-character-enriched prime-coordinate guardrail (P2938/S1888, 2026-06-19)\n\n- P2938 adds the full `U(12)` character table as a genuinely unit-sensitive finite carrier after P2937.\n- The combined candidate `V_p=(gcd(p,12)-1)+#{nontrivial chi with chi(p)=-1}` yields `[1,2,2,2,2]` and has `0` product-additivity defects on the `29` audited products.\n- This repairs the P2937 unit-prime zero rows only at finite-carrier level; it does not export strict nadsoliton provenance, delta/eta source, beta/eta coupling, strict `L_p`, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must prove strict provenance and damping coupling for this exact carrier, or preserve the P2929-P2938 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
