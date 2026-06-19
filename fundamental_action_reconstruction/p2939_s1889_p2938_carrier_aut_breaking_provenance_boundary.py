#!/usr/bin/env python3
"""P2939/S1889: P2938 carrier Aut-breaking and provenance boundary.

P2938 produced a stronger finite carrier [1,2,2,2,2].  P2939 does the next
proof-grade check before any promotion: verify the exact algebraic obligations
that P2935/P2936 demanded for an Aut-breaking prime-coordinate source, then mark
which strict theorem obligations remain missing.

The computation is finite and explicit:
- extend the P2938 prime vector to nodes 1..11 by prime exponents,
- verify product additivity on all audited products d*e<=11,
- compute the U(12) action-defect matrix f(u*d)-f(d),
- extract concrete Aut-breaking witnesses.

Result: the P2938 carrier is nonzero, additive, all-prime nonzero, and
Aut-breaking.  It still is not a strict source theorem because no strict
nadsoliton provenance theorem, delta/eta source, or beta/eta coupling theorem is
exported.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2938 = GEN / "p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate.json"
OUT = GEN / "p2939_s1889_p2938_carrier_aut_breaking_provenance_boundary.json"
MD = GEN / "p2939_s1889_p2938_carrier_aut_breaking_provenance_boundary.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NODES = list(range(1, 12))
PRIMES = [2, 3, 5, 7, 11]
UNITS = [1, 5, 7, 11]


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


def node_values(prime_vector: list[int]) -> dict[int, int]:
    values = {}
    for d in NODES:
        values[d] = sum(exp * coord for exp, coord in zip(factor_vector(d), prime_vector))
    return values


def product_rows(values: dict[int, int]) -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        for e in NODES:
            if d * e <= NODES[-1]:
                defect = values[d * e] - values[d] - values[e]
                rows.append({
                    "d": d,
                    "e": e,
                    "de": d * e,
                    "additive_defect": defect,
                    "passes_additivity": defect == 0,
                })
    return rows


def aut_action_rows(values: dict[int, int]) -> list[dict[str, Any]]:
    rows = []
    for u in UNITS:
        for d in NODES:
            image = (u * d) % MODULUS
            if image == 0:
                raise AssertionError("unit action sent nonzero node to zero")
            defect = values[image] - values[d]
            rows.append({
                "unit": u,
                "node": d,
                "image": image,
                "value_node": values[d],
                "value_image": values[image],
                "aut_invariance_defect": defect,
                "preserves_value": defect == 0,
            })
    return rows


def provenance_acceptance_rows(prime_vector: list[int], prod_rows: list[dict[str, Any]], aut_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "criterion": "nonzero_prime_coordinate_vector",
            "satisfied": any(coord != 0 for coord in prime_vector),
            "evidence": f"prime vector is {prime_vector}",
        },
        {
            "criterion": "all_five_prime_coordinates_nonzero",
            "satisfied": all(coord != 0 for coord in prime_vector),
            "evidence": f"prime vector is {prime_vector}",
        },
        {
            "criterion": "product_additive_on_audited_node_set",
            "satisfied": all(row["passes_additivity"] for row in prod_rows),
            "evidence": "all d*e<=11 product rows have zero additive defect",
        },
        {
            "criterion": "AutZ12_breaking_witness_exists",
            "satisfied": any(not row["preserves_value"] for row in aut_rows if row["unit"] != 1),
            "evidence": "nonidentity unit action changes at least one node value",
        },
        {
            "criterion": "strict_nadsoliton_provenance_theorem",
            "satisfied": False,
            "evidence": "no current artifact proves this carrier is sourced by the strict nadsoliton rather than finite Z12 arithmetic",
        },
        {
            "criterion": "delta_eta_source_law",
            "satisfied": False,
            "evidence": "no current artifact derives delta=4/5 or eta=9/5 from this carrier",
        },
        {
            "criterion": "beta_eta_coupling_theorem",
            "satisfied": False,
            "evidence": "no current artifact couples this carrier to the strict damping denominator",
        },
    ]


def build_payload(p2938: dict[str, Any]) -> dict[str, Any]:
    prime_vector = p2938.get("candidate_certificate", {}).get("prime_coordinate_vector_order_2_3_5_7_11", [1, 2, 2, 2, 2])
    values = node_values(prime_vector)
    prod_rows = product_rows(values)
    aut_rows = aut_action_rows(values)
    nontrivial_defects = [row for row in aut_rows if row["unit"] != 1 and not row["preserves_value"]]
    criteria = provenance_acceptance_rows(prime_vector, prod_rows, aut_rows)
    accepted = all(row["satisfied"] for row in criteria)
    return {
        "status": "P2939_P2938_CARRIER_AUT_BREAKING_PROVENANCE_BOUNDARY_NO_STRICT_SOURCE",
        "input_hashes": {"P2938": hashlib.sha256(P2938.read_bytes()).hexdigest() if P2938.exists() else None},
        "constructed_theoretical_objects": {
            "carrier_under_test": "P2938_UnitCharacter_Enriched_Z12_PrimeCoordinate_Carrier",
            "prime_coordinate_vector_order_2_3_5_7_11": prime_vector,
            "node_values": [{"node": d, "value": values[d], "prime_exponent_vector": factor_vector(d)} for d in NODES],
            "product_additivity_rows": prod_rows,
            "aut_action_defect_rows": aut_rows,
            "provenance_acceptance_rows": criteria,
        },
        "boundary_certificate": {
            "product_pair_count_de_le_11": len(prod_rows),
            "product_additivity_defect_count": sum(1 for row in prod_rows if not row["passes_additivity"]),
            "aut_action_row_count": len(aut_rows),
            "nontrivial_aut_defect_count": len(nontrivial_defects),
            "first_aut_breaking_witness": nontrivial_defects[0] if nontrivial_defects else None,
            "algebraic_readiness_criteria_satisfied_count": sum(1 for row in criteria[:4] if row["satisfied"]),
            "strict_theorem_criteria_satisfied_count": sum(1 for row in criteria[4:] if row["satisfied"]),
            "accepted_strict_source": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "P2938_carrier_nonzero_all_prime": all(coord != 0 for coord in prime_vector),
                "product_additivity_verified": all(row["passes_additivity"] for row in prod_rows),
                "AutZ12_breaking_witness_computed": bool(nontrivial_defects),
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
            "reason": "P2939 verifies that the P2938 carrier satisfies the finite algebraic side of the P2935 target: all five prime coordinates are nonzero, product additivity has zero defects, and a concrete Aut(Z12)-breaking witness exists.  It still fails the strict theorem side: no strict nadsoliton provenance, delta/eta source, or beta/eta coupling theorem is exported.",
            "next_honest_step": "Stop adding finite carrier refinements for this lane.  The next admissible move must be a provenance/coupling theorem for the exact P2938 carrier, or a new typed object outside the prime-coordinate carrier lane; otherwise preserve the P2929-P2939 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["boundary_certificate"]
    lines = [
        "# P2939/S1889 P2938 carrier Aut-breaking provenance boundary",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary certificate",
        f"- product pairs d*e<=11: `{cert['product_pair_count_de_le_11']}`",
        f"- product additivity defects: `{cert['product_additivity_defect_count']}`",
        f"- Aut action rows: `{cert['aut_action_row_count']}`",
        f"- nontrivial Aut defects: `{cert['nontrivial_aut_defect_count']}`",
        f"- first Aut-breaking witness: `{cert['first_aut_breaking_witness']}`",
        f"- algebraic readiness criteria satisfied: `{cert['algebraic_readiness_criteria_satisfied_count']}`",
        f"- strict theorem criteria satisfied: `{cert['strict_theorem_criteria_satisfied_count']}`",
        f"- accepted strict source: `{cert['accepted_strict_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2938))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2939/S1889 P2938 carrier Aut-breaking provenance boundary", "## P2939/S1889 P2938 carrier Aut-breaking provenance boundary\n\n`P2939/S1889` verifies the algebraic side of the exact P2938 carrier `[1,2,2,2,2]`: all five prime coordinates are nonzero, all `29` audited products have zero additivity defect, and the `U(12)` action-defect matrix gives concrete Aut-breaking witnesses such as a nonidentity unit changing a node value.  This still does not export a strict source: strict nadsoliton provenance, delta/eta source, and beta/eta coupling remain absent, so no strict `L_p`, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2939/S1889 P2938 carrier provenance `L_total` guard", "## P2939/S1889 P2938 carrier provenance `L_total` guard\n\n`P2939/S1889` confirms that the P2938 carrier is algebraically ready and Aut-breaking, but not provenance/coupling complete.  Without a strict nadsoliton provenance theorem plus delta/eta and beta/eta coupling, the carrier cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current P2938 carrier provenance-boundary guardrail (P2939/S1889, 2026-06-19)", "## Current P2938 carrier provenance-boundary guardrail (P2939/S1889, 2026-06-19)\n\n- P2939 verifies the P2938 carrier `[1,2,2,2,2]` satisfies the finite algebraic side: all prime coordinates are nonzero, product additivity has zero defects on `29` products, and Aut(Z12)-breaking witnesses exist.\n- P2939 does not export strict nadsoliton provenance, delta/eta source, beta/eta coupling, strict `L_p`, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- Do not add more finite carrier refinements as primary strategy in this lane; the next admissible move is a provenance/coupling theorem for the exact P2938 carrier or a genuinely new typed object outside the prime-coordinate carrier lane.\n- Without such a theorem/object, preserve the P2929-P2939 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
