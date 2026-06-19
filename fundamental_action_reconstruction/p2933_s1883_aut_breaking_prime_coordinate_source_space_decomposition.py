#!/usr/bin/env python3
"""P2933/S1883: Aut-breaking prime-coordinate source-space decomposition.

P2932 proved that Aut(Z12)-invariant additive node laws are zero-only.  P2933
identifies what a future nonzero strict L_p source must therefore be: an explicit
symmetry-breaking vector in the five-dimensional prime-coordinate additive
character space.

The computation is exact.  Product additivity on nodes 1..11 has rank 6 and
nullity 5; a canonical basis is given by the five prime coordinates
Y_2,Y_3,Y_5,Y_7,Y_11.  Adding Aut(Z12)-invariance raises the rank to 11 and
kills the nullity.  Thus every nonzero additive prime-log value source lies in a
five-dimensional Aut-breaking quotient and needs a strict source/selector law;
readiness, naming, or imported numeric logs are not enough.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix import NODES, PRIMES, factor_vector
from p2932_s1882_aut_z12_invariant_additive_value_source_no_go import UNITS, rref

P2932 = GEN / "p2932_s1882_aut_z12_invariant_additive_value_source_no_go.json"
OUT = GEN / "p2933_s1883_aut_breaking_prime_coordinate_source_space_decomposition.json"
MD = GEN / "p2933_s1883_aut_breaking_prime_coordinate_source_space_decomposition.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def zero_row() -> list[int]:
    return [0 for _ in NODES]


def product_equations() -> list[list[int]]:
    rows: list[list[int]] = []
    row = zero_row()
    row[0] = 1
    rows.append(row)
    for d in NODES:
        for e in NODES:
            product = d * e
            if product <= 11:
                row = zero_row()
                row[product - 1] += 1
                row[d - 1] -= 1
                row[e - 1] -= 1
                rows.append(row)
    return rows


def aut_equations() -> list[list[int]]:
    rows: list[list[int]] = []
    for unit in UNITS:
        for d in NODES:
            image = (unit * d) % 12
            if image != 0:
                row = zero_row()
                row[image - 1] += 1
                row[d - 1] -= 1
                rows.append(row)
    return rows


def prime_basis_rows() -> list[dict[str, Any]]:
    rows = []
    for basis_prime in PRIMES:
        values: dict[str, int] = {}
        for node in NODES:
            vector = factor_vector(node)
            values[f"f_{node}"] = vector[PRIMES.index(basis_prime)]
        rows.append({
            "basis_name": f"Y_{basis_prime}",
            "free_prime_coordinate": f"L_{basis_prime}=1, other L_p=0",
            "node_values": values,
            "nonzero": any(value != 0 for value in values.values()),
        })
    return rows


def aut_breaking_witnesses() -> list[dict[str, Any]]:
    basis = prime_basis_rows()
    witnesses = []
    for row in basis:
        values = {int(key.split("_")[1]): value for key, value in row["node_values"].items()}
        defects = []
        for unit in UNITS:
            for node in NODES:
                image = (unit * node) % 12
                if image != 0:
                    defect = values[image] - values[node]
                    if defect != 0:
                        defects.append({"unit": unit, "node": node, "image": image, "defect": defect})
        witnesses.append({
            "basis_name": row["basis_name"],
            "aut_invariance_defect_count": len(defects),
            "first_defect": defects[0] if defects else None,
            "requires_symmetry_breaking_source": len(defects) > 0,
        })
    return witnesses


def candidate_source_rows() -> list[dict[str, Any]]:
    return [
        {
            "candidate": "arbitrary_prime_coordinate_vector",
            "nonzero_additive_vector": True,
            "strict_symmetry_breaking_source_law": False,
            "accepted_as_strict_Lp_source": False,
            "failure": "chooses coordinates without strict source law",
        },
        {
            "candidate": "external_real_log_values",
            "nonzero_additive_vector": True,
            "strict_symmetry_breaking_source_law": False,
            "accepted_as_strict_Lp_source": False,
            "failure": "imports analytic logs rather than deriving values from strict nadsoliton data",
        },
        {
            "candidate": "P2930_residue_distance_vector",
            "nonzero_additive_vector": True,
            "strict_symmetry_breaking_source_law": False,
            "accepted_as_strict_Lp_source": False,
            "failure": "finite labels exist but P2931/P2932 do not export a strict source theorem",
        },
        {
            "candidate": "future_Strict_AutBreaking_PrimeCoordinate_Source_Law",
            "nonzero_additive_vector": None,
            "strict_symmetry_breaking_source_law": None,
            "accepted_as_strict_Lp_source": False,
            "failure": "schema for required future object, not a present artifact",
        },
    ]


def build_payload(p2932: dict[str, Any]) -> dict[str, Any]:
    product = product_equations()
    aut = aut_equations()
    product_rref = rref(product)
    combined_rref = rref(product + aut)
    basis = prime_basis_rows()
    witnesses = aut_breaking_witnesses()
    candidates = candidate_source_rows()
    return {
        "status": "P2933_AUT_BREAKING_PRIME_COORDINATE_SOURCE_SPACE_DECOMPOSITION_NO_SOURCE_LAW",
        "input_hashes": {"P2932": hashlib.sha256(P2932.read_bytes()).hexdigest() if P2932.exists() else None},
        "constructed_theoretical_objects": {
            "source_space": {
                "name": "AutBreaking_PrimeCoordinate_Source_Space",
                "product_additive_rank": product_rref["rank"],
                "product_additive_nullity": product_rref["nullity"],
                "aut_invariant_rank": combined_rref["rank"],
                "aut_invariant_nullity": combined_rref["nullity"],
                "aut_breaking_quotient_dimension": product_rref["nullity"] - combined_rref["nullity"],
            },
            "prime_coordinate_basis_rows": basis,
            "aut_breaking_witnesses": witnesses,
            "candidate_source_rows": candidates,
        },
        "decomposition_certificate": {
            "product_equation_count": len(product),
            "aut_equation_count": len(aut),
            "product_rank": product_rref["rank"],
            "product_nullity": product_rref["nullity"],
            "combined_rank": combined_rref["rank"],
            "combined_nullity": combined_rref["nullity"],
            "prime_coordinate_basis_count": len(basis),
            "all_basis_vectors_break_aut_invariance": all(row["requires_symmetry_breaking_source"] for row in witnesses),
            "accepted_candidate_count": sum(1 for row in candidates if row["accepted_as_strict_Lp_source"]),
        },
        "decision": {
            "positive_witnesses": {
                "five_dimensional_additive_prime_coordinate_space_identified": product_rref["nullity"] == 5,
                "aut_invariant_subspace_zero_confirmed": combined_rref["nullity"] == 0,
                "required_future_object_typed": True,
            },
            "negative_export_flags": {
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
            "reason": "P2933 decomposes the exact additive value-source space after P2932.  Product additivity leaves five prime-coordinate degrees of freedom, while Aut(Z12)-invariance kills all five.  Therefore any nonzero L_p source must be an explicit strict Aut-breaking value law; current coordinate choices, external logs, and residue labels are not such a law.",
            "next_honest_step": "Construct a concrete Strict_AutBreaking_PrimeCoordinate_Source_Law that derives one nonzero vector in the five-dimensional prime-coordinate space from strict nadsoliton data and then run finite additivity, symmetry-breaking provenance, delta/eta, and beta/eta coupling tests.  If no such law is supplied, preserve the P2929-P2933 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["decomposition_certificate"]
    lines = [
        "# P2933/S1883 Aut-breaking prime-coordinate source-space decomposition",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Decomposition certificate",
        f"- product equations: `{cert['product_equation_count']}`",
        f"- product rank/nullity: `{cert['product_rank']}` / `{cert['product_nullity']}`",
        f"- Aut equations: `{cert['aut_equation_count']}`",
        f"- combined rank/nullity: `{cert['combined_rank']}` / `{cert['combined_nullity']}`",
        f"- prime-coordinate basis count: `{cert['prime_coordinate_basis_count']}`",
        f"- all basis vectors break Aut invariance: `{cert['all_basis_vectors_break_aut_invariance']}`",
        f"- accepted candidates: `{cert['accepted_candidate_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2932))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2933/S1883 Aut-breaking prime-coordinate source-space decomposition", "## P2933/S1883 Aut-breaking prime-coordinate source-space decomposition\n\n`P2933/S1883` decomposes the exact additive `L_p` source space after P2932.  Product additivity on the audited node set has rank `6` and nullity `5`, with basis coordinates `Y_2,Y_3,Y_5,Y_7,Y_11`; adding Aut(Z12)-invariance raises the rank to `11` and nullity to `0`.  Thus every nonzero additive prime-log value source is necessarily Aut-breaking and requires an explicit strict source/selector law.  Current arbitrary coordinates, external logs, and residue labels are rejected; no strict `L_p` source, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2933/S1883 Aut-breaking prime-coordinate `L_total` guard", "## P2933/S1883 Aut-breaking prime-coordinate `L_total` guard\n\n`P2933/S1883` identifies the exact remaining shape of a possible nonzero `L_p` source: it must choose a vector in the five-dimensional prime-coordinate space and justify that Aut-breaking choice from strict nadsoliton data.  Until such a strict source law and damping coupling are exported, the prime-coordinate space cannot enter role-bearing nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Aut-breaking prime-coordinate source-space guardrail (P2933/S1883, 2026-06-19)", "## Current Aut-breaking prime-coordinate source-space guardrail (P2933/S1883, 2026-06-19)\n\n- P2933 decomposes the exact additive prime-log value-source space: product additivity leaves five prime-coordinate degrees of freedom, while Aut(Z12)-invariance kills all five.\n- Any nonzero `L_p` source must therefore be an explicit strict Aut-breaking value law, not arbitrary coordinates, external logs, residue labels, or readiness replay.\n- Do not promote the five-dimensional coordinate space to strict damping `beta/eta`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE without a sourced symmetry-breaking law and coupling theorem.\n- The next admissible move is a concrete `Strict_AutBreaking_PrimeCoordinate_Source_Law` with finite additivity, provenance, delta/eta, and beta/eta coupling tests, or a pivot to a different genuinely new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
