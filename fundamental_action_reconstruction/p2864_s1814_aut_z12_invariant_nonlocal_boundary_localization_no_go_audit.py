#!/usr/bin/env python3
"""P2864/S1814: Aut(Z12)-invariant nonlocal boundary localization no-go audit.

P2863 left only a sharper possibility for the P2862 right boundary datum:
perhaps a nonlocal unit-bearing coupling/localization theorem computes
(9/5)*log(11) rather than importing it.  This packet tests the strictest finite
Z12-compatible version of that idea: Aut(Z12)-invariant log-localization weights
on the nonzero Z12 residue/distance classes.

Aut(Z12) has unit orbits {1,5,7,11}, {2,10}, {3,9}, {4,8}, {6}.  Any
Aut-invariant log-localization gives one coefficient to the whole unit orbit
{1,5,7,11}.  Therefore the coefficient of log(11) is tied to coefficients of
log(5) and log(7).  The target boundary vector has prime-11 exponent 9/5 and
zero prime-5/prime-7 exponents, so the exact linear system is inconsistent.
A selector-localized singleton d=11 can represent the target only by importing
both a non-premise selector/localizer and the prime-5 coefficient 9/5.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN, fraction_payload, prime_factors

P2863 = GEN / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.json"
OUT = GEN / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.json"
MD = GEN / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
UNITS = [1, 5, 7, 11]
PRIMES = [2, 3, 5, 7, 11]
TARGET_COEFFICIENT = Fraction(9, 5)
TARGET_VECTOR = {p: Fraction(0, 1) for p in PRIMES} | {11: TARGET_COEFFICIENT}


def factor_vector(n: int) -> dict[int, int]:
    remaining = n
    out = {p: 0 for p in PRIMES}
    for p in PRIMES:
        while remaining % p == 0:
            out[p] += 1
            remaining //= p
    if remaining != 1:
        raise ValueError(f"unexpected factor {remaining} for {n}")
    return out


def aut_orbits() -> list[list[int]]:
    unseen = set(range(1, MODULUS))
    out = []
    while unseen:
        seed = min(unseen)
        orbit = sorted({(unit * seed) % MODULUS for unit in UNITS} - {0})
        out.append(orbit)
        unseen -= set(orbit)
    return out


ORBITS = aut_orbits()


def orbit_prime_vector(orbit: list[int]) -> dict[int, Fraction]:
    vector = {p: Fraction(0, 1) for p in PRIMES}
    for d in orbit:
        for p, exponent in factor_vector(d).items():
            vector[p] += Fraction(exponent, 1)
    return vector


def matrix_rows() -> list[dict[str, Any]]:
    orbit_vectors = [orbit_prime_vector(orbit) for orbit in ORBITS]
    rows = []
    for p in PRIMES:
        rows.append(
            {
                "prime": p,
                "coefficients_by_orbit": [fraction_payload(vector[p]) for vector in orbit_vectors],
                "target": fraction_payload(TARGET_VECTOR[p]),
            }
        )
    return rows


def rank_fraction_matrix(matrix: list[list[Fraction]]) -> int:
    rows = [row[:] for row in matrix]
    rank = 0
    if not rows:
        return 0
    col_count = len(rows[0])
    for col in range(col_count):
        pivot = None
        for r in range(rank, len(rows)):
            if rows[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            continue
        rows[rank], rows[pivot] = rows[pivot], rows[rank]
        pivot_value = rows[rank][col]
        rows[rank] = [value / pivot_value for value in rows[rank]]
        for r in range(len(rows)):
            if r != rank and rows[r][col] != 0:
                factor = rows[r][col]
                rows[r] = [value - factor * pivot_row_value for value, pivot_row_value in zip(rows[r], rows[rank])]
        rank += 1
        if rank == len(rows):
            break
    return rank


def linear_system_certificate() -> dict[str, Any]:
    orbit_vectors = [orbit_prime_vector(orbit) for orbit in ORBITS]
    matrix = [[vector[p] for vector in orbit_vectors] for p in PRIMES]
    target = [TARGET_VECTOR[p] for p in PRIMES]
    augmented = [row + [rhs] for row, rhs in zip(matrix, target)]
    rank_a = rank_fraction_matrix(matrix)
    rank_aug = rank_fraction_matrix(augmented)
    unit_orbit_index = next(index for index, orbit in enumerate(ORBITS) if 11 in orbit)
    unit_orbit_vector = orbit_vectors[unit_orbit_index]
    return {
        "aut_orbits": ORBITS,
        "orbit_prime_vectors": [
            {str(p): fraction_payload(value) for p, value in vector.items()} for vector in orbit_vectors
        ],
        "matrix_rows": matrix_rows(),
        "rank_coefficient_matrix": rank_a,
        "rank_augmented_matrix": rank_aug,
        "linear_system_consistent": rank_a == rank_aug,
        "unit_orbit_index": unit_orbit_index,
        "unit_orbit": ORBITS[unit_orbit_index],
        "unit_orbit_prime5_coefficient": fraction_payload(unit_orbit_vector[5]),
        "unit_orbit_prime7_coefficient": fraction_payload(unit_orbit_vector[7]),
        "unit_orbit_prime11_coefficient": fraction_payload(unit_orbit_vector[11]),
        "contradiction": "Aut-invariance ties log(5), log(7), and log(11) in orbit {1,5,7,11}; target requires coefficient 9/5 on log(11) and 0 on log(5), log(7).",
    }


def selector_localized_singleton_certificate() -> dict[str, Any]:
    return {
        "candidate": "selector_localized_singleton_d11",
        "singleton": 11,
        "required_coefficient": fraction_payload(TARGET_COEFFICIENT),
        "requires_nonpremise_selector_or_localizer": True,
        "requires_prime5_coefficient": 5 in prime_factors(TARGET_COEFFICIENT.denominator),
        "represents_target_if_both_imported": True,
        "exports_boundary_source_law": False,
        "verdict": "A d=11 singleton localizer with coefficient 9/5 represents the target, but imports both the selector/localizer and the prime-5 coefficient.",
    }


def build_payload(p2863: dict[str, Any]) -> dict[str, Any]:
    system = linear_system_certificate()
    singleton = selector_localized_singleton_certificate()
    facts = {
        "p2863_rechecked": p2863.get("status") == "P2863_DIRICHLET_BOUNDARY_DATUM_SOURCE_CLASS_NO_GO_AUDIT_NO_CLOSURE",
        "aut_orbits_computed": system["aut_orbits"] == [[1, 5, 7, 11], [2, 10], [3, 9], [4, 8], [6]],
        "aut_invariant_linear_system_inconsistent": not system["linear_system_consistent"],
        "unit_orbit_ties_5_7_11": system["unit_orbit"] == [1, 5, 7, 11],
        "selector_singleton_requires_imports": singleton["requires_nonpremise_selector_or_localizer"] and singleton["requires_prime5_coefficient"],
        "accepted_count_zero": True,
    }
    return {
        "status": "P2864_AUT_Z12_INVARIANT_NONLOCAL_BOUNDARY_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2863": sha(P2863)},
        "aut_z12_invariant_nonlocal_boundary_localization_no_go_audit": {
            "input_status_rechecked": p2863.get("status"),
            "target_boundary_datum": "log(11^(9/5)) = (9/5)*log(11)",
            "aut_units": UNITS,
            "prime_basis": PRIMES,
            "target_prime_exponent_vector": {str(p): fraction_payload(value) for p, value in TARGET_VECTOR.items()},
            "linear_system_certificate": system,
            "selector_localized_singleton_certificate": singleton,
            "candidate_matrix": [
                {
                    "candidate": "aut_z12_invariant_nonlocal_log_localization",
                    "finite_witness_passes": not system["linear_system_consistent"],
                    "exports_boundary_source_law": False,
                    "exports_eta_source_law": False,
                    "verdict": "exact Aut-invariant localization is inconsistent with isolating only log(11).",
                },
                {
                    "candidate": "selector_localized_singleton_d11",
                    "finite_witness_passes": singleton["represents_target_if_both_imported"],
                    "exports_boundary_source_law": False,
                    "exports_eta_source_law": False,
                    "verdict": singleton["verdict"],
                },
            ],
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "orbit_step": "Aut(Z12) unit orbits on nonzero residues are {1,5,7,11}, {2,10}, {3,9}, {4,8}, and {6}.",
                "rank_step": "The exact prime-exponent linear system has rank(A) < rank([A|b]), so no Aut-invariant orbit weighting solves the target.",
                "selector_step": "Isolating d=11 breaks Aut-invariance and requires a selector/localizer plus the prime-5 coefficient 9/5.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit": all(facts.values()),
            "exports_boundary_source_law": False,
            "exports_eta_source_law": False,
            "exports_selector_closure": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "eta_source_exported": False,
                "prime5_source_exported": False,
                "selector_closure_exported": False,
                "unit_bearing_coupling_localization_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2864 tests the nonlocal localization route left open by P2863.  Aut(Z12)-invariant orbit weights cannot isolate log(11), because the unit orbit ties 11 to 5 and 7 and the exact prime-exponent linear system is inconsistent.  A singleton d=11 localizer represents the target only by importing a selector/localizer and the prime-5 coefficient 9/5.  No boundary source law, selector closure, or eta source is exported.",
            "next_honest_step": "Do not replay Aut-invariant nonlocal log-localization, singleton endpoint localization, selector import, or prime-5 representability as boundary sourcehood.  The next proof-grade move must provide a genuine strict selector/localizer source plus the 9/5 coefficient law, a different unit-bearing coupling theorem not orbit-invariant in this exhausted way, or a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["aut_z12_invariant_nonlocal_boundary_localization_no_go_audit"]
    system = audit["linear_system_certificate"]
    lines = [
        "# P2864/S1814 Aut(Z12)-invariant nonlocal boundary localization no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Nonlocal localization scan",
        f"- Aut(Z12) orbits: `{system['aut_orbits']}`",
        f"- rank(A): `{system['rank_coefficient_matrix']}`",
        f"- rank([A|b]): `{system['rank_augmented_matrix']}`",
        f"- linear system consistent: `{system['linear_system_consistent']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2863))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2864/S1814 Aut(Z12)-invariant nonlocal boundary localization no-go audit",
        "## P2864/S1814 Aut(Z12)-invariant nonlocal boundary localization no-go audit\n\n"
        "`P2864/S1814` tests the nonlocal localization route left after `P2863`: Aut(`Z12`)-invariant log-localization weights over nonzero residue/distance orbits.  The unit orbit `{1,5,7,11}` ties coefficients of `log(5)`, `log(7)`, and `log(11)`, while the target boundary datum requires only `(9/5)log(11)`.  The exact prime-exponent linear system is inconsistent.  A singleton `d=11` localizer represents the target only by importing a selector/localizer and the prime-5 coefficient.  No boundary source law, eta source, selector closure, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2864/S1814 Aut-invariant localization `L_total` guard",
        "## P2864/S1814 Aut-invariant localization `L_total` guard\n\n"
        "`P2864/S1814` adds no strict action term.  Aut-invariant nonlocal log-localization cannot source the endpoint datum, and singleton localization imports a selector plus prime-5 coefficient; neither provides a unit-bearing boundary/source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current Aut(Z12)-invariant boundary localization no-go guardrail (P2864/S1814, 2026-06-18)",
        "## Current Aut(Z12)-invariant boundary localization no-go guardrail (P2864/S1814, 2026-06-18)\n\n"
        "- P2864 tests Aut(`Z12`)-invariant nonlocal log-localization for the P2862/P2863 right boundary datum.\n"
        "- The unit orbit `{1,5,7,11}` prevents isolating `log(11)` without also sourcing unwanted `log(5)` and `log(7)` components; the exact prime-exponent linear system is inconsistent.\n"
        "- Singleton `d=11` localization represents the target only by importing a selector/localizer and the prime-5 coefficient `9/5`.\n"
        "- Do not promote Aut-invariant nonlocal localization, singleton endpoint localization, selector import, prime-5 representability, Dirichlet data, or log-scale harmonicity to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must provide a genuine strict selector/localizer source plus the `9/5` coefficient law, a different unit-bearing coupling theorem outside this exhausted orbit-invariant class, or a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
