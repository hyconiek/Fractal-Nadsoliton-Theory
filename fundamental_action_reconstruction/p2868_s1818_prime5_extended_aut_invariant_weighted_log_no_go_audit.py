#!/usr/bin/env python3
"""P2868/S1818: prime-5-extended Aut-invariant weighted-log no-go audit.

P2867 showed that a coupled Z12-smooth weighted-log functional can match
(9/5)*log(11) only by importing the denominator-prime-5 weight w_11=9/5.
This packet grants the prime-5 coefficient ring as a deliberately stronger
candidate class and asks whether Aut(Z12)-invariant localization can then source
the endpoint without a singleton selector.

The answer is still no.  Aut(Z12) puts 1, 5, 7, and 11 in the same unit orbit.
Any orbit-invariant weighted-log functional assigns one shared coefficient to
log(5), log(7), and log(11).  The target requires coefficient 9/5 on log(11)
and coefficient 0 on log(5) and log(7), so the exact linear system is
inconsistent even after denominator prime 5 is allowed.
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
from p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit import PRIMES, TARGET_VECTOR, factor_vector

P2867 = GEN / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.json"
OUT = GEN / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.json"
MD = GEN / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS_Z12 = [1, 5, 7, 11]
EXTENDED_DENOMINATOR_PRIMES = {2, 3, 5}
BOUNDED_DENOMINATORS = [1, 2, 3, 4, 5, 6, 8, 9, 10, 12, 15]
BOUNDED_NUMERATOR_ABS_MAX = 24


def aut_orbits_nonzero_distances() -> list[tuple[int, ...]]:
    unseen = set(range(1, 12))
    orbits: list[tuple[int, ...]] = []
    while unseen:
        seed = min(unseen)
        orbit = tuple(sorted({(unit * seed) % 12 or 12 for unit in UNITS_Z12 if (unit * seed) % 12 != 0}))
        orbits.append(orbit)
        unseen -= set(orbit)
    return orbits


def orbit_prime_vector(orbit: tuple[int, ...]) -> dict[int, Fraction]:
    out = {p: Fraction(0, 1) for p in PRIMES}
    for d in orbit:
        vector = factor_vector(d)
        for p in PRIMES:
            out[p] += vector[p]
    return out


def matrix_rows(orbits: list[tuple[int, ...]]) -> list[list[Fraction]]:
    return [[orbit_prime_vector(orbit)[p] for orbit in orbits] for p in PRIMES]


def rank(matrix: list[list[Fraction]]) -> int:
    work = [row[:] for row in matrix]
    rows = len(work)
    cols = len(work[0]) if rows else 0
    r = 0
    for c in range(cols):
        pivot = next((i for i in range(r, rows) if work[i][c] != 0), None)
        if pivot is None:
            continue
        work[r], work[pivot] = work[pivot], work[r]
        pivot_value = work[r][c]
        work[r] = [value / pivot_value for value in work[r]]
        for i in range(rows):
            if i != r and work[i][c] != 0:
                factor = work[i][c]
                work[i] = [value - factor * base for value, base in zip(work[i], work[r])]
        r += 1
        if r == rows:
            break
    return r


def linear_system_certificate() -> dict[str, Any]:
    orbits = aut_orbits_nonzero_distances()
    a = matrix_rows(orbits)
    b = [TARGET_VECTOR[p] for p in PRIMES]
    augmented = [row + [rhs] for row, rhs in zip(a, b)]
    return {
        "orbits": [list(orbit) for orbit in orbits],
        "orbit_prime_vectors": {
            str(list(orbit)): {str(p): fraction_payload(value) for p, value in orbit_prime_vector(orbit).items()}
            for orbit in orbits
        },
        "rank_A": rank(a),
        "rank_augmented": rank(augmented),
        "consistent": rank(a) == rank(augmented),
        "unit_orbit_obstruction": "The orbit [1,5,7,11] contributes equally to prime coordinates 5, 7, and 11, but the target asks for 0, 0, and 9/5 respectively.",
    }


def allowed_extended_weights() -> list[Fraction]:
    values = set()
    for denominator in BOUNDED_DENOMINATORS:
        if not set(prime_factors(denominator)).issubset(EXTENDED_DENOMINATOR_PRIMES):
            continue
        for numerator in range(-BOUNDED_NUMERATOR_ABS_MAX, BOUNDED_NUMERATOR_ABS_MAX + 1):
            values.add(Fraction(numerator, denominator))
    return sorted(values)


def unit_orbit_bounded_scan() -> dict[str, Any]:
    values = allowed_extended_weights()
    exact_conditions = [value for value in values if value == 0 and value == Fraction(9, 5)]
    best = min(values, key=lambda value: abs(value - 0) + abs(value - Fraction(9, 5)))
    return {
        "denominators_scanned": BOUNDED_DENOMINATORS,
        "numerator_abs_max": BOUNDED_NUMERATOR_ABS_MAX,
        "candidate_count": len(values),
        "exact_shared_unit_orbit_weight_hits": [fraction_payload(value) for value in exact_conditions],
        "best_shared_unit_orbit_weight": fraction_payload(best),
        "best_l1_error_for_requirements_a_eq_0_and_a_eq_9_over_5": fraction_payload(abs(best) + abs(best - Fraction(9, 5))),
        "exact_shared_weight_absent": not exact_conditions,
    }


def build_payload(p2867: dict[str, Any]) -> dict[str, Any]:
    system = linear_system_certificate()
    bounded = unit_orbit_bounded_scan()
    facts = {
        "p2867_rechecked": p2867.get("status") == "P2867_COUPLED_Z12_SMOOTH_WEIGHTED_LOG_FUNCTIONAL_NO_SOURCE_AUDIT_NO_CLOSURE",
        "prime5_extended_denominators_allowed": 5 in EXTENDED_DENOMINATOR_PRIMES,
        "aut_invariant_linear_system_inconsistent": not system["consistent"] and system["rank_augmented"] > system["rank_A"],
        "bounded_scan_has_no_shared_unit_orbit_weight": bounded["exact_shared_weight_absent"],
        "accepted_count_zero": True,
    }
    return {
        "status": "P2868_PRIME5_EXTENDED_AUT_INVARIANT_WEIGHTED_LOG_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2867": sha(P2867)},
        "prime5_extended_aut_invariant_weighted_log_no_go_audit": {
            "input_status_rechecked": p2867.get("status"),
            "candidate_class": "Aut(Z12)-invariant orbit weighted log functional with denominator support {2,3,5}",
            "target_boundary_datum": "log(11^(9/5)) = (9/5)*log(11)",
            "denominator_primes_allowed": sorted(EXTENDED_DENOMINATOR_PRIMES),
            "linear_system_certificate": system,
            "unit_orbit_bounded_scan": bounded,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "prime5_grant_step": "Unlike P2867, this candidate class allows denominator prime 5 in coefficients.",
                "aut_orbit_step": "Aut(Z12) still ties d=11 to d=1,5,7 in one unit orbit.",
                "rank_step": "The exact orbit-weight linear system has rank(A) < rank([A|b]), so allowing prime-5 denominators does not restore consistency.",
                "sourcehood_step": "The route still needs a non-Aut-invariant strict localizer/selector for d=11, not merely prime-5 coefficient capacity.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_prime5_extended_aut_invariant_weighted_log_no_go_audit": all(facts.values()),
            "exports_boundary_source_law": False,
            "exports_eta_source_law": False,
            "exports_prime5_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "eta_source_exported": False,
                "prime5_source_exported": False,
                "selector_or_localizer_source_exported": False,
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
            "reason": "P2868 grants denominator prime 5 but keeps Aut(Z12)-invariant localization.  The exact orbit-weight system remains inconsistent because the unit orbit ties log(11) to log(5) and log(7), while the target requires only log(11).  Prime-5 coefficient capacity alone is not a d=11 source.",
            "next_honest_step": "Do not replay prime-5 denominator extension under Aut-invariant localization as sourcehood.  A next proof-grade move must introduce a genuine non-premise strict localizer/selector breaking the unit orbit together with a unit-bearing coefficient law, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["prime5_extended_aut_invariant_weighted_log_no_go_audit"]
    lines = [
        "# P2868/S1818 prime-5-extended Aut-invariant weighted-log no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Prime-5-extended Aut-invariant scan",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- denominator primes allowed: `{audit['denominator_primes_allowed']}`",
        f"- rank certificate: `{audit['linear_system_certificate']}`",
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
    payload = build_payload(read_json(P2867))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2868/S1818 prime-5-extended Aut-invariant weighted-log no-go audit",
        "## P2868/S1818 prime-5-extended Aut-invariant weighted-log no-go audit\n\n"
        "`P2868/S1818` grants denominator prime `5` in the weighted-log coefficient ring but keeps Aut(`Z12`)-invariant localization.  The exact orbit-weight system remains inconsistent: the unit orbit `{1,5,7,11}` assigns one shared coefficient to `log(5)`, `log(7)`, and `log(11)`, while the target requires only `(9/5)log(11)`.  Thus prime-5 coefficient capacity alone is not a strict `d=11` localizer/source and does not export a boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2868/S1818 prime-5 Aut-invariant weighted-log `L_total` guard",
        "## P2868/S1818 prime-5 Aut-invariant weighted-log `L_total` guard\n\n"
        "`P2868/S1818` adds no strict action term.  Allowing denominator prime `5` inside an Aut-invariant weighted-log functional does not provide a unit-bearing boundary/source density, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current prime-5-extended Aut-invariant weighted-log no-go guardrail (P2868/S1818, 2026-06-18)",
        "## Current prime-5-extended Aut-invariant weighted-log no-go guardrail (P2868/S1818, 2026-06-18)\n\n"
        "- P2868 grants denominator prime `5` in the weighted-log coefficient ring but retains Aut(`Z12`)-invariant localization.\n"
        "- The exact orbit-weight system is still inconsistent because the unit orbit `{1,5,7,11}` ties `log(11)` to `log(5)` and `log(7)`, while the target requires only `(9/5)log(11)`.\n"
        "- Do not promote prime-5 denominator extension, Aut-invariant orbit weighting, coupled weighted-log functionals, singleton endpoint representation, Dirichlet data, log-scale harmonicity, or pairwise recombination to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce a genuine non-premise strict localizer/selector breaking the unit orbit together with a unit-bearing coefficient law, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
