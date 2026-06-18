#!/usr/bin/env python3
"""P2867/S1817: coupled Z12-smooth weighted log-functional no-source audit.

P2866 ruled out recombining existing localizer/coefficient witnesses as a source
for the right boundary datum (9/5)*log(11).  This packet tests a more coupled
candidate class rather than another product replay: a single weighted log
functional

    F(w) = sum_{d=1}^{11} w_d log(d)

with all weights drawn from the Z12-smooth rational coefficient ring whose
denominators use only primes {2,3}.  This allows localization and coefficient
choice to be solved jointly, so it is a stricter escape hatch than selecting a
singleton localizer and coefficient separately.

Exact prime-exponent bookkeeping gives an immediate obstruction: only log(11)
contains prime 11, so the prime-11 coordinate forces w_11 = 9/5.  That forced
weight has denominator prime 5, hence it is outside the pure Z12-smooth ring.
A bounded finite scan of allowed weights confirms the same obstruction.
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

P2866 = GEN / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.json"
OUT = GEN / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.json"
MD = GEN / "p2867_s1817_coupled_z12_smooth_weighted_log_functional_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

Z12_SMOOTH_DENOMINATOR_PRIMES = {2, 3}
BOUNDED_DENOMINATORS = [1, 2, 3, 4, 6, 8, 9, 12]
BOUNDED_NUMERATOR_ABS_MAX = 24
TARGET_WEIGHT_11 = Fraction(9, 5)


def is_z12_smooth_rational(value: Fraction) -> bool:
    return set(prime_factors(value.denominator)).issubset(Z12_SMOOTH_DENOMINATOR_PRIMES)


def allowed_bounded_weights() -> list[Fraction]:
    values = set()
    for denominator in BOUNDED_DENOMINATORS:
        for numerator in range(-BOUNDED_NUMERATOR_ABS_MAX, BOUNDED_NUMERATOR_ABS_MAX + 1):
            values.add(Fraction(numerator, denominator))
    return sorted(values)


def weighted_prime_vector(weights: dict[int, Fraction]) -> dict[int, Fraction]:
    out = {p: Fraction(0, 1) for p in PRIMES}
    for d, weight in weights.items():
        vector = factor_vector(d)
        for p in PRIMES:
            out[p] += weight * vector[p]
    return out


def coordinate_forcing_certificate() -> dict[str, Any]:
    contributors = [d for d in range(1, 12) if factor_vector(d)[11] != 0]
    forced = TARGET_VECTOR[11]
    return {
        "prime11_contributors": contributors,
        "forced_w_11": fraction_payload(forced),
        "forced_w_11_denominator_primes": sorted(prime_factors(forced.denominator)),
        "forced_w_11_is_z12_smooth": is_z12_smooth_rational(forced),
        "obstruction": "Since only d=11 contributes to the prime-11 coordinate, any exact weighted log functional must set w_11=9/5, which imports denominator prime 5 outside the Z12-smooth denominator support {2,3}.",
    }


def bounded_w11_scan() -> dict[str, Any]:
    allowed = allowed_bounded_weights()
    exact = [value for value in allowed if value == TARGET_WEIGHT_11]
    best = min(allowed, key=lambda value: abs(value - TARGET_WEIGHT_11))
    return {
        "denominators_scanned": BOUNDED_DENOMINATORS,
        "numerator_abs_max": BOUNDED_NUMERATOR_ABS_MAX,
        "candidate_count": len(allowed),
        "exact_w11_hits": [fraction_payload(value) for value in exact],
        "best_w11": fraction_payload(best),
        "best_prime11_coordinate_error": fraction_payload(best - TARGET_WEIGHT_11),
        "exact_w11_absent": not exact,
    }


def build_payload(p2866: dict[str, Any]) -> dict[str, Any]:
    coordinate = coordinate_forcing_certificate()
    bounded = bounded_w11_scan()
    witness_weights = {d: Fraction(0, 1) for d in range(1, 12)}
    witness_weights[11] = TARGET_WEIGHT_11
    witness_vector = weighted_prime_vector(witness_weights)
    facts = {
        "p2866_rechecked": p2866.get("status") == "P2866_PRODUCT_THEOREM_CANDIDATE_PAIR_NO_SOURCE_AUDIT_NO_CLOSURE",
        "prime11_coordinate_forces_w11_9_over_5": coordinate["forced_w_11"]["fraction"] == "9/5",
        "forced_w11_not_z12_smooth": not coordinate["forced_w_11_is_z12_smooth"],
        "bounded_scan_has_no_exact_w11": bounded["exact_w11_absent"],
        "imported_witness_exact_but_not_sourced": witness_vector == TARGET_VECTOR and not is_z12_smooth_rational(TARGET_WEIGHT_11),
    }
    return {
        "status": "P2867_COUPLED_Z12_SMOOTH_WEIGHTED_LOG_FUNCTIONAL_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2866": sha(P2866)},
        "coupled_z12_smooth_weighted_log_functional_no_source_audit": {
            "input_status_rechecked": p2866.get("status"),
            "candidate_class": "F(w)=sum_{d=1}^{11} w_d log(d), with every w_d in Q whose denominator support is contained in {2,3}",
            "target_boundary_datum": "log(11^(9/5)) = (9/5)*log(11)",
            "target_prime_exponent_vector": {str(p): fraction_payload(value) for p, value in TARGET_VECTOR.items()},
            "coordinate_forcing_certificate": coordinate,
            "bounded_w11_scan": bounded,
            "imported_exact_witness": {
                "weights": {str(d): fraction_payload(value) for d, value in witness_weights.items() if value != 0},
                "prime_vector": {str(p): fraction_payload(value) for p, value in witness_vector.items()},
                "exact_representation": witness_vector == TARGET_VECTOR,
                "exports_sourcehood": False,
                "reason": "The exact witness is precisely imported w_11=9/5; it verifies representability but violates the Z12-smooth source class.",
            },
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "coupled_setup": "Localization and coefficient choice are solved jointly as weights w_d, not as a post-hoc product pair.",
                "prime11_coordinate_step": "For d=1..10, the prime-11 exponent is zero; for d=11, it is one. Therefore exact equality forces w_11=9/5.",
                "denominator_step": "The forced value w_11=9/5 has denominator prime 5, outside the Z12-smooth denominator support {2,3}.",
                "finite_scan_step": "A bounded scan over Z12-smooth denominators confirms no exact w_11 hit in the sampled ring fragment.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_coupled_z12_smooth_weighted_log_functional_no_source_audit": all(facts.values()),
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
            "reason": "P2867 tests a genuinely coupled weighted-log functional rather than pairwise recombination.  Exact equality forces w_11=9/5 by the prime-11 coordinate, but that forced weight imports denominator prime 5 and is outside the pure Z12-smooth coefficient/source class.  The imported exact witness remains representability, not sourcehood.",
            "next_honest_step": "Do not replay Z12-smooth weighted log functionals or imported w_11=9/5 as boundary sourcehood.  A next proof-grade move must introduce a new strict mechanism that sources denominator prime 5 and the d=11 endpoint together, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["coupled_z12_smooth_weighted_log_functional_no_source_audit"]
    lines = [
        "# P2867/S1817 coupled Z12-smooth weighted log-functional no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupled weighted-functional scan",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- coordinate forcing: `{audit['coordinate_forcing_certificate']}`",
        f"- bounded scan: `{audit['bounded_w11_scan']}`",
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
    payload = build_payload(read_json(P2866))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2867/S1817 coupled Z12-smooth weighted log-functional no-source audit",
        "## P2867/S1817 coupled Z12-smooth weighted log-functional no-source audit\n\n"
        "`P2867/S1817` tests a genuinely coupled weighted-log functional `F(w)=sum_{d=1}^{11} w_d log(d)` with all coefficients in the pure `Z12`-smooth rational ring whose denominator support is `{2,3}`.  Exact equality to `(9/5)log(11)` forces `w_11=9/5` because only `d=11` carries prime 11.  That forced weight imports denominator prime 5 and is outside the allowed source class; the imported exact witness is representability, not sourcehood.  No boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2867/S1817 coupled weighted-log no-source `L_total` guard",
        "## P2867/S1817 coupled weighted-log no-source `L_total` guard\n\n"
        "`P2867/S1817` adds no strict action term.  The coupled weighted-log obstruction identifies a forced imported coefficient, but it does not supply a unit-bearing boundary/source density, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current coupled Z12-smooth weighted-log no-source guardrail (P2867/S1817, 2026-06-18)",
        "## Current coupled Z12-smooth weighted-log no-source guardrail (P2867/S1817, 2026-06-18)\n\n"
        "- P2867 tests a coupled weighted-log functional rather than pairwise recombination: `F(w)=sum w_d log(d)` with pure `Z12`-smooth rational weights.\n"
        "- Exact endpoint equality forces `w_11=9/5` by the prime-11 coordinate; this imports denominator prime 5 and is outside the `{2,3}` denominator source class.\n"
        "- Do not promote coupled Z12-smooth weighted-log functionals, imported `w_11=9/5`, singleton endpoint representation, Dirichlet data, log-scale harmonicity, or pairwise recombination to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce a new strict mechanism sourcing denominator prime 5 and the `d=11` endpoint together, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
