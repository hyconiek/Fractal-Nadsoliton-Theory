#!/usr/bin/env python3
"""P2855/S1805: Z12 rational phase-lattice source candidate audit.

P2854 says the next proof-grade move needs a genuinely new typed source theorem.
This packet tests one narrow candidate rather than replaying affine transport:
can the strict phase/frequency tuple be sourced by a Z12-compatible rational
phase lattice, whose reduced denominators use only primes already present in
12 = 2^2 * 3?

Result: the strict tuple needs a prime-5 denominator component.  A prime-5
extension can represent the tuple, but that extension is a new unsourced unit,
not a strict source law.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2853 = GEN / "p2853_s1803_phase_frequency_bridge_source_audit.json"
P2854 = GEN / "p2854_s1804_post_p2853_professorial_state_map_no_new_live_frontier.json"
OUT = GEN / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.json"
MD = GEN / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

STRICT_OMEGA = Fraction(743, 4000)
STRICT_PHI = Fraction(13, 80)
LEGACY_OMEGA = math.pi / 4.0
LEGACY_PHI = math.pi / 6.0
DOMAIN = tuple(range(12))
ALLOWED_Z12_PRIMES = {2, 3}
SEARCH_LIMIT = max(STRICT_OMEGA.denominator, STRICT_PHI.denominator)


def prime_factors(n: int) -> dict[int, int]:
    factors: dict[int, int] = {}
    value = n
    candidate = 2
    while candidate * candidate <= value:
        while value % candidate == 0:
            factors[candidate] = factors.get(candidate, 0) + 1
            value //= candidate
        candidate += 1 if candidate == 2 else 2
    if value > 1:
        factors[value] = factors.get(value, 0) + 1
    return factors


def is_z12_compatible_denominator(q: int) -> bool:
    return set(prime_factors(q)).issubset(ALLOWED_Z12_PRIMES)


def fraction_payload(value: Fraction) -> dict[str, Any]:
    return {
        "fraction": f"{value.numerator}/{value.denominator}",
        "numerator": value.numerator,
        "denominator": value.denominator,
        "denominator_prime_factors": prime_factors(value.denominator),
        "z12_compatible_denominator": is_z12_compatible_denominator(value.denominator),
    }


def nearest_fraction_with_denominator(target: Fraction, q: int) -> Fraction:
    scaled = target * q
    floor_n = scaled.numerator // scaled.denominator
    candidates = [Fraction(floor_n, q), Fraction(floor_n + 1, q)]
    if floor_n > 0:
        candidates.append(Fraction(floor_n - 1, q))
    return min(candidates, key=lambda value: (abs(value - target), value.denominator, value.numerator))


def z12_denominators(limit: int) -> list[int]:
    return [q for q in range(1, limit + 1) if is_z12_compatible_denominator(q)]


def best_z12_approximation(target: Fraction, limit: int = SEARCH_LIMIT) -> dict[str, Any]:
    best_value: Fraction | None = None
    best_q: int | None = None
    for q in z12_denominators(limit):
        candidate = nearest_fraction_with_denominator(target, q)
        if best_value is None or (abs(candidate - target), q) < (abs(best_value - target), best_q or q):
            best_value = candidate
            best_q = q
    assert best_value is not None and best_q is not None
    error = abs(best_value - target)
    return {
        "target": fraction_payload(target),
        "search_limit": limit,
        "candidate_denominator_count": len(z12_denominators(limit)),
        "best_value": fraction_payload(best_value),
        "best_denominator_before_reduction": best_q,
        "absolute_error_fraction": f"{error.numerator}/{error.denominator}",
        "absolute_error_float": float(error),
        "exact_match_found": error == 0,
    }


def best_common_z12_pair(limit: int = SEARCH_LIMIT) -> dict[str, Any]:
    best: dict[str, Any] | None = None
    for q in z12_denominators(limit):
        omega = nearest_fraction_with_denominator(STRICT_OMEGA, q)
        phi = nearest_fraction_with_denominator(STRICT_PHI, q)
        omega_error = abs(omega - STRICT_OMEGA)
        phi_error = abs(phi - STRICT_PHI)
        score = (max(omega_error, phi_error), omega_error + phi_error, q)
        if best is None or score < best["score_tuple"]:
            best = {
                "score_tuple": score,
                "common_denominator_before_reduction": q,
                "omega": omega,
                "phi": phi,
                "omega_error": omega_error,
                "phi_error": phi_error,
            }
    assert best is not None
    return {
        "search_limit": limit,
        "common_denominator_before_reduction": best["common_denominator_before_reduction"],
        "omega": fraction_payload(best["omega"]),
        "phi": fraction_payload(best["phi"]),
        "omega_absolute_error_fraction": f"{best['omega_error'].numerator}/{best['omega_error'].denominator}",
        "phi_absolute_error_fraction": f"{best['phi_error'].numerator}/{best['phi_error'].denominator}",
        "omega_absolute_error_float": float(best["omega_error"]),
        "phi_absolute_error_float": float(best["phi_error"]),
        "exact_pair_match_found": best["omega_error"] == 0 and best["phi_error"] == 0,
    }


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def bit_from_sign(value: int) -> int:
    if value == 1:
        return 0
    if value == -1:
        return 1
    raise ValueError("zero sign cannot be converted to a phase bit")


def phase_bits_for(omega: Fraction, phi: Fraction) -> list[int]:
    bits = []
    for d in DOMAIN:
        legacy_cos = math.cos(LEGACY_OMEGA * d + LEGACY_PHI)
        candidate_cos = math.cos(float(omega) * d + float(phi))
        bits.append(bit_from_sign(sign(candidate_cos / legacy_cos)))
    return bits


def candidate_matrix(
    omega_approx: dict[str, Any],
    phi_approx: dict[str, Any],
    common_pair: dict[str, Any],
    exact_bits: list[int],
    p2853_bits: list[int],
) -> list[dict[str, Any]]:
    omega_has_prime5 = 5 in prime_factors(STRICT_OMEGA.denominator)
    phi_has_prime5 = 5 in prime_factors(STRICT_PHI.denominator)
    common_omega = Fraction(
        common_pair["omega"]["numerator"],
        common_pair["omega"]["denominator"],
    )
    common_phi = Fraction(
        common_pair["phi"]["numerator"],
        common_pair["phi"]["denominator"],
    )
    common_bits = phase_bits_for(common_omega, common_phi)
    return [
        {
            "candidate": "pure_z12_denominator_lattice_exact_source",
            "finite_witness_passes": not omega_has_prime5 and not phi_has_prime5,
            "exports_strict_source_law": False,
            "verdict": "blocked: reduced strict denominators contain prime 5, outside the prime support of Z12.",
        },
        {
            "candidate": "bounded_z12_lattice_best_approximation",
            "finite_witness_passes": not omega_approx["exact_match_found"] and not phi_approx["exact_match_found"],
            "exports_strict_source_law": False,
            "verdict": "approximation exists but exact omega/phi sourcing fails; approximation is not a source law.",
        },
        {
            "candidate": "common_z12_lattice_pair_approximation",
            "finite_witness_passes": not common_pair["exact_pair_match_found"],
            "exports_strict_source_law": False,
            "phase_bits_match_p2853": common_bits == p2853_bits,
            "verdict": "common-denominator approximation is nonexact; phase-bit agreement, if present, is only a coarse witness.",
        },
        {
            "candidate": "import_prime5_phase_unit_extension",
            "finite_witness_passes": omega_has_prime5 and phi_has_prime5 and exact_bits == p2853_bits,
            "exports_strict_source_law": False,
            "verdict": "allowing prime 5 represents the strict tuple exactly, but the prime-5 unit is imported rather than sourced from Z12.",
        },
    ]


def build_payload(p2853: dict[str, Any], p2854: dict[str, Any]) -> dict[str, Any]:
    omega_approx = best_z12_approximation(STRICT_OMEGA)
    phi_approx = best_z12_approximation(STRICT_PHI)
    common_pair = best_common_z12_pair()
    p2853_bits = p2853["phase_frequency_bridge_source_audit"]["summary"]["phase_factor_bits"]
    exact_bits = phase_bits_for(STRICT_OMEGA, STRICT_PHI)
    rows = candidate_matrix(omega_approx, phi_approx, common_pair, exact_bits, p2853_bits)
    accepted_count = sum(1 for row in rows if row["exports_strict_source_law"])
    facts = {
        "p2853_rechecked": p2853.get("status") == "P2853_PHASE_FREQUENCY_BRIDGE_SOURCE_AUDIT_NO_CLOSURE",
        "p2854_rechecked": p2854.get("status") == "P2854_POST_P2853_PROFESSORIAL_STATE_MAP_NO_NEW_LIVE_FRONTIER",
        "strict_omega_denominator_has_prime5": 5 in prime_factors(STRICT_OMEGA.denominator),
        "strict_phi_denominator_has_prime5": 5 in prime_factors(STRICT_PHI.denominator),
        "no_exact_z12_lattice_match_for_omega": not omega_approx["exact_match_found"],
        "no_exact_z12_lattice_match_for_phi": not phi_approx["exact_match_found"],
        "exact_strict_bits_match_p2853": exact_bits == p2853_bits,
        "accepted_count_zero": accepted_count == 0,
    }
    return {
        "status": "P2855_Z12_RATIONAL_PHASE_LATTICE_SOURCE_CANDIDATE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2853": sha(P2853), "P2854": sha(P2854)},
        "z12_rational_phase_lattice_source_candidate_audit": {
            "input_statuses_rechecked": {"P2853": p2853.get("status"), "P2854": p2854.get("status")},
            "source_candidate": "Z12-compatible rational phase lattice with denominator prime support subset {2,3}",
            "search_limit": SEARCH_LIMIT,
            "allowed_z12_primes": sorted(ALLOWED_Z12_PRIMES),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
            "omega_best_z12_approximation": omega_approx,
            "phi_best_z12_approximation": phi_approx,
            "best_common_z12_pair": common_pair,
            "exact_strict_phase_bits": exact_bits,
            "p2853_phase_bits": p2853_bits,
            "candidate_matrix": rows,
            "accepted_candidate_count": accepted_count,
            "proof_certificate": {
                "denominator_prime_support_step": "Any reduced rational sourced by a pure Z12 lattice has denominator prime support inside {2,3}.",
                "strict_tuple_step": "omega=743/4000 and phi=13/80 both have reduced denominators containing prime 5.",
                "enumeration_step": f"All Z12-compatible denominators q <= {SEARCH_LIMIT} were enumerated; no exact omega/phi match is found.",
                "extension_step": "A prime-5 extension can represent the tuple, but that is a new source premise, not a theorem exported by Z12 itself.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_z12_rational_phase_lattice_obstruction_audit": all(facts.values()),
            "exports_strict_phase_frequency_source_law": False,
        },
        "decision": {
            "negative_export_flags": {
                "z12_phase_lattice_source_exported": False,
                "prime5_phase_unit_source_exported": False,
                "strict_phase_frequency_source_law_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "selector_closure_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2855 tests one concrete candidate source class: a pure Z12-compatible rational phase lattice.  The strict tuple is exact rational data, but its reduced denominators contain prime 5, outside the Z12 denominator prime support {2,3}.  Enumerated Z12-compatible approximations are nonexact.  Importing a prime-5 unit would represent the tuple, but would be a new unsourced premise rather than a strict source theorem.",
            "next_honest_step": "Do not replay pure Z12 rational-lattice sourcing for omega/phi.  A next proof-grade move must either supply a genuine source for the prime-5 phase unit, or pivot to a genuinely new eta/beta source law.  Without one of those new typed premises, preserve the P2854 no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["z12_rational_phase_lattice_source_candidate_audit"]
    lines = [
        "# P2855/S1805 Z12 rational phase-lattice source candidate audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact denominator audit",
        f"- omega={audit['strict_omega']['fraction']}; factors={audit['strict_omega']['denominator_prime_factors']}; z12-compatible={audit['strict_omega']['z12_compatible_denominator']}",
        f"- phi={audit['strict_phi']['fraction']}; factors={audit['strict_phi']['denominator_prime_factors']}; z12-compatible={audit['strict_phi']['z12_compatible_denominator']}",
        "",
        "## Best Z12-compatible approximations",
        f"- omega best={audit['omega_best_z12_approximation']['best_value']['fraction']}; error={audit['omega_best_z12_approximation']['absolute_error_fraction']}",
        f"- phi best={audit['phi_best_z12_approximation']['best_value']['fraction']}; error={audit['phi_best_z12_approximation']['absolute_error_fraction']}",
        f"- common pair omega={audit['best_common_z12_pair']['omega']['fraction']}; phi={audit['best_common_z12_pair']['phi']['fraction']}; exact_pair={audit['best_common_z12_pair']['exact_pair_match_found']}",
        "",
        "## Candidate matrix",
    ]
    for row in audit["candidate_matrix"]:
        extra = f"; phase_bits_match_p2853={row['phase_bits_match_p2853']}" if "phase_bits_match_p2853" in row else ""
        lines.append(
            f"- {row['candidate']}: finite_witness_passes={row['finite_witness_passes']}; "
            f"exports_strict_source_law={row['exports_strict_source_law']}{extra}; {row['verdict']}"
        )
    lines.extend(
        [
            "",
            "## Boundary",
            payload["decision"]["reason"],
            "",
            "## Recommendation",
            payload["decision"]["next_honest_step"],
        ]
    )
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2853), read_json(P2854))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2855/S1805 Z12 rational phase-lattice source candidate audit",
        "## P2855/S1805 Z12 rational phase-lattice source candidate audit\n\n"
        "`P2855/S1805` tests a pure `Z12`-compatible rational phase-lattice source for strict `omega/phi`.  "
        "The reduced strict values `omega=743/4000` and `phi=13/80` both require denominator prime `5`, outside the pure `Z12` prime support `{2,3}`.  "
        "Finite enumeration of `Z12`-compatible denominators up to `4000` gives only nonexact approximations.  A prime-5 extension would be a new unsourced phase unit, not a closure theorem.  No strict phase/frequency source law, full bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2855/S1805 rational phase lattice `L_total` guard",
        "## P2855/S1805 rational phase lattice `L_total` guard\n\n"
        "`P2855/S1805` adds no action term.  Exact denominator obstruction and rational-lattice approximation do not provide a unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current Z12 rational phase-lattice source candidate guardrail (P2855/S1805, 2026-06-18)",
        "## Current Z12 rational phase-lattice source candidate guardrail (P2855/S1805, 2026-06-18)\n\n"
        "- P2855 tests a pure `Z12`-compatible rational phase-lattice source candidate for strict `omega=743/4000` and `phi=13/80`.\n"
        "- The reduced strict denominators require prime `5`, outside the pure `Z12` denominator prime support `{2,3}`; finite `Z12`-compatible approximations are nonexact, and importing a prime-5 phase unit is a new unsourced premise rather than closure.\n"
        "- Do not promote pure `Z12` rational-lattice approximation, phase-bit agreement, or an imported prime-5 unit to strict phase/frequency source law, selector closure, full bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must source the prime-5 phase unit, supply a different strict phase/frequency source law, provide a genuinely new `eta/beta` source law, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
