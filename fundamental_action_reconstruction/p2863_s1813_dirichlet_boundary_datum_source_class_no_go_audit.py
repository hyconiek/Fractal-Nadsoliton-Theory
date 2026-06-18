#!/usr/bin/env python3
"""P2863/S1813: Dirichlet boundary-datum source-class no-go audit.

P2862 showed that Dirichlet endpoint recovery selects eta=9/5 only after the
right boundary datum log(11^(9/5)) has already been imported.  This packet tests
whether simple strict/Z12-compatible boundary source classes can source that
right endpoint without importing the missing prime-5 coefficient.

The target boundary datum is (9/5)*log(11).  In prime-exponent form this is the
vector with exponent 9/5 on prime 11 and zero on primes 2,3,5,7.  Pure Z12
rational coefficient classes with denominator support {2,3}, and integer/local
moment classes over d=1..11, cannot produce the denominator-5 exponent.  A
prime-5-extended coefficient can represent the target exactly, but that is the
same imported premise already identified in P2855/P2856/P2862, not a source law.
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

P2862 = GEN / "p2862_s1812_log_scale_dirichlet_boundary_eta_recovery_no_source_audit.json"
OUT = GEN / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.json"
MD = GEN / "p2863_s1813_dirichlet_boundary_datum_source_class_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_ETA = Fraction(9, 5)
RIGHT_D = 11
PRIMES = [2, 3, 5, 7, 11]
Z12_PRIME_SUPPORT = {2, 3}
INTEGER_WEIGHT_RANGE = range(-3, 4)
Z12_SMOOTH_DENOMINATOR_LIMIT = 72


def factor_vector(n: int) -> dict[int, int]:
    remaining = n
    out = {p: 0 for p in PRIMES}
    for p in PRIMES:
        while remaining % p == 0:
            out[p] += 1
            remaining //= p
    if remaining != 1:
        raise ValueError(f"unexpected prime factor {remaining} for {n}")
    return out


FACTOR_VECTORS = {d: factor_vector(d) for d in range(1, 12)}
TARGET_VECTOR = {p: Fraction(0, 1) for p in PRIMES} | {11: TARGET_ETA}


def denominator_prime_support(value: Fraction) -> set[int]:
    return set(prime_factors(value.denominator))


def is_z12_smooth_denominator(value: Fraction) -> bool:
    return denominator_prime_support(value).issubset(Z12_PRIME_SUPPORT)


def z12_smooth_fractions(limit: int = Z12_SMOOTH_DENOMINATOR_LIMIT) -> list[Fraction]:
    values: set[Fraction] = set()
    for denominator in range(1, limit + 1):
        if set(prime_factors(denominator)).issubset(Z12_PRIME_SUPPORT):
            for numerator in range(-2 * limit, 2 * limit + 1):
                values.add(Fraction(numerator, denominator))
    return sorted(values)


def coefficient_scan_rows() -> list[dict[str, Any]]:
    values = z12_smooth_fractions()
    best = min(values, key=lambda value: abs(value - TARGET_ETA))
    exact = [value for value in values if value == TARGET_ETA]
    return [
        {
            "candidate": "pure_z12_smooth_coefficient_on_log11",
            "scanned_value_count": len(values),
            "exact_match_count": len(exact),
            "best_value": fraction_payload(best),
            "best_abs_error": float(abs(best - TARGET_ETA)),
            "target_denominator_prime_support": sorted(denominator_prime_support(TARGET_ETA)),
            "allowed_denominator_prime_support": sorted(Z12_PRIME_SUPPORT),
            "exports_boundary_source_law": False,
        },
        {
            "candidate": "prime5_extended_coefficient_on_log11",
            "exact_match_count": 1,
            "exact_value": fraction_payload(TARGET_ETA),
            "imports_prime5": 5 in denominator_prime_support(TARGET_ETA),
            "exports_boundary_source_law": False,
        },
    ]


def integer_moment_scan() -> dict[str, Any]:
    """Bounded integer scan certificate without enumerating the full product space.

    Integer local log moments produce integer prime-exponent vectors.  The
    target vector has only one nonzero component, 9/5 on prime 11, so the
    best bounded integer witness sets the log(11) weight to the nearest
    allowed integer and all other prime components to zero.
    """
    nearest_11_weight = min(INTEGER_WEIGHT_RANGE, key=lambda weight: abs(Fraction(weight, 1) - TARGET_ETA))
    vector = {p: Fraction(0, 1) for p in PRIMES}
    vector[11] = Fraction(nearest_11_weight, 1)
    diff = {p: vector[p] - TARGET_VECTOR[p] for p in PRIMES}
    l1 = sum(abs(value) for value in diff.values())
    exact_possible_by_denominator = TARGET_ETA.denominator == 1
    return {
        "candidate": "bounded_integer_local_log_moment_over_d_2_to_11",
        "weight_range": [min(INTEGER_WEIGHT_RANGE), max(INTEGER_WEIGHT_RANGE)],
        "exact_match_count": 0 if not exact_possible_by_denominator else 1,
        "best_l1_prime_exponent_error": fraction_payload(l1),
        "best_vector": {str(p): fraction_payload(value) for p, value in vector.items()},
        "best_diff": {str(p): fraction_payload(value) for p, value in diff.items()},
        "exports_boundary_source_law": False,
        "certificate": "Integer moment weights yield integer prime exponents; target prime-11 exponent 9/5 is non-integer, so no integer-weight exact match exists in any bounded range.",
    }

def z12_smooth_denominator_obstruction() -> dict[str, Any]:
    return {
        "candidate": "z12_smooth_rational_local_moment_denominator_obstruction",
        "target_prime11_exponent": fraction_payload(TARGET_ETA),
        "target_exponent_denominator_support": sorted(denominator_prime_support(TARGET_ETA)),
        "allowed_denominator_support": sorted(Z12_PRIME_SUPPORT),
        "can_represent_target_exponent_without_prime5": is_z12_smooth_denominator(TARGET_ETA),
        "exports_boundary_source_law": False,
        "verdict": "Any exact local log-moment source for the prime-11 endpoint exponent must produce coefficient 9/5 on log(11); pure Z12-smooth rational denominators cannot produce the denominator prime 5.",
    }


def build_payload(p2862: dict[str, Any]) -> dict[str, Any]:
    coeff_rows = coefficient_scan_rows()
    integer_scan = integer_moment_scan()
    denominator_obstruction = z12_smooth_denominator_obstruction()
    facts = {
        "p2862_rechecked": p2862.get("status") == "P2862_LOG_SCALE_DIRICHLET_BOUNDARY_ETA_RECOVERY_NO_SOURCE_AUDIT_NO_CLOSURE",
        "pure_z12_coefficient_has_no_exact_match": coeff_rows[0]["exact_match_count"] == 0,
        "prime5_extension_exact_but_imported": coeff_rows[1]["exact_match_count"] == 1 and coeff_rows[1]["imports_prime5"],
        "bounded_integer_moment_has_no_exact_match": integer_scan["exact_match_count"] == 0,
        "z12_smooth_denominator_obstruction_holds": not denominator_obstruction["can_represent_target_exponent_without_prime5"],
        "accepted_count_zero": True,
    }
    return {
        "status": "P2863_DIRICHLET_BOUNDARY_DATUM_SOURCE_CLASS_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2862": sha(P2862)},
        "dirichlet_boundary_datum_source_class_no_go_audit": {
            "input_status_rechecked": p2862.get("status"),
            "target_boundary_datum": "log(11^(9/5)) = (9/5)*log(11)",
            "target_eta": fraction_payload(TARGET_ETA),
            "right_endpoint_d": RIGHT_D,
            "prime_basis": PRIMES,
            "target_prime_exponent_vector": {str(p): fraction_payload(value) for p, value in TARGET_VECTOR.items()},
            "coefficient_scan_rows": coeff_rows,
            "integer_moment_scan": integer_scan,
            "denominator_obstruction": denominator_obstruction,
            "candidate_matrix": [
                {
                    "candidate": "pure_z12_smooth_boundary_coefficient",
                    "finite_witness_passes": coeff_rows[0]["exact_match_count"] == 0,
                    "exports_boundary_source_law": False,
                    "verdict": "no exact pure-Z12 coefficient produces 9/5 on log(11).",
                },
                {
                    "candidate": "prime5_extended_boundary_coefficient",
                    "finite_witness_passes": coeff_rows[1]["exact_match_count"] == 1,
                    "exports_boundary_source_law": False,
                    "verdict": "exact representation is possible only by importing the prime-5 coefficient already under audit.",
                },
                {
                    "candidate": "bounded_integer_local_log_moment",
                    "finite_witness_passes": integer_scan["exact_match_count"] == 0,
                    "exports_boundary_source_law": False,
                    "verdict": "bounded integer moments cannot yield the fractional prime-11 exponent 9/5.",
                },
                {
                    "candidate": "z12_smooth_denominator_obstruction",
                    "finite_witness_passes": not denominator_obstruction["can_represent_target_exponent_without_prime5"],
                    "exports_boundary_source_law": False,
                    "verdict": denominator_obstruction["verdict"],
                },
            ],
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "target_step": "The P2862 right boundary is exactly (9/5)*log(11), i.e. prime exponent 9/5 on prime 11.",
                "z12_step": "Pure Z12 denominator support allows only primes {2,3}; it cannot create denominator prime 5 in the coefficient 9/5.",
                "integer_step": "Integer local log moments over d=2..11 produce integer prime-exponent vectors, so they cannot equal exponent 9/5.",
                "extension_step": "Allowing prime 5 represents the target but imports the missing source rather than deriving it.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_dirichlet_boundary_datum_source_class_no_go_audit": all(facts.values()),
            "exports_boundary_source_law": False,
            "exports_eta_source_law": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "eta_source_exported": False,
                "prime5_source_exported": False,
                "target_independent_beta_unit_source_exported": False,
                "unit_bearing_coupling_localization_theorem_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "full_kernel_bridge_exported": False,
                "role_transfer_exported": False,
                "ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2863 audits source classes for the P2862 right boundary datum.  Pure Z12-smooth coefficients and bounded integer local log moments cannot produce the required 9/5 coefficient on log(11); prime-5 extension can represent it exactly but only by importing the missing denominator source.  No boundary source law or eta source is exported.",
            "next_honest_step": "Do not replay Dirichlet endpoint data, pure Z12 coefficient scans, integer log moments, or prime-5 representability as boundary sourcehood.  A next proof-grade move must introduce a genuinely sourced prime-5/9 numerator boundary law, a nonlocal unit-bearing coupling/localization theorem that computes the endpoint datum, or a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["dirichlet_boundary_datum_source_class_no_go_audit"]
    coeff = audit["coefficient_scan_rows"][0]
    lines = [
        "# P2863/S1813 Dirichlet boundary-datum source-class no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Boundary-source scan",
        f"- target boundary datum: `{audit['target_boundary_datum']}`",
        f"- pure Z12 exact coefficient matches: `{coeff['exact_match_count']}`",
        f"- integer moment exact matches: `{audit['integer_moment_scan']['exact_match_count']}`",
        f"- prime-5 extension imports prime 5: `{audit['coefficient_scan_rows'][1]['imports_prime5']}`",
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
    payload = build_payload(read_json(P2862))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2863/S1813 Dirichlet boundary-datum source-class no-go audit",
        "## P2863/S1813 Dirichlet boundary-datum source-class no-go audit\n\n"
        "`P2863/S1813` audits source classes for the `P2862` right Dirichlet datum `log(11^(9/5))=(9/5)log(11)`.  Pure `Z12`-smooth rational coefficients with denominator support `{2,3}` and bounded integer local log-moments over `d=2..11` cannot produce the required denominator-prime-5 coefficient on `log(11)`.  Prime-5 extension represents the datum exactly only by importing the missing source.  No boundary source law, eta source, target-independent beta unit, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2863/S1813 boundary-datum source-class `L_total` guard",
        "## P2863/S1813 boundary-datum source-class `L_total` guard\n\n"
        "`P2863/S1813` adds no strict action term.  Scans over `Z12`-smooth boundary coefficients and integer local log-moments fail to source the endpoint datum; prime-5 representation is imported capacity, not a unit-bearing boundary/source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current Dirichlet boundary-datum source-class no-go guardrail (P2863/S1813, 2026-06-18)",
        "## Current Dirichlet boundary-datum source-class no-go guardrail (P2863/S1813, 2026-06-18)\n\n"
        "- P2863 audits source classes for the P2862 right boundary datum `log(11^(9/5))`.\n"
        "- Pure `Z12`-smooth coefficients and bounded integer local log-moments cannot produce the required `9/5` coefficient on `log(11)`; prime-5 extension is exact representation but imported premise.\n"
        "- Do not promote Dirichlet endpoint data, pure Z12 coefficient scans, integer log moments, prime-5 representability, log-scale harmonicity, or multiplicative covariance to strict damping/compression bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must source the prime-5/9 numerator boundary law, provide a nonlocal unit-bearing coupling/localization theorem computing the endpoint datum, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
