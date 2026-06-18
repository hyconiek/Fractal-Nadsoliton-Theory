#!/usr/bin/env python3
"""P2869/S1819: Aut-character idempotent endpoint-localizer no-source audit.

P2868 showed that Aut(Z12)-invariant orbit weights cannot isolate the endpoint
`d=11`, even after denominator prime 5 is granted.  This packet tests the next
stronger finite class: break the unit orbit by using the full character table of
Aut(Z12)=U(12) on the unit orbit {1,5,7,11}.

The Fourier idempotent for the point 11 exists exactly, so this is a genuine
representability upgrade over P2868.  But it is still not a strict source: it
requires importing a chosen endpoint label / character-polarity projector and,
after multiplying by 9/5, it uses coefficient 9/20 on each character row.  The
current artifacts do not export that projector as a non-premise strict law or a
unit-bearing coupling/localization theorem.
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

P2868 = GEN / "p2868_s1818_prime5_extended_aut_invariant_weighted_log_no_go_audit.json"
OUT = GEN / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.json"
MD = GEN / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

UNITS_Z12 = [1, 5, 7, 11]
TARGET_ENDPOINT = 11
TARGET_COEFFICIENT = Fraction(9, 5)

# Aut(Z12)=V4.  Characters are indexed by their values on 5 and 7; value on 11=5*7.
def character_table() -> dict[str, dict[int, int]]:
    chars: dict[str, dict[int, int]] = {}
    for a in (-1, 1):
        for b in (-1, 1):
            name = f"chi_5={a}_7={b}"
            chars[name] = {1: 1, 5: a, 7: b, 11: a * b}
    return chars


def idempotent_coefficients_for_endpoint(endpoint: int) -> dict[str, Fraction]:
    table = character_table()
    # Delta_endpoint(d) = 1/|G| sum_chi chi(endpoint) chi(d), since all values are real.
    return {name: Fraction(values[endpoint], len(UNITS_Z12)) for name, values in table.items()}


def evaluate_character_linear_combination(coefficients: dict[str, Fraction]) -> dict[int, Fraction]:
    table = character_table()
    return {
        d: sum(coefficients[name] * table[name][d] for name in table)
        for d in UNITS_Z12
    }


def weighted_prime_vector(endpoint_weights: dict[int, Fraction]) -> dict[int, Fraction]:
    out = {p: Fraction(0, 1) for p in PRIMES}
    for d, weight in endpoint_weights.items():
        vector = factor_vector(d)
        for p in PRIMES:
            out[p] += weight * vector[p]
    return out


def denominator_support(values: list[Fraction]) -> list[int]:
    primes: set[int] = set()
    for value in values:
        primes.update(prime_factors(value.denominator))
    return sorted(primes)


def build_payload(p2868: dict[str, Any]) -> dict[str, Any]:
    localizer_coeffs = idempotent_coefficients_for_endpoint(TARGET_ENDPOINT)
    localizer_weights = evaluate_character_linear_combination(localizer_coeffs)
    scaled_coeffs = {name: TARGET_COEFFICIENT * value for name, value in localizer_coeffs.items()}
    scaled_weights = evaluate_character_linear_combination(scaled_coeffs)
    vector = weighted_prime_vector(scaled_weights)
    exact_target = all(vector[p] == TARGET_VECTOR[p] for p in PRIMES)
    imported_requirements = {
        "chosen_endpoint_label_d_11": True,
        "chosen_aut_character_polarity_projector": True,
        "denominator_prime_5_in_scaled_coefficients": 5 in denominator_support(list(scaled_coeffs.values())),
        "unit_bearing_coupling_theorem_missing": True,
    }
    facts = {
        "p2868_rechecked": p2868.get("status") == "P2868_PRIME5_EXTENDED_AUT_INVARIANT_WEIGHTED_LOG_NO_GO_AUDIT_NO_CLOSURE",
        "fourier_idempotent_exactly_isolates_d11": localizer_weights == {1: 0, 5: 0, 7: 0, 11: 1},
        "scaled_idempotent_exactly_matches_target": exact_target,
        "sourcehood_blocked_by_imported_projector_and_unit_link": all(imported_requirements.values()),
        "accepted_count_zero": True,
    }
    return {
        "status": "P2869_AUT_CHARACTER_IDEMPOTENT_ENDPOINT_LOCALIZER_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2868": sha(P2868)},
        "aut_character_idempotent_endpoint_localizer_no_source_audit": {
            "input_status_rechecked": p2868.get("status"),
            "candidate_class": "Aut(Z12)-character Fourier idempotent on the unit orbit {1,5,7,11}, scaled by 9/5",
            "character_table": {name: {str(k): v for k, v in values.items()} for name, values in character_table().items()},
            "endpoint_localizer_coefficients": {name: fraction_payload(value) for name, value in localizer_coeffs.items()},
            "endpoint_localizer_weights_on_unit_orbit": {str(k): fraction_payload(v) for k, v in localizer_weights.items()},
            "scaled_character_coefficients": {name: fraction_payload(value) for name, value in scaled_coeffs.items()},
            "scaled_weights_on_unit_orbit": {str(k): fraction_payload(v) for k, v in scaled_weights.items()},
            "scaled_prime_vector": {str(p): fraction_payload(v) for p, v in vector.items()},
            "target_prime_vector": {str(p): fraction_payload(TARGET_VECTOR[p]) for p in PRIMES},
            "exact_target_representation": exact_target,
            "coefficient_denominator_support": denominator_support(list(scaled_coeffs.values())),
            "imported_requirements": imported_requirements,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "representability_step": "The full Aut-character table gives the point idempotent delta_11 exactly on the unit orbit.",
                "scaling_step": "Multiplying that idempotent by 9/5 exactly reconstructs (9/5)log(11).",
                "sourcehood_step": "The idempotent is a projector onto a pre-chosen endpoint label and character polarity; current artifacts do not export that projector as a non-premise strict law or provide a unit-bearing coupling/localization theorem.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_aut_character_idempotent_source": False,
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
            "reason": "P2869 shows that Aut-character Fourier idempotents can represent the endpoint exactly, but only by importing the endpoint/polarity projector and the scaled 9/20 character coefficients.  This is representability, not strict sourcehood or a unit-bearing coupling theorem.",
            "next_honest_step": "Do not replay Aut-character idempotent endpoint projectors as sourcehood.  A next proof-grade move must supply a strict non-premise law selecting the character projector/polarity and coupling it with a unit-bearing coefficient theorem, or pivot to a different genuinely new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["aut_character_idempotent_endpoint_localizer_no_source_audit"]
    lines = [
        "# P2869/S1819 Aut-character idempotent endpoint-localizer no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact idempotent representation",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- endpoint weights: `{audit['endpoint_localizer_weights_on_unit_orbit']}`",
        f"- scaled prime vector: `{audit['scaled_prime_vector']}`",
        f"- coefficient denominator support: `{audit['coefficient_denominator_support']}`",
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
    payload = build_payload(read_json(P2868))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2869/S1819 Aut-character idempotent endpoint-localizer no-source audit",
        "## P2869/S1819 Aut-character idempotent endpoint-localizer no-source audit\n\n"
        "`P2869/S1819` tests the stronger finite class left after P2868: full Aut(`Z12`)-character Fourier idempotents on the unit orbit `{1,5,7,11}`.  The idempotent `delta_11` exists and, after scaling by `9/5`, exactly represents `(9/5)log(11)`.  This is only representability: it imports a chosen endpoint/polarity projector and scaled `9/20` character coefficients, and it exports no non-premise strict selector/localizer source or unit-bearing coupling/localization theorem.  Therefore no boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2869/S1819 Aut-character idempotent `L_total` guard",
        "## P2869/S1819 Aut-character idempotent `L_total` guard\n\n"
        "`P2869/S1819` adds no strict action term.  The exact Aut-character endpoint projector is a representational idempotent, not an exported unit-bearing boundary/source density, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current Aut-character idempotent endpoint-localizer no-source guardrail (P2869/S1819, 2026-06-18)",
        "## Current Aut-character idempotent endpoint-localizer no-source guardrail (P2869/S1819, 2026-06-18)\n\n"
        "- P2869 tests full Aut(`Z12`)-character Fourier idempotents on the unit orbit `{1,5,7,11}` after the P2868 Aut-invariant obstruction.\n"
        "- The exact endpoint projector `delta_11` exists and scales to `(9/5)log(11)`, but it imports a chosen endpoint/polarity projector and `9/20` character coefficients rather than exporting a strict non-premise source law.\n"
        "- Do not promote Aut-character idempotent projectors, endpoint/polarity imports, prime-5 scaled coefficients, singleton endpoint representation, Dirichlet data, log-scale harmonicity, Aut-invariant orbit weighting, or pairwise recombination to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a strict non-premise law selecting the character projector/polarity and coupling it with a unit-bearing coefficient theorem, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
