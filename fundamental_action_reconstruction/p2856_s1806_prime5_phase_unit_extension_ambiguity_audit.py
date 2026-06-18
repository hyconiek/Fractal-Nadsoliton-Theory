#!/usr/bin/env python3
"""P2856/S1806: prime-5 phase-unit extension ambiguity audit.

P2855 showed that pure Z12 rational phase lattices cannot exactly source
strict omega=743/4000 and phi=13/80 because both require denominator prime 5.
This packet tests the next narrow candidate: allow a prime-5 phase unit in the
rational lattice and ask whether that extension itself selects the strict tuple.

Result: exact representation becomes possible, but selection does not.  The
same Z12 phase-bit profile is realized by many nearby prime-5-extended rational
pairs, so the extension is representational capacity, not a non-premise strict
source law for omega/phi.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import (
    GEN,
    SEARCH_LIMIT,
    STRICT_OMEGA,
    STRICT_PHI,
    fraction_payload,
    is_z12_compatible_denominator,
    phase_bits_for,
    prime_factors,
)

P2855 = GEN / "p2855_s1805_z12_rational_phase_lattice_source_candidate_audit.json"
OUT = GEN / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.json"
MD = GEN / "p2856_s1806_prime5_phase_unit_extension_ambiguity_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

ALLOWED_EXTENDED_PRIMES = {2, 3, 5}
WINDOW = Fraction(1, 200)  # finite local source-neighborhood test around the strict tuple
MAX_WITNESSES = 12


def is_extended_denominator(q: int) -> bool:
    return set(prime_factors(q)).issubset(ALLOWED_EXTENDED_PRIMES)


def extended_denominators(limit: int = SEARCH_LIMIT) -> list[int]:
    return [q for q in range(1, limit + 1) if is_extended_denominator(q)]


def exact_denominator_supports(value: Fraction, denominators: list[int]) -> list[int]:
    return [q for q in denominators if (value * q).denominator == 1]


def local_same_bit_witnesses(target_bits: list[int], limit: int = SEARCH_LIMIT) -> list[dict[str, Any]]:
    """Enumerate a bounded ambiguity witness family.

    We intentionally use a narrow local box around the strict tuple.  The goal is
    not to classify all possible phase-equivalent pairs, but to prove that the
    prime-5 extension does not uniquely select the strict pair even locally.
    """
    witnesses: list[tuple[Fraction, Fraction, Fraction, Fraction, int, int]] = []
    seen_pairs: set[tuple[Fraction, Fraction]] = set()
    denominators = extended_denominators(limit)
    for q_omega in denominators:
        lo_omega = int((STRICT_OMEGA - WINDOW) * q_omega) - 1
        hi_omega = int((STRICT_OMEGA + WINDOW) * q_omega) + 2
        for n_omega in range(max(0, lo_omega), hi_omega + 1):
            omega = Fraction(n_omega, q_omega)
            if abs(omega - STRICT_OMEGA) > WINDOW:
                continue
            for q_phi in denominators:
                lo_phi = int((STRICT_PHI - WINDOW) * q_phi) - 1
                hi_phi = int((STRICT_PHI + WINDOW) * q_phi) + 2
                for n_phi in range(max(0, lo_phi), hi_phi + 1):
                    phi = Fraction(n_phi, q_phi)
                    if omega == STRICT_OMEGA and phi == STRICT_PHI:
                        continue
                    if abs(phi - STRICT_PHI) > WINDOW:
                        continue
                    try:
                        bits = phase_bits_for(omega, phi)
                    except ValueError:
                        continue
                    if bits == target_bits and (omega, phi) not in seen_pairs:
                        seen_pairs.add((omega, phi))
                        distance = abs(omega - STRICT_OMEGA) + abs(phi - STRICT_PHI)
                        witnesses.append((distance, omega, phi, abs(omega - STRICT_OMEGA), abs(phi - STRICT_PHI), q_omega + q_phi))
                        if len(witnesses) >= 500:
                            break
                if len(witnesses) >= 500:
                    break
            if len(witnesses) >= 500:
                break
        if len(witnesses) >= 500:
            break
    witnesses.sort(key=lambda row: (row[0], row[5], row[1].denominator, row[2].denominator, row[1].numerator, row[2].numerator))
    return [
        {
            "omega": fraction_payload(omega),
            "phi": fraction_payload(phi),
            "omega_abs_error_fraction": f"{omega_err.numerator}/{omega_err.denominator}",
            "phi_abs_error_fraction": f"{phi_err.numerator}/{phi_err.denominator}",
            "total_l1_error_fraction": f"{dist.numerator}/{dist.denominator}",
        }
        for dist, omega, phi, omega_err, phi_err, _ in witnesses[:MAX_WITNESSES]
    ]


def build_payload(p2855: dict[str, Any]) -> dict[str, Any]:
    bits = p2855["z12_rational_phase_lattice_source_candidate_audit"]["p2853_phase_bits"]
    denominators = extended_denominators()
    omega_supports = exact_denominator_supports(STRICT_OMEGA, denominators)
    phi_supports = exact_denominator_supports(STRICT_PHI, denominators)
    witnesses = local_same_bit_witnesses(bits)
    facts = {
        "p2855_rechecked": p2855.get("status") == "P2855_Z12_RATIONAL_PHASE_LATTICE_SOURCE_CANDIDATE_AUDIT_NO_CLOSURE",
        "pure_z12_still_rejected": not is_z12_compatible_denominator(STRICT_OMEGA.denominator) and not is_z12_compatible_denominator(STRICT_PHI.denominator),
        "prime5_extension_represents_omega": bool(omega_supports),
        "prime5_extension_represents_phi": bool(phi_supports),
        "strict_bits_reproduced": phase_bits_for(STRICT_OMEGA, STRICT_PHI) == bits,
        "local_ambiguity_witness_exists": len(witnesses) > 0,
    }
    return {
        "status": "P2856_PRIME5_PHASE_UNIT_EXTENSION_AMBIGUITY_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2855": sha(P2855)},
        "prime5_phase_unit_extension_ambiguity_audit": {
            "input_status_rechecked": p2855.get("status"),
            "candidate": "prime-5-extended rational phase lattice with denominator prime support subset {2,3,5}",
            "search_limit": SEARCH_LIMIT,
            "local_window": f"{WINDOW.numerator}/{WINDOW.denominator}",
            "extended_denominator_count": len(denominators),
            "strict_omega": fraction_payload(STRICT_OMEGA),
            "strict_phi": fraction_payload(STRICT_PHI),
            "omega_exact_supporting_denominators_first12": omega_supports[:12],
            "phi_exact_supporting_denominators_first12": phi_supports[:12],
            "target_phase_bits": bits,
            "same_phase_bit_local_witnesses": witnesses,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "representation_step": "Adding prime 5 to the rational lattice prime support makes exact omega/phi representation possible.",
                "ambiguity_step": "The same Z12 phase-bit profile is realized by non-strict nearby extended-lattice pairs in the bounded local enumeration.",
                "source_step": "Representability plus phase-bit preservation is weaker than a source law selecting 743/4000 and 13/80.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_prime5_extension_ambiguity_audit": all(facts.values()),
            "exports_prime5_phase_unit_source_law": False,
            "exports_strict_phase_frequency_source_law": False,
        },
        "decision": {
            "negative_export_flags": {
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
            "reason": "P2856 grants the missing prime-5 denominator support and confirms exact representation, but finite local enumeration finds multiple non-strict extended-lattice pairs with the same Z12 phase-bit profile.  The prime-5 extension therefore supplies coordinate capacity, not a non-premise source law selecting omega=743/4000 and phi=13/80.",
            "next_honest_step": "Do not replay pure Z12 lattices, prime-5 representability, or phase-bit equivalence as a source.  The next proof-grade move requires a genuinely new source-selection law for the prime-5 phase unit and the exact omega/phi numerators, or a genuinely new eta/beta source law.  Without that new typed premise, preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["prime5_phase_unit_extension_ambiguity_audit"]
    lines = [
        "# P2856/S1806 prime-5 phase-unit extension ambiguity audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Exact representation after importing prime 5",
        f"- omega={audit['strict_omega']['fraction']}; first exact denominators={audit['omega_exact_supporting_denominators_first12']}",
        f"- phi={audit['strict_phi']['fraction']}; first exact denominators={audit['phi_exact_supporting_denominators_first12']}",
        "",
        "## Local ambiguity witnesses",
        f"- local window: `{audit['local_window']}` around the strict tuple",
        f"- witness count reported: `{len(audit['same_phase_bit_local_witnesses'])}`",
    ]
    for row in audit["same_phase_bit_local_witnesses"][:5]:
        lines.append(
            f"- omega={row['omega']['fraction']}, phi={row['phi']['fraction']}, "
            f"L1 error={row['total_l1_error_fraction']}"
        )
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2855))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2856/S1806 prime-5 phase-unit extension ambiguity audit",
        "## P2856/S1806 prime-5 phase-unit extension ambiguity audit\n\n"
        "`P2856/S1806` grants the prime-5 rational phase-unit extension blocked in `P2855/S1805`.  Exact representation of `omega=743/4000` and `phi=13/80` becomes possible, but bounded local enumeration finds other prime-5-extended rational pairs with the same `Z12` phase-bit profile.  Thus prime-5 representability is coordinate capacity, not a non-premise strict source-selection law.  No strict phase/frequency source law, full kernel bridge, role transfer, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2856/S1806 prime-5 extension `L_total` guard",
        "## P2856/S1806 prime-5 extension `L_total` guard\n\n"
        "`P2856/S1806` adds no action term.  Prime-5 rational representability and phase-bit ambiguity witnesses do not provide a unit-bearing source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current prime-5 phase-unit extension ambiguity guardrail (P2856/S1806, 2026-06-18)",
        "## Current prime-5 phase-unit extension ambiguity guardrail (P2856/S1806, 2026-06-18)\n\n"
        "- P2856 grants a prime-5-extended rational phase lattice after P2855 and confirms that exact `omega=743/4000`, `phi=13/80` representation is possible.\n"
        "- The same `Z12` phase-bit profile is realized by multiple nearby prime-5-extended rational pairs; the extension is representational capacity, not a source-selection law for the exact strict tuple.\n"
        "- Do not promote prime-5 representability, local phase-bit equivalence, or imported denominator support to strict phase/frequency source law, selector closure, full bridge, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a genuinely new source-selection law for the prime-5 phase unit and exact `omega/phi` numerators, provide a genuinely new `eta/beta` source law, or preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
