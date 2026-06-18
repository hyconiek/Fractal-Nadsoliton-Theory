#!/usr/bin/env python3
"""P2865/S1815: singleton localizer/coefficient product-obligation no-go audit.

P2864 showed that Aut(Z12)-invariant nonlocal localization cannot isolate the
P2862/P2863 right boundary datum, while a singleton d=11 localizer would work
only after importing both a selector/localizer and the prime-5 coefficient 9/5.
This packet makes that product obligation finite and explicit by enumerating
all singleton residue/distance localizers d=1..11 and solving c*log(d) =
(9/5)*log(11) in prime-exponent vectors.

Result: the only exact singleton support is d=11 with coefficient c=9/5.  All
other singleton supports carry the wrong prime support or zero vector.  Thus a
source theorem would need two independent ingredients at once: a strict
selector/localizer for d=11 and a strict coefficient law for 9/5.  Supplying
only one side leaves no boundary/eta source.
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

P2864 = GEN / "p2864_s1814_aut_z12_invariant_nonlocal_boundary_localization_no_go_audit.json"
OUT = GEN / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.json"
MD = GEN / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

PRIMES = [2, 3, 5, 7, 11]
TARGET_COEFFICIENT = Fraction(9, 5)
TARGET_VECTOR = {p: Fraction(0, 1) for p in PRIMES} | {11: TARGET_COEFFICIENT}
Z12_PRIME_SUPPORT = {2, 3}


def factor_vector(n: int) -> dict[int, Fraction]:
    if n == 1:
        return {p: Fraction(0, 1) for p in PRIMES}
    remaining = n
    out = {p: Fraction(0, 1) for p in PRIMES}
    for p in PRIMES:
        while remaining % p == 0:
            out[p] += 1
            remaining //= p
    if remaining != 1:
        raise ValueError(f"unexpected factor {remaining} for {n}")
    return out


def scalar_to_target(base_vector: dict[int, Fraction]) -> Fraction | None:
    nonzero_primes = [p for p, value in base_vector.items() if value != 0]
    if not nonzero_primes:
        return None
    reference = nonzero_primes[0]
    scalar = TARGET_VECTOR[reference] / base_vector[reference]
    for p in PRIMES:
        if scalar * base_vector[p] != TARGET_VECTOR[p]:
            return None
    return scalar


def singleton_rows() -> list[dict[str, Any]]:
    rows = []
    for d in range(1, 12):
        vector = factor_vector(d)
        scalar = scalar_to_target(vector)
        exact = scalar is not None
        rows.append(
            {
                "singleton_d": d,
                "prime_vector": {str(p): fraction_payload(value) for p, value in vector.items()},
                "exact_target_possible_with_scalar": exact,
                "required_scalar": fraction_payload(scalar) if scalar is not None else None,
                "requires_prime5_coefficient": bool(scalar is not None and 5 in prime_factors(scalar.denominator)),
                "requires_nonpremise_singleton_localizer": exact,
                "exports_boundary_source_law": False,
            }
        )
    return rows


def product_obligation_matrix(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    exact_rows = [row for row in rows if row["exact_target_possible_with_scalar"]]
    d11_row = next(row for row in rows if row["singleton_d"] == 11)
    return [
        {
            "candidate": "localizer_without_coefficient_law",
            "finite_witness_passes": d11_row["exact_target_possible_with_scalar"],
            "exports_boundary_source_law": False,
            "verdict": "Choosing d=11 without sourcing coefficient 9/5 does not compute the boundary datum.",
        },
        {
            "candidate": "coefficient_without_localizer_law",
            "finite_witness_passes": TARGET_COEFFICIENT == Fraction(9, 5),
            "exports_boundary_source_law": False,
            "verdict": "Coefficient 9/5 without a strict d=11 localizer does not isolate log(11).",
        },
        {
            "candidate": "singleton_localizer_times_prime5_coefficient",
            "finite_witness_passes": len(exact_rows) == 1 and exact_rows[0]["singleton_d"] == 11,
            "exports_boundary_source_law": False,
            "verdict": "The product represents the target only if both d=11 localizer and coefficient 9/5 are supplied as new premises.",
        },
    ]


def build_payload(p2864: dict[str, Any]) -> dict[str, Any]:
    rows = singleton_rows()
    exact_rows = [row for row in rows if row["exact_target_possible_with_scalar"]]
    matrix = product_obligation_matrix(rows)
    d11 = exact_rows[0] if exact_rows else None
    facts = {
        "p2864_rechecked": p2864.get("status") == "P2864_AUT_Z12_INVARIANT_NONLOCAL_BOUNDARY_LOCALIZATION_NO_GO_AUDIT_NO_CLOSURE",
        "unique_exact_singleton_is_d11": len(exact_rows) == 1 and exact_rows[0]["singleton_d"] == 11,
        "d11_requires_scalar_9_over_5": d11 is not None and d11["required_scalar"]["fraction"] == "9/5",
        "d11_scalar_requires_prime5": d11 is not None and d11["requires_prime5_coefficient"],
        "d11_requires_nonpremise_localizer": d11 is not None and d11["requires_nonpremise_singleton_localizer"],
        "accepted_count_zero": True,
    }
    return {
        "status": "P2865_SINGLETON_LOCALIZER_COEFFICIENT_PRODUCT_OBLIGATION_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2864": sha(P2864)},
        "singleton_localizer_coefficient_product_obligation_no_go_audit": {
            "input_status_rechecked": p2864.get("status"),
            "target_boundary_datum": "log(11^(9/5)) = (9/5)*log(11)",
            "target_prime_exponent_vector": {str(p): fraction_payload(value) for p, value in TARGET_VECTOR.items()},
            "singleton_rows": rows,
            "exact_singleton_rows": exact_rows,
            "product_obligation_matrix": matrix,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "enumeration_step": "All singleton supports d=1..11 were checked by exact prime-exponent vectors.",
                "uniqueness_step": "Only d=11 has prime support compatible with the target vector, and it requires scalar 9/5.",
                "product_obligation_step": "The boundary datum is sourced only if both a strict d=11 localizer and a strict 9/5 coefficient law are exported; neither is supplied by current artifacts.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_singleton_localizer_coefficient_product_obligation_no_go_audit": all(facts.values()),
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
            "reason": "P2865 turns the P2864 singleton escape hatch into an exact product-obligation matrix.  Among all singleton supports d=1..11, only d=11 can match the target, and it requires coefficient 9/5.  Therefore the route needs both a non-premise d=11 localizer and a prime-5 coefficient law; current artifacts provide neither as a strict source.",
            "next_honest_step": "Do not replay singleton localization or coefficient representability separately as boundary sourcehood.  A next proof-grade move must export both sides of the product obligation in one strict theorem, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["singleton_localizer_coefficient_product_obligation_no_go_audit"]
    lines = [
        "# P2865/S1815 singleton localizer/coefficient product-obligation no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Product-obligation scan",
        f"- singleton rows: `{len(audit['singleton_rows'])}`",
        f"- exact singleton rows: `{audit['exact_singleton_rows']}`",
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
    payload = build_payload(read_json(P2864))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2865/S1815 singleton localizer/coefficient product-obligation no-go audit",
        "## P2865/S1815 singleton localizer/coefficient product-obligation no-go audit\n\n"
        "`P2865/S1815` enumerates all singleton residue/distance localizers `d=1..11` after the `P2864` Aut-invariant localization no-go.  Exact prime-exponent matching of `c*log(d)` to `(9/5)log(11)` has a unique support: `d=11` with coefficient `c=9/5`.  Thus the singleton escape hatch requires two new strict inputs at once, a non-premise `d=11` localizer and a prime-5 coefficient law.  Providing either side alone does not source the boundary datum, eta, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2865/S1815 singleton product-obligation `L_total` guard",
        "## P2865/S1815 singleton product-obligation `L_total` guard\n\n"
        "`P2865/S1815` adds no strict action term.  The exact singleton product obligation identifies missing premises rather than supplying a unit-bearing boundary/source density, coupling coefficient, localization/pullback, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current singleton localizer/coefficient product-obligation guardrail (P2865/S1815, 2026-06-18)",
        "## Current singleton localizer/coefficient product-obligation guardrail (P2865/S1815, 2026-06-18)\n\n"
        "- P2865 enumerates all singleton localizers after P2864 and finds a unique exact support: `d=11` with coefficient `9/5`.\n"
        "- This is a product obligation, not closure: current artifacts do not export both a strict non-premise `d=11` localizer and a strict prime-5 coefficient law.\n"
        "- Do not promote singleton localization, coefficient representability, selector import, prime-5 representability, Dirichlet data, Aut-invariant localization, or log-scale harmonicity to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export both sides of the product obligation in one strict theorem, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
