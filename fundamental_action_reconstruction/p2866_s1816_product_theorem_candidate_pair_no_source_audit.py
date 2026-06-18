#!/usr/bin/env python3
"""P2866/S1816: product-theorem candidate-pair no-source audit.

P2865 sharpened the boundary route to a product obligation: an exact source for
(9/5)*log(11) must export both a non-premise d=11 localizer and an independent
prime-5 coefficient law c=9/5.  This packet tests the next honest finite
candidate class: all pairings of currently available localization witnesses
with currently available coefficient witnesses.

The acceptance rule is intentionally strict: a pair represents a boundary source
only if it is exact, the localizer is strict-sourced/non-premise, the coefficient
law is strict-sourced/non-premise, and the pair supplies a unit-bearing coupling
link rather than merely multiplying two imported facts.

Result: exact representation appears only for imported singleton d=11 together
with imported 9/5 coefficient.  All sourced/currently permitted localizer or
coefficient classes either fail exactness or lack sourcehood.  Thus the product
obligation is not discharged by recombining existing artifacts.
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

P2865 = GEN / "p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit.json"
OUT = GEN / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.json"
MD = GEN / "p2866_s1816_product_theorem_candidate_pair_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET = Fraction(9, 5)


def localization_candidates() -> list[dict[str, Any]]:
    return [
        {
            "name": "aut_z12_invariant_orbit_weights",
            "isolates_log11": False,
            "strict_sourced_nonpremise": True,
            "imports_selector_or_localizer": False,
            "reason": "P2864 exact rank certificate is inconsistent: the unit orbit ties 11 to 5 and 7.",
        },
        {
            "name": "singleton_d11_localizer",
            "isolates_log11": True,
            "strict_sourced_nonpremise": False,
            "imports_selector_or_localizer": True,
            "reason": "P2865 shows d=11 is the unique exact singleton support, but the singleton itself is an imported selector/localizer.",
        },
        {
            "name": "all_singletons_except_d11",
            "isolates_log11": False,
            "strict_sourced_nonpremise": False,
            "imports_selector_or_localizer": True,
            "reason": "P2865 enumerates d=1..10 and finds wrong prime support or zero vector.",
        },
    ]


def coefficient_candidates() -> list[dict[str, Any]]:
    return [
        {
            "name": "pure_z12_smooth_denominator_support_2_3",
            "coefficient": None,
            "coefficient_exact_9_over_5": False,
            "strict_sourced_nonpremise": True,
            "imports_prime5_coefficient": False,
            "reason": "P2863 denominator obstruction: support {2,3} cannot generate denominator prime 5.",
        },
        {
            "name": "integer_local_log_moments",
            "coefficient": None,
            "coefficient_exact_9_over_5": False,
            "strict_sourced_nonpremise": True,
            "imports_prime5_coefficient": False,
            "reason": "P2863 integer moment weights give integer prime-exponent vectors, not 9/5.",
        },
        {
            "name": "imported_prime5_coefficient_9_over_5",
            "coefficient": TARGET,
            "coefficient_exact_9_over_5": True,
            "strict_sourced_nonpremise": False,
            "imports_prime5_coefficient": True,
            "reason": "This represents the required coefficient exactly but imports the missing prime-5 law.",
        },
    ]


def pair_matrix() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for loc in localization_candidates():
        for coeff in coefficient_candidates():
            exact_representation = bool(loc["isolates_log11"] and coeff["coefficient_exact_9_over_5"])
            sourced = bool(loc["strict_sourced_nonpremise"] and coeff["strict_sourced_nonpremise"])
            unit_bearing_coupling = False
            rows.append(
                {
                    "localizer": loc["name"],
                    "coefficient_law": coeff["name"],
                    "exact_representation_of_target": exact_representation,
                    "localizer_strict_sourced_nonpremise": loc["strict_sourced_nonpremise"],
                    "coefficient_strict_sourced_nonpremise": coeff["strict_sourced_nonpremise"],
                    "unit_bearing_coupling_link_exported": unit_bearing_coupling,
                    "accepted_as_boundary_source_theorem": bool(exact_representation and sourced and unit_bearing_coupling),
                    "failure_modes": [
                        label
                        for label, failed in [
                            ("not_exact_boundary_datum", not exact_representation),
                            ("localizer_not_strict_sourced", not loc["strict_sourced_nonpremise"]),
                            ("coefficient_not_strict_sourced", not coeff["strict_sourced_nonpremise"]),
                            ("no_unit_bearing_coupling_link", not unit_bearing_coupling),
                        ]
                        if failed
                    ],
                }
            )
    return rows


def build_payload(p2865: dict[str, Any]) -> dict[str, Any]:
    rows = pair_matrix()
    exact_rows = [row for row in rows if row["exact_representation_of_target"]]
    accepted = [row for row in rows if row["accepted_as_boundary_source_theorem"]]
    facts = {
        "p2865_rechecked": p2865.get("status") == "P2865_SINGLETON_LOCALIZER_COEFFICIENT_PRODUCT_OBLIGATION_NO_GO_AUDIT_NO_CLOSURE",
        "candidate_pair_count_is_nine": len(rows) == 9,
        "unique_exact_pair_is_imported_singleton_times_imported_9_over_5": len(exact_rows) == 1
        and exact_rows[0]["localizer"] == "singleton_d11_localizer"
        and exact_rows[0]["coefficient_law"] == "imported_prime5_coefficient_9_over_5",
        "accepted_count_zero": len(accepted) == 0,
    }
    return {
        "status": "P2866_PRODUCT_THEOREM_CANDIDATE_PAIR_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2865": sha(P2865)},
        "product_theorem_candidate_pair_no_source_audit": {
            "input_status_rechecked": p2865.get("status"),
            "target_boundary_datum": "log(11^(9/5)) = (9/5)*log(11)",
            "acceptance_rule": "exact d=11 localization + exact 9/5 coefficient + strict non-premise sourcehood for both sides + unit-bearing coupling link",
            "localization_candidates": localization_candidates(),
            "coefficient_candidates": [
                {**c, "coefficient": fraction_payload(c["coefficient"]) if isinstance(c["coefficient"], Fraction) else None}
                for c in coefficient_candidates()
            ],
            "candidate_pair_matrix": rows,
            "exact_representation_rows": exact_rows,
            "accepted_candidate_count": len(accepted),
            "proof_certificate": {
                "finite_enumeration": "Three localization witnesses times three coefficient witnesses give nine candidate product-theorem pairs.",
                "exactness_step": "The only exact pair is singleton_d11_localizer times imported_prime5_coefficient_9_over_5.",
                "sourcehood_step": "That exact pair fails because both sides are imported/non-sourced and no unit-bearing coupling link is exported; all sourced-side pairs fail exactness.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_product_theorem_candidate_pair_no_source_audit": all(facts.values()),
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
            "reason": "P2866 tests the finite recombination escape hatch after P2865.  Among nine localizer/coefficient candidate pairs, only imported singleton d=11 times imported coefficient 9/5 represents the endpoint exactly, and that pair has no strict sourcehood or unit-bearing coupling theorem.  Existing sourced classes fail exactness.",
            "next_honest_step": "Do not replay pairwise recombination of existing localizer/coefficient witnesses as boundary sourcehood.  A next proof-grade move must introduce one genuinely new coupled theorem that simultaneously sources d=11 localization, the 9/5 coefficient, and a unit-bearing link, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["product_theorem_candidate_pair_no_source_audit"]
    lines = [
        "# P2866/S1816 product-theorem candidate-pair no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate-pair scan",
        f"- candidate pairs: `{len(audit['candidate_pair_matrix'])}`",
        f"- exact representation rows: `{audit['exact_representation_rows']}`",
        f"- accepted source-theorem rows: `{audit['accepted_candidate_count']}`",
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
    payload = build_payload(read_json(P2865))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2866/S1816 product-theorem candidate-pair no-source audit",
        "## P2866/S1816 product-theorem candidate-pair no-source audit\n\n"
        "`P2866/S1816` tests the finite recombination escape hatch left by `P2865`: three currently available localization witnesses paired with three currently available coefficient witnesses.  The only exact representation of `(9/5)log(11)` is imported singleton `d=11` times imported coefficient `9/5`; all sourced/current classes fail exactness, and the exact pair lacks strict sourcehood plus a unit-bearing coupling link.  Therefore no boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2866/S1816 product-pair no-source `L_total` guard",
        "## P2866/S1816 product-pair no-source `L_total` guard\n\n"
        "`P2866/S1816` adds no strict action term.  Pairwise multiplication of existing localization and coefficient witnesses does not provide a unit-bearing boundary/source density, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current product-theorem candidate-pair no-source guardrail (P2866/S1816, 2026-06-18)",
        "## Current product-theorem candidate-pair no-source guardrail (P2866/S1816, 2026-06-18)\n\n"
        "- P2866 scans the finite recombination class after P2865: current localization witnesses paired with current coefficient witnesses.\n"
        "- The only exact endpoint representation is imported singleton `d=11` times imported coefficient `9/5`; sourced/current pairs fail exactness and the exact pair lacks strict sourcehood plus a unit-bearing coupling link.\n"
        "- Do not promote pairwise recombination, imported singleton localization, imported prime-5 coefficient, Aut-invariant localization, Z12-smooth coefficients, integer log moments, Dirichlet data, or log-scale harmonicity to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must introduce one genuinely new coupled theorem sourcing `d=11`, `9/5`, and the unit-bearing link together, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
