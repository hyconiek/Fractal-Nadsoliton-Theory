#!/usr/bin/env python3
"""P2872/S1822: cyclic-predecessor endpoint orientation no-source audit.

P2871 exhausted GL(2,2)-invariant Boolean predicates on the unit-projector layer.
This packet tests a different typed candidate: use the ambient cyclic order of
Z12 and select the predecessor endpoint d=-1 mod 12, i.e. d=11.

The finite calculation confirms that this predecessor localizer represents the
boundary target exactly when multiplied by 9/5.  But the selector is not a
strict source on current artifacts: reflection of the Z12 circle swaps the two
neighbors +1 and -1, so choosing predecessor rather than successor imports an
orientation/boundary-arrow convention.  Reflection-invariant adjacent-endpoint
predicates cannot select singleton 11.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN, fraction_payload
from p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit import PRIMES, TARGET_VECTOR, factor_vector
from p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit import TARGET_COEFFICIENT

P2871 = GEN / "p2871_s1821_exhaustive_gl22_invariant_projector_predicate_no_go_audit.json"
OUT = GEN / "p2872_s1822_cyclic_predecessor_endpoint_orientation_no_source_audit.json"
MD = GEN / "p2872_s1822_cyclic_predecessor_endpoint_orientation_no_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
IDENTITY = 0
SUCCESSOR_ENDPOINT = 1
PREDECESSOR_ENDPOINT = 11
ADJACENT_ENDPOINTS = [SUCCESSOR_ENDPOINT, PREDECESSOR_ENDPOINT]


def reflection(endpoint: int) -> int:
    return (-endpoint) % MODULUS


def adjacent_predicates() -> list[frozenset[int]]:
    return [frozenset(), frozenset({SUCCESSOR_ENDPOINT}), frozenset({PREDECESSOR_ENDPOINT}), frozenset(ADJACENT_ENDPOINTS)]


def predicate_record(predicate: frozenset[int]) -> dict[str, Any]:
    reflected = frozenset(reflection(endpoint) for endpoint in predicate)
    return {
        "selected": sorted(predicate),
        "reflected_selected": sorted(reflected),
        "reflection_invariant": reflected == predicate,
        "selects_predecessor_11": predicate == frozenset({PREDECESSOR_ENDPOINT}),
    }


def singleton_prime_vector(endpoint: int, coefficient: Fraction) -> dict[int, Fraction]:
    vector = factor_vector(endpoint)
    return {p: coefficient * vector[p] for p in PRIMES}


def build_payload(p2871: dict[str, Any]) -> dict[str, Any]:
    predecessor_vector = singleton_prime_vector(PREDECESSOR_ENDPOINT, TARGET_COEFFICIENT)
    successor_vector = singleton_prime_vector(SUCCESSOR_ENDPOINT, TARGET_COEFFICIENT)
    exact_predecessor_target = all(predecessor_vector[p] == TARGET_VECTOR[p] for p in PRIMES)
    records = [predicate_record(predicate) for predicate in adjacent_predicates()]
    invariant_records = [record for record in records if record["reflection_invariant"]]
    accepted = [record for record in invariant_records if record["selects_predecessor_11"]]
    facts = {
        "p2871_rechecked": p2871.get("status") == "P2871_EXHAUSTIVE_GL22_INVARIANT_PROJECTOR_PREDICATE_NO_GO_AUDIT_NO_CLOSURE",
        "predecessor_11_exactly_represents_target_with_9_over_5": exact_predecessor_target,
        "reflection_swaps_successor_1_and_predecessor_11": reflection(SUCCESSOR_ENDPOINT) == PREDECESSOR_ENDPOINT and reflection(PREDECESSOR_ENDPOINT) == SUCCESSOR_ENDPOINT,
        "reflection_invariant_predicates_do_not_select_singleton_11": accepted == [],
        "successor_endpoint_is_target_blind": all(successor_vector[p] == 0 for p in PRIMES),
        "accepted_count_zero": True,
    }
    return {
        "status": "P2872_CYCLIC_PREDECESSOR_ENDPOINT_ORIENTATION_NO_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2871": sha(P2871)},
        "cyclic_predecessor_endpoint_orientation_no_source_audit": {
            "input_status_rechecked": p2871.get("status"),
            "candidate_class": "ambient Z12 cyclic-order predecessor localizer d=-1 mod 12, with reflection-symmetry acceptance test",
            "predecessor_endpoint": PREDECESSOR_ENDPOINT,
            "successor_endpoint": SUCCESSOR_ENDPOINT,
            "predecessor_prime_vector_scaled_by_9_over_5": {str(p): fraction_payload(v) for p, v in predecessor_vector.items()},
            "successor_prime_vector_scaled_by_9_over_5": {str(p): fraction_payload(v) for p, v in successor_vector.items()},
            "exact_predecessor_target_representation": exact_predecessor_target,
            "adjacent_predicate_records": records,
            "reflection_invariant_adjacent_predicates": invariant_records,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "representability_step": "The convention d=-1 mod 12 gives d=11, and (9/5)log(11) matches the target prime vector exactly.",
                "reflection_step": "The circle reflection d -> -d swaps +1 and -1, so predecessor selection is orientation-sensitive.",
                "sourcehood_step": "Current artifacts do not export a strict non-premise orientation/boundary-arrow law plus unit-bearing coefficient/coupling theorem; reflection-invariant adjacent predicates cannot select singleton 11.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_cyclic_predecessor_source": False,
            "exports_boundary_source_law": False,
            "exports_eta_source_law": False,
            "exports_prime5_source_law": False,
            "exports_selector_or_localizer_source": False,
            "exports_orientation_source_law": False,
            "exports_unit_bearing_coupling_localization_theorem": False,
            "exports_strict_damping_compression_bridge": False,
        },
        "decision": {
            "negative_export_flags": {
                "boundary_source_law_exported": False,
                "eta_source_exported": False,
                "prime5_source_exported": False,
                "selector_or_localizer_source_exported": False,
                "orientation_source_law_exported": False,
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
            "reason": "P2872 shows that the cyclic predecessor convention d=-1 mod 12 exactly represents the endpoint target, but it is orientation-sensitive: reflection swaps d=1 and d=11, and current artifacts do not export a strict orientation/boundary-arrow source law or unit-bearing coupling theorem.",
            "next_honest_step": "Do not replay cyclic predecessor, d=-1, or clockwise/counterclockwise endpoint conventions as sourcehood.  A next proof-grade move must export a strict non-premise orientation/boundary-arrow law that fixes predecessor over successor and supplies the unit-bearing 9/5 coefficient/coupling theorem, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["cyclic_predecessor_endpoint_orientation_no_source_audit"]
    lines = [
        "# P2872/S1822 cyclic-predecessor endpoint orientation no-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Cyclic predecessor audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- predecessor endpoint: `{audit['predecessor_endpoint']}`",
        f"- exact predecessor target representation: `{audit['exact_predecessor_target_representation']}`",
        f"- reflection-invariant adjacent predicates: `{audit['reflection_invariant_adjacent_predicates']}`",
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
    payload = build_payload(read_json(P2871))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2872/S1822 cyclic-predecessor endpoint orientation no-source audit",
        "## P2872/S1822 cyclic-predecessor endpoint orientation no-source audit\n\n"
        "`P2872/S1822` tests a different typed escape hatch after P2871: the ambient `Z12` cyclic-order predecessor convention `d=-1 mod 12`, hence `d=11`.  This convention exactly represents `(9/5)log(11)`, but it is orientation-sensitive because reflection swaps successor `d=1` and predecessor `d=11`.  Reflection-invariant adjacent predicates cannot select singleton `11`, and current artifacts export no strict orientation/boundary-arrow source law or unit-bearing `9/5` coupling theorem.  No boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2872/S1822 cyclic predecessor `L_total` guard",
        "## P2872/S1822 cyclic predecessor `L_total` guard\n\n"
        "`P2872/S1822` adds no strict action term.  The cyclic predecessor convention does not export a unit-bearing boundary/source density, orientation law, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current cyclic-predecessor endpoint orientation no-source guardrail (P2872/S1822, 2026-06-18)",
        "## Current cyclic-predecessor endpoint orientation no-source guardrail (P2872/S1822, 2026-06-18)\n\n"
        "- P2872 tests the ambient `Z12` cyclic-order predecessor convention `d=-1 mod 12` as a distinct endpoint-localizer candidate after P2871.\n"
        "- The predecessor `d=11` exactly represents `(9/5)log(11)`, but reflection swaps predecessor `d=11` with successor `d=1`; selecting predecessor imports an orientation/boundary-arrow convention.\n"
        "- Do not promote cyclic predecessor, `d=-1`, clockwise/counterclockwise endpoint conventions, reflection-sensitive adjacency, Boolean projector predicates, Aut-character idempotents, prime-5 scaled coefficients, Dirichlet data, or log-scale harmonicity to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must export a strict non-premise orientation/boundary-arrow law fixing predecessor over successor and a unit-bearing `9/5` coefficient/coupling theorem, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
