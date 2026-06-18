#!/usr/bin/env python3
"""P2870/S1820: intrinsic character-projector selector no-go audit.

P2869 proved that an Aut(Z12)-character Fourier idempotent can represent the
endpoint d=11 exactly, but only after importing the endpoint/polarity projector.
This packet tests the honest next question: can the character-idempotent layer
itself select the d=11 projector without an external endpoint label?

We enumerate the four point idempotents on the unit orbit {1,5,7,11}, the full
GL(2,2) relabeling symmetry of the V4 unit/character square, and a finite family
of intrinsic score features (coefficient moments, sign counts, support size,
norms, and target-prime exactness).  The computation finds that 5, 7, and 11 are
in one intrinsic relabeling orbit.  Exact target-prime matching singles out 11
only by reimporting the already named boundary target log(11), so it is not an
independent strict projector/source law.
"""
from __future__ import annotations

import itertools
import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN, fraction_payload
from p2865_s1815_singleton_localizer_coefficient_product_obligation_no_go_audit import PRIMES, TARGET_VECTOR, factor_vector
from p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit import (
    TARGET_COEFFICIENT,
    UNITS_Z12,
    character_table,
    denominator_support,
    evaluate_character_linear_combination,
    idempotent_coefficients_for_endpoint,
    weighted_prime_vector,
)

P2869 = GEN / "p2869_s1819_aut_character_idempotent_endpoint_localizer_no_source_audit.json"
OUT = GEN / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.json"
MD = GEN / "p2870_s1820_intrinsic_character_projector_selector_no_go_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NONIDENTITY_UNITS = [5, 7, 11]


def all_v4_relabelings() -> list[dict[int, int]]:
    """Return all group automorphisms of V4 fixing 1 and permuting nonidentity units."""
    return [{1: 1, **dict(zip(NONIDENTITY_UNITS, perm))} for perm in itertools.permutations(NONIDENTITY_UNITS)]


def orbit_of_endpoint(endpoint: int) -> list[int]:
    return sorted({mapping[endpoint] for mapping in all_v4_relabelings()})


def projector_payload(endpoint: int) -> dict[str, Any]:
    coeffs = idempotent_coefficients_for_endpoint(endpoint)
    weights = evaluate_character_linear_combination(coeffs)
    scaled = {name: TARGET_COEFFICIENT * value for name, value in coeffs.items()}
    scaled_weights = evaluate_character_linear_combination(scaled)
    vector = weighted_prime_vector(scaled_weights)
    coefficient_values = list(coeffs.values())
    positive = sum(1 for value in coefficient_values if value > 0)
    negative = sum(1 for value in coefficient_values if value < 0)
    zero = sum(1 for value in coefficient_values if value == 0)
    moments = {k: sum(value**k for value in coefficient_values) for k in range(1, 5)}
    exact_target = all(vector[p] == TARGET_VECTOR[p] for p in PRIMES)
    return {
        "endpoint": endpoint,
        "intrinsic_relabeling_orbit": orbit_of_endpoint(endpoint),
        "coefficients": {name: fraction_payload(value) for name, value in coeffs.items()},
        "weights_on_unit_orbit": {str(k): fraction_payload(v) for k, v in weights.items()},
        "scaled_prime_vector": {str(p): fraction_payload(v) for p, v in vector.items()},
        "intrinsic_features": {
            "support_size": sum(1 for value in coefficient_values if value != 0),
            "positive_count": positive,
            "negative_count": negative,
            "zero_count": zero,
            "coefficient_multiset": sorted(fraction_payload(value)["fraction"] for value in coefficient_values),
            "moments": {str(k): fraction_payload(value) for k, value in moments.items()},
            "l2_norm_squared": fraction_payload(sum(value * value for value in coefficient_values)),
            "scaled_denominator_support": denominator_support(list(scaled.values())),
        },
        "exact_target_prime_vector_match": exact_target,
    }


def intrinsic_feature_signature(payload: dict[str, Any]) -> tuple[Any, ...]:
    features = payload["intrinsic_features"]
    return (
        features["support_size"],
        features["positive_count"],
        features["negative_count"],
        features["zero_count"],
        tuple(features["coefficient_multiset"]),
        tuple(features["moments"][str(k)]["fraction"] for k in range(1, 5)),
        features["l2_norm_squared"]["fraction"],
        tuple(features["scaled_denominator_support"]),
    )


def build_payload(p2869: dict[str, Any]) -> dict[str, Any]:
    projectors = {endpoint: projector_payload(endpoint) for endpoint in UNITS_Z12}
    nonidentity_signatures = {endpoint: intrinsic_feature_signature(projectors[endpoint]) for endpoint in NONIDENTITY_UNITS}
    target_hits = [endpoint for endpoint, payload in projectors.items() if payload["exact_target_prime_vector_match"]]
    facts = {
        "p2869_rechecked": p2869.get("status") == "P2869_AUT_CHARACTER_IDEMPOTENT_ENDPOINT_LOCALIZER_NO_SOURCE_AUDIT_NO_CLOSURE",
        "gl22_ties_5_7_11": orbit_of_endpoint(11) == NONIDENTITY_UNITS,
        "intrinsic_features_do_not_single_out_11_among_nonidentity_projectors": len(set(nonidentity_signatures.values())) == 1,
        "target_match_singles_11_only_by_reusing_named_target": target_hits == [11],
        "accepted_count_zero": True,
    }
    return {
        "status": "P2870_INTRINSIC_CHARACTER_PROJECTOR_SELECTOR_NO_GO_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2869": sha(P2869)},
        "intrinsic_character_projector_selector_no_go_audit": {
            "input_status_rechecked": p2869.get("status"),
            "candidate_class": "label-free intrinsic selector over Aut(Z12)-character point idempotents on {1,5,7,11}",
            "v4_relabelings": [{str(k): v for k, v in mapping.items()} for mapping in all_v4_relabelings()],
            "projector_table": {str(endpoint): payload for endpoint, payload in projectors.items()},
            "nonidentity_projector_orbit": orbit_of_endpoint(11),
            "target_prime_vector_hits": target_hits,
            "accepted_candidate_count": 0,
            "proof_certificate": {
                "symmetry_step": "The intrinsic V4/GL(2,2) relabeling symmetry fixes only the identity and permutes 5, 7, and 11.",
                "feature_step": "Coefficient moments, sign counts, support size, norms, and denominator support give identical signatures for the nonidentity projectors.",
                "target_step": "The only way this audit singles out 11 is by asking for the already named target prime vector, which reimports the boundary datum rather than sourcing the projector.",
            },
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_intrinsic_character_projector_selector": False,
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
            "reason": "P2870 finds no label-free intrinsic character-projector selector for d=11.  The nonidentity idempotents for 5, 7, and 11 are tied by V4/GL(2,2) relabeling and by the audited finite intrinsic feature family.  Exact target matching selects 11 only by reusing the named boundary target.",
            "next_honest_step": "Do not replay intrinsic Aut-character projector scoring or target-vector matching as sourcehood.  A next proof-grade move must introduce a genuinely new strict law breaking the nonidentity V4 relabeling symmetry and supplying the unit-bearing coefficient/coupling theorem, or pivot to a different new typed object; otherwise preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["intrinsic_character_projector_selector_no_go_audit"]
    lines = [
        "# P2870/S1820 intrinsic character-projector selector no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite selector audit",
        f"- candidate class: `{audit['candidate_class']}`",
        f"- nonidentity projector orbit: `{audit['nonidentity_projector_orbit']}`",
        f"- target prime-vector hits: `{audit['target_prime_vector_hits']}`",
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
    payload = build_payload(read_json(P2869))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2870/S1820 intrinsic character-projector selector no-go audit",
        "## P2870/S1820 intrinsic character-projector selector no-go audit\n\n"
        "`P2870/S1820` audits whether the P2869 Aut-character idempotent layer can select `d=11` internally.  The finite V4/GL(2,2) relabeling computation ties the nonidentity projectors `5`, `7`, and `11`; coefficient moments, sign counts, support size, norms, and denominator support do not distinguish `11`.  Exact target-prime matching singles out `11` only by reusing the named boundary target `(9/5)log(11)`, so it is not an independent source law.  No boundary source law, eta source, strict damping bridge, `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2870/S1820 intrinsic character-projector selector `L_total` guard",
        "## P2870/S1820 intrinsic character-projector selector `L_total` guard\n\n"
        "`P2870/S1820` adds no strict action term.  Intrinsic character-projector scoring does not export a unit-bearing boundary/source density, coupling coefficient, localization/pullback theorem, variational chain rule, nonproxy `L_total`, EOM, or Hamiltonian source.\n",
    )
    append_once(
        AGENTS,
        "Current intrinsic character-projector selector no-go guardrail (P2870/S1820, 2026-06-18)",
        "## Current intrinsic character-projector selector no-go guardrail (P2870/S1820, 2026-06-18)\n\n"
        "- P2870 audits label-free intrinsic selection over the Aut(`Z12`)-character point idempotents after P2869.\n"
        "- The V4/GL(2,2) relabeling symmetry and the audited finite intrinsic features tie the nonidentity projectors `5`, `7`, and `11`; exact target matching selects `11` only by reimporting the named boundary target.\n"
        "- Do not promote intrinsic character-projector scoring, target-vector matching, Aut-character idempotent projectors, endpoint/polarity imports, prime-5 scaled coefficients, singleton endpoint representation, Dirichlet data, or log-scale harmonicity to strict damping/compression bridge, selector closure, role transfer, unit-bearing `L_total`, EOM, Hamiltonian, or ToE closure.\n"
        "- A next admissible proof-grade move must supply a genuinely new strict law breaking the nonidentity V4 relabeling symmetry and a unit-bearing coefficient/coupling theorem, or use a different new typed object; otherwise preserve no-new-live-frontier.\n",
    )
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
