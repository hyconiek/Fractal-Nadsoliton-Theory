#!/usr/bin/env python3
"""P2952/S1902: torsion-character aggregate weight-source obstruction.

P2951 says the next proof-grade move must supply one concrete missing atom.  This
packet attacks the first atom: strict provenance for the P2938 torsion-character
aggregate.  It does not add a new ratio scan, count alias, beta-scale sample, or
another normal-form lattice.

The P2938 carrier is a sum of two finite ingredients:

    K_p = gcd(p,12)-1                         (multiplication-kernel excess)
    C_p = #{nontrivial chi:U(12)->±1, chi(p)=-1}  (unit-character negativity)

P2938 uses equal weights V_p = 1*K_p + 1*C_p.  P2952 constructs the missing
weight-source theorem interface by auditing the parametric family
V_p(a,b)=a*K_p+b*C_p.  The exact target [1,2,2,2,2] forces a=b=1, but the finite
ingredients alone do not export a strict theorem selecting that equal-weight
normalization.  Therefore the P2938 provenance atom remains undischarged.
"""
from __future__ import annotations

import hashlib
import itertools
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate import OUT as P2938, PRIMES
from p2951_s1901_ratio_package_strict_source_normal_form_lattice import OUT as P2951

OUT = GEN / "p2952_s1902_torsion_character_aggregate_weight_source_obstruction.json"
MD = GEN / "p2952_s1902_torsion_character_aggregate_weight_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_VECTOR = [1, 2, 2, 2, 2]
WEIGHT_DOMAIN = [1, 2, 3, 4]


def ingredient_rows(p2938: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "prime": row["prime"],
            "kernel_excess": row["kernel_excess"],
            "unit_character_negativity_count": row["unit_character_negativity_count"],
            "p2938_equal_weight_coordinate": row["combined_coordinate"],
        }
        for row in p2938["constructed_theoretical_objects"]["prime_coordinate_rows"]
    ]


def vector_for_weights(rows: list[dict[str, Any]], a: int, b: int) -> list[int]:
    return [a * row["kernel_excess"] + b * row["unit_character_negativity_count"] for row in rows]


def weight_family_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    family = []
    for a, b in itertools.product(WEIGHT_DOMAIN, repeat=2):
        vector = vector_for_weights(rows, a, b)
        family.append({
            "kernel_excess_weight_a": a,
            "character_negativity_weight_b": b,
            "prime_vector_order_2_3_5_7_11": vector,
            "vector_sum": sum(vector),
            "matches_p2938_target_vector": vector == TARGET_VECTOR,
            "all_coordinates_positive": all(value > 0 for value in vector),
            "product_additivity_inherited_by_prime_exponent_extension": True,
            "strict_weight_source_exported": False,
        })
    return family


def coefficient_equation_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    equations = []
    for row, target in zip(rows, TARGET_VECTOR):
        equations.append({
            "prime": row["prime"],
            "equation": f"{row['kernel_excess']}*a + {row['unit_character_negativity_count']}*b = {target}",
            "satisfied_by_equal_weights_a1_b1": row["kernel_excess"] + row["unit_character_negativity_count"] == target,
        })
    return equations


def theorem_obligation_rows() -> list[dict[str, Any]]:
    return [
        {
            "obligation": "finite_two_ingredient_decomposition_available",
            "satisfied": True,
            "evidence": "P2938 exposes kernel_excess and unit_character_negativity_count separately for each audited prime",
        },
        {
            "obligation": "target_vector_forces_equal_weights_within_linear_family",
            "satisfied": True,
            "evidence": "nonunit primes force a=1 and unit-prime coordinates force b=1 for [1,2,2,2,2]",
        },
        {
            "obligation": "strict_equal_weight_source_theorem_exported",
            "satisfied": False,
            "evidence": "current artifacts provide the finite sum K_p+C_p but not a nadsoliton theorem selecting equal weights over other positive weight pairs",
        },
        {
            "obligation": "p2951_torsion_character_provenance_atom_discharged",
            "satisfied": False,
            "evidence": "the P2951 provenance atom requires a strict source for the aggregate, including its relative weights",
        },
    ]


def build_payload(p2938: dict[str, Any], p2951: dict[str, Any]) -> dict[str, Any]:
    ingredients = ingredient_rows(p2938)
    family = weight_family_rows(ingredients)
    target_rows = [row for row in family if row["matches_p2938_target_vector"]]
    obligations = theorem_obligation_rows()
    accepted = all(row["satisfied"] for row in obligations)
    return {
        "status": "P2952_TORSION_CHARACTER_AGGREGATE_WEIGHT_SOURCE_OBSTRUCTION_NO_STRICT_PROVENANCE",
        "input_hashes": {
            "P2938": hashlib.sha256(P2938.read_bytes()).hexdigest() if P2938.exists() else None,
            "P2951": hashlib.sha256(P2951.read_bytes()).hexdigest() if P2951.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "TorsionCharacter_Aggregate_EqualWeight_Source_TheoremInterface",
            "ingredient_rows": ingredients,
            "coefficient_equation_rows": coefficient_equation_rows(ingredients),
            "bounded_positive_weight_family_rows": family,
            "theorem_obligation_rows": obligations,
        },
        "weight_source_certificate": {
            "weight_domain": WEIGHT_DOMAIN,
            "family_row_count": len(family),
            "all_family_rows_product_additive_by_construction": all(row["product_additivity_inherited_by_prime_exponent_extension"] for row in family),
            "target_weight_pair_count_in_domain": len(target_rows),
            "target_weight_pairs": [[row["kernel_excess_weight_a"], row["character_negativity_weight_b"]] for row in target_rows],
            "target_vector_forces_equal_weights_in_a_b_family": len(target_rows) == 1 and target_rows[0]["kernel_excess_weight_a"] == 1 and target_rows[0]["character_negativity_weight_b"] == 1,
            "strict_equal_weight_source_theorem_exported": False,
            "p2951_torsion_character_provenance_atom_discharged": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "two_ingredient_decomposition_audited": True,
                "parametric_weight_family_constructed": True,
                "target_equal_weight_pair_identified_finitely": True,
            },
            "negative_export_flags": {
                "strict_equal_weight_source_theorem_exported": False,
                "strict_torsion_character_source_theorem_exported": False,
                "strict_ratio_package_source_theorem_exported": False,
                "strict_delta_eta_source_law_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2952 attacks the P2951 torsion-character provenance atom by splitting the P2938 aggregate into kernel-excess and character-negativity ingredients.  The exact P2938 vector forces the equal-weight pair a=b=1 inside V_p(a,b)=a*K_p+b*C_p, but current artifacts do not export a strict nadsoliton theorem selecting that equal-weight aggregation over the positive weight family.  The provenance atom therefore remains undischarged.",
            "next_honest_step": "Do not add more bounded weight scans or alternative linear aggregate weights.  A next proof-grade move must export an actual strict equal-weight/source theorem for K_p+C_p, or attack a different P2951 atom with a concrete theorem object; otherwise pivot outside the ratio-package lane and preserve the P2929-P2952 no-strict-provenance boundary.",
            "p2951_current_missing_atoms": p2951["normal_form_certificate"]["current_missing_atoms"],
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["weight_source_certificate"]
    lines = [
        "# P2952/S1902 torsion-character aggregate weight-source obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Weight-source certificate",
        f"- weight domain: `{cert['weight_domain']}`",
        f"- family rows: `{cert['family_row_count']}`",
        f"- all rows product-additive by construction: `{cert['all_family_rows_product_additive_by_construction']}`",
        f"- target weight pair count in domain: `{cert['target_weight_pair_count_in_domain']}`",
        f"- target weight pairs: `{cert['target_weight_pairs']}`",
        f"- target vector forces equal weights: `{cert['target_vector_forces_equal_weights_in_a_b_family']}`",
        f"- strict equal-weight source theorem exported: `{cert['strict_equal_weight_source_theorem_exported']}`",
        f"- P2951 torsion-character provenance atom discharged: `{cert['p2951_torsion_character_provenance_atom_discharged']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2938), read_json(P2951))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2952/S1902 torsion-character aggregate weight-source obstruction", "## P2952/S1902 torsion-character aggregate weight-source obstruction\n\n`P2952/S1902` attacks the P2951 torsion-character provenance atom by decomposing the P2938 aggregate into multiplication-kernel excess `K_p=gcd(p,12)-1` and unit-character negativity `C_p`.  The parametric family `V_p(a,b)=a*K_p+b*C_p` remains product-additive by prime-exponent extension; the exact target `[1,2,2,2,2]` forces the equal-weight pair `a=b=1` inside the audited positive family.  But current artifacts do not export a strict nadsoliton theorem selecting that equal-weight aggregation, so no strict provenance, ratio-package source theorem, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2952/S1902 torsion-character weight-source `L_total` guard", "## P2952/S1902 torsion-character weight-source `L_total` guard\n\n`P2952/S1902` shows that the finite P2938 aggregate still contains an unsourced equal-weight premise between kernel-excess and character-negativity ingredients.  Even though `a=b=1` is the unique audited way to recover `[1,2,2,2,2]` in the two-weight family, uniqueness relative to a target is not a strict source theorem.  This cannot enter `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE without an independent equal-weight/source theorem.\n")
    append_once(AGENTS, "Current torsion-character aggregate weight-source guardrail (P2952/S1902, 2026-06-19)", "## Current torsion-character aggregate weight-source guardrail (P2952/S1902, 2026-06-19)\n\n- P2952 attacks one concrete P2951 atom: strict provenance for the P2938 torsion-character aggregate.\n- The decomposition `V_p(a,b)=a*(gcd(p,12)-1)+b*C_p` shows the exact target vector `[1,2,2,2,2]` forces `a=b=1` inside the audited positive two-weight family, but current artifacts do not source the equal-weight aggregation theorem.\n- This does not export strict torsion-character provenance, strict ratio-package source, beta/eta coupling, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- Do not continue bounded weight scans or alternative linear aggregate weights as primary strategy.  A next admissible move must export a strict equal-weight/source theorem, attack a different P2951 atom with a concrete theorem object, or pivot outside this ratio-package lane while preserving the P2929-P2952 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
