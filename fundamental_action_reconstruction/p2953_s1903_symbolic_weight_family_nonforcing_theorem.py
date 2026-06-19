#!/usr/bin/env python3
"""P2953/S1903: symbolic weight-family nonforcing theorem.

P2952 exposed an equal-weight premise in the P2938 torsion-character aggregate
using a bounded positive two-weight audit.  P2953 removes the bounded-scan caveat
and replaces it with an exact symbolic theorem over positive weights.

Let
  K = [1,2,0,0,0]  from multiplication-kernel excess, and
  C = [0,0,2,2,2]  from unit-character negativity.
For any positive rational weights a,b, the carrier
  V(a,b)=a*K+b*C = [a,2a,2b,2b,2b]
remains product-additive by prime-exponent extension and has positive prime
coordinates.  The P2938 target [1,2,2,2,2] is equivalent to a=1 and b=1, but
that is a target equation, not a strict nadsoliton source law selecting equal
weights.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate import OUT as P2938
from p2952_s1902_torsion_character_aggregate_weight_source_obstruction import OUT as P2952

OUT = GEN / "p2953_s1903_symbolic_weight_family_nonforcing_theorem.json"
MD = GEN / "p2953_s1903_symbolic_weight_family_nonforcing_theorem.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_VECTOR = [1, 2, 2, 2, 2]


def frac_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "as_string": f"{value.numerator}/{value.denominator}",
        "as_float": float(value),
    }


def ingredient_vectors(p2938: dict[str, Any]) -> dict[str, list[int]]:
    rows = p2938["constructed_theoretical_objects"]["prime_coordinate_rows"]
    return {
        "kernel_excess_vector_K_order_2_3_5_7_11": [row["kernel_excess"] for row in rows],
        "character_negativity_vector_C_order_2_3_5_7_11": [row["unit_character_negativity_count"] for row in rows],
        "p2938_target_vector_order_2_3_5_7_11": [row["combined_coordinate"] for row in rows],
    }


def symbolic_theorem_rows(vectors: dict[str, list[int]]) -> list[dict[str, Any]]:
    k_vec = vectors["kernel_excess_vector_K_order_2_3_5_7_11"]
    c_vec = vectors["character_negativity_vector_C_order_2_3_5_7_11"]
    return [
        {
            "theorem_row": "symbolic_weight_family_closed_form",
            "statement": "For positive rational a,b, V(a,b)=a*K+b*C=[a,2a,2b,2b,2b].",
            "verified": k_vec == [1, 2, 0, 0, 0] and c_vec == [0, 0, 2, 2, 2],
        },
        {
            "theorem_row": "positive_cone_for_all_positive_weights",
            "statement": "If a>0 and b>0, every coordinate of [a,2a,2b,2b,2b] is positive.",
            "verified": True,
        },
        {
            "theorem_row": "product_additivity_for_all_weights",
            "statement": "Prime-exponent extension makes V(a,b) additive on audited products for all a,b because it is linear in prime exponents.",
            "verified": True,
        },
        {
            "theorem_row": "target_equivalence",
            "statement": "V(a,b)=[1,2,2,2,2] iff a=1 and b=1, using the p=2 and p=5 coordinates.",
            "verified": True,
        },
        {
            "theorem_row": "strict_equal_weight_source_export",
            "statement": "A strict theorem selects a=1 and b=1 independently of target-fitting.",
            "verified": False,
        },
    ]


def rational_witness_rows() -> list[dict[str, Any]]:
    witnesses = [(Fraction(1, 2), Fraction(1, 1)), (Fraction(1, 1), Fraction(3, 2)), (Fraction(2, 1), Fraction(1, 3))]
    rows = []
    for a, b in witnesses:
        vector = [a, 2 * a, 2 * b, 2 * b, 2 * b]
        rows.append({
            "a": frac_payload(a),
            "b": frac_payload(b),
            "vector_order_2_3_5_7_11": [frac_payload(value) for value in vector],
            "positive_coordinates": all(value > 0 for value in vector),
            "matches_p2938_target": vector == [Fraction(value, 1) for value in TARGET_VECTOR],
            "product_additive_by_symbolic_theorem": True,
        })
    return rows


def theorem_obligation_rows() -> list[dict[str, Any]]:
    return [
        {
            "obligation": "unbounded_positive_weight_family_symbolically_closed",
            "satisfied": True,
            "evidence": "V(a,b)=[a,2a,2b,2b,2b] for all positive rational a,b",
        },
        {
            "obligation": "p2938_target_equivalent_to_equal_weights",
            "satisfied": True,
            "evidence": "target coordinates force a=1 from p=2 and b=1 from p=5",
        },
        {
            "obligation": "strict_equal_weight_source_theorem_exported",
            "satisfied": False,
            "evidence": "target equivalence does not supply a nadsoliton source law selecting equal weights",
        },
        {
            "obligation": "p2951_torsion_character_provenance_atom_discharged",
            "satisfied": False,
            "evidence": "the provenance atom still needs a strict source for K+C rather than target-fit equal weights",
        },
    ]


def build_payload(p2938: dict[str, Any], p2952: dict[str, Any]) -> dict[str, Any]:
    vectors = ingredient_vectors(p2938)
    theorem_rows = symbolic_theorem_rows(vectors)
    witness_rows = rational_witness_rows()
    obligations = theorem_obligation_rows()
    accepted = all(row["satisfied"] for row in obligations)
    return {
        "status": "P2953_SYMBOLIC_WEIGHT_FAMILY_NONFORCING_THEOREM_NO_STRICT_SOURCE",
        "input_hashes": {
            "P2938": hashlib.sha256(P2938.read_bytes()).hexdigest() if P2938.exists() else None,
            "P2952": hashlib.sha256(P2952.read_bytes()).hexdigest() if P2952.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "Symbolic_TorsionCharacter_WeightFamily_Nonforcing_Theorem",
            "ingredient_vectors": vectors,
            "symbolic_theorem_rows": theorem_rows,
            "positive_rational_witness_rows": witness_rows,
            "theorem_obligation_rows": obligations,
        },
        "symbolic_weight_certificate": {
            "closed_form": "V(a,b)=[a,2a,2b,2b,2b]",
            "unbounded_positive_rational_family_constructed": True,
            "sample_non_target_positive_witness_count": sum(1 for row in witness_rows if row["positive_coordinates"] and not row["matches_p2938_target"]),
            "target_vector_equivalent_to_a_equals_1_b_equals_1": True,
            "p2952_bounded_scan_caveat_removed": True,
            "strict_equal_weight_source_theorem_exported": False,
            "p2951_torsion_character_provenance_atom_discharged": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "symbolic_unbounded_family_constructed": True,
                "target_equal_weight_equivalence_proved": True,
                "bounded_scan_caveat_removed": True,
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
            "reason": "P2953 removes the bounded-weight caveat from P2952 by proving the exact symbolic family V(a,b)=[a,2a,2b,2b,2b] for all positive rational weights.  The target P2938 vector is equivalent to a=b=1, but this is still only a target-equation equivalence; no strict nadsoliton theorem selects the equal-weight aggregation K+C as a source law.",
            "next_honest_step": "Do not continue weight-family variants, bounded scans, or target-equivalence restatements.  A next proof-grade move must either supply an independent strict equal-weight/source theorem for K+C, attack a different P2951 atom with new source data, or pivot outside the ratio-package lane while preserving the P2929-P2953 no-strict-source boundary.",
            "p2952_target_pairs": p2952["weight_source_certificate"]["target_weight_pairs"],
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["symbolic_weight_certificate"]
    lines = [
        "# P2953/S1903 symbolic weight-family nonforcing theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Symbolic weight certificate",
        f"- closed form: `{cert['closed_form']}`",
        f"- unbounded positive rational family constructed: `{cert['unbounded_positive_rational_family_constructed']}`",
        f"- sample non-target positive witnesses: `{cert['sample_non_target_positive_witness_count']}`",
        f"- target vector equivalent to a=1,b=1: `{cert['target_vector_equivalent_to_a_equals_1_b_equals_1']}`",
        f"- P2952 bounded-scan caveat removed: `{cert['p2952_bounded_scan_caveat_removed']}`",
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
    payload = build_payload(read_json(P2938), read_json(P2952))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2953/S1903 symbolic weight-family nonforcing theorem", "## P2953/S1903 symbolic weight-family nonforcing theorem\n\n`P2953/S1903` removes the bounded-weight caveat from P2952.  For the P2938 ingredients `K=[1,2,0,0,0]` and `C=[0,0,2,2,2]`, the exact positive rational family is `V(a,b)=a*K+b*C=[a,2a,2b,2b,2b]`; it is product-additive by prime-exponent extension for all positive `a,b`.  The target vector `[1,2,2,2,2]` is equivalent to `a=b=1`, but that equivalence is target-fitting, not an exported strict equal-weight/source theorem.  Therefore no strict torsion-character provenance, ratio-package source theorem, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2953/S1903 symbolic weight-family `L_total` guard", "## P2953/S1903 symbolic weight-family `L_total` guard\n\n`P2953/S1903` proves the unbounded positive rational weight-family obstruction exactly: `V(a,b)=[a,2a,2b,2b,2b]` stays product-additive and positive for all positive weights, while the desired P2938 vector merely imposes the target equations `a=b=1`.  Since no independent equal-weight/source theorem is exported, this cannot enter `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current symbolic weight-family nonforcing guardrail (P2953/S1903, 2026-06-19)", "## Current symbolic weight-family nonforcing guardrail (P2953/S1903, 2026-06-19)\n\n- P2953 removes the bounded-scan caveat from P2952: the P2938 ingredient family is exactly `V(a,b)=[a,2a,2b,2b,2b]` for all positive rational weights.\n- The target vector `[1,2,2,2,2]` is equivalent to `a=b=1`, but this is target-equation equivalence rather than a strict nadsoliton source law selecting equal weights.\n- No strict torsion-character provenance, ratio-package source, beta/eta coupling, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n- Do not continue weight-family variants, bounded scans, or target-equivalence restatements as primary strategy.  A next admissible move must export an independent strict equal-weight/source theorem, attack a different P2951 atom with new source data, or pivot outside this ratio-package lane while preserving the P2929-P2953 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
