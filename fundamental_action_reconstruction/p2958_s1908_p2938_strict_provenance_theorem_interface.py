#!/usr/bin/env python3
"""P2958/S1908: P2938 strict provenance theorem interface.

After P2957, the honest ratio-package frontier is no longer another beta-scale
normalization or Euler insertion.  P2958 attacks the other recommended route:
strict provenance for the P2938 torsion-character aggregate.  It does not
repeat P2952/P2953 weight-family nonforcing.  Instead it constructs the missing
provenance theorem interface whose job would be to turn the already-computed
finite aggregate V=[1,2,2,2,2] into a strict nadsoliton-sourced object.

The interface separates four typed layers:
  1. finite aggregate arithmetic (already positive),
  2. decomposition provenance for kernel-excess K and character-negativity C,
  3. an exported strict functor/localizer from nadsoliton structure to K+C,
  4. downstream coupling to the exact ratio package and beta/unit data.

The finite checks verify the target decomposition V=K+C and the existing
product-additive carrier.  The obstruction is not arithmetic: current artifacts
still do not export a strict nadsoliton functor/localizer selecting this U(12)
aggregate as a source law, do not prove a nonconventional equal-summand law as
strict provenance, and do not provide the downstream beta/unit coupling blocked
by P2957.
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
from p2953_s1903_symbolic_weight_family_nonforcing_theorem import OUT as P2953
from p2957_s1907_positive_beta_scale_unit_source_obstruction import OUT as P2957

OUT = GEN / "p2958_s1908_p2938_strict_provenance_theorem_interface.json"
MD = GEN / "p2958_s1908_p2938_strict_provenance_theorem_interface.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
PRIMES = [2, 3, 5, 7, 11]


def vector_add(a: list[int], b: list[int]) -> list[int]:
    return [x + y for x, y in zip(a, b)]


def dot(a: list[int], b: list[int]) -> int:
    return sum(x * y for x, y in zip(a, b))


def decomposition_rows(p2938: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2938["constructed_theoretical_objects"]["prime_coordinate_rows"]:
        k = row["kernel_excess"]
        c = row["unit_character_negativity_count"]
        rows.append({
            "prime": row["prime"],
            "kernel_excess_component": k,
            "unit_character_negativity_component": c,
            "combined_coordinate": row["combined_coordinate"],
            "decomposition_passes": k + c == row["combined_coordinate"],
            "strict_provenance_exported_for_prime_row": False,
        })
    return rows


def finite_provenance_square(rows: list[dict[str, Any]]) -> dict[str, Any]:
    k = [row["kernel_excess_component"] for row in rows]
    c = [row["unit_character_negativity_component"] for row in rows]
    v = [row["combined_coordinate"] for row in rows]
    return {
        "primes": PRIMES,
        "kernel_excess_vector_K": k,
        "character_negativity_vector_C": c,
        "target_vector_V": v,
        "K_plus_C": vector_add(k, c),
        "K_plus_C_equals_V": vector_add(k, c) == v,
        "sum_V": sum(v),
        "eta_from_sum_over_five_primes": {"numerator": sum(v), "denominator": 5, "as_string": f"{sum(v)}/5"},
        "K_dot_C": dot(k, c),
        "support_overlap_count": sum(1 for x, y in zip(k, c) if x and y),
        "interpretation": "K and C live on disjoint prime supports, so K+C is arithmetically sharp; this still does not make the sum a strict nadsoliton source law.",
    }


def provenance_obligation_rows(p2938: dict[str, Any], p2953: dict[str, Any], p2957: dict[str, Any]) -> list[dict[str, Any]]:
    cert2938 = p2938["candidate_certificate"]
    cert2957 = p2957["positive_beta_scale_unit_certificate"]
    return [
        {
            "obligation": "finite_U12_character_aggregate_constructed",
            "satisfied": cert2938["prime_coordinate_vector_order_2_3_5_7_11"] == [1, 2, 2, 2, 2],
            "evidence": "P2938 constructs the unit-character enriched aggregate vector [1,2,2,2,2]",
        },
        {
            "obligation": "product_additive_carrier_available",
            "satisfied": cert2938["product_additivity_defect_count"] == 0,
            "evidence": "P2938 product-additivity defects vanish on audited d*e<=11 products",
        },
        {
            "obligation": "strict_nadsoliton_functor_to_U12_aggregate_exported",
            "satisfied": False,
            "evidence": "current artifacts do not export a strict map from nadsoliton ontology to this U(12) character aggregate",
        },
        {
            "obligation": "nonconventional_equal_summand_law_exported",
            "satisfied": p2953["symbolic_weight_certificate"]["strict_equal_weight_source_theorem_exported"],
            "evidence": "P2953 shows a=b=1 is target-equation equivalence, not an exported equal-weight source theorem",
        },
        {
            "obligation": "aggregate_localizer_not_just_finite_carrier",
            "satisfied": False,
            "evidence": "no theorem localizes why exactly K+C, rather than another positive family member, is nadsoliton-sourced",
        },
        {
            "obligation": "ratio_package_beta_unit_coupling_exported",
            "satisfied": cert2957["p2951_positive_beta_scale_unit_atom_discharged"],
            "evidence": "P2957 leaves beta/unit source and unit-bearing nonproxy coupling unexported",
        },
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = [
        "finite_U12_character_aggregate_constructed",
        "product_additive_carrier_available",
        "strict_nadsoliton_functor_to_U12_aggregate_exported",
        "nonconventional_equal_summand_law_exported",
        "aggregate_localizer_not_just_finite_carrier",
        "ratio_package_beta_unit_coupling_exported",
    ]
    matrix = []
    for mask in range(1 << len(names)):
        present = {name: bool(mask & (1 << i)) for i, name in enumerate(names)}
        matrix.append({
            "mask": mask,
            "present": present,
            "missing": [name for name, value in present.items() if not value],
            "accepts_strict_p2938_provenance_theorem": all(present.values()),
        })
    return matrix


def build_payload(p2938: dict[str, Any], p2953: dict[str, Any], p2957: dict[str, Any]) -> dict[str, Any]:
    decomp = decomposition_rows(p2938)
    square = finite_provenance_square(decomp)
    obligations = provenance_obligation_rows(p2938, p2953, p2957)
    matrix = acceptance_matrix()
    current_accepts = all(row["satisfied"] for row in obligations)
    return {
        "status": "P2958_P2938_STRICT_PROVENANCE_THEOREM_INTERFACE_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2938": hashlib.sha256(P2938.read_bytes()).hexdigest() if P2938.exists() else None,
            "P2953": hashlib.sha256(P2953.read_bytes()).hexdigest() if P2953.exists() else None,
            "P2957": hashlib.sha256(P2957.read_bytes()).hexdigest() if P2957.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "P2938_StrictNadsolitonProvenance_TheoremInterface",
            "decomposition_rows": decomp,
            "finite_provenance_square": square,
            "provenance_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "provenance_certificate": {
            "target_vector_V": square["target_vector_V"],
            "K_plus_C_equals_V": square["K_plus_C_equals_V"],
            "sum_V": square["sum_V"],
            "eta_from_sum_over_five_primes": square["eta_from_sum_over_five_primes"],
            "product_additive_carrier_available": obligations[1]["satisfied"],
            "strict_nadsoliton_functor_exported": False,
            "nonconventional_equal_summand_law_exported": obligations[3]["satisfied"],
            "aggregate_localizer_exported": False,
            "ratio_package_beta_unit_coupling_exported": obligations[5]["satisfied"],
            "p2951_p2938_strict_provenance_atom_discharged": current_accepts,
            "acceptance_matrix_row_count": len(matrix),
            "acceptance_matrix_accepted_row_count": sum(1 for row in matrix if row["accepts_strict_p2938_provenance_theorem"]),
        },
        "lay_summary": {
            "positive_progress": "There is real positive progress: the finite aggregate is arithmetically rigid enough that K+C exactly gives [1,2,2,2,2], sum 9, hence eta=9/5 at the package level.",
            "not_a_breakthrough": "It is not yet a breakthrough/closure, because this proves a clean candidate interface rather than a strict nadsoliton source law.  The missing step is a theorem explaining why the nadsoliton itself must generate this aggregate and coupling, not merely that the aggregate works once chosen.",
        },
        "decision": {
            "positive_witnesses": {
                "finite_decomposition_K_plus_C_verified": True,
                "product_additive_carrier_reused": True,
                "strict_provenance_acceptance_matrix_constructed": True,
            },
            "negative_export_flags": {
                "strict_p2938_torsion_character_provenance_exported": False,
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
            "reason": "P2958 attacks strict provenance for the P2938 torsion-character aggregate.  The finite arithmetic is positive: K+C exactly recovers [1,2,2,2,2], sum 9, and eta=9/5, with product-additive carrier support.  But current artifacts still lack the strict nadsoliton functor/localizer selecting this U(12) aggregate as a source law, the nonconventional equal-summand provenance theorem blocked by P2953, and the beta/unit coupling blocked by P2957.",
            "next_honest_step": "Do not continue finite K+C decompositions, symbolic weight-family variants, beta-scale normalization, P2601 replay, scalar Euler insertion, or count/role aliases.  A next proof-grade move must either export the actual strict nadsoliton functor/localizer from ontology to the P2938 U(12) aggregate, construct a genuinely unit-bearing nonproxy coupling that can receive the aggregate, or pivot outside the ratio-package lane while preserving the P2929-P2958 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["provenance_certificate"]
    lines = [
        "# P2958/S1908 P2938 strict provenance theorem interface",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Provenance certificate",
        f"- target vector V: `{cert['target_vector_V']}`",
        f"- K+C equals V: `{cert['K_plus_C_equals_V']}`",
        f"- sum V: `{cert['sum_V']}`",
        f"- eta from sum over five primes: `{cert['eta_from_sum_over_five_primes']['as_string']}`",
        f"- product-additive carrier available: `{cert['product_additive_carrier_available']}`",
        f"- strict nadsoliton functor exported: `{cert['strict_nadsoliton_functor_exported']}`",
        f"- nonconventional equal-summand law exported: `{cert['nonconventional_equal_summand_law_exported']}`",
        f"- aggregate localizer exported: `{cert['aggregate_localizer_exported']}`",
        f"- ratio-package beta/unit coupling exported: `{cert['ratio_package_beta_unit_coupling_exported']}`",
        f"- P2951 P2938 strict provenance atom discharged: `{cert['p2951_p2938_strict_provenance_atom_discharged']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_row_count']}/{cert['acceptance_matrix_accepted_row_count']}`",
        "",
        "## Lay summary",
        payload["lay_summary"]["positive_progress"],
        "",
        payload["lay_summary"]["not_a_breakthrough"],
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2938), read_json(P2953), read_json(P2957))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2958/S1908 P2938 strict provenance theorem interface", "## P2958/S1908 P2938 strict provenance theorem interface\n\n`P2958/S1908` attacks the remaining P2951 strict-provenance atom for the P2938 torsion-character aggregate.  The positive finite side is real: the decomposition `K=[1,2,0,0,0]`, `C=[0,0,2,2,2]` satisfies `K+C=[1,2,2,2,2]`, sum `9`, and package-level `eta=9/5`, while the P2938 product-additive carrier remains available.  The obstruction is source-theoretic rather than arithmetic: no strict nadsoliton functor/localizer selects the U(12) aggregate as a source law, no nonconventional equal-summand provenance theorem is exported, and P2957 still blocks beta/unit coupling.  Thus no strict ratio-package source theorem, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2958/S1908 P2938 provenance `L_total` guard", "## P2958/S1908 P2938 provenance `L_total` guard\n\n`P2958/S1908` verifies the finite `K+C` provenance interface for the P2938 aggregate, but it does not export the strict nadsoliton functor/localizer, nonconventional equal-summand law, or beta/unit coupling required for a sourced damping coefficient.  Therefore the aggregate still cannot enter `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE as a strict source.\n")
    append_once(AGENTS, "Current P2938 strict provenance theorem-interface guardrail (P2958/S1908, 2026-06-20)", "## Current P2938 strict provenance theorem-interface guardrail (P2958/S1908, 2026-06-20)\n\n- P2958 attacks strict provenance for the P2938 torsion-character aggregate rather than replaying finite K+C decompositions, symbolic weight-family variants, beta-scale normalization, P2601 prose, scalar Euler insertion, count aliases, or role-signature routes.\n- The finite arithmetic is positive: `K+C=[1,2,2,2,2]`, sum `9`, and package-level `eta=9/5` are verified with the P2938 product-additive carrier.\n- This is meaningful proof progress but not closure: no strict nadsoliton functor/localizer selects the U(12) aggregate as a source law, no nonconventional equal-summand provenance theorem is exported, and P2957 still blocks beta/unit coupling.\n- Do not promote P2958 to strict ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.  A next admissible move must export the actual strict nadsoliton functor/localizer, construct a genuinely unit-bearing nonproxy coupling that can receive the aggregate, or pivot outside the ratio-package lane while preserving the P2929-P2958 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
