#!/usr/bin/env python3
"""P2959/S1909: P2938 U(12) aggregate localizer acceptance no-go.

P2958 made the missing strict nadsoliton functor/localizer explicit.  P2959
attacks that object directly and computationally.  It does not rerun the P2958
K+C certificate as a conclusion, and it does not add another weight-family
variant.  Instead it enumerates a finite localizer target space around the
already-known P2938 decomposition:

    V(a,b) = a*K + b*C, with K=[1,2,0,0,0], C=[0,0,2,2,2]

for small nonnegative provenance-coefficients a,b.  The point is to test which
localizer predicates would actually select the desired aggregate K+C and
whether those predicates are strict-source data or target-coded conventions.

The computation finds that sum=9 is target-coded and still not unique in the bounded lattice, while K+C can be selected by an extra primitive equal-summand predicate.  But current artifacts do not export
those predicates as a strict nadsoliton functor/localizer.  Thus P2959 records a
bounded no-go for the currently available localizer predicates: they can certify
readiness, not provenance.
"""
from __future__ import annotations

import hashlib
import json
import math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2953_s1903_symbolic_weight_family_nonforcing_theorem import OUT as P2953
from p2958_s1908_p2938_strict_provenance_theorem_interface import OUT as P2958

OUT = GEN / "p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go.json"
MD = GEN / "p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

K = [1, 2, 0, 0, 0]
C = [0, 0, 2, 2, 2]
TARGET = [1, 2, 2, 2, 2]


def combine(a: int, b: int) -> list[int]:
    return [a * x + b * y for x, y in zip(K, C)]


def candidate_rows(max_coeff: int = 3) -> list[dict[str, Any]]:
    rows = []
    for a in range(max_coeff + 1):
        for b in range(max_coeff + 1):
            if a == 0 and b == 0:
                continue
            vector = combine(a, b)
            rows.append({
                "a_kernel_excess_weight": a,
                "b_character_negativity_weight": b,
                "vector": vector,
                "sum": sum(vector),
                "all_prime_coordinates_positive": all(x > 0 for x in vector),
                "equal_summand": a == b,
                "primitive_pair": math.gcd(a, b) == 1,
                "sum_9_target_predicate": sum(vector) == 9,
                "equals_target_K_plus_C": vector == TARGET,
            })
    return rows


def localizer_predicate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    predicates = [
        ("positive_coordinates", lambda r: r["all_prime_coordinates_positive"], False, "finite readiness only; selects every positive a,b pair"),
        ("sum_9", lambda r: r["sum_9_target_predicate"], False, "selects K+C here, but imports the target package sum instead of sourcing it"),
        ("equal_summand", lambda r: r["equal_summand"], False, "selects multiple diagonal multiples, not unique without a primitive/unit premise"),
        ("primitive_equal_summand", lambda r: r["equal_summand"] and r["primitive_pair"], False, "selects K+C in this bounded lattice, but current artifacts do not export this as a strict nadsoliton law"),
        ("positive_and_primitive_equal_summand", lambda r: r["all_prime_coordinates_positive"] and r["equal_summand"] and r["primitive_pair"], False, "same target as primitive_equal_summand; still an unsourced localizer premise"),
    ]
    out = []
    for name, pred, strict_exported, boundary in predicates:
        selected = [r for r in rows if pred(r)]
        out.append({
            "predicate": name,
            "selected_count": len(selected),
            "selected_pairs": [[r["a_kernel_excess_weight"], r["b_character_negativity_weight"]] for r in selected],
            "selects_exact_target_pair_only": len(selected) == 1 and selected[0]["equals_target_K_plus_C"],
            "strict_nadsoliton_localizer_exported": strict_exported,
            "boundary": boundary,
        })
    return out


def acceptance_obligation_rows(p2953: dict[str, Any], p2958: dict[str, Any], predicate_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    cert2958 = p2958["provenance_certificate"]
    primitive_equal_row = next(row for row in predicate_rows if row["predicate"] == "primitive_equal_summand")
    return [
        {
            "obligation": "finite_localizer_lattice_constructed",
            "satisfied": True,
            "evidence": "P2959 enumerates V(a,b)=a*K+b*C for 0<=a,b<=3, excluding (0,0)",
        },
        {
            "obligation": "target_pair_selectable_by_some_predicate",
            "satisfied": primitive_equal_row["selects_exact_target_pair_only"],
            "evidence": "primitive equal-summand selects K+C in the bounded localizer lattice",
        },
        {
            "obligation": "predicate_not_target_coded_or_conventional",
            "satisfied": False,
            "evidence": "sum=9 imports the package target; primitive equal-summand is an extra premise unless sourced",
        },
        {
            "obligation": "strict_nadsoliton_functor_exports_predicate",
            "satisfied": cert2958["strict_nadsoliton_functor_exported"],
            "evidence": "P2958 leaves the nadsoliton-to-U12 functor/localizer unexported",
        },
        {
            "obligation": "strict_equal_weight_source_theorem_exported",
            "satisfied": p2953["symbolic_weight_certificate"]["strict_equal_weight_source_theorem_exported"],
            "evidence": "P2953 proves target equivalence a=b=1 but no strict equal-weight source theorem",
        },
        {
            "obligation": "downstream_beta_unit_coupling_exported",
            "satisfied": cert2958["ratio_package_beta_unit_coupling_exported"],
            "evidence": "P2958/P2957 leave the beta/unit coupling absent",
        },
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = [
        "finite_localizer_lattice_constructed",
        "target_pair_selectable_by_some_predicate",
        "predicate_not_target_coded_or_conventional",
        "strict_nadsoliton_functor_exports_predicate",
        "strict_equal_weight_source_theorem_exported",
        "downstream_beta_unit_coupling_exported",
    ]
    matrix = []
    for mask in range(1 << len(names)):
        present = {name: bool(mask & (1 << i)) for i, name in enumerate(names)}
        matrix.append({
            "mask": mask,
            "present": present,
            "missing": [name for name, value in present.items() if not value],
            "accepts_strict_u12_aggregate_localizer": all(present.values()),
        })
    return matrix


def build_payload(p2953: dict[str, Any], p2958: dict[str, Any]) -> dict[str, Any]:
    candidates = candidate_rows()
    predicates = localizer_predicate_rows(candidates)
    obligations = acceptance_obligation_rows(p2953, p2958, predicates)
    matrix = acceptance_matrix()
    accepted = all(row["satisfied"] for row in obligations)
    exact_predicates = [row["predicate"] for row in predicates if row["selects_exact_target_pair_only"]]
    return {
        "status": "P2959_P2938_U12_AGGREGATE_LOCALIZER_ACCEPTANCE_NO_GO_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2953": hashlib.sha256(P2953.read_bytes()).hexdigest() if P2953.exists() else None,
            "P2958": hashlib.sha256(P2958.read_bytes()).hexdigest() if P2958.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "P2938_U12Aggregate_StrictLocalizer_AcceptanceNoGo",
            "localizer_candidate_rows": candidates,
            "localizer_predicate_rows": predicates,
            "acceptance_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "localizer_certificate": {
            "candidate_lattice_row_count": len(candidates),
            "exact_target_pair": [1, 1],
            "exact_target_vector": TARGET,
            "predicates_selecting_exact_target_only": exact_predicates,
            "all_exact_selecting_predicates_strict_exported": all(next(row for row in predicates if row["predicate"] == pred)["strict_nadsoliton_localizer_exported"] for pred in exact_predicates),
            "finite_localizer_lattice_constructed": obligations[0]["satisfied"],
            "target_pair_selectable_by_some_predicate": obligations[1]["satisfied"],
            "predicate_not_target_coded_or_conventional": obligations[2]["satisfied"],
            "strict_nadsoliton_functor_exports_predicate": obligations[3]["satisfied"],
            "strict_equal_weight_source_theorem_exported": obligations[4]["satisfied"],
            "downstream_beta_unit_coupling_exported": obligations[5]["satisfied"],
            "p2958_functor_localizer_obligation_discharged": accepted,
            "acceptance_matrix_row_count": len(matrix),
            "acceptance_matrix_accepted_row_count": sum(1 for row in matrix if row["accepts_strict_u12_aggregate_localizer"]),
        },
        "decision": {
            "positive_witnesses": {
                "bounded_localizer_lattice_constructed": True,
                "target_pair_selectable_by_formal_predicate": bool(exact_predicates),
                "localizer_predicate_boundaries_classified": True,
            },
            "negative_export_flags": {
                "strict_nadsoliton_to_u12_functor_localizer_exported": False,
                "strict_p2938_torsion_character_provenance_exported": False,
                "strict_equal_weight_source_theorem_exported": False,
                "strict_ratio_package_source_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2959 attacks the missing P2958 nadsoliton-to-U12 functor/localizer.  In the bounded localizer lattice, sum=9 is target-coded and nonunique, while the target pair (a,b)=(1,1) is selectable by primitive equal-summand predicates; none is strict unless a nadsoliton functor exports it.  Current artifacts export no such functor, no strict equal-weight theorem, and no beta/unit coupling.",
            "next_honest_step": "Do not continue bounded localizer-predicate variants, K+C decompositions, symbolic weight-family variants, beta-scale normalization, P2601 replay, scalar Euler insertion, or count/role aliases.  A next proof-grade move must either introduce a genuinely new strict nadsoliton localizer law not equivalent to target-coded sum=9 cuts or primitive equal-weight convention, construct a unit-bearing nonproxy coupling outside this localizer lane, or pivot outside the ratio-package lane while preserving the P2929-P2959 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["localizer_certificate"]
    lines = [
        "# P2959/S1909 P2938 U(12) aggregate localizer acceptance no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Localizer certificate",
        f"- candidate lattice rows: `{cert['candidate_lattice_row_count']}`",
        f"- exact target pair: `{cert['exact_target_pair']}`",
        f"- exact target vector: `{cert['exact_target_vector']}`",
        f"- predicates selecting exact target only: `{cert['predicates_selecting_exact_target_only']}`",
        f"- all exact-selecting predicates strict-exported: `{cert['all_exact_selecting_predicates_strict_exported']}`",
        f"- predicate not target-coded or conventional: `{cert['predicate_not_target_coded_or_conventional']}`",
        f"- strict nadsoliton functor exports predicate: `{cert['strict_nadsoliton_functor_exports_predicate']}`",
        f"- strict equal-weight source theorem exported: `{cert['strict_equal_weight_source_theorem_exported']}`",
        f"- downstream beta/unit coupling exported: `{cert['downstream_beta_unit_coupling_exported']}`",
        f"- P2958 functor/localizer obligation discharged: `{cert['p2958_functor_localizer_obligation_discharged']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_row_count']}/{cert['acceptance_matrix_accepted_row_count']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2953), read_json(P2958))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2959/S1909 P2938 U(12) aggregate localizer acceptance no-go", "## P2959/S1909 P2938 U(12) aggregate localizer acceptance no-go\n\n`P2959/S1909` attacks the missing P2958 nadsoliton-to-U(12) functor/localizer directly.  It enumerates the bounded localizer lattice `V(a,b)=a*K+b*C` for `0<=a,b<=3` and classifies predicates that select the target pair `(1,1)`.  Formal selection is possible through primitive equal-summand predicates, while `sum=9` is target-coded and nonunique in the bounded lattice; none of these predicates is strict unless a nadsoliton localizer exports it.  Current artifacts still export no such functor/localizer, no strict equal-weight theorem, and no beta/unit coupling.  Thus no strict P2938 provenance, ratio-package source theorem, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2959/S1909 U(12) localizer `L_total` guard", "## P2959/S1909 U(12) localizer `L_total` guard\n\n`P2959/S1909` shows that bounded formal localizer predicates either import the target-coded but nonunique `sum=9` condition or pick the P2938 target through an unsourced primitive equal-summand convention.  Without a strict nadsoliton functor/localizer and beta/unit coupling, the aggregate still cannot enter `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE as a sourced damping coefficient.\n")
    append_once(AGENTS, "Current P2938 U(12) aggregate localizer no-go guardrail (P2959/S1909, 2026-06-20)", "## Current P2938 U(12) aggregate localizer no-go guardrail (P2959/S1909, 2026-06-20)\n\n- P2959 attacks the missing P2958 nadsoliton-to-U(12) functor/localizer by enumerating bounded localizer predicates over `V(a,b)=a*K+b*C`, rather than replaying K+C decompositions, symbolic weight-family variants, beta-scale normalization, P2601 prose, scalar Euler insertion, count aliases, or role-signature routes.\n- Formal predicates either make target-coded/nonunique `sum=9` cuts or select `(a,b)=(1,1)` through unsourced primitive equal-summand assumptions; current artifacts do not export these as strict nadsoliton localizer laws.\n- Do not promote P2959 to strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce a genuinely new strict nadsoliton localizer law not equivalent to target-coded `sum=9` cuts or primitive equal-weight convention, construct a unit-bearing nonproxy coupling outside this localizer lane, or pivot outside the ratio-package lane while preserving the P2929-P2959 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
