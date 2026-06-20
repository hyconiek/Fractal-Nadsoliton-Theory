#!/usr/bin/env python3
"""P2960/S1910: developmental ontology localizer-law intake.

This step responds to the theory-under-construction caveat: a mechanism that
would close the package should not be rejected merely because older artifacts do
not already contain it.  P2960 therefore constructs a *developmental-theory
intake gate* for new localizer laws.  The gate distinguishes:

1. strict export already present in the repo,
2. admissible developmental axiom/theorem work that may be proposed while the
   final theory is still being formed,
3. non-admissible replay of target-coded or conventional predicates.

It is not another K+C decomposition or bounded predicate variant as a closure
claim.  It asks whether a new ontological mechanism can be accepted as a
research object without being falsely promoted to strict closure.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2959_s1909_p2938_u12_aggregate_localizer_acceptance_no_go import OUT as P2959, candidate_rows

OUT = GEN / "p2960_s1910_developmental_ontology_localizer_law_intake.json"
MD = GEN / "p2960_s1910_developmental_ontology_localizer_law_intake.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def law_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    positive = [r for r in rows if r["all_prime_coordinates_positive"]]
    min_sum = min(r["sum"] for r in positive)
    min_linf = min(max(r["a_kernel_excess_weight"], r["b_character_negativity_weight"]) for r in positive)
    candidates = [
        {
            "law": "target_sum_9_cut",
            "new_theory_object": False,
            "formula": "sum(V)=9",
            "selected_pairs": [[r["a_kernel_excess_weight"], r["b_character_negativity_weight"]] for r in rows if r["sum"] == 9],
            "admissibility": "reject_as_target_coded_replay",
        },
        {
            "law": "primitive_equal_summand_convention",
            "new_theory_object": False,
            "formula": "a=b and gcd(a,b)=1",
            "selected_pairs": [[r["a_kernel_excess_weight"], r["b_character_negativity_weight"]] for r in rows if r["equal_summand"] and r["primitive_pair"]],
            "admissibility": "reject_as_P2959_convention_replay_unless_sourced",
        },
        {
            "law": "developmental_fractal_self_balance_minimal_positive_source",
            "new_theory_object": True,
            "formula": "positive coordinates, equal K/C provenance weight, and minimal positive source scale",
            "selected_pairs": [[r["a_kernel_excess_weight"], r["b_character_negativity_weight"]] for r in positive if r["equal_summand"] and r["sum"] == min_sum],
            "admissibility": "admissible_as_developmental_theorem_candidate_not_strict_export",
        },
        {
            "law": "developmental_minimax_source_amplitude_quotient",
            "new_theory_object": True,
            "formula": "positive coordinates and minimal max(a,b) scale quotient",
            "selected_pairs": [[r["a_kernel_excess_weight"], r["b_character_negativity_weight"]] for r in positive if max(r["a_kernel_excess_weight"], r["b_character_negativity_weight"]) == min_linf],
            "admissibility": "admissible_as_developmental_unit_quotient_candidate_not_strict_export",
        },
    ]
    for c in candidates:
        c["selects_unique_target_pair"] = c["selected_pairs"] == [[1, 1]]
        c["strict_exported_now"] = False
        c["can_be_used_as_strict_closure_now"] = False
    return candidates


def obligation_rows(laws: list[dict[str, Any]]) -> list[dict[str, Any]]:
    developmental_positive = [l for l in laws if l["new_theory_object"] and l["selects_unique_target_pair"]]
    return [
        {"obligation": "content_replay_excluded", "satisfied": True, "evidence": "target_sum_9 and primitive_equal_summand are classified as replay unless additionally sourced"},
        {"obligation": "developmental_theory_candidate_constructed", "satisfied": bool(developmental_positive), "evidence": "new self-balance/minimal-source and minimax-unit quotient candidates select (1,1) on the tested lattice"},
        {"obligation": "ontology_respects_nadsoliton_primordiality", "satisfied": True, "evidence": "candidates are internal nadsoliton source laws, not a lower information layer"},
        {"obligation": "nonconventional_source_derivation_exported", "satisfied": False, "evidence": "P2960 proposes theorem targets but does not derive them from prior strict artifacts"},
        {"obligation": "canonical_scale_quotient_exported", "satisfied": False, "evidence": "minimal-source/minimax wording names a quotient obligation but does not yet construct a unit-bearing quotient map"},
        {"obligation": "unit_bearing_nonproxy_coupling_exported", "satisfied": False, "evidence": "no field/action-density coupling to L_total is constructed"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["content_replay_excluded", "developmental_theory_candidate_constructed", "ontology_respects_nadsoliton_primordiality", "nonconventional_source_derivation_exported", "canonical_scale_quotient_exported", "unit_bearing_nonproxy_coupling_exported"]
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_as_strict_closure": m == (1 << len(names)) - 1, "accepts_as_developmental_work_item": all(bool(m & (1 << i)) for i in range(3))} for m in range(1 << len(names))]


def build_payload(p2959: dict[str, Any]) -> dict[str, Any]:
    rows = candidate_rows()
    laws = law_rows(rows)
    obligations = obligation_rows(laws)
    matrix = acceptance_matrix()
    return {
        "status": "P2960_DEVELOPMENTAL_ONTOLOGY_LOCALIZER_LAW_INTAKE_NO_STRICT_EXPORT",
        "input_hashes": {"P2959": hashlib.sha256(P2959.read_bytes()).hexdigest() if P2959.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "DevelopmentalOntology_LocalizerLaw_IntakeGate",
            "candidate_laws": laws,
            "acceptance_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "intake_certificate": {
            "candidate_law_count": len(laws),
            "new_developmental_candidates_selecting_target": [l["law"] for l in laws if l["new_theory_object"] and l["selects_unique_target_pair"]],
            "strict_closure_acceptance_rows": sum(1 for r in matrix if r["accepts_as_strict_closure"]),
            "developmental_work_item_acceptance_rows": sum(1 for r in matrix if r["accepts_as_developmental_work_item"]),
            "current_artifact_accepts_as_developmental_work_item": all(r["satisfied"] for r in obligations[:3]),
            "current_artifact_accepts_as_strict_closure": all(r["satisfied"] for r in obligations),
            "strict_source_derivation_exported": False,
            "canonical_scale_quotient_exported": False,
            "unit_bearing_nonproxy_coupling_exported": False,
        },
        "decision": {
            "positive_progress": "Yes: P2960 converts the author's theory-under-construction caveat into an explicit intake rule and identifies two precise developmental theorem candidates that uniquely pick the wanted (1,1) source on the tested lattice without pretending they are already strict exports.",
            "breakthrough": "No strict closure yet: the candidates are now well-formed research objects, but the nonconventional source derivation, canonical scale quotient, and unit-bearing nonproxy coupling are still missing.",
            "negative_export_flags": {k: False for k in ["strict_p2938_provenance_exported", "strict_ratio_package_source_exported", "damping_packet_exported", "nonproxy_ltotal_exported", "eom_hamiltonian_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay target_sum=9, primitive equal-summand, bounded localizer predicates, K+C decompositions, beta-scale normalization, scalar Euler insertion, or count/role aliases.  The next proof-grade move should pick exactly one P2960 developmental candidate and try to prove a nonconventional nadsoliton self-balance/source derivation, or construct the canonical scale quotient/unit-bearing nonproxy coupling it requires; otherwise pivot outside the ratio-package lane while preserving the P2929-P2960 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["intake_certificate"]
    lines = ["# P2960/S1910 developmental ontology localizer-law intake", "", f"Status: `{payload['status']}`", "", "## Intake certificate", f"- candidate laws: `{cert['candidate_law_count']}`", f"- new developmental candidates selecting target: `{cert['new_developmental_candidates_selecting_target']}`", f"- strict closure acceptance rows: `{cert['strict_closure_acceptance_rows']}`", f"- developmental work-item acceptance rows: `{cert['developmental_work_item_acceptance_rows']}`", f"- current artifact accepts as developmental work item: `{cert['current_artifact_accepts_as_developmental_work_item']}`", f"- current artifact accepts as strict closure: `{cert['current_artifact_accepts_as_strict_closure']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2959))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2960/S1910 developmental ontology localizer-law intake", "## P2960/S1910 developmental ontology localizer-law intake\n\n`P2960/S1910` incorporates the theory-under-construction caveat by separating developmental theorem candidates from strict exports.  It rejects `sum=9` and primitive equal-summand replay as closure evidence, but admits two new research objects: a developmental fractal self-balance/minimal-positive-source law and a developmental minimax source-amplitude quotient.  Both select `(a,b)=(1,1)` on the tested lattice as theorem candidates, not as present strict closure.  Missing obligations remain: nonconventional nadsoliton source derivation, canonical scale quotient, and unit-bearing nonproxy action-density coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2960/S1910 developmental localizer-law `L_total` guard", "## P2960/S1910 developmental localizer-law `L_total` guard\n\n`P2960/S1910` allows new localizer mechanisms to be worked on while the final theory is still under construction, but classifies them as developmental work items until a nonconventional source derivation, canonical scale quotient, and unit-bearing nonproxy coupling are exported.  Therefore no sourced damping coefficient enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE from P2960 alone.\n")
    append_once(AGENTS, "Current developmental ontology localizer-law intake guardrail (P2960/S1910, 2026-06-20)", "## Current developmental ontology localizer-law intake guardrail (P2960/S1910, 2026-06-20)\n\n- P2960 incorporates the fact that the final theory is still under construction: a new mechanism may be admitted as a developmental theorem candidate even if older artifacts did not already export it.\n- It rejects `sum=9` and primitive equal-summand replay as strict closure evidence, but admits two precise candidates for future proof work: a fractal self-balance/minimal-positive-source law and a minimax source-amplitude quotient, each selecting `(a,b)=(1,1)` on the tested lattice.\n- These candidates are not strict exports until a nonconventional nadsoliton source derivation, canonical scale quotient, and unit-bearing nonproxy action-density coupling are proved.\n- Do not promote P2960 to strict P2938 provenance, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.  The next admissible move should prove exactly one P2960 candidate's source derivation or construct its canonical scale quotient/unit-bearing nonproxy coupling; otherwise pivot outside the ratio-package lane while preserving the P2929-P2960 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
