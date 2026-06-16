#!/usr/bin/env python3
"""P2810/S1760: spectral-only source functional obstruction.

P2809 found no exported source law/coupling theorem.  P2810 tests the most
concrete next candidate class without inventing dynamics: any source functional
that factors only through the P2804 adjacency/complement characteristic-polynomial
pair.  P2804/P2805/P2808 already give enough finite data for a theorem-style
obstruction: the spectral-pair quotient has fewer classes than the canonical
isomorphism quotient, and P2805 proves the residual spectral collisions are not
isomorphism duplicates.
"""
from __future__ import annotations

import json
from math import comb
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import EXPECTED_GRAPH_COUNT, sha

GEN = ROOT / "generated"
P2804 = GEN / "p2804_s1754_girth4_spectral_complement_quotient_audit.json"
P2805 = GEN / "p2805_s1755_girth4_spectral_collision_isomorphism_refinement.json"
P2808 = GEN / "p2808_s1758_pynauty_canonical_digest_manifest.json"
P2809 = GEN / "p2809_s1759_canonical_quotient_source_law_coupling_audit.json"
OUT = GEN / "p2810_s1760_spectral_only_source_functional_obstruction.json"
MD = GEN / "p2810_s1760_spectral_only_source_functional_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def obstruction_matrix(p2804: dict[str, Any], p2805: dict[str, Any], p2808: dict[str, Any], p2809: dict[str, Any]) -> dict[str, Any]:
    q2804 = p2804["girth4_spectral_complement_quotient_audit"]["adjacency_complement_pair_quotient"]
    r2805 = p2805["spectral_collision_isomorphism_refinement"]
    d2808 = p2808["canonical_digest_manifest"]
    canonical_count = d2808["canonical_certificate_hash_class_count"]
    spectral_pair_count = q2804["class_count"]
    collision_graph_count = r2805["spectral_pair_collision_graph_count"]
    collision_class_count = r2805["spectral_pair_collision_class_count"]
    singleton_count = r2805["spectral_pair_singleton_class_count"]
    quotient_defect = canonical_count - spectral_pair_count
    indistinguishable_pair_count = sum(comb(row["input_size"], 2) for row in r2805["top_refined_collision_classes"])
    # P2805 stores the exact global pairwise checks separately; use it as the
    # exact count of checked non-isomorphism pairs inside all collision classes.
    exact_checked_pair_count = r2805["negative_isomorphism_rejections_inside_collisions"]
    return {
        "candidate_class": "F(G)=Phi(chi_A_G, chi_A_complement_G), i.e. any spectral source/action functional factoring only through the P2804 adjacency/complement characteristic-polynomial pair",
        "canonical_certificate_class_count_from_p2808": canonical_count,
        "spectral_pair_class_count_from_p2804": spectral_pair_count,
        "spectral_pair_singleton_class_count": singleton_count,
        "spectral_pair_collision_class_count": collision_class_count,
        "spectral_pair_collision_graph_count": collision_graph_count,
        "quotient_defect_canonical_minus_spectral_pair": quotient_defect,
        "p2805_positive_isomorphism_matches_inside_collisions": r2805["positive_isomorphism_matches_inside_collisions"],
        "p2805_negative_isomorphism_rejections_inside_collisions": exact_checked_pair_count,
        "top20_indistinguishable_pair_count_sample": indistinguishable_pair_count,
        "component_size_histogram_inside_collisions": r2805["component_size_histogram_inside_collisions"],
        "formal_obstruction": {
            "if_source_functional_factors_through_spectral_pair": "then it is constant on each P2804 spectral-pair collision class",
            "but_p2805_refinement": "finds no isomorphic duplicates inside those collision classes",
            "therefore": "spectral-pair-only functionals cannot injectively address the P2808 canonical quotient and cannot by themselves be the missing strict graph-source law/coupling",
        },
        "sample_collision_witnesses": r2805["top_refined_collision_classes"][:5],
        "p2809_no_export_input_status": p2809["status"],
    }


def acceptance_matrix(matrix: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "canonical_quotient_count_is_16828": matrix["canonical_certificate_class_count_from_p2808"] == EXPECTED_GRAPH_COUNT,
        "spectral_pair_quotient_is_smaller_than_canonical": matrix["spectral_pair_class_count_from_p2804"] < matrix["canonical_certificate_class_count_from_p2808"],
        "spectral_collision_graphs_are_nonisomorphic_under_p2805": matrix["p2805_positive_isomorphism_matches_inside_collisions"] == 0 and matrix["p2805_negative_isomorphism_rejections_inside_collisions"] > 0,
        "quotient_defect_positive": matrix["quotient_defect_canonical_minus_spectral_pair"] > 0,
        "strict_spectral_source_law_exported": False,
        "typed_coupling_to_K_or_Ltotal_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_spectral_pair_only_source_functional_obstruction": all(facts[key] for key in [
            "canonical_quotient_count_is_16828",
            "spectral_pair_quotient_is_smaller_than_canonical",
            "spectral_collision_graphs_are_nonisomorphic_under_p2805",
            "quotient_defect_positive",
        ]),
        "accepted_as_strict_source_law_or_coupling": False,
        "missing_for_promotion": [key for key, value in facts.items() if not value],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["spectral_only_obstruction_matrix"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2810/S1760 spectral-only source functional obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate class",
        matrix["candidate_class"],
        "",
        "## Finite obstruction counts",
        f"- canonical_certificate_class_count_from_p2808={matrix['canonical_certificate_class_count_from_p2808']}",
        f"- spectral_pair_class_count_from_p2804={matrix['spectral_pair_class_count_from_p2804']}",
        f"- quotient_defect_canonical_minus_spectral_pair={matrix['quotient_defect_canonical_minus_spectral_pair']}",
        f"- spectral_pair_collision_class_count={matrix['spectral_pair_collision_class_count']}",
        f"- spectral_pair_collision_graph_count={matrix['spectral_pair_collision_graph_count']}",
        f"- p2805_positive_isomorphism_matches_inside_collisions={matrix['p2805_positive_isomorphism_matches_inside_collisions']}",
        f"- p2805_negative_isomorphism_rejections_inside_collisions={matrix['p2805_negative_isomorphism_rejections_inside_collisions']}",
        "",
        "## Decision",
        f"- accepted_as_spectral_pair_only_source_functional_obstruction={acceptance['accepted_as_spectral_pair_only_source_functional_obstruction']}",
        f"- accepted_as_strict_source_law_or_coupling={acceptance['accepted_as_strict_source_law_or_coupling']}",
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
    p2804 = read_json(P2804)
    p2805 = read_json(P2805)
    p2808 = read_json(P2808)
    p2809 = read_json(P2809)
    matrix = obstruction_matrix(p2804, p2805, p2808, p2809)
    payload: dict[str, Any] = {
        "status": "P2810_SPECTRAL_PAIR_ONLY_SOURCE_FUNCTIONAL_OBSTRUCTION_NO_STRICT_SOURCE_LAW_NO_CLOSURE",
        "input_hashes": {"P2804": sha(P2804), "P2805": sha(P2805), "P2808": sha(P2808), "P2809": sha(P2809)},
        "input_statuses": {"P2804": p2804.get("status"), "P2805": p2805.get("status"), "P2808": p2808.get("status"), "P2809": p2809.get("status")},
        "spectral_only_obstruction_matrix": matrix,
        "decision": {
            "negative_export_flags": {
                "strict_spectral_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "canonical_geometry_closure_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2810 proves a bounded finite obstruction for the first explicit post-P2809 candidate class: any source/action functional depending only on the P2804 adjacency/complement characteristic-polynomial pair collapses 16,828 P2808 canonical classes to 16,211 spectral-pair classes.  P2805 verifies the residual collision members are non-isomorphic, so the collapse is not duplicate bookkeeping.  Therefore spectral-pair-only data cannot be the missing strict source law or K/L_total coupling.",
            "next_honest_step": "If continuing this lane, add exactly one strictly richer invariant/coupling ingredient beyond the adjacency/complement spectral-pair data, such as an explicit automorphism-weighted local motif functional or a typed graph-to-K/L_total variational term, and rerun the obstruction/acceptance matrix.  Do not replay spectral-pair-only functionals as source-law candidates.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(matrix)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2810/S1760 spectral-only source functional obstruction", "## P2810/S1760 spectral-only source functional obstruction\n\n`P2810/S1760` tests the first explicit post-P2809 candidate class: source/action functionals that factor only through the P2804 adjacency/complement characteristic-polynomial pair.  The finite obstruction is exact on current artifacts: P2808 has `16,828` canonical classes, while the spectral-pair quotient has only `16,211` classes; P2805 shows the residual collision members are non-isomorphic, not duplicate records.  Thus spectral-pair-only functionals cannot provide the missing strict spectral source law or `K`/`L_total` coupling.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2810/S1760 spectral-only source Ltotal guard", "## P2810/S1760 spectral-only source Ltotal guard\n\n`P2810/S1760` adds no `L_total` term.  It rules out a whole candidate family: any variational/source term depending only on the adjacency/complement spectral-pair data is too coarse for the P2808 canonical quotient because P2805-certified non-isomorphic graphs remain collapsed.  A future Lagrangian candidate must add a strictly richer typed coupling ingredient, not reuse spectral-pair data alone.\n")
    append_once(AGENTS, "Current spectral-pair-only source obstruction guardrail (P2810/S1760, 2026-06-16)", "## Current spectral-pair-only source obstruction guardrail (P2810/S1760, 2026-06-16)\n\n- P2810 proves a finite obstruction for any source/action functional factoring only through the P2804 adjacency/complement characteristic-polynomial pair: `16,828` P2808 canonical classes collapse to `16,211` spectral-pair classes, and P2805 verifies the collapsed collision members are non-isomorphic.\n- Do not replay spectral-pair-only functionals as strict source-law or `K`/`L_total` coupling candidates.  This obstruction is not selector closure, bridge closure, role transfer, role-bearing `L_total`, or ToE closure.\n- A next admissible move in this lane must add exactly one strictly richer typed ingredient beyond spectral-pair data, such as an automorphism-weighted local motif functional or an explicit graph-to-`K`/`L_total` variational coupling, and test it with a falsifiable obstruction/acceptance matrix.\n")
    return payload


if __name__ == "__main__":
    main()
