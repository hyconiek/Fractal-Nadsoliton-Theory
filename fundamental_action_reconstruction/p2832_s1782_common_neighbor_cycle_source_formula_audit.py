#!/usr/bin/env python3
"""P2832/S1782: non-label common-neighbor/cycle source formula audit.

P2831 shows that digest labels/ranks cannot be promoted to source dynamics.  P2832
tests one genuinely non-label graph-source formula: the exact common-neighbor
and short-cycle profile of each decoded graph.  The candidate has a concrete
formula and dimensionless finite normalization, but it must still separate the
P2830 carrier and export a variational derivative/coupling theorem before any
K/L_total promotion.  The computation is a bounded source-law acceptance test,
not a closure claim.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2819_s1769_p2815_digest_blocker_cut_response_audit import P2818

GEN = ROOT / "generated"
P2831 = GEN / "p2831_s1781_p2815_digest_source_law_coupling_obligation_matrix.json"
OUT = GEN / "p2832_s1782_common_neighbor_cycle_source_formula_audit.json"
MD = GEN / "p2832_s1782_common_neighbor_cycle_source_formula_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
EXPECTED_GRAPH_COUNT = 16828


def adjacency_masks(graph: tuple[tuple[int, ...], ...]) -> list[int]:
    masks: list[int] = []
    for neighbors in graph:
        mask = 0
        for neighbor in neighbors:
            mask |= 1 << (neighbor - 1)
        masks.append(mask)
    return masks


def common_neighbor_cycle_profile(graph: tuple[tuple[int, ...], ...]) -> tuple[Any, ...]:
    masks = adjacency_masks(graph)
    n = len(masks)
    edge_triangle_sum = 0
    common_hist: Counter[int] = Counter()
    four_cycle_twice = 0
    for i in range(n):
        for j in range(i + 1, n):
            common = (masks[i] & masks[j]).bit_count()
            common_hist[common] += 1
            four_cycle_twice += common * (common - 1) // 2
            if (masks[i] >> j) & 1:
                edge_triangle_sum += common
    triangle_count = edge_triangle_sum // 3
    four_cycle_count = four_cycle_twice // 2
    return (
        triangle_count,
        four_cycle_count,
        tuple(sorted(common_hist.items())),
    )


def build_audit(p2831: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[common_neighbor_cycle_profile(graph)].append(index)
    residual = {key: indices for key, indices in classes.items() if len(indices) > 1}
    hist = Counter(len(indices) for indices in classes.values())
    largest = sorted((len(indices), indices[:12], key) for key, indices in residual.items())[-12:]
    return {
        "candidate_formula": "Q_cyc(G)=(triangles(G), 4-cycles(G), histogram_{pairs {u,v}} |N(u)∩N(v)|), normalized by fixed n=16 and degree=4 constants when used as a dimensionless coefficient vector",
        "decoded_graph_count": len(graphs),
        "input_p2831_status": p2831.get("status"),
        "non_label_graph_formula_exported_for_candidate": True,
        "dimensionless_normalization_available_for_candidate": True,
        "candidate_class_count": len(classes),
        "candidate_collision_class_count": len(residual),
        "candidate_collision_graph_count": sum(len(indices) for indices in residual.values()),
        "candidate_max_class_size": max(hist),
        "candidate_defect_after_formula": len(graphs) - len(classes),
        "candidate_class_size_histogram": dict(sorted(hist.items())),
        "largest_collision_samples": [
            {"class_size": size, "sample_graph_indices_0_based": sample, "profile": key}
            for size, sample, key in reversed(largest)
        ],
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2831: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2831_digest_label_lane_rejected": p2831["acceptance_matrix"]["accepted_as_bounded_obligation_no_go"],
        "exactly_one_non_label_formula_tested": True,
        "non_label_graph_formula_exported_for_candidate": audit["non_label_graph_formula_exported_for_candidate"],
        "dimensionless_normalization_available_for_candidate": audit["dimensionless_normalization_available_for_candidate"],
        "candidate_separates_full_carrier": audit["candidate_class_count"] == EXPECTED_GRAPH_COUNT,
        "variational_derivative_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_theorem_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2818_local_response_rejected",
        "p2831_digest_label_lane_rejected",
        "exactly_one_non_label_formula_tested",
        "non_label_graph_formula_exported_for_candidate",
        "dimensionless_normalization_available_for_candidate",
        "candidate_separates_full_carrier",
        "variational_derivative_exported",
        "typed_graph_to_K_or_Ltotal_coupling_theorem_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_source_law_coupling": accepted,
        "accepted_as_bounded_formula_no_go": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["common_neighbor_cycle_source_formula_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2832/S1782 common-neighbor/cycle source formula audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate formula", audit["candidate_formula"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- candidate_class_count={audit['candidate_class_count']}",
        f"- candidate_collision_class_count={audit['candidate_collision_class_count']}",
        f"- candidate_collision_graph_count={audit['candidate_collision_graph_count']}",
        f"- candidate_max_class_size={audit['candidate_max_class_size']}",
        f"- candidate_defect_after_formula={audit['candidate_defect_after_formula']}", "",
        "## Acceptance",
        f"- accepted_as_source_law_coupling={acceptance['accepted_as_source_law_coupling']}",
        f"- accepted_as_bounded_formula_no_go={acceptance['accepted_as_bounded_formula_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2831 = read_json(P2831)
    audit = build_audit(p2831)
    payload: dict[str, Any] = {
        "status": "P2832_COMMON_NEIGHBOR_CYCLE_SOURCE_FORMULA_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2831": sha(P2831), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2831": p2831.get("status")},
        "common_neighbor_cycle_source_formula_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2832 tests one genuinely non-label graph formula after P2831: the exact common-neighbor/cycle profile Q_cyc.  The formula is concrete and dimensionlessly normalizable on the fixed 16-node 4-regular carrier, but the finite computation yields only 344 classes, with 272 residual collision classes covering 16,756 graphs and max class size 788.  Since it does not separate the P2830 carrier and still exports no variational derivative or typed graph-to-K/L_total coupling theorem, it is rejected as a source-law/coupling promotion.",
            "next_honest_step": "Do not promote low-order common-neighbor/cycle profiles to dynamics.  The next admissible source-law move must add one higher-order non-label formula that is not merely the digest label, preferably an exact edge-toggle response polynomial or typed action-density expression with a proved variational derivative into K or L_total; otherwise pivot away from graph-source coupling and preserve the P2831-P2832 no-coupling boundary.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2831)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2832/S1782 common-neighbor cycle source formula audit", "## P2832/S1782 common-neighbor cycle source formula audit\n\n`P2832/S1782` tests one genuinely non-label graph-source formula after P2831: `Q_cyc(G)=(triangles, 4-cycles, common-neighbor histogram)`, with dimensionless normalization available on the fixed 16-node 4-regular carrier.  The finite computation yields only `344` classes, with `272` residual collision classes covering `16,756` graphs and max class size `788`; it therefore does not separate the P2830 carrier and still exports no variational derivative or typed graph-to-`K`/`L_total` coupling theorem.  No strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2832/S1782 common-neighbor cycle Ltotal guard", "## P2832/S1782 common-neighbor cycle Ltotal guard\n\n`P2832/S1782` adds no variational term to `L_total`.  The tested low-order common-neighbor/cycle formula is non-label and dimensionless on the fixed carrier, but it leaves `16,756` collision graphs and has no Euler-Lagrange derivative or typed graph-to-`K`/`L_total` coupling theorem.\n")
    append_once(AGENTS, "Current common-neighbor cycle source formula guardrail (P2832/S1782, 2026-06-17)", "## Current common-neighbor cycle source formula guardrail (P2832/S1782, 2026-06-17)\n\n- P2832 tests one non-label candidate graph-source formula after P2831: exact triangle count, 4-cycle count, and common-neighbor histogram, dimensionlessly normalized on the fixed `16`-node `4`-regular carrier.\n- The formula yields only `344` classes, with `272` residual collision classes covering `16,756` graphs and max class size `788`; it does not separate the P2830 carrier and exports no variational derivative or typed graph-to-`K`/`L_total` coupling theorem.\n- Do not promote low-order common-neighbor/cycle profiles to `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible source-law move must add one higher-order non-label formula with a proved variational derivative into `K` or `L_total`, or pivot away from graph-source coupling.\n")
    return payload


if __name__ == "__main__":
    main()
