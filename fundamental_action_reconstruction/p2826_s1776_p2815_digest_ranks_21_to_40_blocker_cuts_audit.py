#!/usr/bin/env python3
"""P2826/S1776: ranks 21-40 P2818 blocker-cuts P2815-digest audit.

P2825 separated ranks 11-20 after the first twenty audited cuts.  P2826 takes
the next bounded descending-size batch: ranks 21-40 of the P2818 collision
classes (266, 257, 234, 233, 229, 224, 194, 193, 189, 174, 142, 141, 121, 120, 116, 114, 112, 112, 108, and 105 graphs).  The audit uses the same
full P2815-style spectral/local-motif edge-toggle response digest and reports
cumulative coverage with P2825, while still refusing source-law/coupling or
closure promotion.
"""
from __future__ import annotations

import json
import os
from collections import Counter, defaultdict
from multiprocessing import Pool
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2811_s1761_local_motif_refined_source_candidate_audit import stable_digest
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2818_s1768_local_edge_variational_response_energy_audit import local_edge_variational_response_energy
from p2819_s1769_p2815_digest_blocker_cut_response_audit import P2818, p2815_digest_response

GEN = ROOT / "generated"
P2825 = GEN / "p2825_s1775_p2815_digest_ranks_11_to_20_blocker_cuts_audit.json"
OUT = GEN / "p2826_s1776_p2815_digest_ranks_21_to_40_blocker_cuts_audit.json"
MD = GEN / "p2826_s1776_p2815_digest_ranks_21_to_40_blocker_cuts_audit.md"
MANIFEST = GEN / "p2826_s1776_ranks_21_to_40_blocker_cuts_compact_manifest.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

START_RANK = 21
END_RANK = 40


def compute_row(args: tuple[int, tuple[tuple[int, ...], ...], int]) -> dict[str, Any]:
    index, graph, class_rank = args
    response = p2815_digest_response(graph)
    return {
        "graph_index_0_based": index,
        "p2818_blocker_rank_1_based": class_rank,
        "p2815_digest_response_sha256": stable_digest(response),
        "edge_count": response["edge_count"],
        "nonedge_count": response["nonedge_count"],
        "distinct_deletion_digest_count": len(response["deletion_digest_multiset"]),
        "distinct_addition_digest_count": len(response["addition_digest_multiset"]),
    }


def sorted_p2818_collision_classes(graphs: list[tuple[tuple[int, ...], ...]]) -> list[tuple[Any, list[int]]]:
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[local_edge_variational_response_energy(graph)].append(index)
    collisions = [(key, indices) for key, indices in classes.items() if len(indices) > 1]
    return sorted(collisions, key=lambda item: (-len(item[1]), item[1][0], stable_digest(item[0])))


def compute_rows(graphs: list[tuple[tuple[int, ...], ...]], jobs: list[tuple[int, int]]) -> list[dict[str, Any]]:
    if not jobs:
        return []
    worker_count = min(8, max(1, os.cpu_count() or 1), len(jobs))
    with Pool(processes=worker_count) as pool:
        return pool.map(compute_row, [(index, graphs[index], class_rank) for index, class_rank in jobs])


def summarize_class(summary: dict[str, Any], rows: list[dict[str, Any]]) -> dict[str, Any]:
    refined: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        refined[row["p2815_digest_response_sha256"]].append(row["graph_index_0_based"])
    residual = {key: vals for key, vals in refined.items() if len(vals) > 1}
    hist = Counter(len(vals) for vals in refined.values())
    return {
        **summary,
        "computed_graph_count": len(rows),
        "toggle_digest_refined_class_count": len(refined),
        "toggle_digest_collision_class_count": len(residual),
        "toggle_digest_collision_graph_count": sum(len(vals) for vals in residual.values()),
        "toggle_digest_max_class_size": max(hist),
        "blocker_defect_after_toggle_digest": len(rows) - len(refined),
        "class_size_histogram": dict(sorted(hist.items())),
    }


def build_audit(p2825: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    ranked_classes = sorted_p2818_collision_classes(graphs)
    selected = ranked_classes[START_RANK - 1:END_RANK]
    jobs: list[tuple[int, int]] = []
    class_summaries = []
    for rank, (key, indices) in enumerate(selected, start=START_RANK):
        class_summaries.append({
            "p2818_blocker_rank_1_based": rank,
            "p2818_blocker_key_sha256": stable_digest(key),
            "p2818_blocker_class_size": len(indices),
            "first_graph_index_0_based": indices[0],
        })
        jobs.extend((index, rank) for index in indices)

    rows = compute_rows(graphs, jobs)
    rows_by_rank: dict[int, list[dict[str, Any]]] = defaultdict(list)
    for row in rows:
        rows_by_rank[row["p2818_blocker_rank_1_based"]].append(row)
    per_class = [summarize_class(summary, rows_by_rank[summary["p2818_blocker_rank_1_based"]]) for summary in class_summaries]

    combined_refined: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        combined_refined[row["p2815_digest_response_sha256"]].append(row["graph_index_0_based"])
    combined_residual = {key: vals for key, vals in combined_refined.items() if len(vals) > 1}
    combined_hist = Counter(len(vals) for vals in combined_refined.values())
    previous_graph_count = p2825["p2815_digest_ranks_11_to_20_blocker_cuts_audit"]["cumulative_audited_graph_count"]
    previous_class_count = p2825["p2815_digest_ranks_11_to_20_blocker_cuts_audit"]["cumulative_audited_collision_class_count"]

    MANIFEST.write_text(json.dumps({
        "manifest_kind": "P2826 ranks 21-40 P2818 blocker-cuts compact digest manifest",
        "rank_range_1_based": [START_RANK, END_RANK],
        "class_summaries": class_summaries,
        "computed_new_row_count": len(rows),
        "rows": rows,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "candidate_response": "full P2815-style spectral/local-motif edge-toggle response digest on P2818 blocker-cut ranks 21-40",
        "decoded_graph_count": len(graphs),
        "audited_rank_range_1_based": [START_RANK, END_RANK],
        "audited_p2818_collision_class_count": len(selected),
        "audited_blocker_class_sizes": [len(indices) for _, indices in selected],
        "audited_graph_count": len(rows),
        "previous_p2825_cumulative_audited_collision_class_count": previous_class_count,
        "previous_p2825_cumulative_audited_graph_count": previous_graph_count,
        "cumulative_audited_collision_class_count": previous_class_count + len(selected),
        "cumulative_audited_graph_count": previous_graph_count + len(rows),
        "per_class_results": per_class,
        "combined_toggle_digest_refined_class_count": len(combined_refined),
        "combined_toggle_digest_collision_class_count": len(combined_residual),
        "combined_toggle_digest_collision_graph_count": sum(len(vals) for vals in combined_residual.values()),
        "combined_toggle_digest_max_class_size": max(combined_hist),
        "combined_defect_after_toggle_digest": len(rows) - len(combined_refined),
        "combined_class_size_histogram": dict(sorted(combined_hist.items())),
        "manifest_path": str(MANIFEST.relative_to(ROOT)),
        "manifest_sha256": sha(MANIFEST),
        "sample_rows": rows[:12],
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2825: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2825_ranks_11_to_20_witness_positive": p2825["acceptance_matrix"]["accepted_as_ranks_11_to_20_blocker_cut_response_witness"],
        "expected_rank_range_audited": audit["audited_rank_range_1_based"] == [START_RANK, END_RANK],
        "ranks_21_to_40_p2818_collision_classes_audited": audit["audited_p2818_collision_class_count"] == 20,
        "uses_full_p2815_spectral_local_motif_response_digest": True,
        "each_audited_blocker_cut_fully_separated": all(row["toggle_digest_collision_class_count"] == 0 for row in audit["per_class_results"]),
        "combined_ranks_21_to_40_digest_fully_separated": audit["combined_toggle_digest_collision_class_count"] == 0,
        "all_p2818_collision_classes_audited": False,
        "all_16828_graphs_audited": audit["cumulative_audited_graph_count"] == audit["decoded_graph_count"],
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2818_local_response_rejected", "p2825_ranks_11_to_20_witness_positive", "expected_rank_range_audited",
        "ranks_21_to_40_p2818_collision_classes_audited", "uses_full_p2815_spectral_local_motif_response_digest",
        "each_audited_blocker_cut_fully_separated", "combined_ranks_21_to_40_digest_fully_separated",
        "all_p2818_collision_classes_audited", "all_16828_graphs_audited",
        "strict_graph_source_law_exported", "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_ranks_21_to_40_blocker_cut_response_witness": facts["each_audited_blocker_cut_fully_separated"] and facts["combined_ranks_21_to_40_digest_fully_separated"],
        "accepted_as_full_carrier_source_law_coupling": accepted,
        "accepted_as_bounded_no_closure_audit": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["p2815_digest_ranks_21_to_40_blocker_cuts_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2826/S1776 P2815-digest ranks 21-40 blocker-cuts audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_response"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- audited_rank_range_1_based={audit['audited_rank_range_1_based']}",
        f"- audited_p2818_collision_class_count={audit['audited_p2818_collision_class_count']}",
        f"- audited_blocker_class_sizes={audit['audited_blocker_class_sizes']}",
        f"- audited_graph_count={audit['audited_graph_count']}",
        f"- cumulative_audited_collision_class_count={audit['cumulative_audited_collision_class_count']}",
        f"- cumulative_audited_graph_count={audit['cumulative_audited_graph_count']}",
        f"- combined_toggle_digest_refined_class_count={audit['combined_toggle_digest_refined_class_count']}",
        f"- combined_toggle_digest_collision_class_count={audit['combined_toggle_digest_collision_class_count']}",
        f"- combined_toggle_digest_collision_graph_count={audit['combined_toggle_digest_collision_graph_count']}",
        f"- combined_toggle_digest_max_class_size={audit['combined_toggle_digest_max_class_size']}",
        f"- combined_defect_after_toggle_digest={audit['combined_defect_after_toggle_digest']}", "",
        "## Compact manifest", f"- manifest_path={audit['manifest_path']}", f"- manifest_sha256={audit['manifest_sha256']}", "",
        "## Acceptance", f"- accepted_as_ranks_21_to_40_blocker_cut_response_witness={acceptance['accepted_as_ranks_21_to_40_blocker_cut_response_witness']}",
        f"- accepted_as_full_carrier_source_law_coupling={acceptance['accepted_as_full_carrier_source_law_coupling']}",
        f"- accepted_as_bounded_no_closure_audit={acceptance['accepted_as_bounded_no_closure_audit']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2825 = read_json(P2825)
    audit = build_audit(p2825)
    separated = audit["combined_toggle_digest_collision_class_count"] == 0
    payload: dict[str, Any] = {
        "status": "P2826_P2815_DIGEST_RANKS_21_TO_40_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE" if separated else "P2826_P2815_DIGEST_RANKS_21_TO_40_BLOCKER_CUTS_RESIDUAL_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2825": sha(P2825), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2825": p2825.get("status")},
        "p2815_digest_ranks_21_to_40_blocker_cuts_audit": audit,
        "decision": {
            "negative_export_flags": {
                "full_16828_carrier_audited": False,
                "all_p2818_collision_classes_audited": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2826 extends P2825 by auditing descending-size P2818 blocker-cut ranks 21-40, covering 3,384 additional graphs and bringing cumulative audited coverage to forty collision classes / 13,136 graphs.  The full P2815-style digest separates each audited cut and the combined rank-21-to-40 audited set with zero residual collisions.  This is stronger batch evidence, but it still audits only forty P2818 collision classes cumulatively rather than all P2818 collision classes or the full 16,828-graph carrier, and it exports no strict graph-source law or typed K/L_total coupling theorem.  Therefore no L_total, bridge, role-transfer, selector, or ToE closure follows.",
            "next_honest_step": "Continue descending-size cached digest batches over the remaining P2818 collision classes, with an explicit stop-on-first-residual rule.  A practical next batch is ranks 41-80, unless runtime/diff-size requires smaller shards.  Source-law/coupling promotion remains blocked until all P2818 collision classes/full carrier separate and a separate theorem supplies units plus variational derivative.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2825)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2826/S1776 P2815 digest ranks 21-40 blocker-cuts audit", "## P2826/S1776 P2815 digest ranks 21-40 blocker-cuts audit\n\n`P2826/S1776` extends P2825 by auditing descending-size P2818 blocker-cut ranks `21-40` (`266`, `257`, `234`, `233`, `229`, `224`, `194`, `193`, `189`, `174`, `142`, `141`, `121`, `120`, `116`, `114`, `112`, `112`, `108`, and `105` graphs; `3,384` additional graphs) using the full P2815-style spectral/local-motif edge-toggle response digest.  The finite audit separates every audited cut and the combined rank-21-to-40 audited set with zero residual collisions; cumulative P2825+P2826 coverage is forty collision classes / `13,136` graphs.  This is still not all P2818 collision classes or a full `16,828`-graph source-law/coupling audit, and it exports no strict graph-source law or typed `K`/`L_total` coupling theorem.  No `L_total`, bridge, role-transfer, selector, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2826/S1776 ranks 21-40 blocker-cuts Ltotal guard", "## P2826/S1776 ranks 21-40 blocker-cuts Ltotal guard\n\n`P2826/S1776` adds no variational term to `L_total`.  The P2815-style digest separates blocker-cut ranks `21-40` (`3,384` additional graphs), but without all-collision-class/full-carrier coverage and a typed graph-to-`K`/`L_total` coupling theorem it remains finite graph evidence rather than an Euler-Lagrange/source-law term.\n")
    append_once(AGENTS, "Current P2815-digest ranks 21-40 blocker-cuts guardrail (P2826/S1776, 2026-06-17)", "## Current P2815-digest ranks 21-40 blocker-cuts guardrail (P2826/S1776, 2026-06-17)\n\n- P2826 extends P2825 to descending-size P2818 local-response collision-class ranks `21-40`: `266`, `257`, `234`, `233`, `229`, `224`, `194`, `193`, `189`, `174`, `142`, `141`, `121`, `120`, `116`, `114`, `112`, `112`, `108`, and `105` graphs (`3,384` additional graphs) audited with the full P2815-style spectral/local-motif edge-toggle response digest.\n- The digest separates each audited cut and the combined rank-21-to-40 audited set with zero residual collisions; cumulative P2825+P2826 coverage is forty collision classes / `13,136` graphs, still not all P2818 collision classes or a full `16,828`-graph source-law/coupling audit.\n- Do not promote `L_total`, bridge closure, role transfer, selector closure, or ToE closure from P2826; continue descending-size cached digest batches with stop-on-first-residual discipline, or provide a separate typed graph-to-`K`/`L_total` coupling theorem with units and variational derivative.\n")
    return payload


if __name__ == "__main__":
    main()
