#!/usr/bin/env python3
"""P2822/S1772: top-two P2818 blocker-cuts P2815-digest response audit.

P2821 closed the complete largest P2818 local-response blocker cut under the
full P2815-style spectral/local-motif edge-toggle response digest.  P2822 takes
the next bounded computational step: audit the two largest P2818 collision
classes together, reusing the P2821 cached rows for the 788-graph class and
computing the second 776-graph class.  This is still a blocker-cut escalation,
not full-carrier/source-law/coupling closure.
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
P2821 = GEN / "p2821_s1771_p2815_digest_full_largest_blocker_cut_audit.json"
P2821_MANIFEST = GEN / "p2821_s1771_full_largest_blocker_cut_compact_manifest.json"
OUT = GEN / "p2822_s1772_p2815_digest_top_two_blocker_cuts_audit.json"
MD = GEN / "p2822_s1772_p2815_digest_top_two_blocker_cuts_audit.md"
MANIFEST = GEN / "p2822_s1772_top_two_blocker_cuts_compact_manifest.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


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


def top_p2818_collision_classes(graphs: list[tuple[tuple[int, ...], ...]], limit: int = 2) -> list[tuple[Any, list[int]]]:
    classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        classes[local_edge_variational_response_energy(graph)].append(index)
    collisions = [(key, indices) for key, indices in classes.items() if len(indices) > 1]
    return sorted(collisions, key=lambda item: (-len(item[1]), item[1][0], stable_digest(item[0])))[:limit]


def cached_p2821_rows() -> dict[int, dict[str, Any]]:
    if not P2821_MANIFEST.exists():
        return {}
    manifest = read_json(P2821_MANIFEST)
    rows: dict[int, dict[str, Any]] = {}
    for row in manifest.get("rows", []):
        if "graph_index_0_based" in row and "p2815_digest_response_sha256" in row:
            rows[row["graph_index_0_based"]] = dict(row)
    return rows


def compute_missing_rows(graphs: list[tuple[tuple[int, ...], ...]], jobs: list[tuple[int, int]]) -> list[dict[str, Any]]:
    if not jobs:
        return []
    worker_count = min(8, max(1, os.cpu_count() or 1), len(jobs))
    with Pool(processes=worker_count) as pool:
        return pool.map(compute_row, [(index, graphs[index], class_rank) for index, class_rank in jobs])


def build_audit() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    top_classes = top_p2818_collision_classes(graphs, 2)
    cached_rows = cached_p2821_rows()
    rows_by_index: dict[int, dict[str, Any]] = {}
    missing_jobs: list[tuple[int, int]] = []
    class_summaries = []
    for rank, (key, indices) in enumerate(top_classes, start=1):
        for index in indices:
            if index in cached_rows:
                row = dict(cached_rows[index])
                row["p2818_blocker_rank_1_based"] = rank
                rows_by_index[index] = row
            else:
                missing_jobs.append((index, rank))
        class_summaries.append({
            "p2818_blocker_rank_1_based": rank,
            "p2818_blocker_key_sha256": stable_digest(key),
            "p2818_blocker_class_size": len(indices),
            "first_graph_index_0_based": indices[0],
        })
    for row in compute_missing_rows(graphs, missing_jobs):
        rows_by_index[row["graph_index_0_based"]] = row

    rows = []
    per_class: list[dict[str, Any]] = []
    for summary, (_, indices) in zip(class_summaries, top_classes):
        class_rows = [rows_by_index[index] for index in indices]
        rows.extend(class_rows)
        refined: dict[str, list[int]] = defaultdict(list)
        for row in class_rows:
            refined[row["p2815_digest_response_sha256"]].append(row["graph_index_0_based"])
        residual = {key: vals for key, vals in refined.items() if len(vals) > 1}
        hist = Counter(len(vals) for vals in refined.values())
        per_class.append({
            **summary,
            "computed_graph_count": len(class_rows),
            "toggle_digest_refined_class_count": len(refined),
            "toggle_digest_collision_class_count": len(residual),
            "toggle_digest_collision_graph_count": sum(len(vals) for vals in residual.values()),
            "toggle_digest_max_class_size": max(hist),
            "blocker_defect_after_toggle_digest": len(class_rows) - len(refined),
            "class_size_histogram": dict(sorted(hist.items())),
        })

    combined_refined: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        combined_refined[row["p2815_digest_response_sha256"]].append(row["graph_index_0_based"])
    combined_residual = {key: vals for key, vals in combined_refined.items() if len(vals) > 1}
    combined_hist = Counter(len(vals) for vals in combined_refined.values())
    reused = len([row for row in rows if row["graph_index_0_based"] in cached_rows])

    MANIFEST.write_text(json.dumps({
        "manifest_kind": "P2822 top-two P2818 blocker-cuts compact digest manifest",
        "class_summaries": class_summaries,
        "reused_p2821_cached_row_count": reused,
        "computed_new_row_count": len(missing_jobs),
        "rows": rows,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "candidate_response": "full P2815-style spectral/local-motif edge-toggle response digest on the two largest P2818 blocker cuts",
        "decoded_graph_count": len(graphs),
        "audited_p2818_collision_class_count": len(top_classes),
        "audited_graph_count": len(rows),
        "audited_blocker_class_sizes": [len(indices) for _, indices in top_classes],
        "reused_p2821_cached_row_count": reused,
        "computed_new_row_count": len(missing_jobs),
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


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2821: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2821_full_largest_blocker_witness_positive": p2821["acceptance_matrix"]["accepted_as_full_largest_blocker_cut_response_witness"],
        "top_two_p2818_collision_classes_audited": audit["audited_p2818_collision_class_count"] == 2,
        "uses_full_p2815_spectral_local_motif_response_digest": True,
        "each_audited_blocker_cut_fully_separated": all(row["toggle_digest_collision_class_count"] == 0 for row in audit["per_class_results"]),
        "combined_top_two_digest_fully_separated": audit["combined_toggle_digest_collision_class_count"] == 0,
        "all_p2818_collision_classes_audited": False,
        "all_16828_graphs_audited": audit["audited_graph_count"] == audit["decoded_graph_count"],
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2818_local_response_rejected", "p2821_full_largest_blocker_witness_positive",
        "top_two_p2818_collision_classes_audited", "uses_full_p2815_spectral_local_motif_response_digest",
        "each_audited_blocker_cut_fully_separated", "combined_top_two_digest_fully_separated",
        "all_p2818_collision_classes_audited", "all_16828_graphs_audited",
        "strict_graph_source_law_exported", "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_top_two_blocker_cut_response_witness": facts["each_audited_blocker_cut_fully_separated"] and facts["combined_top_two_digest_fully_separated"],
        "accepted_as_full_carrier_source_law_coupling": accepted,
        "accepted_as_bounded_no_closure_audit": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["p2815_digest_top_two_blocker_cuts_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2822/S1772 P2815-digest top-two blocker-cuts audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_response"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- audited_p2818_collision_class_count={audit['audited_p2818_collision_class_count']}",
        f"- audited_blocker_class_sizes={audit['audited_blocker_class_sizes']}",
        f"- audited_graph_count={audit['audited_graph_count']}",
        f"- reused_p2821_cached_row_count={audit['reused_p2821_cached_row_count']}",
        f"- computed_new_row_count={audit['computed_new_row_count']}",
        f"- combined_toggle_digest_refined_class_count={audit['combined_toggle_digest_refined_class_count']}",
        f"- combined_toggle_digest_collision_class_count={audit['combined_toggle_digest_collision_class_count']}",
        f"- combined_toggle_digest_collision_graph_count={audit['combined_toggle_digest_collision_graph_count']}",
        f"- combined_toggle_digest_max_class_size={audit['combined_toggle_digest_max_class_size']}",
        f"- combined_defect_after_toggle_digest={audit['combined_defect_after_toggle_digest']}", "",
        "## Compact manifest", f"- manifest_path={audit['manifest_path']}", f"- manifest_sha256={audit['manifest_sha256']}", "",
        "## Acceptance", f"- accepted_as_top_two_blocker_cut_response_witness={acceptance['accepted_as_top_two_blocker_cut_response_witness']}",
        f"- accepted_as_full_carrier_source_law_coupling={acceptance['accepted_as_full_carrier_source_law_coupling']}",
        f"- accepted_as_bounded_no_closure_audit={acceptance['accepted_as_bounded_no_closure_audit']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2821 = read_json(P2821)
    audit = build_audit()
    separated = audit["combined_toggle_digest_collision_class_count"] == 0
    payload: dict[str, Any] = {
        "status": "P2822_P2815_DIGEST_TOP_TWO_BLOCKER_CUTS_WITNESS_NO_FULL_COUPLING_NO_CLOSURE" if separated else "P2822_P2815_DIGEST_TOP_TWO_BLOCKER_CUTS_RESIDUAL_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2821": sha(P2821), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2821": p2821.get("status")},
        "p2815_digest_top_two_blocker_cuts_audit": audit,
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
            "reason": "P2822 extends P2821 from the single largest P2818 blocker cut to the two largest P2818 blocker cuts, covering 1,564 graphs total.  The full P2815-style digest separates each audited cut and the combined top-two audited set with zero residual collisions.  This is stronger blocker-cut evidence, but it still audits only two P2818 collision classes rather than all P2818 collision classes or the full 16,828-graph carrier, and it exports no strict graph-source law or typed K/L_total coupling theorem.  Therefore no L_total, bridge, role-transfer, selector, or ToE closure follows.",
            "next_honest_step": "Continue the same cached digest audit over the remaining P2818 collision classes in descending blocker-size batches, recording compact per-class residual counts.  If any class produces a collision, stop and localize that residual; if all classes separate, only then attempt a separate source-law/coupling theorem with units and variational derivative.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2821)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2822/S1772 P2815 digest top-two blocker-cuts audit", "## P2822/S1772 P2815 digest top-two blocker-cuts audit\n\n`P2822/S1772` extends P2821 from the complete largest P2818 local-response blocker cut to the top two P2818 blocker cuts (`788` and `776` graphs, `1,564` total) using the full P2815-style spectral/local-motif edge-toggle response digest.  The finite audit separates both audited cuts and their combined audited set with zero residual collisions, but it is still not all P2818 collision classes or a full `16,828`-graph source-law/coupling audit, and it exports no strict graph-source law or typed `K`/`L_total` coupling theorem.  No `L_total`, bridge, role-transfer, selector, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2822/S1772 top-two blocker-cuts Ltotal guard", "## P2822/S1772 top-two blocker-cuts Ltotal guard\n\n`P2822/S1772` adds no variational term to `L_total`.  The P2815-style digest separates the top two P2818 blocker cuts (`1,564` graphs total), but without all-collision-class/full-carrier coverage and a typed graph-to-`K`/`L_total` coupling theorem it remains finite graph evidence rather than an Euler-Lagrange/source-law term.\n")
    append_once(AGENTS, "Current P2815-digest top-two blocker-cuts guardrail (P2822/S1772, 2026-06-17)", "## Current P2815-digest top-two blocker-cuts guardrail (P2822/S1772, 2026-06-17)\n\n- P2822 extends P2821 to the top two P2818 local-response collision classes: `788` and `776` graphs (`1,564` total) audited with the full P2815-style spectral/local-motif edge-toggle response digest.\n- The digest separates both audited cuts and the combined top-two audited set with zero residual collisions, but this is still not all P2818 collision classes or a full `16,828`-graph source-law/coupling audit.\n- Do not promote `L_total`, bridge closure, role transfer, selector closure, or ToE closure from P2822; the next admissible computational move is descending-size cached digest batches over every remaining P2818 collision class, stopping on any residual collision, or a separate typed graph-to-`K`/`L_total` coupling theorem with units and variational derivative.\n")
    return payload


if __name__ == "__main__":
    main()
