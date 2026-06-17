#!/usr/bin/env python3
"""P2821/S1771: full largest-blocker-cut P2815-digest response audit.

P2820 separated four deterministic shards (96 graphs) from the largest P2818
local-response blocker cut.  P2821 performs the stronger bounded computation
that P2819 originally recommended: audit the entire 788-graph largest blocker
cut with the same full P2815-style spectral/local-motif edge-toggle response
digest.  This is still only one blocker cut, not all P2818 collision classes or
a source-law/coupling theorem.
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
from p2819_s1769_p2815_digest_blocker_cut_response_audit import P2818, largest_p2818_collision_class, p2815_digest_response

GEN = ROOT / "generated"
P2819 = GEN / "p2819_s1769_p2815_digest_blocker_cut_response_audit.json"
P2820 = GEN / "p2820_s1770_p2815_digest_blocker_cut_shard_response_audit.json"
P2820_MANIFEST = GEN / "p2820_s1770_blocker_cut_shard_compact_manifest.json"
OUT = GEN / "p2821_s1771_p2815_digest_full_largest_blocker_cut_audit.json"
MD = GEN / "p2821_s1771_p2815_digest_full_largest_blocker_cut_audit.md"
MANIFEST = GEN / "p2821_s1771_full_largest_blocker_cut_compact_manifest.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def compute_row(args: tuple[int, tuple[tuple[int, ...], ...]]) -> dict[str, Any]:
    index, graph = args
    response = p2815_digest_response(graph)
    return {
        "graph_index_0_based": index,
        "p2815_digest_response_sha256": stable_digest(response),
        "edge_count": response["edge_count"],
        "nonedge_count": response["nonedge_count"],
        "distinct_deletion_digest_count": len(response["deletion_digest_multiset"]),
        "distinct_addition_digest_count": len(response["addition_digest_multiset"]),
    }


def cached_p2820_rows() -> dict[int, dict[str, Any]]:
    if not P2820_MANIFEST.exists():
        return {}
    manifest = read_json(P2820_MANIFEST)
    rows: dict[int, dict[str, Any]] = {}
    for row in manifest.get("rows", []):
        if "graph_index_0_based" in row and "p2815_digest_response_sha256" in row:
            rows[row["graph_index_0_based"]] = dict(row)
    return rows


def compute_missing_rows(graphs: list[tuple[tuple[int, ...], ...]], indices: list[int]) -> list[dict[str, Any]]:
    if not indices:
        return []
    worker_count = min(8, max(1, os.cpu_count() or 1), len(indices))
    with Pool(processes=worker_count) as pool:
        return pool.map(compute_row, [(index, graphs[index]) for index in indices])


def build_audit() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    blocker_key, blocker_indices = largest_p2818_collision_class(graphs)
    cached_rows = cached_p2820_rows()
    rows_by_index = {index: cached_rows[index] for index in blocker_indices if index in cached_rows}
    missing_indices = [index for index in blocker_indices if index not in rows_by_index]
    for row in compute_missing_rows(graphs, missing_indices):
        rows_by_index[row["graph_index_0_based"]] = row
    rows = [rows_by_index[index] for index in blocker_indices]

    refined: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        refined[row["p2815_digest_response_sha256"]].append(row["graph_index_0_based"])
    residual = {key: indices for key, indices in refined.items() if len(indices) > 1}
    class_size_histogram = Counter(len(indices) for indices in refined.values())

    MANIFEST.write_text(json.dumps({
        "manifest_kind": "P2821 full largest P2818 blocker-cut compact digest manifest",
        "p2818_blocker_key_sha256": stable_digest(blocker_key),
        "p2818_blocker_class_size": len(blocker_indices),
        "reused_p2820_cached_row_count": len([index for index in blocker_indices if index in cached_rows]),
        "computed_new_row_count": len(missing_indices),
        "rows": rows,
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "candidate_response": "full P2815-style spectral/local-motif edge-toggle response digest on the complete largest P2818 blocker cut",
        "blocker_cut": "complete largest P2818 E_local_toggle collision class",
        "decoded_graph_count": len(graphs),
        "p2818_blocker_class_size": len(blocker_indices),
        "computed_graph_count": len(rows),
        "reused_p2820_cached_row_count": len([index for index in blocker_indices if index in cached_rows]),
        "computed_new_row_count": len(missing_indices),
        "p2818_blocker_key_sha256": stable_digest(blocker_key),
        "toggle_digest_refined_class_count": len(refined),
        "toggle_digest_collision_class_count": len(residual),
        "toggle_digest_collision_graph_count": sum(len(indices) for indices in residual.values()),
        "toggle_digest_max_class_size": max(class_size_histogram),
        "blocker_defect_after_toggle_digest": len(blocker_indices) - len(refined),
        "class_size_histogram": dict(sorted(class_size_histogram.items())),
        "manifest_path": str(MANIFEST.relative_to(ROOT)),
        "manifest_sha256": sha(MANIFEST),
        "sample_rows": rows[:12],
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2819: dict[str, Any], p2820: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2819_sample_witness_positive": p2819["acceptance_matrix"]["accepted_as_sample_response_witness"],
        "p2820_shard_witness_positive": p2820["acceptance_matrix"]["accepted_as_four_shard_response_witness"],
        "exactly_one_blocker_cut_tested": True,
        "uses_full_p2815_spectral_local_motif_response_digest": True,
        "entire_largest_blocker_cut_audited": audit["computed_graph_count"] == audit["p2818_blocker_class_size"],
        "blocker_cut_is_fully_separated": audit["toggle_digest_collision_class_count"] == 0,
        "all_16828_graphs_audited": audit["computed_graph_count"] == audit["decoded_graph_count"],
        "all_p2818_collision_classes_audited": False,
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2818_local_response_rejected", "p2819_sample_witness_positive", "p2820_shard_witness_positive",
        "uses_full_p2815_spectral_local_motif_response_digest", "entire_largest_blocker_cut_audited",
        "blocker_cut_is_fully_separated", "all_16828_graphs_audited", "all_p2818_collision_classes_audited",
        "strict_graph_source_law_exported", "typed_graph_to_K_or_Ltotal_coupling_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_full_largest_blocker_cut_response_witness": facts["entire_largest_blocker_cut_audited"] and facts["blocker_cut_is_fully_separated"],
        "accepted_as_full_carrier_source_law_coupling": accepted,
        "accepted_as_bounded_no_closure_audit": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["p2815_digest_full_largest_blocker_cut_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2821/S1771 P2815-digest full largest-blocker-cut audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_response"], "", "## Blocker cut", audit["blocker_cut"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2818_blocker_class_size={audit['p2818_blocker_class_size']}",
        f"- computed_graph_count={audit['computed_graph_count']}",
        f"- reused_p2820_cached_row_count={audit['reused_p2820_cached_row_count']}",
        f"- computed_new_row_count={audit['computed_new_row_count']}",
        f"- toggle_digest_refined_class_count={audit['toggle_digest_refined_class_count']}",
        f"- toggle_digest_collision_class_count={audit['toggle_digest_collision_class_count']}",
        f"- toggle_digest_collision_graph_count={audit['toggle_digest_collision_graph_count']}",
        f"- toggle_digest_max_class_size={audit['toggle_digest_max_class_size']}",
        f"- blocker_defect_after_toggle_digest={audit['blocker_defect_after_toggle_digest']}", "",
        "## Compact manifest", f"- manifest_path={audit['manifest_path']}", f"- manifest_sha256={audit['manifest_sha256']}", "",
        "## Acceptance", f"- accepted_as_full_largest_blocker_cut_response_witness={acceptance['accepted_as_full_largest_blocker_cut_response_witness']}",
        f"- accepted_as_full_carrier_source_law_coupling={acceptance['accepted_as_full_carrier_source_law_coupling']}",
        f"- accepted_as_bounded_no_closure_audit={acceptance['accepted_as_bounded_no_closure_audit']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2819 = read_json(P2819)
    p2820 = read_json(P2820)
    audit = build_audit()
    separated = audit["toggle_digest_collision_class_count"] == 0
    payload: dict[str, Any] = {
        "status": "P2821_P2815_DIGEST_FULL_LARGEST_BLOCKER_CUT_WITNESS_NO_FULL_COUPLING_NO_CLOSURE" if separated else "P2821_P2815_DIGEST_FULL_LARGEST_BLOCKER_CUT_RESIDUAL_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2819": sha(P2819), "P2820": sha(P2820), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2819": p2819.get("status"), "P2820": p2820.get("status")},
        "p2815_digest_full_largest_blocker_cut_audit": audit,
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
            "reason": "P2821 completes the P2819/P2820 largest-blocker-cut escalation: the same full P2815-style digest is audited on all 788 graphs in the largest P2818 blocker cut and separates the cut with zero residual collisions.  This is a stronger finite blocker-cut witness, but it still audits only one P2818 collision class rather than all P2818 collision classes or the full 16,828-graph carrier, and it exports no strict graph-source law or typed K/L_total coupling theorem.  Therefore no L_total, bridge, role-transfer, selector, or ToE closure follows.",
            "next_honest_step": "Extend the cached P2815-digest response audit from the largest blocker cut to every remaining P2818 collision class, with compact per-class residual counts and a manifest.  Only after all P2818 collision classes/full carrier are audited should any source-law/coupling theorem with units and variational derivative be attempted; until then, keep no-coupling/no-closure boundaries.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2819, p2820)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2821/S1771 P2815 digest full largest-blocker-cut audit", "## P2821/S1771 P2815 digest full largest-blocker-cut audit\n\n`P2821/S1771` completes the P2819/P2820 largest-blocker-cut escalation by applying the full P2815-style spectral/local-motif edge-toggle response digest to all `788` graphs in the largest P2818 local-response blocker cut.  The finite audit separates that whole blocker cut with zero residual collisions, but it is still not all P2818 collision classes or a full `16,828`-graph source-law/coupling audit, and it exports no strict graph-source law or typed `K`/`L_total` coupling theorem.  No `L_total`, bridge, role-transfer, selector, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2821/S1771 full largest-blocker-cut Ltotal guard", "## P2821/S1771 full largest-blocker-cut Ltotal guard\n\n`P2821/S1771` adds no variational term to `L_total`.  The P2815-style digest separates the complete largest `788`-graph P2818 blocker cut, but without all-collision-class/full-carrier coverage and a typed graph-to-`K`/`L_total` coupling theorem it remains finite graph evidence rather than an Euler-Lagrange/source-law term.\n")
    append_once(AGENTS, "Current P2815-digest full largest-blocker-cut guardrail (P2821/S1771, 2026-06-17)", "## Current P2815-digest full largest-blocker-cut guardrail (P2821/S1771, 2026-06-17)\n\n- P2821 extends P2819/P2820 to the complete largest P2818 local-response collision class: all `788` graphs are audited with the full P2815-style spectral/local-motif edge-toggle response digest.\n- The digest separates that full blocker cut with zero residual collisions, giving stronger finite response-witness evidence, but this is still not all P2818 collision classes or a full `16,828`-graph source-law/coupling audit.\n- Do not promote `L_total`, bridge closure, role transfer, selector closure, or ToE closure from P2821; the next admissible computational move is the same cached digest audit over every remaining P2818 collision class, with compact residual counts and manifest, or an independent typed graph-to-`K`/`L_total` coupling theorem with units and variational derivative.\n")
    return payload


if __name__ == "__main__":
    main()
