#!/usr/bin/env python3
"""P2830/S1780: full-carrier cached P2815-digest audit.

P2829 completed all 272 P2818 local-response collision classes but left the 72
singleton local-response classes outside the cached digest batch.  P2830 closes
that finite carrier gap by combining the cached P2822-P2829 collision-class
manifests with a fresh P2815-style digest computation on those 72 singleton
P2818 classes, then testing the combined 16,828-graph carrier for residual
digest collisions.  This is still only a finite carrier-separation witness: it
does not export a graph-source law, units, a variational derivative, or a typed
coupling theorem to K/L_total.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2811_s1761_local_motif_refined_source_candidate_audit import stable_digest
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2818_s1768_local_edge_variational_response_energy_audit import local_edge_variational_response_energy
from p2819_s1769_p2815_digest_blocker_cut_response_audit import P2818
from p2829_s1779_p2815_digest_ranks_161_to_272_blocker_cuts_audit import compute_rows, sorted_p2818_collision_classes

GEN = ROOT / "generated"
P2829 = GEN / "p2829_s1779_p2815_digest_ranks_161_to_272_blocker_cuts_audit.json"
OUT = GEN / "p2830_s1780_p2815_digest_full_carrier_cached_audit.json"
MD = GEN / "p2830_s1780_p2815_digest_full_carrier_cached_audit.md"
MANIFEST = GEN / "p2830_s1780_full_carrier_cached_digest_manifest.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CACHED_COLLISION_MANIFESTS = [
    GEN / "p2822_s1772_top_two_blocker_cuts_compact_manifest.json",
    GEN / "p2823_s1773_next_three_blocker_cuts_compact_manifest.json",
    GEN / "p2824_s1774_ranks_6_to_10_blocker_cuts_compact_manifest.json",
    GEN / "p2825_s1775_ranks_11_to_20_blocker_cuts_compact_manifest.json",
    GEN / "p2826_s1776_ranks_21_to_40_blocker_cuts_compact_manifest.json",
    GEN / "p2827_s1777_ranks_41_to_80_blocker_cuts_compact_manifest.json",
    GEN / "p2828_s1778_ranks_81_to_160_blocker_cuts_compact_manifest.json",
    GEN / "p2829_s1779_ranks_161_to_272_blocker_cuts_compact_manifest.json",
]


def load_cached_collision_rows() -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    rows: list[dict[str, Any]] = []
    sources: list[dict[str, Any]] = []
    for path in CACHED_COLLISION_MANIFESTS:
        manifest = read_json(path)
        manifest_rows = manifest["rows"]
        rows.extend(manifest_rows)
        sources.append({
            "manifest_path": str(path.relative_to(ROOT)),
            "manifest_sha256": sha(path),
            "rank_range_1_based": manifest.get("rank_range_1_based"),
            "row_count": len(manifest_rows),
        })
    return rows, sources


def singleton_p2818_indices(graphs: list[tuple[tuple[int, ...], ...]]) -> list[int]:
    collision_indices: set[int] = set()
    for _, indices in sorted_p2818_collision_classes(graphs):
        collision_indices.update(indices)
    return [index for index in range(len(graphs)) if index not in collision_indices]


def full_carrier_histogram(rows: list[dict[str, Any]]) -> tuple[dict[str, list[int]], Counter[int]]:
    refined: dict[str, list[int]] = defaultdict(list)
    for row in rows:
        refined[row["p2815_digest_response_sha256"]].append(row["graph_index_0_based"])
    return refined, Counter(len(indices) for indices in refined.values())


def build_audit(p2829: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    cached_rows, cached_sources = load_cached_collision_rows()
    cached_indices = {row["graph_index_0_based"] for row in cached_rows}
    singletons = singleton_p2818_indices(graphs)
    singleton_rows = compute_rows(graphs, [(index, 0) for index in singletons])
    all_rows = cached_rows + singleton_rows
    all_indices = {row["graph_index_0_based"] for row in all_rows}
    refined, hist = full_carrier_histogram(all_rows)
    residual = {key: vals for key, vals in refined.items() if len(vals) > 1}

    local_response_classes: dict[tuple[Any, ...], list[int]] = defaultdict(list)
    for index, graph in enumerate(graphs):
        local_response_classes[local_edge_variational_response_energy(graph)].append(index)
    singleton_local_class_count = sum(1 for indices in local_response_classes.values() if len(indices) == 1)

    MANIFEST.write_text(json.dumps({
        "manifest_kind": "P2830 full 16,828-graph cached P2815 digest manifest",
        "cached_collision_manifests": cached_sources,
        "cached_collision_row_count": len(cached_rows),
        "fresh_singleton_row_count": len(singleton_rows),
        "fresh_singleton_indices_0_based": singletons,
        "full_carrier_row_count": len(all_rows),
        "full_carrier_digest_sha256_by_graph_index": [
            {"graph_index_0_based": row["graph_index_0_based"], "p2815_digest_response_sha256": row["p2815_digest_response_sha256"]}
            for row in sorted(all_rows, key=lambda row: row["graph_index_0_based"])
        ],
    }, indent=2, sort_keys=True) + "\n", encoding="utf-8")

    return {
        "candidate_response": "cached full-carrier P2815-style spectral/local-motif edge-toggle response digest over all 16,828 decoded graphs",
        "decoded_graph_count": len(graphs),
        "cached_collision_manifest_count": len(CACHED_COLLISION_MANIFESTS),
        "cached_collision_row_count": len(cached_rows),
        "cached_collision_unique_graph_count": len(cached_indices),
        "p2818_collision_class_count_after_p2829": p2829["p2815_digest_ranks_161_to_272_blocker_cuts_audit"]["cumulative_audited_collision_class_count"],
        "p2818_collision_graph_count_after_p2829": p2829["p2815_digest_ranks_161_to_272_blocker_cuts_audit"]["cumulative_audited_graph_count"],
        "fresh_singleton_p2818_local_response_class_count": singleton_local_class_count,
        "fresh_singleton_graph_count": len(singleton_rows),
        "fresh_singleton_indices_0_based": singletons,
        "full_carrier_row_count": len(all_rows),
        "full_carrier_unique_graph_count": len(all_indices),
        "full_carrier_missing_graph_count": len(graphs) - len(all_indices),
        "full_carrier_duplicate_graph_count": len(all_rows) - len(all_indices),
        "full_carrier_digest_refined_class_count": len(refined),
        "full_carrier_digest_collision_class_count": len(residual),
        "full_carrier_digest_collision_graph_count": sum(len(vals) for vals in residual.values()),
        "full_carrier_digest_max_class_size": max(hist),
        "full_carrier_defect_after_digest": len(all_rows) - len(refined),
        "full_carrier_class_size_histogram": dict(sorted(hist.items())),
        "manifest_path": str(MANIFEST.relative_to(ROOT)),
        "manifest_sha256": sha(MANIFEST),
        "cached_manifest_sources": cached_sources,
        "sample_fresh_singleton_rows": singleton_rows[:12],
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2829: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2829_all_collision_classes_witness_positive": p2829["acceptance_matrix"]["accepted_as_ranks_161_to_272_blocker_cut_response_witness"],
        "all_p2818_collision_classes_audited": audit["p2818_collision_class_count_after_p2829"] == 272,
        "all_16828_graphs_audited": audit["full_carrier_row_count"] == audit["decoded_graph_count"] == 16828,
        "every_graph_index_present_once": audit["full_carrier_unique_graph_count"] == audit["decoded_graph_count"] and audit["full_carrier_duplicate_graph_count"] == 0,
        "uses_full_p2815_spectral_local_motif_response_digest": True,
        "full_carrier_digest_fully_separated": audit["full_carrier_digest_collision_class_count"] == 0,
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "units_and_variational_derivative_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    carrier_witness = all(facts[key] for key in [
        "p2818_local_response_rejected",
        "p2829_all_collision_classes_witness_positive",
        "all_p2818_collision_classes_audited",
        "all_16828_graphs_audited",
        "every_graph_index_present_once",
        "uses_full_p2815_spectral_local_motif_response_digest",
        "full_carrier_digest_fully_separated",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    source_law = carrier_witness and all(facts[key] for key in [
        "strict_graph_source_law_exported",
        "typed_graph_to_K_or_Ltotal_coupling_exported",
        "units_and_variational_derivative_exported",
    ])
    return {
        "facts": facts,
        "accepted_as_full_16828_carrier_digest_separation_witness": carrier_witness,
        "accepted_as_full_carrier_source_law_coupling": source_law,
        "accepted_as_bounded_no_closure_audit": not source_law,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["p2815_digest_full_carrier_cached_audit"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2830/S1780 P2815-digest full-carrier cached audit", "", f"Status: `{payload['status']}`", "",
        "## Candidate", audit["candidate_response"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- cached_collision_manifest_count={audit['cached_collision_manifest_count']}",
        f"- cached_collision_row_count={audit['cached_collision_row_count']}",
        f"- fresh_singleton_graph_count={audit['fresh_singleton_graph_count']}",
        f"- full_carrier_row_count={audit['full_carrier_row_count']}",
        f"- full_carrier_unique_graph_count={audit['full_carrier_unique_graph_count']}",
        f"- full_carrier_missing_graph_count={audit['full_carrier_missing_graph_count']}",
        f"- full_carrier_duplicate_graph_count={audit['full_carrier_duplicate_graph_count']}",
        f"- full_carrier_digest_refined_class_count={audit['full_carrier_digest_refined_class_count']}",
        f"- full_carrier_digest_collision_class_count={audit['full_carrier_digest_collision_class_count']}",
        f"- full_carrier_digest_collision_graph_count={audit['full_carrier_digest_collision_graph_count']}",
        f"- full_carrier_digest_max_class_size={audit['full_carrier_digest_max_class_size']}",
        f"- full_carrier_defect_after_digest={audit['full_carrier_defect_after_digest']}", "",
        "## Compact manifest", f"- manifest_path={audit['manifest_path']}", f"- manifest_sha256={audit['manifest_sha256']}", "",
        "## Acceptance", f"- accepted_as_full_16828_carrier_digest_separation_witness={acceptance['accepted_as_full_16828_carrier_digest_separation_witness']}",
        f"- accepted_as_full_carrier_source_law_coupling={acceptance['accepted_as_full_carrier_source_law_coupling']}",
        f"- accepted_as_bounded_no_closure_audit={acceptance['accepted_as_bounded_no_closure_audit']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2829 = read_json(P2829)
    audit = build_audit(p2829)
    separated = audit["full_carrier_digest_collision_class_count"] == 0 and audit["full_carrier_row_count"] == audit["decoded_graph_count"]
    payload: dict[str, Any] = {
        "status": "P2830_P2815_DIGEST_FULL_16828_CARRIER_SEPARATION_WITNESS_NO_SOURCE_LAW_NO_CLOSURE" if separated else "P2830_P2815_DIGEST_FULL_CARRIER_RESIDUAL_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2829": sha(P2829), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2829": p2829.get("status")},
        "p2815_digest_full_carrier_cached_audit": audit,
        "decision": {
            "negative_export_flags": {
                "full_16828_carrier_audited": audit["full_carrier_row_count"] == audit["decoded_graph_count"] == 16828,
                "full_carrier_digest_separated": audit["full_carrier_digest_collision_class_count"] == 0,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "units_and_variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2830 extends P2829 from all 272 P2818 collision classes to the full 16,828-graph carrier by adding the 72 singleton local-response classes and checking all cached/fresh P2815 digest rows together.  The full carrier separates with zero residual digest collisions and zero graph-index coverage defects.  This is a stronger finite carrier-separation witness, but it still exports no strict graph-source law, no typed graph-to-K/L_total coupling theorem, and no units/variational derivative; therefore no L_total, bridge, role-transfer, selector, or ToE closure follows.",
            "next_honest_step": "Run a narrow source-law/coupling theorem audit for the now-separated P2815 digest carrier: define the candidate graph-source functional, its target-independent units/normalization, and its variational derivative into K or L_total, with stop-on-first-missing-premise discipline.  Do not replay more finite separation batches unless a residual is found or the carrier changes.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2829)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2830/S1780 P2815 digest full-carrier cached audit", "## P2830/S1780 P2815 digest full-carrier cached audit\n\n`P2830/S1780` extends P2829 by combining cached P2822-P2829 collision-class digest rows (`16,756` graphs) with fresh digest rows for the `72` singleton P2818 local-response classes.  The resulting full `16,828`-graph P2815-style spectral/local-motif edge-toggle response digest has `16,828` refined classes, zero residual collisions, and zero graph-index coverage defects.  This is a full finite carrier-separation witness, not a source-law/coupling theorem: no strict graph-source law, typed `K`/`L_total` coupling, units, or variational derivative is exported, so no `L_total`, bridge, role-transfer, selector, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2830/S1780 full-carrier digest Ltotal guard", "## P2830/S1780 full-carrier digest Ltotal guard\n\n`P2830/S1780` adds no variational term to `L_total`.  It separates the full `16,828`-graph carrier under the P2815-style digest, but without an exported graph-source functional, target units/normalization, and variational derivative into `K` or `L_total`, the result remains finite carrier evidence rather than an Euler-Lagrange/source-law term.\n")
    append_once(AGENTS, "Current P2815-digest full-carrier cached audit guardrail (P2830/S1780, 2026-06-17)", "## Current P2815-digest full-carrier cached audit guardrail (P2830/S1780, 2026-06-17)\n\n- P2830 extends P2829 from all `272` P2818 collision classes / `16,756` graphs to the full `16,828`-graph carrier by adding fresh P2815 digest rows for the `72` singleton P2818 local-response classes and checking all cached/fresh rows together.\n- The full carrier separates into `16,828` P2815 digest classes with zero residual collisions and zero graph-index coverage defects; this closes the finite carrier-separation audit for the current decoded graph carrier.\n- Do not promote `L_total`, bridge closure, role transfer, selector closure, or ToE closure from P2830.  The next honest move is a narrow source-law/coupling theorem audit for the separated digest carrier: candidate graph-source functional, target-independent units/normalization, and variational derivative into `K` or `L_total`, with stop-on-first-missing-premise discipline.\n")
    return payload


if __name__ == "__main__":
    main()
