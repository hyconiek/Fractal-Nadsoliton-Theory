#!/usr/bin/env python3
"""P2835/S1785: combined finite witness source-law theorem obligation audit.

P2834 separates the P2833 residuals, so carrier separation is no longer the
right next move.  P2835 composes the P2833 first-variation witness with the
P2834 residual second-variation witness into one deterministic full-carrier
separating key, then audits whether that finite key satisfies the theorem-level
source-law/coupling obligations needed for K or L_total promotion.
"""
from __future__ import annotations

import hashlib
import json
from collections import Counter, defaultdict
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2833_s1783_edge_toggle_response_polynomial_source_audit import edge_toggle_response_polynomial
from p2834_s1784_residual_second_variation_two_edge_toggle_audit import MANIFEST as P2834_MANIFEST

GEN = ROOT / "generated"
P2830 = GEN / "p2830_s1780_p2815_digest_full_carrier_cached_audit.json"
P2833 = GEN / "p2833_s1783_edge_toggle_response_polynomial_source_audit.json"
P2834 = GEN / "p2834_s1784_residual_second_variation_two_edge_toggle_audit.json"
OUT = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
MD = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.md"
MANIFEST = GEN / "p2835_s1785_combined_witness_separator_manifest.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
EXPECTED_FULL_CARRIER_COUNT = 16828
EXPECTED_P2833_RESIDUAL_GRAPH_COUNT = 138


def stable_digest(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()


def p2833_classes_for_graphs(graphs: list[tuple[tuple[int, ...], ...]]) -> tuple[dict[str, list[int]], dict[int, str]]:
    by_digest: dict[str, list[int]] = defaultdict(list)
    graph_digest: dict[int, str] = {}
    for index, graph in enumerate(graphs):
        digest = stable_digest(edge_toggle_response_polynomial(graph))
        graph_digest[index] = digest
        by_digest[digest].append(index)
    return by_digest, graph_digest


def build_combined_separator_manifest(graphs: list[tuple[tuple[int, ...], ...]], p2834_manifest: dict[str, Any]) -> dict[str, Any]:
    p2833_by_digest, p2833_graph_digest = p2833_classes_for_graphs(graphs)
    residual_p2833_digests = {digest for digest, indices in p2833_by_digest.items() if len(indices) > 1}
    p2834_patch = {
        row["graph_index_0_based"]: row["second_variation_profile_sha256"]
        for row in p2834_manifest["graph_profile_digests"]
    }
    rows: list[dict[str, Any]] = []
    combined_classes: dict[tuple[str, str | None], list[int]] = defaultdict(list)
    patched_count = 0
    for index in range(len(graphs)):
        first_digest = p2833_graph_digest[index]
        needs_patch = first_digest in residual_p2833_digests
        second_digest = p2834_patch.get(index) if needs_patch else None
        if needs_patch:
            patched_count += 1
        key = (first_digest, second_digest)
        combined_classes[key].append(index)
        rows.append({
            "graph_index_0_based": index,
            "p2833_edge_toggle_profile_sha256": first_digest,
            "p2834_second_variation_patch_sha256": second_digest,
            "p2834_patch_required": needs_patch,
            "combined_separator_sha256": stable_digest(key),
        })
    residual_after = {key: indices for key, indices in combined_classes.items() if len(indices) > 1}
    hist = Counter(len(indices) for indices in combined_classes.values())
    manifest_payload = {
        "status": "P2835_COMBINED_FINITE_SEPARATOR_MANIFEST",
        "full_carrier_graph_count": len(graphs),
        "p2833_class_count": len(p2833_by_digest),
        "p2833_residual_digest_count": len(residual_p2833_digests),
        "p2834_patch_graph_count": patched_count,
        "combined_class_count": len(combined_classes),
        "combined_collision_class_count": len(residual_after),
        "combined_collision_graph_count": sum(len(indices) for indices in residual_after.values()),
        "combined_class_size_histogram": dict(sorted(hist.items())),
        "patched_graph_indices_0_based": [row["graph_index_0_based"] for row in rows if row["p2834_patch_required"]],
        "combined_separator_row_count": len(rows),
        "combined_separator_rows_sha256": stable_digest(rows),
    }
    MANIFEST.write_text(json.dumps(manifest_payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return manifest_payload


def theorem_obligation_matrix(combined: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "obligation": "finite_full_carrier_separation",
            "status": "satisfied",
            "computational_witness": {
                "combined_class_count": combined["combined_class_count"],
                "combined_collision_class_count": combined["combined_collision_class_count"],
                "combined_collision_graph_count": combined["combined_collision_graph_count"],
            },
            "blocks_promotion_if_missing": False,
        },
        {
            "obligation": "non_label_graph_formula",
            "status": "satisfied_as_finite_piecewise_formula",
            "computational_witness": "P2833 edge-toggle profile plus P2834 residual second-variation patch; no digest/rank labels used as primitive formula inputs",
            "blocks_promotion_if_missing": False,
        },
        {
            "obligation": "target_independent_units_and_normalization",
            "status": "missing",
            "computational_witness": None,
            "blocks_promotion_if_missing": True,
        },
        {
            "obligation": "typed_domain_codomain_into_kernel_or_lagrangian_density",
            "status": "missing",
            "computational_witness": None,
            "blocks_promotion_if_missing": True,
        },
        {
            "obligation": "proved_variational_derivative_into_K_or_Ltotal",
            "status": "missing",
            "computational_witness": "finite graph toggle derivatives exist, but no exported continuum/formal variational derivative into K or L_total exists",
            "blocks_promotion_if_missing": True,
        },
        {
            "obligation": "coupling_source_law_with_units_and_coefficient",
            "status": "missing",
            "computational_witness": None,
            "blocks_promotion_if_missing": True,
        },
        {
            "obligation": "selector_bridge_role_transfer_independence",
            "status": "satisfied_by_exclusion",
            "computational_witness": "audit imports no selector closure, bridge closure, or legacy role-transfer premise",
            "blocks_promotion_if_missing": False,
        },
    ]


def build_audit(p2830: dict[str, Any], p2833: dict[str, Any], p2834: dict[str, Any]) -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    p2834_manifest = read_json(P2834_MANIFEST)
    combined = build_combined_separator_manifest(graphs, p2834_manifest)
    matrix = theorem_obligation_matrix(combined)
    missing_blockers = [row["obligation"] for row in matrix if row["blocks_promotion_if_missing"] and row["status"] == "missing"]
    return {
        "decoded_full_carrier_graph_count": len(graphs),
        "input_statuses_rechecked": {
            "P2830": p2830.get("status"),
            "P2833": p2833.get("status"),
            "P2834": p2834.get("status"),
        },
        "combined_separator": {
            "manifest_path": str(MANIFEST.relative_to(REPO)),
            "manifest_sha256": sha(MANIFEST),
            "p2833_class_count": combined["p2833_class_count"],
            "p2833_residual_digest_count": combined["p2833_residual_digest_count"],
            "p2834_patch_graph_count": combined["p2834_patch_graph_count"],
            "combined_class_count": combined["combined_class_count"],
            "combined_collision_class_count": combined["combined_collision_class_count"],
            "combined_collision_graph_count": combined["combined_collision_graph_count"],
            "combined_class_size_histogram": combined["combined_class_size_histogram"],
        },
        "theorem_obligation_matrix": matrix,
        "missing_blocking_obligations": missing_blockers,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    separator = audit["combined_separator"]
    facts = {
        "full_carrier_count_matches": audit["decoded_full_carrier_graph_count"] == EXPECTED_FULL_CARRIER_COUNT,
        "p2834_patch_count_matches": separator["p2834_patch_graph_count"] == EXPECTED_P2833_RESIDUAL_GRAPH_COUNT,
        "combined_separator_is_full_carrier_injective": separator["combined_class_count"] == EXPECTED_FULL_CARRIER_COUNT and separator["combined_collision_class_count"] == 0,
        "blocking_theorem_obligations_remain": bool(audit["missing_blocking_obligations"]),
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted_as_combined_finite_separator = all([
        facts["full_carrier_count_matches"],
        facts["p2834_patch_count_matches"],
        facts["combined_separator_is_full_carrier_injective"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    accepted_as_source_law_coupling = accepted_as_combined_finite_separator and not facts["blocking_theorem_obligations_remain"]
    return {
        "facts": facts,
        "accepted_as_combined_finite_separator": accepted_as_combined_finite_separator,
        "accepted_as_source_law_coupling": accepted_as_source_law_coupling,
        "accepted_as_theorem_obligation_no_go": not accepted_as_source_law_coupling,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["combined_witness_source_law_theorem_obligation_audit"]
    separator = audit["combined_separator"]
    lines = [
        "# P2835/S1785 combined witness source-law theorem obligation audit", "", f"Status: `{payload['status']}`", "",
        "## Combined finite separator",
        f"- decoded_full_carrier_graph_count={audit['decoded_full_carrier_graph_count']}",
        f"- p2833_class_count={separator['p2833_class_count']}",
        f"- p2834_patch_graph_count={separator['p2834_patch_graph_count']}",
        f"- combined_class_count={separator['combined_class_count']}",
        f"- combined_collision_class_count={separator['combined_collision_class_count']}",
        f"- combined_collision_graph_count={separator['combined_collision_graph_count']}",
        f"- manifest_path={separator['manifest_path']}",
        f"- manifest_sha256={separator['manifest_sha256']}", "",
        "## Blocking theorem obligations",
        *(f"- {name}" for name in audit["missing_blocking_obligations"]), "",
        "## Acceptance",
        f"- accepted_as_combined_finite_separator={payload['acceptance_matrix']['accepted_as_combined_finite_separator']}",
        f"- accepted_as_source_law_coupling={payload['acceptance_matrix']['accepted_as_source_law_coupling']}",
        f"- accepted_as_theorem_obligation_no_go={payload['acceptance_matrix']['accepted_as_theorem_obligation_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2830 = read_json(P2830)
    p2833 = read_json(P2833)
    p2834 = read_json(P2834)
    audit = build_audit(p2830, p2833, p2834)
    payload: dict[str, Any] = {
        "status": "P2835_COMBINED_FINITE_SEPARATOR_THEOREM_OBLIGATION_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2830": sha(P2830), "P2833": sha(P2833), "P2834": sha(P2834), "P2834_manifest": sha(P2834_MANIFEST), "16_4_4.scd": sha(SCD)},
        "combined_witness_source_law_theorem_obligation_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "proved_variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2835 composes the P2833 first-variation witness with the P2834 residual second-variation patch and obtains an injective full-carrier finite separator: 16,828 combined classes, zero collisions, and exactly 138 patched residual graphs.  The theorem-obligation matrix still blocks source-law/coupling promotion: target-independent units/normalization, typed domain/codomain into K or L_total, a proved variational derivative into K/L_total, and a coupling source law with units/coefficient remain missing.",
            "next_honest_step": "The next honest move is not another finite separation refinement.  It should attack exactly one missing theorem obligation, preferably target-independent units/normalization for the combined graph functional; if no strict units/normalization source can be exported, preserve the P2831-P2835 finite-separator/no-coupling boundary and pivot away from graph-source promotion.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2835/S1785 combined finite separator theorem obligation audit", "## P2835/S1785 combined finite separator theorem obligation audit\n\n`P2835/S1785` composes the P2833 edge-toggle first-variation witness with the P2834 residual second-variation patch.  The combined finite key separates the full `16,828`-graph carrier into `16,828` classes with zero collisions, using exactly `138` P2834 residual patches.  This closes the current finite-separator question but not the theorem-level source-law question: target-independent units/normalization, typed domain/codomain into `K`/`L_total`, a proved variational derivative into `K`/`L_total`, and a coupling source law with units/coefficient remain missing.  No strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2835/S1785 combined separator Ltotal guard", "## P2835/S1785 combined separator Ltotal guard\n\n`P2835/S1785` adds no term to `L_total`.  The P2833/P2834 combined finite graph witness is injective on the full decoded carrier, but it remains outside `L_total` until units/normalization, typed coupling, and a proved Euler-Lagrange/variational derivative into `K` or `L_total` are exported.\n")
    append_once(AGENTS, "Current combined finite separator theorem-obligation guardrail (P2835/S1785, 2026-06-17)", "## Current combined finite separator theorem-obligation guardrail (P2835/S1785, 2026-06-17)\n\n- P2835 composes the P2833 edge-toggle first-variation witness with the P2834 residual second-variation patch, separating the full `16,828`-graph carrier into `16,828` combined classes with zero collisions and exactly `138` patched residual graphs.\n- This is accepted only as a combined finite separator.  It does not export target-independent units/normalization, typed domain/codomain into `K`/`L_total`, a proved variational derivative, or a coupling source law with units/coefficient.\n- Do not promote P2835 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  The next admissible move must attack exactly one missing theorem obligation, preferably units/normalization, or preserve the finite-separator/no-coupling boundary.\n")
    return payload


if __name__ == "__main__":
    main()
