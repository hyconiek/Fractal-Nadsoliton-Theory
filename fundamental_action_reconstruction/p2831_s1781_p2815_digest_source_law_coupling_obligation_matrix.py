#!/usr/bin/env python3
"""P2831/S1781: source-law/coupling obligation matrix after P2830.

P2830 gives a full finite 16,828-graph P2815 digest separation witness.  P2831
performs the next narrow theorem-audit step: test whether the currently exported
carrier/digest labels themselves satisfy the obligations for promotion to a
strict graph-source law coupled to K or L_total.  The matrix is deliberately
finite and stop-on-first-missing-premise: separation is verified from the P2830
manifest, then each available carrier-derived candidate is rejected unless it
has a non-label formula, units/normalization, a variational derivative, and a
typed coupling theorem.
"""
from __future__ import annotations

import json
from collections import Counter
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2819_s1769_p2815_digest_blocker_cut_response_audit import P2818

GEN = ROOT / "generated"
P2830 = GEN / "p2830_s1780_p2815_digest_full_carrier_cached_audit.json"
P2830_MANIFEST = GEN / "p2830_s1780_full_carrier_cached_digest_manifest.json"
OUT = GEN / "p2831_s1781_p2815_digest_source_law_coupling_obligation_matrix.json"
MD = GEN / "p2831_s1781_p2815_digest_source_law_coupling_obligation_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
EXPECTED_GRAPH_COUNT = 16828


def digest_rows(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    return manifest["full_carrier_digest_sha256_by_graph_index"]


def candidate_obligation_row(name: str, description: str, values: list[Any], *, non_label_formula: bool) -> dict[str, Any]:
    counts = Counter(values)
    residual = {str(key): count for key, count in counts.items() if count > 1}
    unique_count = len(counts)
    facts = {
        "deterministic_on_p2830_manifest": len(values) == EXPECTED_GRAPH_COUNT,
        "separates_full_carrier": unique_count == EXPECTED_GRAPH_COUNT,
        "non_label_graph_formula_exported": non_label_formula,
        "target_independent_units_or_normalization_exported": False,
        "variational_derivative_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_theorem_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "deterministic_on_p2830_manifest",
        "separates_full_carrier",
        "non_label_graph_formula_exported",
        "target_independent_units_or_normalization_exported",
        "variational_derivative_exported",
        "typed_graph_to_K_or_Ltotal_coupling_theorem_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "candidate": name,
        "description": description,
        "value_count": len(values),
        "unique_value_count": unique_count,
        "collision_class_count": len(residual),
        "collision_graph_count": sum(residual.values()),
        "max_class_size": max(counts.values()),
        "sample_residual_collisions": dict(list(residual.items())[:5]),
        "facts": facts,
        "accepted_as_source_law_coupling": accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def build_audit(p2830: dict[str, Any], manifest: dict[str, Any]) -> dict[str, Any]:
    rows = sorted(digest_rows(manifest), key=lambda row: row["graph_index_0_based"])
    graph_indices = [row["graph_index_0_based"] for row in rows]
    digest_values = [row["p2815_digest_response_sha256"] for row in rows]
    lex_rank = {digest: rank for rank, digest in enumerate(sorted(digest_values))}
    digest_int_mod_4096 = [int(digest[:12], 16) % 4096 for digest in digest_values]
    digest_prefix_16 = [digest[:16] for digest in digest_values]
    candidate_rows = [
        candidate_obligation_row(
            "Q_digest_sha256_label",
            "Use the P2830 P2815 digest SHA-256 label itself as a scalar/source label.",
            digest_values,
            non_label_formula=False,
        ),
        candidate_obligation_row(
            "Q_digest_lex_rank",
            "Use the lexicographic rank of the P2830 digest labels as a normalized ordinal coefficient.",
            [lex_rank[digest] for digest in digest_values],
            non_label_formula=False,
        ),
        candidate_obligation_row(
            "Q_digest_prefix_16",
            "Use the first 16 hex characters of the digest label as a shortened source label.",
            digest_prefix_16,
            non_label_formula=False,
        ),
        candidate_obligation_row(
            "Q_digest_int_mod_4096",
            "Use a finite integer projection of the digest label modulo 4096 as a compact coefficient.",
            digest_int_mod_4096,
            non_label_formula=False,
        ),
    ]
    return {
        "audited_question": "Can P2830 full-carrier digest separation alone be promoted to a strict graph-source law/coupling theorem?",
        "decoded_graph_count": p2830["p2815_digest_full_carrier_cached_audit"]["decoded_graph_count"],
        "p2830_full_carrier_rows": len(rows),
        "p2830_full_carrier_digest_classes": p2830["p2815_digest_full_carrier_cached_audit"]["full_carrier_digest_refined_class_count"],
        "p2830_full_carrier_digest_collisions": p2830["p2815_digest_full_carrier_cached_audit"]["full_carrier_digest_collision_class_count"],
        "graph_index_min": min(graph_indices),
        "graph_index_max": max(graph_indices),
        "graph_index_coverage_is_exact": graph_indices == list(range(EXPECTED_GRAPH_COUNT)),
        "candidate_obligation_rows": candidate_rows,
        "accepted_candidate_count": sum(1 for row in candidate_rows if row["accepted_as_source_law_coupling"]),
        "rejected_candidate_count": sum(1 for row in candidate_rows if not row["accepted_as_source_law_coupling"]),
    }


def acceptance_matrix(audit: dict[str, Any], p2818: dict[str, Any], p2830: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "p2818_local_response_rejected": p2818["acceptance_matrix"]["accepted_as_bounded_candidate_rejection"],
        "p2830_full_carrier_separation_witness": p2830["acceptance_matrix"]["accepted_as_full_16828_carrier_digest_separation_witness"],
        "graph_index_coverage_exact": audit["graph_index_coverage_is_exact"],
        "obligation_matrix_executed": len(audit["candidate_obligation_rows"]) == 4,
        "some_candidate_accepted": audit["accepted_candidate_count"] > 0,
        "strict_graph_source_law_exported": False,
        "typed_graph_to_K_or_Ltotal_coupling_exported": False,
        "units_and_variational_derivative_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted = all(facts[key] for key in [
        "p2818_local_response_rejected",
        "p2830_full_carrier_separation_witness",
        "graph_index_coverage_exact",
        "obligation_matrix_executed",
        "some_candidate_accepted",
        "strict_graph_source_law_exported",
        "typed_graph_to_K_or_Ltotal_coupling_exported",
        "units_and_variational_derivative_exported",
    ]) and not facts["selector_bridge_or_role_transfer_imported"]
    return {
        "facts": facts,
        "accepted_as_source_law_coupling_theorem": accepted,
        "accepted_as_bounded_obligation_no_go": not accepted,
        "missing_for_promotion": [key for key, value in facts.items() if not value and key != "selector_bridge_or_role_transfer_imported"],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["p2815_digest_source_law_coupling_obligation_matrix"]
    acceptance = payload["acceptance_matrix"]
    lines = [
        "# P2831/S1781 P2815 digest source-law/coupling obligation matrix", "", f"Status: `{payload['status']}`", "",
        "## Audited question", audit["audited_question"], "", "## Finite counts",
        f"- decoded_graph_count={audit['decoded_graph_count']}",
        f"- p2830_full_carrier_rows={audit['p2830_full_carrier_rows']}",
        f"- p2830_full_carrier_digest_classes={audit['p2830_full_carrier_digest_classes']}",
        f"- p2830_full_carrier_digest_collisions={audit['p2830_full_carrier_digest_collisions']}",
        f"- graph_index_coverage_is_exact={audit['graph_index_coverage_is_exact']}",
        f"- accepted_candidate_count={audit['accepted_candidate_count']}",
        f"- rejected_candidate_count={audit['rejected_candidate_count']}", "",
        "## Candidate rows",
    ]
    for row in audit["candidate_obligation_rows"]:
        lines.extend([
            f"- {row['candidate']}: unique={row['unique_value_count']}, collisions={row['collision_class_count']}, accepted={row['accepted_as_source_law_coupling']}, missing={row['missing_for_promotion']}",
        ])
    lines.extend([
        "", "## Acceptance",
        f"- accepted_as_source_law_coupling_theorem={acceptance['accepted_as_source_law_coupling_theorem']}",
        f"- accepted_as_bounded_obligation_no_go={acceptance['accepted_as_bounded_obligation_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2818 = read_json(P2818)
    p2830 = read_json(P2830)
    manifest = read_json(P2830_MANIFEST)
    audit = build_audit(p2830, manifest)
    payload: dict[str, Any] = {
        "status": "P2831_P2815_DIGEST_SOURCE_LAW_COUPLING_OBLIGATION_NO_GO_NO_CLOSURE",
        "input_hashes": {"P2818": sha(P2818), "P2830": sha(P2830), "P2830_MANIFEST": sha(P2830_MANIFEST), "16_4_4.scd": sha(SCD)},
        "input_statuses": {"P2818": p2818.get("status"), "P2830": p2830.get("status")},
        "p2815_digest_source_law_coupling_obligation_matrix": audit,
        "decision": {
            "negative_export_flags": {
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "units_and_variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2831 audits the exact theorem gap left by P2830: whether full finite P2815 digest separation by itself supplies a strict graph-source law/coupling.  Four carrier-derived candidates are tested.  Full digest labels, lexicographic ranks, and digest prefixes are injective but remain labels without non-label graph formulas, target-independent units, variational derivatives, or typed graph-to-K/L_total coupling theorems; the compact modulo projection is not even separating.  Therefore the source-law/coupling promotion is rejected without importing selector, bridge, or role-transfer premises.",
            "next_honest_step": "Do not replay digest labels, hashes, or ordinal ranks as dynamics.  The next admissible move requires one genuinely new non-label graph-source formula with target-independent units/normalization and an explicit variational derivative into K or L_total; if no such formula is supplied, pivot to a different new typed object or preserve the P2830-P2831 no-coupling boundary.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit, p2818, p2830)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2831/S1781 P2815 digest source-law coupling obligation matrix", "## P2831/S1781 P2815 digest source-law coupling obligation matrix\n\n`P2831/S1781` audits the theorem gap left by P2830: whether full finite P2815 digest separation alone supplies a strict graph-source law/coupling.  Four carrier-derived candidates are tested: full digest labels, lexicographic digest ranks, digest prefixes, and a modulo integer projection.  The label/rank/prefix candidates separate but remain label functionals with no non-label graph formula, target-independent units, variational derivative, or typed graph-to-`K`/`L_total` coupling theorem; the modulo projection is not separating.  Thus no strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2831/S1781 digest source-law coupling Ltotal guard", "## P2831/S1781 digest source-law coupling Ltotal guard\n\n`P2831/S1781` adds no variational term to `L_total`.  It proves that current P2830 carrier-derived digest labels/ranks/prefixes/projections do not satisfy the source-law obligations: no non-label formula with target-independent units, no Euler-Lagrange derivative, and no typed graph-to-`K`/`L_total` coupling theorem is exported.\n")
    append_once(AGENTS, "Current P2815 digest source-law coupling obligation guardrail (P2831/S1781, 2026-06-17)", "## Current P2815 digest source-law coupling obligation guardrail (P2831/S1781, 2026-06-17)\n\n- P2831 audits whether P2830 full-carrier digest separation can itself be promoted to a strict graph-source law/coupling theorem.\n- Four carrier-derived candidates are tested: digest SHA labels, lexicographic digest ranks, digest prefixes, and a modulo integer projection.  The separating candidates are still labels/ranks with no non-label graph formula, no target-independent units, no variational derivative, and no typed graph-to-`K`/`L_total` coupling theorem; the modulo projection is non-separating.\n- Do not promote digest labels, hashes, or ordinal ranks to `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must supply one genuinely new non-label graph-source formula with units/normalization and an explicit variational derivative into `K` or `L_total`, or pivot to a different new typed object.\n")
    return payload


if __name__ == "__main__":
    main()
