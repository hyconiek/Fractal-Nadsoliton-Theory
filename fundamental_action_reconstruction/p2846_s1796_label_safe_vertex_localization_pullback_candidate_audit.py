#!/usr/bin/env python3
"""P2846/S1796: label-safe vertex localization/pullback candidate audit.

P2845 recommends isolating one missing source premise: a strict localization /
pullback object x_G or rho_G(x).  P2846 attacks that single premise with a
finite, label-safe graph calculation on the full 16-node 4-regular carrier.  It
tests vertex-density candidates that are invariant under relabeling and asks
whether any one can be promoted to a field/spacetime pullback object.
"""
from __future__ import annotations

import json
from collections import Counter, defaultdict
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2845 = GEN / "p2845_s1795_unit_bearing_typed_source_coupling_dimensional_obstruction_audit.json"
OUT = GEN / "p2846_s1796_label_safe_vertex_localization_pullback_candidate_audit.json"
MD = GEN / "p2846_s1796_label_safe_vertex_localization_pullback_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
EXPECTED_GRAPH_COUNT = 16828

REQUIRED_LOCALIZATION_PREMISES = (
    "label_gauge_invariant",
    "normalized_finite_vertex_measure",
    "nonconstant_on_full_carrier",
    "canonical_vertex_to_field_support",
    "spacetime_pullback_formula",
    "target_independent_units_or_volume_form",
    "locality_covariance",
    "variational_chain_rule",
    "coupling_coefficient_rule",
)


def adjacency_sets(graph: tuple[tuple[int, ...], ...]) -> tuple[frozenset[int], ...]:
    return tuple(frozenset(neighbor - 1 for neighbor in neighbors) for neighbors in graph)


def triangle_counts(adj: tuple[frozenset[int], ...]) -> tuple[int, ...]:
    counts = []
    for v, neighbors in enumerate(adj):
        total = 0
        for a, b in combinations(sorted(neighbors), 2):
            if b in adj[a]:
                total += 1
        counts.append(total)
    return tuple(counts)


def four_cycle_counts(adj: tuple[frozenset[int], ...]) -> tuple[int, ...]:
    counts = []
    vertices = range(len(adj))
    for v in vertices:
        total = 0
        for opposite in vertices:
            if opposite == v or opposite in adj[v]:
                continue
            common = len(adj[v] & adj[opposite])
            total += common * (common - 1) // 2
        counts.append(total)
    return tuple(counts)


def local_motif_wl_profile(adj: tuple[frozenset[int], ...], tri: tuple[int, ...], sq: tuple[int, ...], rounds: int = 3) -> tuple[tuple[int, int], ...]:
    colors: tuple[Any, ...] = tuple((tri[v], sq[v]) for v in range(len(adj)))
    for _ in range(rounds):
        raw = tuple((colors[v], tuple(sorted(colors[n] for n in adj[v]))) for v in range(len(adj)))
        color_ids = {color: idx for idx, color in enumerate(sorted(set(raw), key=repr))}
        colors = tuple(color_ids[color] for color in raw)
    return tuple(sorted(Counter(colors).items()))

def profile_payload(profiles: list[Any]) -> dict[str, Any]:
    histogram = Counter(profiles)
    class_sizes = sorted(histogram.values(), reverse=True)
    return {
        "distinct_profile_count": len(histogram),
        "singleton_profile_count": sum(1 for size in class_sizes if size == 1),
        "max_profile_class_size": class_sizes[0] if class_sizes else 0,
        "top_profile_class_sizes": class_sizes[:10],
        "nonconstant_on_full_carrier": len(histogram) > 1,
    }


def finite_candidate_profiles() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    profiles: dict[str, list[Any]] = defaultdict(list)
    for graph in graphs:
        adj = adjacency_sets(graph)
        tri = triangle_counts(adj)
        sq = four_cycle_counts(adj)
        profiles["uniform_vertex_measure"].append(((VERTEX_COUNT, 1),))
        profiles["triangle_count_vertex_measure"].append(tuple(sorted(Counter(tri).items())))
        profiles["four_cycle_count_vertex_measure"].append(tuple(sorted(Counter(sq).items())))
        profiles["triangle_square_joint_vertex_measure"].append(tuple(sorted(Counter(zip(tri, sq)).items())))
        profiles["local_motif_wl_vertex_measure"].append(local_motif_wl_profile(adj, tri, sq))
    return {
        "decoded_graph_count": len(graphs),
        "expected_graph_count": EXPECTED_GRAPH_COUNT,
        "candidate_profile_stats": {name: profile_payload(rows) for name, rows in profiles.items()},
    }


def premise_matrix(stats: dict[str, Any]) -> dict[str, Any]:
    rows = {}
    for name, profile in stats["candidate_profile_stats"].items():
        premises = {
            "label_gauge_invariant": True,
            "normalized_finite_vertex_measure": True,
            "nonconstant_on_full_carrier": profile["nonconstant_on_full_carrier"],
            "canonical_vertex_to_field_support": False,
            "spacetime_pullback_formula": False,
            "target_independent_units_or_volume_form": False,
            "locality_covariance": False,
            "variational_chain_rule": False,
            "coupling_coefficient_rule": False,
        }
        rows[name] = {
            "profile_stats": profile,
            "premises": premises,
            "missing_premises": [key for key in REQUIRED_LOCALIZATION_PREMISES if not premises[key]],
            "accepted_as_strict_localization_pullback": all(premises.values()),
        }
    return rows


def build_audit(p2845: dict[str, Any]) -> dict[str, Any]:
    stats = finite_candidate_profiles()
    rows = premise_matrix(stats)
    return {
        "input_statuses_rechecked": {"P2845": p2845.get("status")},
        "carrier_check": {
            "decoded_graph_count": stats["decoded_graph_count"],
            "expected_graph_count": stats["expected_graph_count"],
            "coverage_ok": stats["decoded_graph_count"] == stats["expected_graph_count"],
        },
        "required_localization_premises": list(REQUIRED_LOCALIZATION_PREMISES),
        "candidate_rows": rows,
        "accepted_candidate_count": sum(1 for row in rows.values() if row["accepted_as_strict_localization_pullback"]),
        "nonconstant_candidate_count": sum(1 for row in rows.values() if row["premises"]["nonconstant_on_full_carrier"]),
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    rows = audit["candidate_rows"]
    facts = {
        "p2845_rechecked": audit["input_statuses_rechecked"]["P2845"] == "P2845_UNIT_BEARING_TYPED_SOURCE_COUPLING_DIMENSIONAL_OBSTRUCTION_NO_CLOSURE",
        "full_carrier_covered": audit["carrier_check"]["coverage_ok"],
        "at_least_one_nonconstant_label_safe_measure_found": audit["nonconstant_candidate_count"] > 0,
        "no_candidate_exports_strict_pullback": audit["accepted_candidate_count"] == 0,
        "all_candidates_lack_field_support": all(not row["premises"]["canonical_vertex_to_field_support"] for row in rows.values()),
        "all_candidates_lack_variational_chain_rule": all(not row["premises"]["variational_chain_rule"] for row in rows.values()),
        "no_closure_promoted": True,
    }
    return {
        "facts": facts,
        "accepted_as_label_safe_localization_obstruction_audit": all(facts.values()),
        "exports_strict_localization_pullback": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["label_safe_vertex_localization_pullback_audit"]
    lines = [
        "# P2846/S1796 label-safe vertex localization/pullback candidate audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Carrier check",
        f"- decoded_graph_count={audit['carrier_check']['decoded_graph_count']}",
        f"- coverage_ok={audit['carrier_check']['coverage_ok']}",
        "",
        "## Candidate rows",
    ]
    for name, row in audit["candidate_rows"].items():
        stats = row["profile_stats"]
        lines.extend([
            f"### {name}",
            f"- distinct_profile_count={stats['distinct_profile_count']}",
            f"- singleton_profile_count={stats['singleton_profile_count']}",
            f"- max_profile_class_size={stats['max_profile_class_size']}",
            f"- nonconstant_on_full_carrier={stats['nonconstant_on_full_carrier']}",
            f"- accepted_as_strict_localization_pullback={row['accepted_as_strict_localization_pullback']}",
            f"- missing_premises={row['missing_premises']}",
        ])
    lines.extend([
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2845 = read_json(P2845)
    audit = build_audit(p2845)
    payload: dict[str, Any] = {
        "status": "P2846_LABEL_SAFE_VERTEX_LOCALIZATION_PULLBACK_CANDIDATE_NO_GO_NO_CLOSURE",
        "input_hashes": {"P2845": sha(P2845), "16_4_4.scd": sha(SCD)},
        "label_safe_vertex_localization_pullback_audit": audit,
        "decision": {
            "negative_export_flags": {
                "strict_localization_pullback_exported": False,
                "field_support_exported": False,
                "target_independent_units_or_volume_form_exported": False,
                "coupling_coefficient_rule_exported": False,
                "variational_chain_rule_exported": False,
                "nonproxy_ltotal_term_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2846 finds real nonconstant label-safe vertex-density profiles on the full 16,828 graph carrier, but every candidate remains an anonymous finite vertex measure.  None exports a canonical vertex-to-field support, spacetime pullback formula, target-independent volume/unit form, locality/covariance theorem, coupling coefficient rule, or variational chain rule.  Thus the P2845 localization premise remains blocked.",
            "next_honest_step": "Do not replay anonymous vertex-profile localization.  The next admissible move must add one genuinely new bridge from anonymous vertices to field support: either a strict coordinate/support theorem for x_G/rho_G(x), or a target-independent volume form tied to one named density.  If no such object is supplied, pivot to the alternative P2845 route: one target-independent coupling coefficient/unit source or one concrete kernel bridge atom with source premise.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2846/S1796 label-safe vertex localization/pullback candidate audit", "## P2846/S1796 label-safe vertex localization/pullback candidate audit\n\n`P2846/S1796` isolates the P2845 localization/pullback premise and tests label-safe finite vertex-density candidates on the full `16,828` graph carrier: uniform, triangle-count, four-cycle-count, triangle/square joint, and local motif-WL vertex measures.  Several candidates are nonconstant and label-safe, but all remain anonymous vertex measures with no canonical vertex-to-field support, spacetime pullback formula, target-independent volume/unit form, locality/covariance theorem, coupling coefficient rule, or variational chain rule.  No strict localization/pullback object, `L_total`, EOM, Hamiltonian, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2846/S1796 vertex localization Ltotal guard", "## P2846/S1796 vertex localization Ltotal guard\n\n`P2846/S1796` adds no action term or pullback.  Nonconstant label-safe vertex profiles exist on the finite graph carrier, but anonymous vertex measures are not field support and do not provide target-independent volume units, coupling coefficients, or variational chain rules for `L_total`.\n")
    append_once(AGENTS, "Current label-safe vertex localization/pullback guardrail (P2846/S1796, 2026-06-18)", "## Current label-safe vertex localization/pullback guardrail (P2846/S1796, 2026-06-18)\n\n- P2846 isolates the P2845 localization/pullback premise and tests label-safe finite vertex-density candidates across the full `16,828` graph carrier.\n- Nonconstant anonymous vertex profiles exist, but no candidate exports canonical vertex-to-field support, spacetime pullback, target-independent volume/unit form, locality/covariance, coupling coefficient rule, or variational chain rule.\n- Do not promote anonymous vertex profiles to strict localization/pullback, unit-bearing `L_total` source/coupling, EOM, Hamiltonian, or ToE closure.\n- A next admissible move must introduce a genuine vertex-to-field support theorem `x_G`/`rho_G(x)`, a target-independent volume form tied to one named density, a coupling coefficient/unit source, or one concrete kernel bridge atom with source premise.\n")
    return payload


if __name__ == "__main__":
    main()
