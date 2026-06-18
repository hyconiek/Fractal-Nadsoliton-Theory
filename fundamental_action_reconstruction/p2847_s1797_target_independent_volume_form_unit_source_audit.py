#!/usr/bin/env python3
"""P2847/S1797: target-independent volume-form/unit-source audit.

After P2846, nonconstant label-safe anonymous vertex profiles exist, but no
field support/pullback theorem exists.  P2847 attacks the adjacent single
premise recommended by P2846: can one named finite vertex density export a
target-independent volume form or unit source for an L_total coupling?

The calculation is intentionally finite and conservative.  It reuses the full
16,828 graph carrier, builds raw and normalized mass totals for the P2846
profile families, and separates formal finite probability normalization from a
unit-bearing spacetime volume/coupling source law.
"""
from __future__ import annotations

import json
from collections import Counter
from fractions import Fraction
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, decode_scd_bytes, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2846 = GEN / "p2846_s1796_label_safe_vertex_localization_pullback_candidate_audit.json"
OUT = GEN / "p2847_s1797_target_independent_volume_form_unit_source_audit.json"
MD = GEN / "p2847_s1797_target_independent_volume_form_unit_source_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
EXPECTED_GRAPH_COUNT = 16828
REQUIRED_VOLUME_PREMISES = (
    "label_gauge_invariant_density",
    "finite_normalization_exists",
    "nonconstant_density_available",
    "target_independent_total_volume",
    "unit_dimension_source",
    "canonical_vertex_to_field_support",
    "spacetime_pullback_formula",
    "locality_covariance",
    "coupling_coefficient_rule",
    "variational_chain_rule",
)


def adjacency_sets(graph: tuple[tuple[int, ...], ...]) -> tuple[frozenset[int], ...]:
    return tuple(frozenset(neighbor - 1 for neighbor in neighbors) for neighbors in graph)


def triangle_counts(adj: tuple[frozenset[int], ...]) -> tuple[int, ...]:
    return tuple(sum(1 for a, b in combinations(sorted(adj[v]), 2) if b in adj[a]) for v in range(len(adj)))


def four_cycle_counts(adj: tuple[frozenset[int], ...]) -> tuple[int, ...]:
    counts = []
    for v in range(len(adj)):
        total = 0
        for opposite in range(len(adj)):
            if opposite == v or opposite in adj[v]:
                continue
            common = len(adj[v] & adj[opposite])
            total += common * (common - 1) // 2
        counts.append(total)
    return tuple(counts)


def density_vectors(adj: tuple[frozenset[int], ...]) -> dict[str, tuple[int, ...]]:
    tri = triangle_counts(adj)
    sq = four_cycle_counts(adj)
    return {
        "uniform_vertex_density": tuple(1 for _ in range(VERTEX_COUNT)),
        "triangle_plus_one_density": tuple(t + 1 for t in tri),
        "four_cycle_plus_one_density": tuple(s + 1 for s in sq),
        "triangle_square_plus_one_density": tuple((t + 1) * (s + 1) for t, s in zip(tri, sq)),
    }


def finite_volume_stats() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    totals: dict[str, list[int]] = {}
    profile_counts: dict[str, Counter[Any]] = {}
    for graph in graphs:
        adj = adjacency_sets(graph)
        for name, vec in density_vectors(adj).items():
            totals.setdefault(name, []).append(sum(vec))
            # A normalized probability profile is represented exactly by numerator multiset over total mass.
            total = sum(vec)
            normalized_profile = tuple(sorted(Counter(Fraction(v, total) for v in vec).items()))
            profile_counts.setdefault(name, Counter())[normalized_profile] += 1
    rows = {}
    for name, masses in totals.items():
        mass_hist = Counter(masses)
        rows[name] = {
            "distinct_raw_total_mass_count": len(mass_hist),
            "raw_total_mass_range": [min(masses), max(masses)],
            "top_raw_total_masses": mass_hist.most_common(10),
            "distinct_normalized_probability_profile_count": len(profile_counts[name]),
            "finite_normalization_exists": all(m > 0 for m in masses),
            "raw_total_mass_constant_on_carrier": len(mass_hist) == 1,
            "nonconstant_density_available": len(profile_counts[name]) > 1,
        }
    return {"decoded_graph_count": len(graphs), "candidate_rows": rows}


def premise_matrix(stats: dict[str, Any]) -> dict[str, Any]:
    rows = {}
    for name, row in stats["candidate_rows"].items():
        premises = {
            "label_gauge_invariant_density": True,
            "finite_normalization_exists": row["finite_normalization_exists"],
            "nonconstant_density_available": row["nonconstant_density_available"],
            "target_independent_total_volume": row["raw_total_mass_constant_on_carrier"],
            "unit_dimension_source": False,
            "canonical_vertex_to_field_support": False,
            "spacetime_pullback_formula": False,
            "locality_covariance": False,
            "coupling_coefficient_rule": False,
            "variational_chain_rule": False,
        }
        # Constant uniform mass is target-independent only as a finite counting measure;
        # it still lacks units, field support, and a coupling law.
        rows[name] = {
            "finite_stats": row,
            "premises": premises,
            "missing_premises": [key for key in REQUIRED_VOLUME_PREMISES if not premises[key]],
            "accepted_as_target_independent_unit_volume_source": all(premises.values()),
        }
    return rows


def build_payload(p2846: dict[str, Any]) -> dict[str, Any]:
    stats = finite_volume_stats()
    rows = premise_matrix(stats)
    accepted_count = sum(1 for row in rows.values() if row["accepted_as_target_independent_unit_volume_source"])
    facts = {
        "p2846_rechecked": p2846.get("status") == "P2846_LABEL_SAFE_VERTEX_LOCALIZATION_PULLBACK_CANDIDATE_NO_GO_NO_CLOSURE",
        "full_carrier_covered": stats["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "finite_probability_normalizations_exist": all(row["premises"]["finite_normalization_exists"] for row in rows.values()),
        "at_least_one_nonconstant_density_profile_exists": any(row["premises"]["nonconstant_density_available"] for row in rows.values()),
        "no_row_exports_unit_dimension_source": all(not row["premises"]["unit_dimension_source"] for row in rows.values()),
        "no_row_exports_field_support_or_pullback": all((not row["premises"]["canonical_vertex_to_field_support"] and not row["premises"]["spacetime_pullback_formula"]) for row in rows.values()),
        "accepted_count_zero": accepted_count == 0,
    }
    payload: dict[str, Any] = {
        "status": "P2847_TARGET_INDEPENDENT_VOLUME_FORM_UNIT_SOURCE_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2846": sha(P2846), "16_4_4.scd": sha(SCD)},
        "volume_form_unit_source_audit": {
            "input_statuses_rechecked": {"P2846": p2846.get("status")},
            "carrier_check": {"decoded_graph_count": stats["decoded_graph_count"], "expected_graph_count": EXPECTED_GRAPH_COUNT, "coverage_ok": stats["decoded_graph_count"] == EXPECTED_GRAPH_COUNT},
            "required_volume_premises": list(REQUIRED_VOLUME_PREMISES),
            "candidate_rows": rows,
            "accepted_candidate_count": accepted_count,
        },
        "acceptance_matrix": {"facts": facts, "accepted_as_volume_form_unit_source_obstruction_audit": all(facts.values()), "exports_target_independent_unit_volume_source": False},
        "decision": {
            "negative_export_flags": {
                "target_independent_volume_form_exported": False,
                "unit_dimension_source_exported": False,
                "canonical_field_support_exported": False,
                "coupling_coefficient_rule_exported": False,
                "variational_chain_rule_exported": False,
                "nonproxy_ltotal_term_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2847 finds that finite probability normalization is available for the tested label-safe densities, and nonconstant normalized profiles exist.  However, variable raw masses make nonconstant density-weighted total volumes carrier-dependent, while the uniform count is only a dimensionless counting measure.  No candidate exports a unit dimension source, canonical vertex-to-field support, spacetime pullback, coupling coefficient rule, or variational chain rule.",
            "next_honest_step": "Do not replay finite probability normalization as a unit-bearing volume form.  The next admissible move should isolate the alternative P2845 premise: one target-independent coupling coefficient/unit source for a single named density, with an explicit source law and units.  If that source law is not supplied, pivot to exactly one concrete kernel bridge atom with exported source premise, or preserve the no-new-live-frontier certificate.",
        },
    }
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["volume_form_unit_source_audit"]
    lines = [
        "# P2847/S1797 target-independent volume-form/unit-source audit",
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
        stats = row["finite_stats"]
        lines.extend([
            f"### {name}",
            f"- distinct_raw_total_mass_count={stats['distinct_raw_total_mass_count']}",
            f"- raw_total_mass_range={stats['raw_total_mass_range']}",
            f"- distinct_normalized_probability_profile_count={stats['distinct_normalized_probability_profile_count']}",
            f"- accepted_as_target_independent_unit_volume_source={row['accepted_as_target_independent_unit_volume_source']}",
            f"- missing_premises={row['missing_premises']}",
        ])
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2846))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2847/S1797 target-independent volume-form/unit-source audit", "## P2847/S1797 target-independent volume-form/unit-source audit\n\n`P2847/S1797` attacks the P2846-admissible volume/unit premise for named label-safe densities on the full `16,828` graph carrier.  Finite probability normalizations exist, and several normalized density profiles are nonconstant, but nonconstant raw density totals are carrier-dependent while the uniform total is only a dimensionless counting measure.  No candidate exports a target-independent unit-bearing volume form, canonical vertex-to-field support, spacetime pullback, coupling coefficient rule, variational chain rule, nonproxy `L_total`, EOM, Hamiltonian, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2847/S1797 volume unit Ltotal guard", "## P2847/S1797 volume/unit `L_total` guard\n\n`P2847/S1797` adds no action term.  Finite probability normalization of anonymous vertex densities is not a unit-bearing spacetime volume form and does not provide field support, pullback, coupling coefficients, or variational chain rules for `L_total`.\n")
    append_once(AGENTS, "Current target-independent volume/unit source guardrail (P2847/S1797, 2026-06-18)", "## Current target-independent volume/unit source guardrail (P2847/S1797, 2026-06-18)\n\n- P2847 tests the P2846-admissible target-independent volume/unit premise for named label-safe densities on the full `16,828` graph carrier.\n- Finite probability normalization exists and nonconstant normalized density profiles exist, but nonconstant raw total volumes are carrier-dependent and the uniform count is dimensionless; no unit source, field support, spacetime pullback, coupling coefficient rule, or variational chain rule is exported.\n- Do not promote finite probability normalization to unit-bearing volume form, strict localization/pullback, unit-bearing `L_total` source/coupling, EOM, Hamiltonian, bridge, role-transfer, or ToE closure.\n- The next admissible move should isolate one target-independent coupling coefficient/unit source for a single named density with explicit source law and units, or else one concrete kernel bridge atom with exported source premise; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
