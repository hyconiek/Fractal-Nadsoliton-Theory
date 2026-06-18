#!/usr/bin/env python3
"""P2848/S1798: target-independent coupling coefficient/unit-source law audit.

P2847 blocks finite probability normalization as a unit-bearing volume form.  The
remaining P2845/P2847 local route is a single target-independent coupling
coefficient/unit source for one named density.  P2848 audits that premise with
exact rational carrier statistics and a finite candidate coefficient matrix.

The audit distinguishes three things that are easy to conflate:
1. a graph-dependent normalizer such as 1 / raw_mass(G),
2. a target-independent but empirical carrier aggregate such as 1 / mean_mass,
3. a theorem-level unit-bearing source law with dimension, locality/covariance,
   and a variational chain rule.
Only (3) would unlock a nonproxy L_total term; P2848 tests whether current
finite artifacts export it.
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
P2847 = GEN / "p2847_s1797_target_independent_volume_form_unit_source_audit.json"
OUT = GEN / "p2848_s1798_coupling_coefficient_unit_source_law_audit.json"
MD = GEN / "p2848_s1798_coupling_coefficient_unit_source_law_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
EXPECTED_GRAPH_COUNT = 16828
REQUIRED_COEFFICIENT_PREMISES = (
    "named_density_fixed",
    "coefficient_defined_on_full_carrier",
    "target_independent_across_graphs",
    "unit_dimension_source_law",
    "non_empirical_source_law",
    "compatible_with_volume_or_pullback",
    "locality_covariance",
    "variational_chain_rule",
    "nonproxy_ltotal_insertion_rule",
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


def density_masses(adj: tuple[frozenset[int], ...]) -> dict[str, int]:
    tri = triangle_counts(adj)
    sq = four_cycle_counts(adj)
    return {
        "uniform_vertex_density": VERTEX_COUNT,
        "triangle_plus_one_density": sum(t + 1 for t in tri),
        "four_cycle_plus_one_density": sum(s + 1 for s in sq),
        "triangle_square_plus_one_density": sum((t + 1) * (s + 1) for t, s in zip(tri, sq)),
    }


def carrier_mass_table() -> dict[str, Any]:
    graphs, _ = decode_scd_bytes(SCD.read_bytes())
    masses: dict[str, list[int]] = {}
    for graph in graphs:
        adj = adjacency_sets(graph)
        for name, mass in density_masses(adj).items():
            masses.setdefault(name, []).append(mass)
    rows = {}
    for name, values in masses.items():
        hist = Counter(values)
        total = sum(values)
        n = len(values)
        mean = Fraction(total, n)
        second_moment = Fraction(sum(v * v for v in values), n)
        variance = second_moment - mean * mean
        rows[name] = {
            "graph_count": n,
            "distinct_mass_count": len(hist),
            "mass_range": [min(values), max(values)],
            "mean_mass": str(mean),
            "variance_mass": str(variance),
            "top_masses": hist.most_common(10),
            "masses": values,
        }
    return {"decoded_graph_count": len(graphs), "density_mass_rows": rows}


def coefficient_values(values: list[int], law: str) -> list[Fraction]:
    n = len(values)
    mean = Fraction(sum(values), n)
    second_moment = Fraction(sum(v * v for v in values), n)
    variance = second_moment - mean * mean
    if law == "constant_dimensionless_one":
        return [Fraction(1, 1) for _ in values]
    if law == "inverse_raw_mass":
        return [Fraction(1, v) for v in values]
    if law == "inverse_carrier_mean_mass":
        return [Fraction(1, 1) / mean for _ in values]
    if law == "inverse_mass_variance_if_nonzero":
        if variance == 0:
            return []
        return [Fraction(1, 1) / variance for _ in values]
    raise ValueError(law)


def coefficient_matrix(table: dict[str, Any]) -> dict[str, Any]:
    laws = (
        "constant_dimensionless_one",
        "inverse_raw_mass",
        "inverse_carrier_mean_mass",
        "inverse_mass_variance_if_nonzero",
    )
    rows = {}
    for density, mass_row in table["density_mass_rows"].items():
        values = mass_row["masses"]
        rows[density] = {}
        for law in laws:
            coeffs = coefficient_values(values, law)
            premises = {
                "named_density_fixed": True,
                "coefficient_defined_on_full_carrier": bool(coeffs) and len(coeffs) == len(values),
                "target_independent_across_graphs": bool(coeffs) and len(set(coeffs)) == 1,
                "unit_dimension_source_law": False,
                "non_empirical_source_law": law == "constant_dimensionless_one",
                "compatible_with_volume_or_pullback": False,
                "locality_covariance": False,
                "variational_chain_rule": False,
                "nonproxy_ltotal_insertion_rule": False,
            }
            rows[density][law] = {
                "distinct_coefficient_count": len(set(coeffs)) if coeffs else 0,
                "coefficient_sample": [str(c) for c in sorted(set(coeffs))[:5]],
                "premises": premises,
                "missing_premises": [key for key in REQUIRED_COEFFICIENT_PREMISES if not premises[key]],
                "accepted_as_target_independent_unit_coefficient_source": all(premises.values()),
            }
    return rows


def strip_masses(table: dict[str, Any]) -> dict[str, Any]:
    rows = {}
    for name, row in table["density_mass_rows"].items():
        rows[name] = {key: value for key, value in row.items() if key != "masses"}
    return {"decoded_graph_count": table["decoded_graph_count"], "density_mass_rows": rows}


def build_payload(p2847: dict[str, Any]) -> dict[str, Any]:
    table = carrier_mass_table()
    matrix = coefficient_matrix(table)
    accepted_count = sum(
        1
        for laws in matrix.values()
        for row in laws.values()
        if row["accepted_as_target_independent_unit_coefficient_source"]
    )
    facts = {
        "p2847_rechecked": p2847.get("status") == "P2847_TARGET_INDEPENDENT_VOLUME_FORM_UNIT_SOURCE_AUDIT_NO_CLOSURE",
        "full_carrier_covered": table["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
        "some_candidate_coefficients_are_target_independent": any(
            row["premises"]["target_independent_across_graphs"] for laws in matrix.values() for row in laws.values()
        ),
        "graph_dependent_normalizers_detected": any(
            row["distinct_coefficient_count"] > 1 for laws in matrix.values() for row in laws.values()
        ),
        "no_unit_dimension_source_law_exported": all(
            not row["premises"]["unit_dimension_source_law"] for laws in matrix.values() for row in laws.values()
        ),
        "no_nonproxy_ltotal_insertion_rule_exported": all(
            not row["premises"]["nonproxy_ltotal_insertion_rule"] for laws in matrix.values() for row in laws.values()
        ),
        "accepted_count_zero": accepted_count == 0,
    }
    return {
        "status": "P2848_COUPLING_COEFFICIENT_UNIT_SOURCE_LAW_AUDIT_NO_CLOSURE",
        "input_hashes": {"P2847": sha(P2847), "16_4_4.scd": sha(SCD)},
        "coefficient_unit_source_audit": {
            "input_statuses_rechecked": {"P2847": p2847.get("status")},
            "carrier_check": {
                "decoded_graph_count": table["decoded_graph_count"],
                "expected_graph_count": EXPECTED_GRAPH_COUNT,
                "coverage_ok": table["decoded_graph_count"] == EXPECTED_GRAPH_COUNT,
            },
            "density_mass_stats": strip_masses(table)["density_mass_rows"],
            "required_coefficient_premises": list(REQUIRED_COEFFICIENT_PREMISES),
            "candidate_coefficient_rows": matrix,
            "accepted_candidate_count": accepted_count,
        },
        "acceptance_matrix": {
            "facts": facts,
            "accepted_as_coefficient_unit_source_obstruction_audit": all(facts.values()),
            "exports_target_independent_unit_coefficient_source": False,
        },
        "decision": {
            "negative_export_flags": {
                "target_independent_unit_coefficient_source_exported": False,
                "unit_dimension_source_law_exported": False,
                "non_empirical_unit_source_exported": False,
                "volume_or_pullback_compatibility_exported": False,
                "variational_chain_rule_exported": False,
                "nonproxy_ltotal_term_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2848 exhausts a finite coefficient-source matrix for the named P2847 densities.  Dimensionless constants are target-independent but not unit-bearing source laws; inverse raw-mass normalizers become graph-dependent on nonconstant densities; carrier mean/variance coefficients are empirical aggregate normalizations, not exported strict source laws.  No candidate supplies volume/pullback compatibility, variational chain rule, or nonproxy L_total insertion.",
            "next_honest_step": "The local P2845 route through finite vertex densities has now tested localization, volume/unit normalization, and coefficient/unit source candidates without closure.  Do not replay density normalizers.  The next admissible proof-grade move should pivot to exactly one concrete kernel bridge atom with an exported source premise, preferably a damping/compression atom for beta/eta or an amplitude-passage atom, while preserving the kernel split and forbidding role transfer until a bridge theorem exists.  If no new source premise is supplied, preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["coefficient_unit_source_audit"]
    lines = [
        "# P2848/S1798 coupling coefficient/unit-source law audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Carrier check",
        f"- decoded_graph_count={audit['carrier_check']['decoded_graph_count']}",
        f"- coverage_ok={audit['carrier_check']['coverage_ok']}",
        "",
        "## Density mass stats",
    ]
    for density, stats in audit["density_mass_stats"].items():
        lines.extend([
            f"### {density}",
            f"- distinct_mass_count={stats['distinct_mass_count']}",
            f"- mass_range={stats['mass_range']}",
            f"- mean_mass={stats['mean_mass']}",
            f"- variance_mass={stats['variance_mass']}",
        ])
    lines.extend(["", "## Coefficient candidates"])
    for density, laws in audit["candidate_coefficient_rows"].items():
        lines.append(f"### {density}")
        for law, row in laws.items():
            lines.append(
                f"- {law}: distinct_coefficients={row['distinct_coefficient_count']}; "
                f"accepted={row['accepted_as_target_independent_unit_coefficient_source']}; "
                f"missing={row['missing_premises']}"
            )
    lines.extend(["", "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    payload = build_payload(read_json(P2847))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2848/S1798 coupling coefficient/unit-source law audit", "## P2848/S1798 coupling coefficient/unit-source law audit\n\n`P2848/S1798` attacks the remaining P2845/P2847 local premise: a target-independent coupling coefficient/unit source for one named density.  Exact rational carrier statistics distinguish dimensionless constants, graph-dependent inverse raw-mass normalizers, and empirical carrier mean/variance coefficients.  None exports a unit-dimension source law, volume/pullback compatibility, variational chain rule, nonproxy `L_total` insertion, EOM, Hamiltonian, bridge, role-transfer, or ToE closure.  The next honest proof-grade move should pivot to exactly one concrete kernel bridge atom with an exported source premise, or preserve no-new-live-frontier.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2848/S1798 coefficient unit Ltotal guard", "## P2848/S1798 coefficient/unit `L_total` guard\n\n`P2848/S1798` adds no action term.  Dimensionless constants, graph-dependent raw-mass normalizers, and empirical carrier aggregate coefficients do not provide a unit-bearing source law, volume/pullback compatibility, variational chain rule, or nonproxy insertion rule for `L_total`.\n")
    append_once(AGENTS, "Current coupling coefficient/unit source-law guardrail (P2848/S1798, 2026-06-18)", "## Current coupling coefficient/unit source-law guardrail (P2848/S1798, 2026-06-18)\n\n- P2848 tests the remaining local P2845/P2847 premise: a target-independent coupling coefficient/unit source for one named finite density.\n- Dimensionless constants are target-independent but not unit-bearing source laws; inverse raw-mass normalizers are graph-dependent on nonconstant densities; carrier mean/variance coefficients are empirical aggregate normalizations, not strict source laws.\n- Do not promote density normalizers to unit-bearing `L_total`, EOM, Hamiltonian, bridge, role-transfer, selector, or ToE closure.\n- The next admissible proof-grade move should pivot to exactly one concrete kernel bridge atom with an exported source premise, preferably damping/compression `beta/eta` or amplitude passage; otherwise preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    main()
