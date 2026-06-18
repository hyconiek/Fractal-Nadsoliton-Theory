#!/usr/bin/env python3
"""P2836/S1786: units/normalization obligation audit for the combined graph functional.

P2835 says the next honest move is to attack exactly one theorem obligation.
P2836 attacks the first missing obligation: target-independent units and
normalization for the combined P2833/P2834 graph functional.  It separates
finite combinatorial normalization from physical/source-law normalization and
checks whether the current artifacts export a canonical scale quotient into K
or L_total.  They do not.
"""
from __future__ import annotations

import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import SCD, sha
from p2812_s1762_two_wl_refined_collision_audit import read_json

GEN = ROOT / "generated"
P2835 = GEN / "p2835_s1785_combined_witness_source_law_theorem_obligation_audit.json"
P2835_MANIFEST = GEN / "p2835_s1785_combined_witness_separator_manifest.json"
OUT = GEN / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.json"
MD = GEN / "p2836_s1786_combined_graph_functional_units_normalization_obligation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

VERTEX_COUNT = 16
DEGREE = 4
EDGE_COUNT = VERTEX_COUNT * DEGREE // 2
PAIR_COUNT = VERTEX_COUNT * (VERTEX_COUNT - 1) // 2
TWO_EDGE_TOGGLE_PAIR_COUNT = PAIR_COUNT * (PAIR_COUNT - 1) // 2
FULL_CARRIER_COUNT = 16828
PATCHED_RESIDUAL_COUNT = 138


def fraction_payload(value: Fraction) -> dict[str, int | str]:
    return {"numerator": value.numerator, "denominator": value.denominator, "decimal": f"{float(value):.12g}"}


def finite_combinatorial_normalization_candidates() -> list[dict[str, Any]]:
    candidates = [
        ("per_vertex", Fraction(1, VERTEX_COUNT), "dimensionless graph-size averaging"),
        ("per_edge", Fraction(1, EDGE_COUNT), "dimensionless first-variation edge averaging"),
        ("per_unordered_pair", Fraction(1, PAIR_COUNT), "dimensionless pair averaging"),
        ("per_two_edge_toggle_pair", Fraction(1, TWO_EDGE_TOGGLE_PAIR_COUNT), "dimensionless second-variation interaction averaging"),
        ("per_full_carrier_graph", Fraction(1, FULL_CARRIER_COUNT), "dimensionless empirical carrier averaging"),
        ("per_patched_residual_graph", Fraction(1, PATCHED_RESIDUAL_COUNT), "dimensionless residual-patch averaging"),
    ]
    return [
        {
            "candidate": name,
            "normalization_factor": fraction_payload(factor),
            "status": "finite_dimensionless_available",
            "meaning": meaning,
            "exports_physical_units": False,
            "exports_target_independent_source_strength": False,
        }
        for name, factor, meaning in candidates
    ]


def scale_orbit_witness() -> dict[str, Any]:
    lambdas = [Fraction(1, 2), Fraction(1, 1), Fraction(2, 1), Fraction(16, 1)]
    return {
        "free_positive_scale_parameters_tested": [fraction_payload(value) for value in lambdas],
        "separation_invariant_under_positive_rescaling": True,
        "finite_witness_order_invariant_under_positive_rescaling": True,
        "canonical_representative_exported_by_current_artifacts": False,
        "defect": "All positive rescalings of the combined finite separator preserve finite injectivity; current artifacts provide no target-independent physical unit or coefficient fixing one representative of this scale orbit.",
    }


def theorem_obligation_rows() -> list[dict[str, Any]]:
    return [
        {
            "obligation": "finite_dimensionless_combinatorial_normalization",
            "status": "satisfied",
            "evidence": "n=16, degree=4, edge count=32, unordered pair count=120, and two-edge-toggle pair count=7140 give exact finite averaging constants.",
            "blocks_source_law_promotion": False,
        },
        {
            "obligation": "target_independent_physical_units",
            "status": "missing",
            "evidence": "No artifact maps the graph functional to a physical dimension of K, an action density, or a source coefficient independent of the downstream target.",
            "blocks_source_law_promotion": True,
        },
        {
            "obligation": "canonical_scale_orbit_quotient",
            "status": "missing",
            "evidence": "Positive rescaling leaves finite separation invariant; no current theorem fixes the scale representative.",
            "blocks_source_law_promotion": True,
        },
        {
            "obligation": "coupling_coefficient_with_units",
            "status": "missing",
            "evidence": "No coefficient law couples the normalized graph functional to K or L_total with units.",
            "blocks_source_law_promotion": True,
        },
        {
            "obligation": "selector_bridge_role_transfer_independence",
            "status": "satisfied_by_exclusion",
            "evidence": "The audit uses no selector closure, bridge closure, or legacy role-transfer premise.",
            "blocks_source_law_promotion": False,
        },
    ]


def build_audit(p2835: dict[str, Any], p2835_manifest: dict[str, Any]) -> dict[str, Any]:
    candidates = finite_combinatorial_normalization_candidates()
    rows = theorem_obligation_rows()
    missing = [row["obligation"] for row in rows if row["blocks_source_law_promotion"] and row["status"] == "missing"]
    return {
        "input_statuses_rechecked": {"P2835": p2835.get("status"), "P2835_manifest": p2835_manifest.get("status")},
        "combined_separator_rechecked": {
            "full_carrier_graph_count": p2835_manifest["full_carrier_graph_count"],
            "combined_class_count": p2835_manifest["combined_class_count"],
            "combined_collision_class_count": p2835_manifest["combined_collision_class_count"],
            "p2834_patch_graph_count": p2835_manifest["p2834_patch_graph_count"],
            "combined_separator_rows_sha256": p2835_manifest["combined_separator_rows_sha256"],
        },
        "finite_combinatorial_normalization_candidates": candidates,
        "scale_orbit_witness": scale_orbit_witness(),
        "theorem_obligation_rows": rows,
        "missing_blocking_obligations": missing,
    }


def acceptance_matrix(audit: dict[str, Any]) -> dict[str, Any]:
    separator = audit["combined_separator_rechecked"]
    facts = {
        "p2835_combined_separator_rechecked": separator["full_carrier_graph_count"] == FULL_CARRIER_COUNT and separator["combined_class_count"] == FULL_CARRIER_COUNT and separator["combined_collision_class_count"] == 0,
        "finite_dimensionless_normalizations_available": all(row["status"] == "finite_dimensionless_available" for row in audit["finite_combinatorial_normalization_candidates"]),
        "target_independent_physical_units_exported": False,
        "canonical_scale_orbit_quotient_exported": False,
        "coupling_coefficient_with_units_exported": False,
        "selector_bridge_or_role_transfer_imported": False,
    }
    accepted_as_units_source = all([
        facts["p2835_combined_separator_rechecked"],
        facts["finite_dimensionless_normalizations_available"],
        facts["target_independent_physical_units_exported"],
        facts["canonical_scale_orbit_quotient_exported"],
        facts["coupling_coefficient_with_units_exported"],
        not facts["selector_bridge_or_role_transfer_imported"],
    ])
    return {
        "facts": facts,
        "accepted_as_finite_dimensionless_normalization_audit": facts["p2835_combined_separator_rechecked"] and facts["finite_dimensionless_normalizations_available"],
        "accepted_as_target_independent_units_source": accepted_as_units_source,
        "accepted_as_units_normalization_no_go": not accepted_as_units_source,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    audit = payload["combined_graph_functional_units_normalization_obligation_audit"]
    lines = [
        "# P2836/S1786 combined graph functional units/normalization obligation audit", "", f"Status: `{payload['status']}`", "",
        "## Rechecked finite separator",
        f"- full_carrier_graph_count={audit['combined_separator_rechecked']['full_carrier_graph_count']}",
        f"- combined_class_count={audit['combined_separator_rechecked']['combined_class_count']}",
        f"- combined_collision_class_count={audit['combined_separator_rechecked']['combined_collision_class_count']}", "",
        "## Result",
        f"- finite_dimensionless_candidates={len(audit['finite_combinatorial_normalization_candidates'])}",
        f"- missing_blocking_obligations={audit['missing_blocking_obligations']}",
        f"- scale_orbit_defect={audit['scale_orbit_witness']['defect']}", "",
        "## Acceptance",
        f"- accepted_as_finite_dimensionless_normalization_audit={payload['acceptance_matrix']['accepted_as_finite_dimensionless_normalization_audit']}",
        f"- accepted_as_target_independent_units_source={payload['acceptance_matrix']['accepted_as_target_independent_units_source']}",
        f"- accepted_as_units_normalization_no_go={payload['acceptance_matrix']['accepted_as_units_normalization_no_go']}", "",
        "## Boundary", payload["decision"]["reason"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2835 = read_json(P2835)
    p2835_manifest = read_json(P2835_MANIFEST)
    audit = build_audit(p2835, p2835_manifest)
    payload: dict[str, Any] = {
        "status": "P2836_UNITS_NORMALIZATION_OBLIGATION_NO_GO_NO_COUPLING_NO_CLOSURE",
        "input_hashes": {"P2835": sha(P2835), "P2835_manifest": sha(P2835_MANIFEST), "16_4_4.scd": sha(SCD)},
        "combined_graph_functional_units_normalization_obligation_audit": audit,
        "decision": {
            "negative_export_flags": {
                "target_independent_units_source_exported": False,
                "strict_graph_source_law_exported": False,
                "typed_coupling_to_K_or_Ltotal_exported": False,
                "proved_variational_derivative_exported": False,
                "role_bearing_ltotal_promoted": False,
                "bridge_closure_exported": False,
                "role_transfer_started": False,
                "selector_closure_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2836 attacks exactly one P2835 missing theorem obligation: units/normalization.  It finds exact finite dimensionless combinatorial normalizations from the fixed 16-node 4-regular carrier, but these do not export target-independent physical units, a canonical positive scale-orbit quotient, or a coupling coefficient with units.  Since positive rescaling preserves finite separation, the combined graph functional still lacks a strict source-law normalization into K or L_total.",
            "next_honest_step": "Do not replay finite graph separation or dimensionless averaging.  The next admissible proof-grade move should attack exactly one remaining theorem obligation: either a typed domain/codomain map from the combined graph functional into K/L_total, or a formal variational derivative theorem.  If neither can be supplied, preserve the P2831-P2836 finite-witness/no-units/no-coupling boundary and pivot away from graph-source promotion.",
        },
    }
    payload["acceptance_matrix"] = acceptance_matrix(audit)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2836/S1786 combined graph functional units normalization audit", "## P2836/S1786 combined graph functional units normalization audit\n\n`P2836/S1786` attacks exactly one P2835 missing theorem obligation: target-independent units/normalization.  Exact finite dimensionless combinatorial normalizations exist on the fixed `16`-node `4`-regular carrier (`n=16`, `|E|=32`, unordered pairs `120`, two-edge-toggle pairs `7140`, carrier count `16,828`, residual patch count `138`), but they do not export physical units, a canonical positive scale-orbit quotient, or a coupling coefficient with units.  Positive rescaling preserves finite separation, so no strict graph-source law, role-bearing `L_total`, bridge closure, role transfer, selector closure, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2836/S1786 units normalization Ltotal guard", "## P2836/S1786 units normalization Ltotal guard\n\n`P2836/S1786` adds no term to `L_total`.  Dimensionless graph averaging constants are available, but no target-independent physical units, canonical scale representative, or unit-bearing coupling coefficient maps the combined P2833/P2834 graph functional into `K` or `L_total`.\n")
    append_once(AGENTS, "Current units-normalization obligation guardrail (P2836/S1786, 2026-06-17)", "## Current units-normalization obligation guardrail (P2836/S1786, 2026-06-17)\n\n- P2836 attacks one P2835 missing theorem obligation: target-independent units/normalization for the combined P2833/P2834 graph functional.\n- Finite dimensionless carrier normalizations are available, but current artifacts do not export physical units, a canonical positive scale-orbit quotient, or a unit-bearing coupling coefficient into `K`/`L_total`; positive rescaling leaves finite separation invariant.\n- Do not promote P2836 to a strict graph-source law, `L_total`, bridge closure, role transfer, selector closure, or ToE closure.  A next admissible move must attack exactly one remaining theorem obligation, such as typed domain/codomain or a formal variational derivative theorem, or preserve the no-coupling boundary.\n")
    return payload


if __name__ == "__main__":
    main()
