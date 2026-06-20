#!/usr/bin/env python3
"""P2971/S1921: typed support-provenance incidence complex.

P2970 closed the current ratio-package k-installation lane.  The remaining
honest move is to introduce exactly one typed structural object outside that
lane.  This module constructs a support/provenance incidence complex for the
P2938/P2961 aggregate that is richer than the old K/C split: it keeps component
provenance, slot support, weights, and aggregate-coordinate reception as typed
incidence data instead of collapsing to scalar ratios or K/C exchange.

The finite audit is positive as a structural construction, but negative as a
strict closure theorem: the object preserves the K/C support mismatch and has no
current source-localizer, unit-bearing coupling, variational chain rule, or
nonproxy L_total installation.
"""
from __future__ import annotations

import hashlib, itertools, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2970_s1920_nonconventional_nonproxy_k_installation_law_obstruction import OUT as P2970

OUT = GEN / "p2971_s1921_typed_support_provenance_incidence_complex.json"
MD = GEN / "p2971_s1921_typed_support_provenance_incidence_complex.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SLOTS = [
    {"slot": "K0", "component": "K", "aggregate_index": 0, "weight": 1, "provenance": "torsion_character_K_linear_seed"},
    {"slot": "K1", "component": "K", "aggregate_index": 1, "weight": 2, "provenance": "torsion_character_K_double_seed"},
    {"slot": "C0", "component": "C", "aggregate_index": 2, "weight": 2, "provenance": "torsion_character_C_equal_summand_a"},
    {"slot": "C1", "component": "C", "aggregate_index": 3, "weight": 2, "provenance": "torsion_character_C_equal_summand_b"},
    {"slot": "C2", "component": "C", "aggregate_index": 4, "weight": 2, "provenance": "torsion_character_C_equal_summand_c"},
]


def automorphism_count(slots: list[dict[str, Any]], preserve_provenance: bool) -> int:
    count = 0
    for perm in itertools.permutations(range(len(slots))):
        ok = True
        for i, j in enumerate(perm):
            a, b = slots[i], slots[j]
            if a["component"] != b["component"] or a["weight"] != b["weight"]:
                ok = False
                break
            if preserve_provenance and a["provenance"] != b["provenance"]:
                ok = False
                break
        if ok:
            count += 1
    return count


def component_certificate(slots: list[dict[str, Any]]) -> dict[str, Any]:
    components = sorted({s["component"] for s in slots})
    rows = []
    for component in components:
        sub = [s for s in slots if s["component"] == component]
        rows.append({
            "component": component,
            "support_size": len(sub),
            "weight_sum": sum(s["weight"] for s in sub),
            "weight_multiset": sorted(s["weight"] for s in sub),
            "provenance_labels": [s["provenance"] for s in sub],
        })
    return {
        "component_rows": rows,
        "aggregate_vector": [s["weight"] for s in sorted(slots, key=lambda r: r["aggregate_index"])],
        "aggregate_sum": sum(s["weight"] for s in slots),
        "primitive_mean": f"{sum(s['weight'] for s in slots)}/{len(slots)}",
        "automorphisms_preserving_component_and_weight": automorphism_count(slots, preserve_provenance=False),
        "automorphisms_preserving_component_weight_and_provenance": automorphism_count(slots, preserve_provenance=True),
    }


def structural_rows(cert: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "object": "typed_support_provenance_incidence_complex",
            "new_typed_structural_object": True,
            "outside_ratio_package_arithmetic": True,
            "preserves_KC_support_mismatch": True,
            "preserves_provenance_labels": True,
            "strict_source_localizer_exported": False,
            "unit_bearing_coupling_exported": False,
            "nonproxy_variational_chain_rule_exported": False,
            "accepted_current_strict_object": False,
            "witness": f"component automorphisms drop from {cert['automorphisms_preserving_component_and_weight']} to {cert['automorphisms_preserving_component_weight_and_provenance']} when provenance labels are kept",
        },
        {
            "object": "source_localizer_enriched_incidence_complex_schema",
            "new_typed_structural_object": True,
            "outside_ratio_package_arithmetic": True,
            "preserves_KC_support_mismatch": True,
            "preserves_provenance_labels": True,
            "strict_source_localizer_exported": True,
            "unit_bearing_coupling_exported": False,
            "nonproxy_variational_chain_rule_exported": False,
            "accepted_current_strict_object": False,
            "witness": "schema only; no strict nadsoliton localizer selects this incidence object",
        },
        {
            "object": "unit_bearing_action_coupled_incidence_schema",
            "new_typed_structural_object": True,
            "outside_ratio_package_arithmetic": True,
            "preserves_KC_support_mismatch": True,
            "preserves_provenance_labels": True,
            "strict_source_localizer_exported": True,
            "unit_bearing_coupling_exported": True,
            "nonproxy_variational_chain_rule_exported": True,
            "accepted_current_strict_object": True,
            "witness": "completed theorem schema; unavailable in current artifacts",
        },
    ]


def obligation_rows(rows: list[dict[str, Any]], cert: dict[str, Any]) -> list[dict[str, Any]]:
    current = rows[0]
    return [
        {"obligation": "new_typed_structural_object_constructed", "satisfied": current["new_typed_structural_object"], "evidence": "typed slots keep component, support, weight, aggregate coordinate, and provenance"},
        {"obligation": "outside_ratio_package_arithmetic", "satisfied": current["outside_ratio_package_arithmetic"], "evidence": "object is incidence/provenance data, not a scalar quotient, unit convention, or k-selection rule"},
        {"obligation": "KC_mismatch_preserved", "satisfied": current["preserves_KC_support_mismatch"], "evidence": "K support size 2 and C support size 3 remain distinct"},
        {"obligation": "provenance_noncollapse_checked", "satisfied": cert["automorphisms_preserving_component_weight_and_provenance"] < cert["automorphisms_preserving_component_and_weight"], "evidence": "provenance labels cut the component/weight automorphism group"},
        {"obligation": "strict_source_localizer_exported", "satisfied": current["strict_source_localizer_exported"], "evidence": "no current localizer selects this incidence complex as a strict source"},
        {"obligation": "unit_bearing_coupling_exported", "satisfied": current["unit_bearing_coupling_exported"], "evidence": "no unit-bearing coupling from the incidence complex into P2964/L_total is exported"},
        {"obligation": "nonproxy_variational_chain_rule_exported", "satisfied": current["nonproxy_variational_chain_rule_exported"], "evidence": "no variational chain rule is attached to the incidence object"},
        {"obligation": "accepted_current_strict_object", "satisfied": any(r["accepted_current_strict_object"] and r["object"] == "typed_support_provenance_incidence_complex" for r in rows), "evidence": "current object is developmental, not strict closure"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["typed_incidence_object", "outside_ratio_lane", "KC_mismatch_preserved", "provenance_preserved", "strict_source_localizer", "unit_bearing_coupling", "variational_chain_rule"]
    full = (1 << len(names)) - 1
    developmental_mask = sum(1 << i for i in range(4))
    return [{"mask": m, "present": {n: bool(m & (1 << i)) for i, n in enumerate(names)}, "accepts_developmental_structural_object": (m & developmental_mask) == developmental_mask, "accepts_strict_nonproxy_source_object": m == full} for m in range(1 << len(names))]


def build_payload(p2970_path: Any) -> dict[str, Any]:
    cert = component_certificate(SLOTS)
    rows = structural_rows(cert)
    obligations = obligation_rows(rows, cert)
    matrix = acceptance_matrix()
    return {
        "status": "P2971_TYPED_SUPPORT_PROVENANCE_INCIDENCE_COMPLEX_DEVELOPMENTAL_NO_STRICT_EXPORT",
        "input_hashes": {"P2970": hashlib.sha256(p2970_path.read_bytes()).hexdigest() if p2970_path.exists() else None},
        "constructed_theoretical_objects": {"candidate_object": "TypedSupportProvenanceIncidenceComplex", "slots": SLOTS, "component_certificate": cert, "structural_rows": rows, "proof_obligation_rows": obligations, "finite_acceptance_matrix": matrix},
        "structural_certificate": {"aggregate_vector": cert["aggregate_vector"], "aggregate_sum": cert["aggregate_sum"], "primitive_mean": cert["primitive_mean"], "component_weight_automorphisms": cert["automorphisms_preserving_component_and_weight"], "provenance_preserving_automorphisms": cert["automorphisms_preserving_component_weight_and_provenance"], "acceptance_matrix_rows": len(matrix), "developmental_rows": sum(1 for r in matrix if r["accepts_developmental_structural_object"]), "strict_rows": sum(1 for r in matrix if r["accepts_strict_nonproxy_source_object"])},
        "decision": {
            "positive_progress": "P2971 constructs a genuinely new typed support/provenance incidence object outside scalar ratio-package arithmetic and preserves the K/C support and provenance mismatch instead of erasing it.",
            "breakthrough": "No strict source or L_total closure is exported: the incidence complex lacks a strict source-localizer, unit-bearing coupling, and nonproxy variational chain rule.",
            "negative_export_flags": {k: False for k in ["strict_source_localizer_exported", "unit_bearing_coupling_exported", "nonproxy_variational_chain_rule_exported", "strict_ratio_package_source_exported", "nonproxy_ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay ratio-package arithmetic, K/C exchange, scalar mediator collapse, unit conventions, k-selection predicates, or formal Euler placeholders.  The next proof-grade move must add exactly one missing theorem to this incidence object: a strict nadsoliton source-localizer selecting it, a unit-bearing coupling into a named action density, or a nonproxy variational chain rule; otherwise preserve the P2929-P2971 developmental/no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["structural_certificate"]
    lines = ["# P2971/S1921 typed support-provenance incidence complex", "", f"Status: `{payload['status']}`", "", "## Structural certificate", f"- aggregate vector: `{cert['aggregate_vector']}`", f"- aggregate sum / primitive mean: `{cert['aggregate_sum']}` / `{cert['primitive_mean']}`", f"- component-weight automorphisms: `{cert['component_weight_automorphisms']}`", f"- provenance-preserving automorphisms: `{cert['provenance_preserving_automorphisms']}`", f"- acceptance matrix rows / developmental / strict: `{cert['acceptance_matrix_rows']}` / `{cert['developmental_rows']}` / `{cert['strict_rows']}`", "", "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2970)
    payload = build_payload(P2970)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2971/S1921 typed support-provenance incidence complex", "## P2971/S1921 typed support-provenance incidence complex\n\n`P2971/S1921` executes the P2970 pivot by constructing one typed structural object outside scalar ratio-package arithmetic: a support/provenance incidence complex for the P2938/P2961 aggregate.  It preserves component provenance, K support size `2`, C support size `3`, aggregate vector `[1,2,2,2,2]`, and shows provenance labels reduce the component/weight automorphism count from `6` to `1`.  This is developmental structural progress, not strict closure: no strict source-localizer, unit-bearing coupling, nonproxy variational chain rule, nonproxy `L_total`, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2971/S1921 incidence-complex `L_total` guard", "## P2971/S1921 incidence-complex `L_total` guard\n\n`P2971/S1921` adds a typed support/provenance incidence object rather than a scalar coefficient or unit convention.  Because no strict source-localizer, unit-bearing action-density coupling, or nonproxy variational chain rule is attached to the incidence complex, no sourced term enters `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current typed support-provenance incidence complex guardrail (P2971/S1921, 2026-06-20)", "## Current typed support-provenance incidence complex guardrail (P2971/S1921, 2026-06-20)\n\n- P2971 constructs a new typed support/provenance incidence complex outside scalar ratio-package arithmetic: it preserves K/C support mismatch, provenance labels, and aggregate-coordinate reception.\n- The finite audit is developmental only: provenance labels reduce component/weight automorphisms, but no strict nadsoliton source-localizer, unit-bearing coupling, or nonproxy variational chain rule is exported.\n- Do not promote the incidence complex to strict ratio-package source, nonproxy `L_total`, bridge closure, role transfer, or ToE by replaying K/C exchange, scalar mediator collapse, unit conventions, k-selection predicates, or formal Euler placeholders.\n- A next admissible move must add exactly one missing theorem to this incidence object: a strict source-localizer, a unit-bearing coupling into a named action density, or a nonproxy variational chain rule; otherwise preserve the P2929-P2971 boundary.\n")
    return payload

if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
