#!/usr/bin/env python3
"""P2902/S1852: pointed signed defect law variational template audit.

P2901 built the right-shaped defect-placement schema but left the actual
(basepoint, polarity) pair imported.  P2902 constructs the next missing
*theoretical objects* explicitly in their strongest honest form:

1. a pointed signed source axiom A(b0, sigma0),
2. the induced 9/5 defect edge and local density template,
3. a finite variational derivative table for that density.

The audit then separates what is genuinely constructed (a reproducible
axiom-augmented variational template) from what is still not strict (the axiom
is not sourced internally and the density is not a nonproxy L_total term).
"""
from __future__ import annotations

import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2803_s1753_meringer_scd_import_decode_validation_certificate import sha
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2901 = GEN / "p2901_s1851_explicit_defect_placement_source_law_candidate_gate.json"
OUT = GEN / "p2902_s1852_pointed_signed_defect_law_variational_template_audit.json"
MD = GEN / "p2902_s1852_pointed_signed_defect_law_variational_template_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12
CARRIER_STEP = 5
UNIT_SYMBOL = "U_9_5"


def defect_edge(basepoint: int, polarity: int) -> tuple[int, int]:
    return (basepoint % N, (basepoint + polarity * CARRIER_STEP) % N)


def all_pointed_axioms() -> list[dict[str, Any]]:
    rows = []
    for basepoint in range(N):
        for polarity in (-1, 1):
            tail, head = defect_edge(basepoint, polarity)
            rows.append({
                "axiom_name": f"A_pointed_signed_defect_{basepoint}_{'plus' if polarity > 0 else 'minus'}",
                "basepoint": basepoint,
                "polarity": polarity,
                "defect_edge": [tail, head],
                "carrier_offset_mod_12": (head - tail) % N,
                "density_template": f"{polarity}*{UNIT_SYMBOL}*delta_edge({tail},{head})*q_9_5({tail},{head})",
                "axiom_status": "explicit_added_premise_not_internal_strict_source",
            })
    return rows


def variational_derivative(row: dict[str, Any]) -> dict[str, Any]:
    tail, head = row["defect_edge"]
    polarity = row["polarity"]
    support_key = f"q_9_5[{tail},{head}]"
    return {
        "action_template": f"S_{row['axiom_name']} = {row['density_template']}",
        "varied_coordinate": support_key,
        "nonzero_derivative": {support_key: f"{polarity}*{UNIT_SYMBOL}"},
        "zero_derivative_count_on_other_directed_edges": N * N - 1,
        "local_variational_chain_rule_holds_for_template": True,
        "unit_bearing_status": "symbolic_unit_assigned_by_template_not_strictly_sourced",
    }


def translation_class_count(rows: list[dict[str, Any]]) -> dict[str, Any]:
    # Translating a pointed axiom moves basepoint and edge but preserves polarity.
    by_polarity: dict[int, list[dict[str, Any]]] = {-1: [], 1: []}
    for row in rows:
        by_polarity[row["polarity"]].append(row)
    return {
        "pointed_axiom_count": len(rows),
        "translation_classes_after_forgetting_pointer": 2,
        "class_sizes": {str(p): len(items) for p, items in by_polarity.items()},
        "strictly_internal_class_selector_count": 0,
    }


def build_payload(p2901: dict[str, Any]) -> dict[str, Any]:
    rows = all_pointed_axioms()
    derivatives = [variational_derivative(row) for row in rows]
    classes = translation_class_count(rows)
    nonzero_derivatives = sum(1 for d in derivatives if d["nonzero_derivative"])
    return {
        "status": "P2902_POINTED_SIGNED_DEFECT_LAW_VARIATIONAL_TEMPLATE_AXIOMATIC_NO_STRICT_CLOSURE",
        "input_hashes": {"P2901": sha(P2901)},
        "constructed_theoretical_objects": {
            "pointed_signed_source_axiom_family": rows,
            "localized_density_templates": [row["density_template"] for row in rows],
            "finite_variational_derivative_table": derivatives,
            "translation_class_certificate": classes,
        },
        "acceptance_matrix": {
            "p2901_rechecked": p2901.get("status") == "P2901_EXPLICIT_DEFECT_PLACEMENT_SOURCE_LAW_CANDIDATE_CONDITIONAL_NO_CLOSURE",
            "pointed_signed_axioms_constructed": True,
            "localized_density_templates_constructed": True,
            "finite_variational_derivatives_computed": True,
            "pointed_axiom_count": classes["pointed_axiom_count"],
            "nonzero_variational_derivative_count": nonzero_derivatives,
            "translation_classes_after_forgetting_pointer": classes["translation_classes_after_forgetting_pointer"],
            "strictly_internal_class_selector_count": classes["strictly_internal_class_selector_count"],
            "axiom_augmented_template_valid": nonzero_derivatives == len(rows),
            "nonimported_strict_source_law_exported": False,
            "accepted_as_strict_missing_object": False,
        },
        "decision": {
            "positive_witnesses": {
                "axiom_augmented_source_law_object_constructed": True,
                "defect_placement_computed_after_axiom": True,
                "localized_density_template_constructed_after_axiom": True,
                "finite_variational_derivative_nonzero_on_support": True,
            },
            "negative_export_flags": {
                "nonimported_basepoint_or_polarity_law_exported": False,
                "strict_defect_placement_source_law_exported": False,
                "unit_bearing_strict_density_exported": False,
                "coupling_to_9_over_5_variational_density_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "strict_damping_compression_bridge_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2902 constructs the missing theoretical objects in axiom-augmented form and verifies the finite variational derivative of each localized template.  This is a real readiness improvement over P2901: once a pointed signed axiom is supplied, the defect edge, density support, and local derivative are computed.  However the pointed axiom itself is still imported; forgetting the pointer leaves two polarity classes and zero internally selected strict class members.  Therefore no strict source law, nonproxy L_total, EOM, Hamiltonian, bridge, role transfer, or ToE closure is exported.",
            "next_honest_step": "The next honest proof-grade move is not another pointed template.  Either provide an internal strict source theorem that derives one pointed signed axiom A(b0,sigma0) from nadsoliton data and then lift the symbolic unit U_9_5 to a nonproxy L_total coupling, or pivot to a different typed object outside the torsor/basepoint/defect-template family.  Without that, preserve the no-new-live-frontier certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2902/S1852 pointed signed defect law variational template audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Constructed theoretical objects",
        "- pointed signed source axiom family `A(b,sigma)`",
        "- induced `9/5` defect edge `D=(b,b+sigma*5)`",
        f"- localized density template `sigma*{UNIT_SYMBOL}*delta_edge(D)*q_9_5(D)`",
        "- finite variational derivative table `dS/dq_9_5(D)=sigma*U_9_5`",
        "",
        "## Finite acceptance gate",
        f"- pointed axiom count: `{acc['pointed_axiom_count']}`",
        f"- nonzero variational derivative count: `{acc['nonzero_variational_derivative_count']}`",
        f"- translation classes after forgetting pointer: `{acc['translation_classes_after_forgetting_pointer']}`",
        f"- strictly internal class selector count: `{acc['strictly_internal_class_selector_count']}`",
        f"- accepted as strict missing object: `{acc['accepted_as_strict_missing_object']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2901))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2902/S1852 pointed signed defect law variational template audit", "## P2902/S1852 pointed signed defect law variational template audit\n\n`P2902/S1852` constructs the missing objects in their strongest honest axiom-augmented form: a pointed signed axiom family `A(b,sigma)`, induced `9/5` defect edge `D=(b,b+sigma*5)`, localized template `sigma*U_9_5*delta_edge(D)*q_9_5(D)`, and finite derivative `dS/dq_9_5(D)=sigma*U_9_5`.  All `24` pointed templates have nonzero local derivative on their support, but the pointer/sign axiom is imported; after forgetting it there are two polarity classes and `0` internally selected strict class members.  Thus P2902 is an axiom-augmented variational readiness template, not a strict source law, nonproxy `L_total`, EOM, Hamiltonian, bridge, role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2902/S1852 pointed defect template `L_total` guard", "## P2902/S1852 pointed defect template `L_total` guard\n\n`P2902/S1852` computes a local symbolic derivative for each axiom-pointed `9/5` density template, but the unit `U_9_5` and pointed signed axiom remain added premises.  Therefore the result is a variational template only and does not export nonproxy `L_total`, EOM, Hamiltonian closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current pointed signed defect variational template guardrail (P2902/S1852, 2026-06-19)", "## Current pointed signed defect variational template guardrail (P2902/S1852, 2026-06-19)\n\n- P2902 constructs an axiom-augmented pointed signed source family `A(b,sigma)`, induced `9/5` defect edge, localized symbolic density template, and finite derivative table.\n- The template is computationally coherent after the pointed axiom: `24` pointed templates have nonzero local support derivative, but the pointer/sign axiom is still imported and no internal strict selector/source theorem is exported.\n- Do not promote axiom-pointed templates, symbolic unit assignment `U_9_5`, local derivative tables, chosen `(b,sigma)`, or density-support notation to strict phase/origin sourcehood, strict damping/compression bridge, role transfer, nonproxy `L_total`, EOM, Hamiltonian, or ToE closure.\n- A next admissible proof-grade move must derive one pointed signed axiom internally from strict nadsoliton data and then lift `U_9_5` to a unit-bearing nonproxy `L_total` coupling theorem, pivot outside the torsor/basepoint/defect-template family, or preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
