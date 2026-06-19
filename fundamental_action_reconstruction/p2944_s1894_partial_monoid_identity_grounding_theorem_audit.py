#!/usr/bin/env python3
"""P2944/S1894: partial-monoid identity-grounding theorem audit.

P2943 showed that the P2938 carrier has a finite identity-grounded positive
cone: f(1)=0 and every nonidentity audited node has positive value.  The next
missing theorem object is identity grounding itself.  P2944 constructs the
finite algebraic part of that object by auditing the partial multiplication
monoid on nodes 1..11 with products retained exactly when d*e<=11.

The finite result is positive but bounded: node 1 is the unique two-sided total
identity on the audited partial multiplication structure, and it is the unique
zero of the P2938 positive cone.  This is an algebraic identity-grounding
witness, not yet a strict nadsoliton identity-grounding theorem and not a
delta/eta or beta/eta coupling theorem.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2940_s1890_p2938_carrier_aut_orbit_selector_burden import NODES, carrier_value
from p2943_s1893_identity_grounded_positive_cone_orientation_theorem_audit import OUT as P2943

OUT = GEN / "p2944_s1894_partial_monoid_identity_grounding_theorem_audit.json"
MD = GEN / "p2944_s1894_partial_monoid_identity_grounding_theorem_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def product_defined(a: int, b: int) -> bool:
    return a * b <= NODES[-1]


def identity_candidate_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for candidate in NODES:
        left_defined_all = all(product_defined(candidate, node) for node in NODES)
        right_defined_all = all(product_defined(node, candidate) for node in NODES)
        left_identity_all = left_defined_all and all(candidate * node == node for node in NODES)
        right_identity_all = right_defined_all and all(node * candidate == node for node in NODES)
        rows.append({
            "candidate": candidate,
            "left_products_defined_for_all_nodes": left_defined_all,
            "right_products_defined_for_all_nodes": right_defined_all,
            "left_identity_for_all_nodes": left_identity_all,
            "right_identity_for_all_nodes": right_identity_all,
            "two_sided_total_identity": left_identity_all and right_identity_all,
            "carrier_value": carrier_value(candidate),
            "zero_carrier_value": carrier_value(candidate) == 0,
        })
    return rows


def partial_product_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for a in NODES:
        for b in NODES:
            if product_defined(a, b):
                rows.append({
                    "a": a,
                    "b": b,
                    "product": a * b,
                    "carrier_additive": carrier_value(a * b) == carrier_value(a) + carrier_value(b),
                    "identity_left_stable": a != 1 or a * b == b,
                    "identity_right_stable": b != 1 or a * b == a,
                })
    return rows


def theorem_rows(identity_rows: list[dict[str, Any]], product_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    identities = [row for row in identity_rows if row["two_sided_total_identity"]]
    zero_nodes = [row for row in identity_rows if row["zero_carrier_value"]]
    return [
        {
            "premise": "unique_two_sided_total_identity_in_audited_partial_monoid",
            "satisfied": [row["candidate"] for row in identities] == [1],
            "evidence": f"identity candidates: {[row['candidate'] for row in identities]}",
        },
        {
            "premise": "unique_zero_carrier_node_matches_identity",
            "satisfied": [row["candidate"] for row in zero_nodes] == [1],
            "evidence": f"zero carrier nodes: {[row['candidate'] for row in zero_nodes]}",
        },
        {
            "premise": "identity_products_stable",
            "satisfied": all(row["identity_left_stable"] and row["identity_right_stable"] for row in product_rows),
            "evidence": "all retained products with identity factor are stable",
        },
        {
            "premise": "carrier_additive_on_partial_products",
            "satisfied": all(row["carrier_additive"] for row in product_rows),
            "evidence": f"{len(product_rows)} retained partial products audited",
        },
        {
            "premise": "finite_identity_grounding_witness_constructed",
            "satisfied": True,
            "evidence": "the unique partial-monoid identity and the unique positive-cone zero coincide at node 1",
        },
    ]


def strict_acceptance_rows(algebraic_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "criterion": "finite_partial_monoid_identity_grounding_verified",
            "satisfied": all(row["satisfied"] for row in algebraic_rows),
            "evidence": "unique two-sided total identity and unique carrier zero coincide",
        },
        {
            "criterion": "strict_nadsoliton_identity_grounding_theorem_exported",
            "satisfied": False,
            "evidence": "finite partial-monoid identity is not yet proved to be a strict nadsoliton source theorem",
        },
        {
            "criterion": "strict_positive_orientation_source_theorem_exported",
            "satisfied": False,
            "evidence": "identity grounding remains algebraic; strict positive orientation provenance is not exported",
        },
        {
            "criterion": "delta_eta_source_law_exported",
            "satisfied": False,
            "evidence": "no delta/eta source law is derived from the partial-monoid identity witness",
        },
        {
            "criterion": "beta_eta_coupling_theorem_exported",
            "satisfied": False,
            "evidence": "no beta/eta coupling theorem is derived from the partial-monoid identity witness",
        },
    ]


def build_payload(_: dict[str, Any]) -> dict[str, Any]:
    identities = identity_candidate_rows()
    products = partial_product_rows()
    algebraic = theorem_rows(identities, products)
    strict_rows = strict_acceptance_rows(algebraic)
    accepted = all(row["satisfied"] for row in strict_rows)
    identity_nodes = [row["candidate"] for row in identities if row["two_sided_total_identity"]]
    zero_nodes = [row["candidate"] for row in identities if row["zero_carrier_value"]]
    return {
        "status": "P2944_PARTIAL_MONOID_IDENTITY_GROUNDING_THEOREM_AUDIT_NO_STRICT_SOURCE",
        "input_hashes": {"P2943": hashlib.sha256(P2943.read_bytes()).hexdigest() if P2943.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "Audited_Partial_Multiplication_Monoid_Identity_Grounding_Witness",
            "identity_candidate_rows": identities,
            "partial_product_rows": products,
            "algebraic_theorem_rows": algebraic,
            "strict_acceptance_rows": strict_rows,
        },
        "identity_grounding_certificate": {
            "node_count": len(NODES),
            "partial_product_count": len(products),
            "two_sided_total_identity_nodes": identity_nodes,
            "zero_carrier_nodes": zero_nodes,
            "identity_equals_unique_zero": identity_nodes == [1] and zero_nodes == [1],
            "finite_identity_grounding_verified": strict_rows[0]["satisfied"],
            "strict_identity_grounding_theorem_exported": False,
            "accepted_strict_source": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "finite_partial_monoid_identity_grounding_verified": strict_rows[0]["satisfied"],
                "identity_node_matches_unique_positive_cone_zero": identity_nodes == [1] and zero_nodes == [1],
            },
            "negative_export_flags": {
                "strict_identity_grounding_theorem_exported": False,
                "strict_positive_orientation_source_theorem_exported": False,
                "strict_selector_source_exported": False,
                "strict_aut_breaking_prime_coordinate_source_law_exported": False,
                "strict_prime_log_value_source_exported": False,
                "strict_delta_eta_source_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2944 proves the finite partial-monoid identity witness demanded after P2943: node 1 is the unique two-sided total identity on the audited retained products and the unique zero of the P2938 positive cone.  This is still an algebraic witness, not a strict nadsoliton identity-grounding theorem and not a delta/eta or beta/eta coupling theorem.",
            "next_honest_step": "A next admissible move must lift this finite partial-monoid identity witness to a strict nadsoliton identity-grounding theorem and derive delta/eta plus beta/eta coupling from it.  Without that lift, pivot to a genuinely new typed object outside the P2938 identity/positive-cone lane or preserve the P2929-P2944 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["identity_grounding_certificate"]
    lines = [
        "# P2944/S1894 partial-monoid identity-grounding theorem audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Identity-grounding certificate",
        f"- node count: `{cert['node_count']}`",
        f"- partial product count: `{cert['partial_product_count']}`",
        f"- two-sided total identity nodes: `{cert['two_sided_total_identity_nodes']}`",
        f"- zero carrier nodes: `{cert['zero_carrier_nodes']}`",
        f"- identity equals unique zero: `{cert['identity_equals_unique_zero']}`",
        f"- finite identity grounding verified: `{cert['finite_identity_grounding_verified']}`",
        f"- strict identity-grounding theorem exported: `{cert['strict_identity_grounding_theorem_exported']}`",
        f"- accepted strict source: `{cert['accepted_strict_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2943))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2944/S1894 partial-monoid identity-grounding theorem audit", "## P2944/S1894 partial-monoid identity-grounding theorem audit\n\n`P2944/S1894` audits the finite partial multiplication monoid on nodes `1..11` with retained products `d*e<=11`.  Node `1` is the unique two-sided total identity and the unique zero of the P2938 positive cone, so the finite algebraic identity-grounding witness demanded after P2943 is constructed.  This is not yet a strict nadsoliton identity-grounding theorem and exports no delta/eta source law, beta/eta coupling theorem, strict `L_p`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2944/S1894 partial-monoid identity `L_total` guard", "## P2944/S1894 partial-monoid identity `L_total` guard\n\n`P2944/S1894` verifies the finite partial-monoid identity witness for the P2938 positive cone, but no strict nadsoliton identity-grounding theorem or damping coupling is exported.  The witness cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current partial-monoid identity-grounding guardrail (P2944/S1894, 2026-06-19)", "## Current partial-monoid identity-grounding guardrail (P2944/S1894, 2026-06-19)\n\n- P2944 verifies the finite partial-monoid identity witness for the exact P2938 positive cone: node `1` is the unique two-sided total identity on retained products `d*e<=11` and the unique zero carrier node.\n- This is an algebraic identity-grounding witness only; it is not a strict nadsoliton identity-grounding theorem and does not export delta/eta or beta/eta coupling.\n- Do not promote P2944 to strict provenance, strict `L_p`, damping source, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must lift this witness to a strict nadsoliton identity-grounding theorem plus damping coupling, or pivot outside the P2938 identity/positive-cone lane; otherwise preserve the P2929-P2944 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
