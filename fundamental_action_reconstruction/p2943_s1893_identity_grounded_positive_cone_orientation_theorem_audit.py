#!/usr/bin/env python3
"""P2943/S1893: identity-grounded positive-cone orientation theorem audit.

P2942 isolated the missing sign atom for the P2941 extremal selector as a
positive value-order orientation `a>0`.  P2943 does not rescan affine orders.
It constructs the next theorem-object candidate: an identity-grounded positive
cone for the exact P2938 carrier.

The finite theorem side is real: the P2938 vector has f(1)=0, all prime
coordinates are positive, every non-identity audited node has positive value,
and product additivity has no defects.  These algebraic facts conditionally
orient the value order away from the identity and therefore justify the positive
orientation sign within the finite carrier algebra.

The strict side is still not exported: current artifacts do not prove that this
identity-grounded positive cone is a strict nadsoliton source, nor do they couple
it to delta/eta and beta/eta.  Thus P2943 is an algebraic orientation theorem
candidate, not strict provenance or a damping source.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate import PRIMES, product_rows
from p2940_s1890_p2938_carrier_aut_orbit_selector_burden import NODES, PRIME_VECTOR, carrier_value
from p2942_s1892_extremal_selector_value_order_orientation_source_gate import OUT as P2942

OUT = GEN / "p2943_s1893_identity_grounded_positive_cone_orientation_theorem_audit.json"
MD = GEN / "p2943_s1893_identity_grounded_positive_cone_orientation_theorem_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def node_value_rows() -> list[dict[str, Any]]:
    return [
        {
            "node": node,
            "carrier_value": carrier_value(node),
            "is_identity": node == 1,
            "satisfies_positive_cone": carrier_value(node) == 0 if node == 1 else carrier_value(node) > 0,
        }
        for node in NODES
    ]


def theorem_premise_rows() -> list[dict[str, Any]]:
    rows = node_value_rows()
    node_values = {node: carrier_value(node) for node in NODES}
    products = product_rows(node_values)
    product_defects = [row for row in products if row["additive_defect"] != 0]
    return [
        {
            "premise": "identity_grounded_zero",
            "satisfied": carrier_value(1) == 0,
            "evidence": f"f(1)={carrier_value(1)}",
        },
        {
            "premise": "all_prime_coordinates_positive",
            "satisfied": all(coord > 0 for coord in PRIME_VECTOR),
            "evidence": f"prime vector {PRIME_VECTOR} for primes {PRIMES}",
        },
        {
            "premise": "all_nonidentity_nodes_positive",
            "satisfied": all(row["carrier_value"] > 0 for row in rows if not row["is_identity"]),
            "evidence": "every node 2..11 has positive P2938 carrier value",
        },
        {
            "premise": "product_additivity_zero_defects",
            "satisfied": not product_defects,
            "evidence": f"{len(products)} product rows, {len(product_defects)} defects",
        },
        {
            "premise": "positive_orientation_sign_conditionally_fixed",
            "satisfied": True,
            "evidence": "identity zero plus positive nonidentity cone orients increasing carrier value away from identity, matching the P2942 a>0 branch",
        },
    ]


def strict_acceptance_rows(algebraic_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "criterion": "finite_positive_cone_theorem_candidate_verified",
            "satisfied": all(row["satisfied"] for row in algebraic_rows),
            "evidence": "identity zero, positive primes, positive nonidentity nodes, and product additivity all hold",
        },
        {
            "criterion": "strict_nadsoliton_identity_grounding_theorem_exported",
            "satisfied": False,
            "evidence": "current artifacts do not export a strict theorem identifying this positive cone as nadsoliton-sourced rather than algebraic carrier structure",
        },
        {
            "criterion": "nonconventional_positive_orientation_source_exported",
            "satisfied": False,
            "evidence": "the positive orientation is conditionally fixed inside the carrier algebra, but its strict source status remains unproved",
        },
        {
            "criterion": "delta_eta_source_law_exported",
            "satisfied": False,
            "evidence": "no delta/eta source law is derived from the positive cone",
        },
        {
            "criterion": "beta_eta_coupling_theorem_exported",
            "satisfied": False,
            "evidence": "no beta/eta coupling theorem is derived from the positive cone",
        },
    ]


def build_payload(_: dict[str, Any]) -> dict[str, Any]:
    node_rows = node_value_rows()
    algebraic = theorem_premise_rows()
    strict_rows = strict_acceptance_rows(algebraic)
    finite_verified = strict_rows[0]["satisfied"]
    accepted = all(row["satisfied"] for row in strict_rows)
    return {
        "status": "P2943_IDENTITY_GROUNDED_POSITIVE_CONE_ORIENTATION_THEOREM_CANDIDATE_NO_STRICT_SOURCE",
        "input_hashes": {"P2942": hashlib.sha256(P2942.read_bytes()).hexdigest() if P2942.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "Identity_Grounded_Positive_Cone_Orientation_Theorem_For_P2938_Carrier",
            "prime_vector_order_2_3_5_7_11": PRIME_VECTOR,
            "node_value_rows": node_rows,
            "algebraic_theorem_premise_rows": algebraic,
            "strict_acceptance_rows": strict_rows,
        },
        "positive_cone_certificate": {
            "node_count": len(node_rows),
            "identity_value": carrier_value(1),
            "positive_prime_coordinate_count": sum(1 for coord in PRIME_VECTOR if coord > 0),
            "positive_nonidentity_node_count": sum(1 for row in node_rows if not row["is_identity"] and row["carrier_value"] > 0),
            "product_additivity_defect_count": sum(1 for row in product_rows({node: carrier_value(node) for node in NODES}) if row["additive_defect"] != 0),
            "finite_positive_cone_theorem_candidate_verified": finite_verified,
            "strict_identity_grounding_theorem_exported": False,
            "strict_provenance_theorem_exported": False,
            "accepted_strict_source": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "algebraic_positive_orientation_sign_fixed": finite_verified,
                "P2942_a_greater_than_zero_branch_conditionally_sourced_inside_carrier_algebra": finite_verified,
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
            "reason": "P2943 proves the finite carrier-algebra side of the P2942 positive orientation sign: f(1)=0, all prime coordinates are positive, all nonidentity nodes are positive, and product additivity has zero defects.  This conditionally orients value growth away from the identity and selects the a>0 branch.  It still does not export a strict nadsoliton identity-grounding theorem or delta/eta plus beta/eta coupling.",
            "next_honest_step": "A next admissible move must prove strict nadsoliton identity-grounding for this positive cone and then derive delta/eta plus beta/eta coupling from it.  If that cannot be exported, pivot to a genuinely new typed object outside the P2938 positive-cone/polarity lane or preserve the P2929-P2943 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["positive_cone_certificate"]
    lines = [
        "# P2943/S1893 identity-grounded positive-cone orientation theorem audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Positive-cone certificate",
        f"- node count: `{cert['node_count']}`",
        f"- identity value: `{cert['identity_value']}`",
        f"- positive prime coordinates: `{cert['positive_prime_coordinate_count']}`",
        f"- positive nonidentity nodes: `{cert['positive_nonidentity_node_count']}`",
        f"- product-additivity defects: `{cert['product_additivity_defect_count']}`",
        f"- finite positive-cone theorem candidate verified: `{cert['finite_positive_cone_theorem_candidate_verified']}`",
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
    payload = build_payload(read_json(P2942))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2943/S1893 identity-grounded positive-cone orientation theorem audit", "## P2943/S1893 identity-grounded positive-cone orientation theorem audit\n\n`P2943/S1893` constructs an identity-grounded positive-cone theorem candidate for the exact P2938 carrier.  The finite carrier algebra verifies `f(1)=0`, all five prime coordinates are positive, all nonidentity audited nodes have positive value, and product additivity has zero defects; this conditionally sources the P2942 positive-orientation branch `a>0` inside the carrier algebra.  The strict side is still not exported: no nadsoliton identity-grounding theorem, delta/eta source law, beta/eta coupling theorem, strict `L_p`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2943/S1893 positive-cone `L_total` guard", "## P2943/S1893 positive-cone `L_total` guard\n\n`P2943/S1893` verifies a finite identity-grounded positive cone for the P2938 carrier and thereby conditionally fixes the P2942 `a>0` orientation within carrier algebra.  Because no strict nadsoliton identity-grounding theorem or damping coupling is exported, the positive cone cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current identity-grounded positive-cone orientation guardrail (P2943/S1893, 2026-06-19)", "## Current identity-grounded positive-cone orientation guardrail (P2943/S1893, 2026-06-19)\n\n- P2943 verifies the finite carrier-algebra positive cone for the exact P2938 carrier: `f(1)=0`, all prime coordinates are positive, all nonidentity audited nodes have positive values, and product additivity has zero defects.\n- This conditionally fixes the P2942 positive orientation branch `a>0` inside the carrier algebra, but it is not yet a strict nadsoliton identity-grounding theorem.\n- Do not promote P2943 to strict provenance, strict `L_p`, delta/eta source, beta/eta coupling, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- The next admissible move must export strict nadsoliton identity-grounding for this positive cone plus damping coupling, or pivot outside the P2938 positive-cone/polarity lane; otherwise preserve the P2929-P2943 no-new-live-frontier boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
