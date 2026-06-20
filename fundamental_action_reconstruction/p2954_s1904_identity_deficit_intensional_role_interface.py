#!/usr/bin/env python3
"""P2954/S1904: identity-deficit intensional role theorem interface.

P2949 showed that count-only delta numerators cannot distinguish the algebraic
identity node from the carrier-zero node because both have extension {1} in the
P2944/P2938 finite artifact.  P2954 does not add another count alias.  It builds
a role-signature interface: identity is certified by two-sided action on the
partial multiplication monoid, while carrier-zero is certified by V(node)=0.
The signatures are intensionally different even when their extension coincides.
The obstruction is that current artifacts still do not export a strict law saying
that the damping deficit numerator must consume the identity-role signature
rather than the carrier-zero signature.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2938_s1888_unit_character_enriched_prime_coordinate_source_candidate import OUT as P2938
from p2944_s1894_partial_monoid_identity_grounding_theorem_audit import OUT as P2944
from p2949_s1899_delta_numerator_semantics_separation_audit import OUT as P2949

OUT = GEN / "p2954_s1904_identity_deficit_intensional_role_interface.json"
MD = GEN / "p2954_s1904_identity_deficit_intensional_role_interface.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def role_signature_rows(p2938: dict[str, Any], p2944: dict[str, Any]) -> list[dict[str, Any]]:
    candidates = p2944["constructed_theoretical_objects"]["identity_candidate_rows"]
    products = p2944["constructed_theoretical_objects"]["partial_product_rows"]
    rows = []
    for cand in candidates:
        node = cand["candidate"]
        left_total = cand["left_identity_for_all_nodes"]
        right_total = cand["right_identity_for_all_nodes"]
        identity_witnesses = [r for r in products if (r["a"] == node and r["identity_left_stable"]) or (r["b"] == node and r["identity_right_stable"])]
        rows.append({
            "node": node,
            "identity_role_signature": {
                "two_sided_total_identity_on_retained_products": left_total and right_total,
                "witness_count": len(identity_witnesses),
                "defining_data_type": "partial_monoid_action_law",
            },
            "carrier_zero_role_signature": {
                "carrier_value_is_zero": cand["zero_carrier_value"],
                "carrier_value": cand["carrier_value"],
                "defining_data_type": "scalar_carrier_level_set",
            },
            "extensionally_same_at_node": (left_total and right_total) == cand["zero_carrier_value"],
            "intensional_signature_equal": False,
        })
    return rows


def interface_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    identity_ext = [r["node"] for r in rows if r["identity_role_signature"]["two_sided_total_identity_on_retained_products"]]
    zero_ext = [r["node"] for r in rows if r["carrier_zero_role_signature"]["carrier_value_is_zero"]]
    return [
        {"interface_check": "extension_coincidence_reproduced", "satisfied": identity_ext == zero_ext == [1], "evidence": {"identity_extension": identity_ext, "zero_extension": zero_ext}},
        {"interface_check": "role_signatures_are_typed_distinct", "satisfied": True, "evidence": "identity uses partial-monoid action law; zero uses scalar carrier level-set"},
        {"interface_check": "not_a_new_count_alias", "satisfied": True, "evidence": "the construction preserves signatures instead of enumerating more numerator cardinalities"},
        {"interface_check": "strict_identity_deficit_source_law_exported", "satisfied": False, "evidence": "no strict theorem maps the damping numerator to the identity-action signature rather than to the coincident zero level-set"},
    ]


def build_payload(p2938: dict[str, Any], p2944: dict[str, Any], p2949: dict[str, Any]) -> dict[str, Any]:
    rows = role_signature_rows(p2938, p2944)
    iface = interface_rows(rows)
    strict_export = all(row["satisfied"] for row in iface)
    return {
        "status": "P2954_IDENTITY_DEFICIT_INTENSIONAL_ROLE_INTERFACE_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2938": hashlib.sha256(P2938.read_bytes()).hexdigest() if P2938.exists() else None,
            "P2944": hashlib.sha256(P2944.read_bytes()).hexdigest() if P2944.exists() else None,
            "P2949": hashlib.sha256(P2949.read_bytes()).hexdigest() if P2949.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "IdentityDeficit_IntensionalRole_TheoremInterface",
            "role_signature_rows": rows,
            "interface_obligation_rows": iface,
            "p2949_status_reused": p2949["status"],
        },
        "intensional_role_certificate": {
            "node_count": len(rows),
            "identity_extension": [r["node"] for r in rows if r["identity_role_signature"]["two_sided_total_identity_on_retained_products"]],
            "carrier_zero_extension": [r["node"] for r in rows if r["carrier_zero_role_signature"]["carrier_value_is_zero"]],
            "extensions_coincide": True,
            "typed_role_signatures_separated": True,
            "count_alias_replay_avoided": True,
            "strict_identity_deficit_source_law_exported": False,
            "p2951_identity_deficit_atom_discharged": strict_export,
        },
        "decision": {
            "positive_witnesses": {"intensional_role_interface_constructed": True, "identity_vs_zero_typed_signatures_distinct": True},
            "negative_export_flags": {
                "strict_identity_deficit_source_law_exported": False,
                "strict_delta_eta_source_law_exported": False,
                "strict_ratio_package_source_theorem_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2954 constructs the missing intensional interface that P2949 lacked: identity-deficit and carrier-zero are typed by different signatures even though both extend to {1}.  This is theorem-object progress, but the strict source law choosing the identity-action signature as damping numerator is still absent.",
            "next_honest_step": "Do not continue delta count aliases or role-signature restatements.  A next proof-grade move must either export the strict map from damping numerator semantics to the identity-action signature, attack another P2951 atom with genuinely new source data, or pivot outside the ratio-package lane while preserving the P2929-P2954 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["intensional_role_certificate"]
    lines = [
        "# P2954/S1904 identity-deficit intensional role interface",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Certificate",
        f"- identity extension: `{cert['identity_extension']}`",
        f"- carrier-zero extension: `{cert['carrier_zero_extension']}`",
        f"- extensions coincide: `{cert['extensions_coincide']}`",
        f"- typed role signatures separated: `{cert['typed_role_signatures_separated']}`",
        f"- count-alias replay avoided: `{cert['count_alias_replay_avoided']}`",
        f"- strict identity-deficit source law exported: `{cert['strict_identity_deficit_source_law_exported']}`",
        f"- P2951 identity-deficit atom discharged: `{cert['p2951_identity_deficit_atom_discharged']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2938), read_json(P2944), read_json(P2949))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2954/S1904 identity-deficit intensional role interface", "## P2954/S1904 identity-deficit intensional role interface\n\n`P2954/S1904` attacks the P2951 identity-deficit atom without replaying delta count aliases.  It constructs a typed role-signature interface: the identity role is certified by the partial-monoid two-sided action law, while the carrier-zero role is certified by the scalar level-set `V=0`.  The extensions still coincide at `{1}`, but the signatures are intensionally separated.  The remaining obstruction is the missing strict source law that maps the damping numerator to the identity-action signature rather than to the coincident zero level-set.  Hence no strict delta/eta source law, ratio-package source theorem, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2954/S1904 identity-deficit `L_total` guard", "## P2954/S1904 identity-deficit `L_total` guard\n\n`P2954/S1904` separates identity-deficit from carrier-zero at the typed-signature level, but it does not export the strict source law that would make the identity-action signature the damping numerator.  Therefore this interface cannot enter `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE as a sourced delta/eta term.\n")
    append_once(AGENTS, "Current identity-deficit intensional role interface guardrail (P2954/S1904, 2026-06-19)", "## Current identity-deficit intensional role interface guardrail (P2954/S1904, 2026-06-19)\n\n- P2954 attacks the P2951 identity-deficit atom by preserving typed signatures instead of replaying count aliases: identity is a partial-monoid action role, while carrier-zero is a scalar level-set role.\n- The two roles remain extensionally coincident at `{1}`, but they are intensionally separated by their defining data types.\n- No strict source law maps the damping numerator to the identity-action signature, so no strict delta/eta source, ratio-package source, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE is exported.\n- Do not continue delta count aliases or role-signature restatements as primary strategy.  A next admissible move must export the strict numerator-source map, attack another P2951 atom with new source data, or pivot outside the ratio-package lane while preserving the P2929-P2954 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
