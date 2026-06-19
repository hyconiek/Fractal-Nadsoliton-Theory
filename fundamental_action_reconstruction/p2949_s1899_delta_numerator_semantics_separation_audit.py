#!/usr/bin/env python3
"""P2949/S1899: delta numerator semantics separation audit.

P2948 identified one remaining theorem premise as a strict semantic choice of the
P2945 delta numerator.  P2949 attacks exactly that premise.  It does not add a
new ratio formula.  It asks whether the current finite artifacts distinguish the
algebraic identity role from the carrier-zero role strongly enough to select

    delta=(prime_count-identity_count)/prime_count

rather than a zero-deficit or alias numerator.

The result is a bounded theorem obstruction: in P2944 the two role predicates
have the same extension on nodes 1..11, namely {1}.  Therefore every count-only
or extension-only functional sees the identity role and zero role identically.
A strict delta numerator theorem would need new intensional/source semantics for
"identity deficit"; current artifacts only provide extensional equality.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2944_s1894_partial_monoid_identity_grounding_theorem_audit import OUT as P2944
from p2948_s1898_torsion_character_ratio_package_theorem_skeleton import OUT as P2948

OUT = GEN / "p2949_s1899_delta_numerator_semantics_separation_audit.json"
MD = GEN / "p2949_s1899_delta_numerator_semantics_separation_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_DELTA = Fraction(4, 5)


def frac_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "as_string": f"{value.numerator}/{value.denominator}",
        "as_float": float(value),
    }


def role_extension_rows(p2944: dict[str, Any]) -> list[dict[str, Any]]:
    rows = []
    for row in p2944["constructed_theoretical_objects"]["identity_candidate_rows"]:
        rows.append({
            "node": row["candidate"],
            "identity_role": row["two_sided_total_identity"],
            "zero_role": row["zero_carrier_value"],
            "roles_agree_on_node": row["two_sided_total_identity"] == row["zero_carrier_value"],
            "carrier_value": row["carrier_value"],
        })
    return rows


def numerator_functional_rows(prime_count: int, identity_nodes: list[int], zero_nodes: list[int]) -> list[dict[str, Any]]:
    identity_count = len(identity_nodes)
    zero_count = len(zero_nodes)
    candidates = [
        ("identity_deficit", Fraction(prime_count - identity_count, prime_count), "uses algebraic two-sided-total-identity role"),
        ("zero_deficit", Fraction(prime_count - zero_count, prime_count), "uses carrier-zero role"),
        ("intersection_deficit", Fraction(prime_count - len(set(identity_nodes) & set(zero_nodes)), prime_count), "uses identity∩zero extension"),
        ("union_deficit", Fraction(prime_count - len(set(identity_nodes) | set(zero_nodes)), prime_count), "uses identity∪zero extension"),
    ]
    return [
        {
            "functional": name,
            "delta": frac_payload(value),
            "matches_delta_4_5": value == TARGET_DELTA,
            "semantic_basis": basis,
            "strict_semantics_selected": False,
        }
        for name, value, basis in candidates
    ]


def separation_theorem_rows(role_rows: list[dict[str, Any]], numerator_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    identity_extension = [row["node"] for row in role_rows if row["identity_role"]]
    zero_extension = [row["node"] for row in role_rows if row["zero_role"]]
    return [
        {
            "theorem_obligation": "identity_and_zero_extensions_separated",
            "satisfied": identity_extension != zero_extension,
            "evidence": f"identity extension={identity_extension}; zero extension={zero_extension}",
        },
        {
            "theorem_obligation": "count_functionals_distinguish_identity_from_zero",
            "satisfied": len({row["delta"]["as_string"] for row in numerator_rows}) > 1,
            "evidence": f"delta outputs={sorted({row['delta']['as_string'] for row in numerator_rows})}",
        },
        {
            "theorem_obligation": "strict_intensional_identity_deficit_semantics_exported",
            "satisfied": False,
            "evidence": "current artifacts verify algebraic identity=carrier-zero extensionally but do not export a strict source law preferring identity-deficit semantics",
        },
        {
            "theorem_obligation": "p2948_delta_numerator_premise_discharged",
            "satisfied": False,
            "evidence": "P2948's finite package still needs a strict theorem choosing the delta numerator semantics",
        },
    ]


def acceptance_rows(separation_rows: list[dict[str, Any]], p2948: dict[str, Any]) -> list[dict[str, Any]]:
    p2948_delta_missing = p2948["package_certificate"]["strict_delta_numerator_semantics_exported"] is False
    return [
        {
            "criterion": "role_extension_audit_constructed",
            "satisfied": True,
            "evidence": "identity and zero predicates are compared node-by-node on P2944 rows",
        },
        {
            "criterion": "p2948_delta_gap_targeted",
            "satisfied": p2948_delta_missing,
            "evidence": "P2948 explicitly leaves strict delta numerator semantics unexported",
        },
        {
            "criterion": "identity_deficit_semantics_strictly_selected",
            "satisfied": False,
            "evidence": separation_rows[2]["evidence"],
        },
        {
            "criterion": "p2948_delta_numerator_premise_discharged",
            "satisfied": False,
            "evidence": separation_rows[3]["evidence"],
        },
    ]


def build_payload(p2944: dict[str, Any], p2948: dict[str, Any]) -> dict[str, Any]:
    cert = p2944["identity_grounding_certificate"]
    prime_count = p2948["constructed_theoretical_objects"]["finite_spine_rows"][2]["value"]["denominator"]
    role_rows = role_extension_rows(p2944)
    numerator_rows = numerator_functional_rows(prime_count, cert["two_sided_total_identity_nodes"], cert["zero_carrier_nodes"])
    theorem_rows = separation_theorem_rows(role_rows, numerator_rows)
    acceptance = acceptance_rows(theorem_rows, p2948)
    accepted = all(row["satisfied"] for row in acceptance)
    return {
        "status": "P2949_DELTA_NUMERATOR_SEMANTICS_SEPARATION_AUDIT_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2944": hashlib.sha256(P2944.read_bytes()).hexdigest() if P2944.exists() else None,
            "P2948": hashlib.sha256(P2948.read_bytes()).hexdigest() if P2948.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "IdentityDeficit_DeltaNumerator_Semantics_Separation_Audit",
            "role_extension_rows": role_rows,
            "numerator_functional_rows": numerator_rows,
            "separation_theorem_rows": theorem_rows,
            "acceptance_rows": acceptance,
        },
        "delta_numerator_semantics_certificate": {
            "identity_extension": cert["two_sided_total_identity_nodes"],
            "zero_extension": cert["zero_carrier_nodes"],
            "identity_zero_extensions_equal": cert["two_sided_total_identity_nodes"] == cert["zero_carrier_nodes"],
            "all_numerator_functionals_match_delta_4_5": all(row["matches_delta_4_5"] for row in numerator_rows),
            "strict_intensional_identity_deficit_semantics_exported": False,
            "p2948_delta_numerator_premise_discharged": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "role_extension_audit_constructed": True,
                "p2948_delta_semantics_gap_targeted": True,
                "finite_identity_zero_coincidence_verified": cert["identity_equals_unique_zero"],
            },
            "negative_export_flags": {
                "strict_delta_numerator_semantics_exported": False,
                "strict_delta_eta_source_law_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2949 targets exactly one P2948 skeleton premise: strict delta numerator semantics.  The current finite artifacts identify the algebraic identity node and the carrier-zero node extensionally as the same singleton {1}.  Consequently count-only numerator functionals using identity, zero, intersection, or union all return 4/5, and no strict intensional theorem selecting identity-deficit semantics is exported.",
            "next_honest_step": "Do not add more count aliases for delta.  A next proof-grade move must either export a genuine intensional/source theorem making identity-deficit the strict numerator, attack a different P2948 skeleton premise such as torsion-character provenance or beta/eta coupling, or pivot outside the ratio-package lane while preserving the P2929-P2949 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["delta_numerator_semantics_certificate"]
    lines = [
        "# P2949/S1899 delta numerator semantics separation audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Delta numerator semantics certificate",
        f"- identity extension: `{cert['identity_extension']}`",
        f"- zero extension: `{cert['zero_extension']}`",
        f"- identity/zero extensions equal: `{cert['identity_zero_extensions_equal']}`",
        f"- all numerator functionals match delta=4/5: `{cert['all_numerator_functionals_match_delta_4_5']}`",
        f"- strict intensional identity-deficit semantics exported: `{cert['strict_intensional_identity_deficit_semantics_exported']}`",
        f"- P2948 delta numerator premise discharged: `{cert['p2948_delta_numerator_premise_discharged']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2944), read_json(P2948))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2949/S1899 delta numerator semantics separation audit", "## P2949/S1899 delta numerator semantics separation audit\n\n`P2949/S1899` targets exactly one P2948 skeleton premise: strict semantics for the P2945 delta numerator.  The audit compares the P2944 algebraic identity role and carrier-zero role node-by-node and finds the same extension `{1}`.  Therefore identity-deficit, zero-deficit, intersection-deficit, and union-deficit count functionals all return `delta=4/5`; current artifacts do not export a strict intensional/source theorem selecting identity-deficit semantics.  No strict delta/eta source law, beta/eta coupling theorem, damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2949/S1899 delta numerator semantics `L_total` guard", "## P2949/S1899 delta numerator semantics `L_total` guard\n\n`P2949/S1899` shows that the P2948 delta numerator semantics premise is not discharged: identity and carrier-zero roles coincide extensionally on the current finite artifact, so count-only delta numerator choices are alias-equivalent.  Without an intensional strict source theorem and beta/eta coupling, this cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current delta numerator semantics separation guardrail (P2949/S1899, 2026-06-19)", "## Current delta numerator semantics separation guardrail (P2949/S1899, 2026-06-19)\n\n- P2949 attacks exactly one P2948 skeleton premise: strict semantic selection of the P2945 delta numerator.\n- The P2944 algebraic identity role and carrier-zero role have the same finite extension `{1}`, so identity-deficit, zero-deficit, intersection-deficit, and union-deficit count functionals all return `delta=4/5`.\n- This extensional coincidence does not export a strict intensional/source theorem preferring identity-deficit semantics, and it does not export strict delta/eta, beta/eta coupling, damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- Do not continue count-alias variants as primary strategy.  A next admissible move must export an intensional identity-deficit source theorem, attack another P2948 premise such as torsion-character provenance or beta/eta coupling, or pivot outside this ratio-package lane while preserving the P2929-P2949 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
