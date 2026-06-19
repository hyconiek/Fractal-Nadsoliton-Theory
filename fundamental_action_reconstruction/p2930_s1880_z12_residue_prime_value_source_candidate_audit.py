#!/usr/bin/env python3
"""P2930/S1880: Z12 residue prime-value source candidate audit.

P2929 says the next admissible move must supply one genuinely new strict typed
object and run the intake gate.  P2930 supplies exactly one bounded candidate
for one missing P2927/P2929 obligation: a possible strict source for the five
prime-log values L_2,L_3,L_5,L_7,L_11.  The candidate uses only the Z12 residue
circle to assign a finite cyclic-distance value to each audited prime.

The result is intentionally conservative.  The candidate computes nonzero
finite labels and extends to the audited multiplicative carrier, but it fails as
a strict source law because the value scale is conventional, it is not an
exported nadsoliton source theorem, and it does not supply the delta/eta source
or beta/eta coupling theorem required by P2927/P2928.
"""
from __future__ import annotations

import hashlib
import json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix import NODES, PRIMES, factor_vector

P2929 = GEN / "p2929_s1879_post_damping_state_map_no_new_live_frontier_certificate.json"
OUT = GEN / "p2930_s1880_z12_residue_prime_value_source_candidate_audit.json"
MD = GEN / "p2930_s1880_z12_residue_prime_value_source_candidate_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def z12_distance(residue: int) -> int:
    r = residue % 12
    return min(r, 12 - r)


def candidate_prime_values() -> dict[str, int]:
    return {f"L_{p}": z12_distance(p) for p in PRIMES}


def additive_extension_rows(values: dict[str, int]) -> list[dict[str, Any]]:
    rows = []
    for d in NODES:
        vector = factor_vector(d)
        y_d = sum(exp * values[f"L_{p}"] for p, exp in zip(PRIMES, vector))
        rows.append({
            "node": d,
            "prime_exponent_vector": vector,
            "candidate_y_d": y_d,
            "formula": "sum_p v_p(d)*cyclic_distance_Z12(p,0)",
        })
    return rows


def product_defect_rows(values: dict[str, int]) -> list[dict[str, Any]]:
    node_values = {row["node"]: row["candidate_y_d"] for row in additive_extension_rows(values)}
    rows = []
    for d in NODES:
        for e in NODES:
            product = d * e
            if product <= 11:
                defect = node_values[product] - node_values[d] - node_values[e]
                rows.append({"d": d, "e": e, "product": product, "additive_defect": defect})
    return rows


def aut_sensitivity_rows(values: dict[str, int]) -> list[dict[str, Any]]:
    units = [1, 5, 7, 11]
    rows = []
    for unit in units:
        image_values = {f"L_{p}": z12_distance(unit * p) for p in PRIMES}
        rows.append({
            "aut_unit": unit,
            "image_prime_residue_values": image_values,
            "preserves_candidate_prime_values": image_values == values,
        })
    return rows


def acceptance_matrix(values: dict[str, int]) -> dict[str, Any]:
    product_rows = product_defect_rows(values)
    aut_rows = aut_sensitivity_rows(values)
    return {
        "candidate_name": "Z12_Residue_Prime_Value_Source_Candidate",
        "targets_one_missing_obligation": "strict prime-log value source L_2,L_3,L_5,L_7,L_11",
        "computes_five_nonzero_prime_labels": all(v != 0 for v in values.values()),
        "audited_product_pair_count": len(product_rows),
        "formal_additive_defect_count": sum(1 for row in product_rows if row["additive_defect"] != 0),
        "aut_unit_rows_preserving_candidate_count": sum(1 for row in aut_rows if row["preserves_candidate_prime_values"]),
        "strict_nadsoliton_source_theorem_exported": False,
        "intrinsic_nonconventional_scale_exported": False,
        "delta_eta_source_exported": False,
        "beta_eta_coupling_theorem_exported": False,
        "accepted_as_strict_prime_log_value_source": False,
    }


def build_payload(p2929: dict[str, Any]) -> dict[str, Any]:
    values = candidate_prime_values()
    matrix = acceptance_matrix(values)
    return {
        "status": "P2930_Z12_RESIDUE_PRIME_VALUE_SOURCE_CANDIDATE_AUDIT_REJECTED_AS_STRICT_SOURCE",
        "input_hashes": {"P2929": hashlib.sha256(P2929.read_bytes()).hexdigest() if P2929.exists() else None},
        "constructed_theoretical_objects": {
            "candidate": {
                "name": "Z12_Residue_Prime_Value_Source_Candidate",
                "definition": "L_p := min(p mod 12, 12-(p mod 12)) for p in {2,3,5,7,11}",
                "prime_values": values,
                "intake_gate_result": "genuinely new bounded candidate, finite test available, but not accepted as strict source",
            },
            "additive_extension_rows": additive_extension_rows(values),
            "product_defect_rows": product_defect_rows(values),
            "aut_sensitivity_rows": aut_sensitivity_rows(values),
        },
        "acceptance_matrix": matrix,
        "decision": {
            "positive_witnesses": {
                "new_typed_candidate_supplied_to_p2929_intake_gate": True,
                "five_nonzero_prime_labels_computed": matrix["computes_five_nonzero_prime_labels"],
                "finite_additive_carrier_has_zero_product_defects": matrix["formal_additive_defect_count"] == 0,
            },
            "negative_export_flags": {
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
            "reason": "P2930 supplies one genuinely new finite candidate for the missing L_p value-source obligation.  It computes nonzero Z12 residue-distance labels and its multiplicative extension has zero formal product defects on the audited node set.  However, the labels are conventional residue-distance numbers, not an exported strict nadsoliton source theorem with intrinsic scale, and they do not provide delta/eta or beta/eta coupling.  The candidate is therefore rejected as a strict prime-log value source.",
            "next_honest_step": "Either provide an explicit strict nadsoliton source theorem that turns the Z12 residue labels into intrinsically scaled L_p values and couples them to the P2927/P2928 damping packet, or pivot to a different genuinely new typed object.  Without that, preserve the P2929/P2930 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    matrix = payload["acceptance_matrix"]
    lines = [
        "# P2930/S1880 Z12 residue prime-value source candidate audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Candidate",
        f"- object: `{payload['constructed_theoretical_objects']['candidate']['name']}`",
        f"- prime values: `{payload['constructed_theoretical_objects']['candidate']['prime_values']}`",
        "",
        "## Finite certificate",
        f"- audited product pairs: `{matrix['audited_product_pair_count']}`",
        f"- formal additive defects: `{matrix['formal_additive_defect_count']}`",
        f"- accepted as strict value source: `{matrix['accepted_as_strict_prime_log_value_source']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2929))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2930/S1880 Z12 residue prime-value source candidate audit", "## P2930/S1880 Z12 residue prime-value source candidate audit\n\n`P2930/S1880` runs the P2929 fresh-object intake gate on one new bounded candidate for the missing strict prime-log value source: `Z12_Residue_Prime_Value_Source_Candidate`, with `L_p := min(p mod 12, 12-(p mod 12))` for `p in {2,3,5,7,11}`.  The candidate computes five nonzero labels and its additive extension has zero formal product defects on the `29` audited products `d*e<=11`.  It is rejected as a strict value source because the residue-distance scale is conventional, no strict nadsoliton source theorem is exported, and no delta/eta or beta/eta coupling theorem is supplied.  No strict damping source packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2930/S1880 Z12 residue prime-value candidate `L_total` guard", "## P2930/S1880 Z12 residue prime-value candidate `L_total` guard\n\n`P2930/S1880` provides a real finite candidate for the missing `L_p` source obligation, but rejects it as role-bearing input: Z12 residue distances are finite labels with zero formal carrier defects, not an intrinsic strict source theorem with scale and damping coupling.  Therefore this candidate cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current Z12 residue prime-value source candidate guardrail (P2930/S1880, 2026-06-19)", "## Current Z12 residue prime-value source candidate guardrail (P2930/S1880, 2026-06-19)\n\n- P2930 supplies one genuinely new bounded candidate for the P2927/P2929 missing strict `L_p` value-source obligation: `Z12_Residue_Prime_Value_Source_Candidate`.\n- The candidate computes nonzero finite residue-distance labels and has zero formal additive product defects on the audited products, so it is a real carrier witness.\n- It is not accepted as a strict source theorem: the value scale is conventional, no strict nadsoliton source theorem is exported, and no `delta/eta` source or beta/eta coupling theorem is supplied.\n- Do not promote P2930 to strict damping `beta/eta`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE; a next move must add an intrinsic strict source theorem/coupling for these labels or pivot to a different new typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
