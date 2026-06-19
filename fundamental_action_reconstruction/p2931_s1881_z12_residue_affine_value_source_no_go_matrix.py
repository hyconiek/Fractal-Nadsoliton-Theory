#!/usr/bin/env python3
"""P2931/S1881: Z12 residue affine value-source no-go matrix.

P2930 supplied one finite Z12 residue-distance candidate for the missing strict
prime-log value source.  P2931 strengthens that follow-up into a proof-grade
obstruction: instead of testing only the literal distance labels, it audits the
whole one-parameter affine cyclic-distance source family

    F_a(d) = a * min(d mod 12, 12-(d mod 12)), with F_a(1)=0.

A prime-log value source cannot be only a per-prime label: it must extend to a
multiplicative/additive character on audited nodes.  P2931 therefore compares
F_a(d*e) with F_a(d)+F_a(e) for every product d*e<=11.  The symbolic defect row
for 2*3=6 is already a, so all nonzero scales fail; the only additive member is
the zero source, which cannot provide the required nonzero prime-log values.
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
from p2928_s1878_beta_eta_coupling_theorem_obstruction_matrix import NODES, PRIMES
from p2930_s1880_z12_residue_prime_value_source_candidate_audit import z12_distance

P2930 = GEN / "p2930_s1880_z12_residue_prime_value_source_candidate_audit.json"
OUT = GEN / "p2931_s1881_z12_residue_affine_value_source_no_go_matrix.json"
MD = GEN / "p2931_s1881_z12_residue_affine_value_source_no_go_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def product_pairs() -> list[tuple[int, int, int]]:
    return [(d, e, d * e) for d in NODES for e in NODES if d * e <= 11]


def symbolic_defect_rows() -> list[dict[str, Any]]:
    rows = []
    for d, e, product in product_pairs():
        coefficient = z12_distance(product) - z12_distance(d) - z12_distance(e)
        rows.append({
            "d": d,
            "e": e,
            "product": product,
            "defect_formula": f"{coefficient}*a",
            "defect_coefficient": coefficient,
            "is_obstruction_row_for_nonzero_a": coefficient != 0,
        })
    return rows


def solve_affine_scale_constraints(rows: list[dict[str, Any]]) -> dict[str, Any]:
    nonzero_coefficients = sorted({abs(row["defect_coefficient"]) for row in rows if row["defect_coefficient"] != 0})
    symbolic_solution = "a=0" if nonzero_coefficients else "a free"
    return {
        "equation_family": "(dist_Z12(d*e)-dist_Z12(d)-dist_Z12(e))*a = 0 for all audited d*e<=11",
        "audited_product_pair_count": len(rows),
        "nonzero_symbolic_defect_row_count": sum(1 for row in rows if row["defect_coefficient"] != 0),
        "nonzero_defect_coefficients_abs": nonzero_coefficients,
        "minimal_witness_row": next(row for row in rows if row["d"] == 2 and row["e"] == 3),
        "symbolic_solution_space": symbolic_solution,
        "solution_dimension": 0 if symbolic_solution == "a=0" else 1,
    }


def finite_scale_scan(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    scan = []
    for numerator in range(-12, 13):
        a = Fraction(numerator, 12)
        defects = [a * row["defect_coefficient"] for row in rows]
        additive = all(defect == 0 for defect in defects)
        prime_values = {f"L_{p}": a * z12_distance(p) for p in PRIMES}
        scan.append({
            "a": str(a),
            "all_product_defects_zero": additive,
            "all_prime_values_nonzero": all(value != 0 for value in prime_values.values()),
            "accepted_strict_value_source_candidate": False,
            "failure": "zero source" if additive else "multiplicativity/additivity defect",
        })
    return scan


def build_payload(p2930: dict[str, Any]) -> dict[str, Any]:
    rows = symbolic_defect_rows()
    solution = solve_affine_scale_constraints(rows)
    scan = finite_scale_scan(rows)
    additive_rows = [row for row in scan if row["all_product_defects_zero"]]
    accepted_rows = [row for row in scan if row["accepted_strict_value_source_candidate"]]
    return {
        "status": "P2931_Z12_RESIDUE_AFFINE_VALUE_SOURCE_NO_GO_MATRIX_ZERO_ONLY",
        "input_hashes": {"P2930": hashlib.sha256(P2930.read_bytes()).hexdigest() if P2930.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_family": {
                "name": "Z12_Residue_Affine_Prime_Value_Source_Family",
                "definition": "F_a(d)=a*min(d mod 12, 12-(d mod 12)) with F_a(1)=0",
                "target": "strict source values for L_2,L_3,L_5,L_7,L_11 compatible with additive product character",
            },
            "symbolic_product_defect_rows": rows,
            "affine_scale_solution_certificate": solution,
            "finite_scale_scan": scan,
        },
        "no_go_certificate": {
            "audited_product_pair_count": solution["audited_product_pair_count"],
            "nonzero_symbolic_defect_row_count": solution["nonzero_symbolic_defect_row_count"],
            "symbolic_solution_space": solution["symbolic_solution_space"],
            "finite_scan_row_count": len(scan),
            "finite_scan_additive_row_count": len(additive_rows),
            "finite_scan_accepted_row_count": len(accepted_rows),
            "zero_source_is_only_additive_member": len(additive_rows) == 1 and additive_rows[0]["a"] == "0",
            "nonzero_strict_prime_log_value_source_exported": False,
        },
        "decision": {
            "positive_witnesses": {
                "p2930_candidate_family_generalized": True,
                "symbolic_defect_system_constructed": True,
                "finite_scale_scan_executed": True,
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
            "reason": "P2931 strengthens P2930 by auditing the full affine Z12 residue-distance source family.  The symbolic product equation for 2*3=6 already forces a=0, and the finite scan confirms that the only additive member is the zero source.  Therefore no nonzero strict prime-log value source is exported from this residue-affine family.",
            "next_honest_step": "Do not continue residue-distance variants as a primary strategy.  A next admissible move must either provide a non-affine strict nadsoliton value-source theorem with its own finite additivity/coupling test, or pivot to a different genuinely new typed object outside the residue-distance family.  Otherwise preserve the P2929-P2931 no-new-live-frontier boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["no_go_certificate"]
    witness = payload["constructed_theoretical_objects"]["affine_scale_solution_certificate"]["minimal_witness_row"]
    lines = [
        "# P2931/S1881 Z12 residue affine value-source no-go matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Symbolic no-go certificate",
        f"- audited product pairs: `{cert['audited_product_pair_count']}`",
        f"- nonzero symbolic defect rows: `{cert['nonzero_symbolic_defect_row_count']}`",
        f"- minimal witness row: `{witness['d']}*{witness['e']}={witness['product']}` with defect `{witness['defect_formula']}`",
        f"- symbolic solution space: `{cert['symbolic_solution_space']}`",
        f"- zero source is only additive member: `{cert['zero_source_is_only_additive_member']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2930))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2931/S1881 Z12 residue affine value-source no-go matrix", "## P2931/S1881 Z12 residue affine value-source no-go matrix\n\n`P2931/S1881` strengthens P2930 by auditing the full affine residue-distance family `F_a(d)=a*min(d mod 12,12-(d mod 12))` as a possible source for `L_2,L_3,L_5,L_7,L_11`.  The symbolic additive-product equations on the `29` audited products contain the witness `2*3=6` with defect `a`, so all nonzero scales fail and the only additive member is the zero source.  The finite scan over `25` rational scales confirms zero accepted rows.  No nonzero strict prime-log value source, damping source packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2931/S1881 residue-affine value-source `L_total` guard", "## P2931/S1881 residue-affine value-source `L_total` guard\n\n`P2931/S1881` blocks the residue-affine continuation of P2930: multiplicativity/additivity forces the affine residue-distance scale to `a=0`, while strict damping needs nonzero prime-log values plus delta/eta and coupling sources.  Therefore this family cannot enter role-bearing nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current residue-affine prime-value source no-go guardrail (P2931/S1881, 2026-06-19)", "## Current residue-affine prime-value source no-go guardrail (P2931/S1881, 2026-06-19)\n\n- P2931 generalizes the P2930 residue-distance candidate to the affine family `F_a(d)=a*dist_Z12(d,0)` and tests additive compatibility on all audited products `d*e<=11`.\n- The symbolic witness `2*3=6` has defect `a`, so multiplicativity/additivity forces `a=0`; the only additive member is the zero source and cannot supply nonzero `L_p` values.\n- Do not continue residue-distance variants as primary strict `L_p` source candidates without a genuinely non-affine strict nadsoliton source theorem and a new finite acceptance test.\n- No strict damping `beta/eta`, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure is exported from P2931.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
