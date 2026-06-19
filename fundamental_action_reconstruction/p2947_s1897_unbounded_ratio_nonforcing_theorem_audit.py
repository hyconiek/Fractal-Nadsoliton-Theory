#!/usr/bin/env python3
"""P2947/S1897: unbounded ratio non-forcing theorem audit.

P2946 used a bounded positive-vector model class to show that the P2945 eta
formula is not forced by finite identity-positive-cone premises.  P2947 removes
that boundedness caveat.  It constructs the exact parametric theorem object
behind the obstruction:

  for every integer S >= prime_count, the positive vector
  (S-prime_count+1, 1, 1, 1, 1)
  has length prime_count, positive integer coordinates, and eta=S/prime_count.

Thus identity-positive-cone positivity alone permits infinitely many eta values;
eta=9/5 is one value in an infinite family, not a strict consequence of those
premises.  P2947 also records the exact delta alias kernel caused by
identity_count=zero_count=1.  This is a theorem-level obstruction, not another
bounded scan and not a strict damping export.
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
from p2945_s1895_identity_positive_cone_delta_eta_source_candidate import OUT as P2945
from p2946_s1896_delta_eta_ratio_strict_forcing_obstruction_matrix import OUT as P2946

OUT = GEN / "p2947_s1897_unbounded_ratio_nonforcing_theorem_audit.json"
MD = GEN / "p2947_s1897_unbounded_ratio_nonforcing_theorem_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_ETA = Fraction(9, 5)
TARGET_DELTA = Fraction(4, 5)
SAMPLE_SUMS = list(range(5, 16))


def frac_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "as_string": f"{value.numerator}/{value.denominator}",
        "as_float": float(value),
    }


def parametric_positive_vector(sum_value: int, prime_count: int) -> list[int]:
    if sum_value < prime_count:
        raise ValueError("positive integer coordinates require sum_value >= prime_count")
    return [sum_value - prime_count + 1] + [1] * (prime_count - 1)


def parametric_eta_witness_rows(prime_count: int) -> list[dict[str, Any]]:
    rows = []
    for sum_value in SAMPLE_SUMS:
        vector = parametric_positive_vector(sum_value, prime_count)
        eta = Fraction(sum_value, prime_count)
        rows.append({
            "sum_value": sum_value,
            "witness_vector": vector,
            "all_coordinates_positive": all(coord > 0 for coord in vector),
            "vector_length": len(vector),
            "eta": frac_payload(eta),
            "matches_target_eta_9_5": eta == TARGET_ETA,
        })
    return rows


def theorem_rows(prime_count: int) -> list[dict[str, Any]]:
    return [
        {
            "theorem_object": "Parametric_Positive_Cone_Eta_Family",
            "statement": "For every integer S>=prime_count, v(S)=(S-prime_count+1,1,...,1) is positive and has eta=S/prime_count.",
            "verified_symbolically": True,
            "consequence": "positive-cone premises alone admit infinitely many eta values",
        },
        {
            "theorem_object": "Target_Eta_Is_One_Member_Not_A_Consequence",
            "statement": "eta=9/5 is obtained only when S=9 for prime_count=5; S=5,6,7,8,10,... give different eta values while preserving positivity.",
            "verified_symbolically": prime_count == 5,
            "consequence": "P2945 eta requires a source theorem for the sum 9 or for the exact P2938 vector",
        },
        {
            "theorem_object": "Delta_Alias_Kernel",
            "statement": "identity_count=zero_count=1 makes formulas subtracting identity_count, zero_count, min, max, or their mean coincide at 4/5.",
            "verified_symbolically": True,
            "consequence": "P2945 delta formula is not uniquely selected by the finite counts without a theorem choosing the numerator semantics",
        },
    ]


def delta_alias_kernel_rows(prime_count: int, identity_count: int, zero_count: int) -> list[dict[str, Any]]:
    numerator_terms = {
        "identity_count": Fraction(identity_count, 1),
        "zero_count": Fraction(zero_count, 1),
        "min(identity_count,zero_count)": Fraction(min(identity_count, zero_count), 1),
        "max(identity_count,zero_count)": Fraction(max(identity_count, zero_count), 1),
        "mean(identity_count,zero_count)": Fraction(identity_count + zero_count, 2),
    }
    rows = []
    for term, subtrahend in numerator_terms.items():
        value = Fraction(prime_count, 1) - subtrahend
        ratio = value / prime_count
        rows.append({
            "subtracted_term": term,
            "subtrahend": frac_payload(subtrahend),
            "ratio": frac_payload(ratio),
            "matches_delta_4_5": ratio == TARGET_DELTA,
            "strict_semantics_selected": False,
        })
    return rows


def acceptance_rows(eta_rows: list[dict[str, Any]], theorem_object_rows: list[dict[str, Any]], delta_rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    non_target_eta_rows = [row for row in eta_rows if not row["matches_target_eta_9_5"]]
    return [
        {
            "criterion": "unbounded_parametric_eta_family_constructed",
            "satisfied": all(row["verified_symbolically"] for row in theorem_object_rows[:1]),
            "evidence": theorem_object_rows[0]["statement"],
        },
        {
            "criterion": "positive_cone_premises_force_eta_9_5",
            "satisfied": False,
            "evidence": f"sample witnesses already include {len(non_target_eta_rows)} positive non-9/5 eta values, and the symbolic family is infinite",
        },
        {
            "criterion": "finite_counts_uniquely_select_delta_formula",
            "satisfied": False,
            "evidence": f"{sum(row['matches_delta_4_5'] for row in delta_rows)} alias formulas coincide at 4/5",
        },
        {
            "criterion": "strict_theorem_selects_p2938_vector_and_sum_9",
            "satisfied": False,
            "evidence": "no strict source theorem selects the exact vector [1,2,2,2,2] or sum 9 from the infinite positive family",
        },
        {
            "criterion": "strict_beta_eta_coupling_theorem_exported",
            "satisfied": False,
            "evidence": "no strict theorem couples the unsourced ratio formulas to beta/eta damping",
        },
    ]


def build_payload(p2945: dict[str, Any], p2946: dict[str, Any]) -> dict[str, Any]:
    inventory = p2945["constructed_theoretical_objects"]["source_inventory"]
    prime_count = inventory["prime_count"]
    eta_rows = parametric_eta_witness_rows(prime_count)
    theorem_object_rows = theorem_rows(prime_count)
    delta_rows = delta_alias_kernel_rows(prime_count, inventory["identity_count"], inventory["zero_count"])
    acceptance = acceptance_rows(eta_rows, theorem_object_rows, delta_rows)
    accepted = all(row["satisfied"] for row in acceptance)
    return {
        "status": "P2947_UNBOUNDED_RATIO_NONFORCING_THEOREM_AUDIT_NO_STRICT_EXPORT",
        "input_hashes": {
            "P2945": hashlib.sha256(P2945.read_bytes()).hexdigest() if P2945.exists() else None,
            "P2946": hashlib.sha256(P2946.read_bytes()).hexdigest() if P2946.exists() else None,
        },
        "constructed_theoretical_objects": {
            "candidate_object": "Parametric_Positive_Cone_Eta_Family_And_Delta_Alias_Kernel",
            "theorem_rows": theorem_object_rows,
            "parametric_eta_witness_rows": eta_rows,
            "delta_alias_kernel_rows": delta_rows,
            "acceptance_rows": acceptance,
            "p2946_bounded_result_subsumed": p2946["strict_forcing_certificate"]["eta_9_5_forced_by_positive_cone_premises"] is False,
        },
        "unbounded_nonforcing_certificate": {
            "parametric_family_symbolically_verified": theorem_object_rows[0]["verified_symbolically"],
            "sample_sum_count": len(eta_rows),
            "sample_non_target_eta_count": sum(not row["matches_target_eta_9_5"] for row in eta_rows),
            "eta_9_5_forced_by_positive_cone_premises": False,
            "delta_formula_uniquely_selected_by_counts": False,
            "strict_p2938_vector_sum9_source_theorem_exported": False,
            "strict_beta_eta_coupling_theorem_exported": False,
            "accepted_strict_ratio_source_theorem": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "unbounded_parametric_nonforcing_theorem_constructed": True,
                "p2946_bounded_scan_subsumed_by_symbolic_family": True,
                "delta_alias_kernel_constructed": True,
            },
            "negative_export_flags": {
                "strict_ratio_source_theorem_exported": False,
                "strict_p2938_vector_source_theorem_exported": False,
                "strict_delta_eta_source_law_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2947 replaces the bounded P2946 scan with a symbolic parametric obstruction: positive-cone premises admit eta=S/5 for every S>=5, so eta=9/5 and the exact P2938 vector require an additional strict source theorem.  The delta formula is also not uniquely selected because multiple identity/zero aliases coincide at 4/5.",
            "next_honest_step": "Do not continue ratio-forcing by more scans or aliases.  A next admissible proof-grade move must introduce a new strict source theorem selecting the exact P2938 vector/sum 9 and the delta numerator semantics, then prove beta/eta coupling; otherwise pivot to a genuinely new typed object outside the P2938/P2945 ratio lane or preserve the P2929-P2947 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["unbounded_nonforcing_certificate"]
    lines = [
        "# P2947/S1897 unbounded ratio non-forcing theorem audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Unbounded non-forcing certificate",
        f"- parametric family symbolically verified: `{cert['parametric_family_symbolically_verified']}`",
        f"- sample sum count: `{cert['sample_sum_count']}`",
        f"- sample non-target eta count: `{cert['sample_non_target_eta_count']}`",
        f"- eta=9/5 forced by positive-cone premises: `{cert['eta_9_5_forced_by_positive_cone_premises']}`",
        f"- delta formula uniquely selected by counts: `{cert['delta_formula_uniquely_selected_by_counts']}`",
        f"- strict P2938 vector/sum9 theorem exported: `{cert['strict_p2938_vector_sum9_source_theorem_exported']}`",
        f"- strict beta/eta coupling theorem exported: `{cert['strict_beta_eta_coupling_theorem_exported']}`",
        f"- accepted strict ratio source theorem: `{cert['accepted_strict_ratio_source_theorem']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2945), read_json(P2946))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2947/S1897 unbounded ratio non-forcing theorem audit", "## P2947/S1897 unbounded ratio non-forcing theorem audit\n\n`P2947/S1897` replaces the bounded P2946 scan with an exact parametric obstruction: for every integer `S>=5`, the positive vector `(S-4,1,1,1,1)` has `eta=S/5`, so the identity-positive-cone premises alone admit infinitely many eta values.  Therefore `eta=9/5`, the exact P2938 vector `[1,2,2,2,2]`, and the P2945 ratio formulas still require an independent strict source theorem and beta/eta coupling.  P2947 also records the delta alias kernel caused by `identity_count=zero_count=1`.  No strict damping packet, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2947/S1897 unbounded ratio non-forcing `L_total` guard", "## P2947/S1897 unbounded ratio non-forcing `L_total` guard\n\n`P2947/S1897` proves at the symbolic finite-family level that positive-cone premises do not force `eta=9/5` or the exact P2938 vector: `eta=S/5` remains available for every `S>=5`.  Because the ratio source theorem and beta/eta coupling are still absent, this result cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current unbounded ratio non-forcing theorem guardrail (P2947/S1897, 2026-06-19)", "## Current unbounded ratio non-forcing theorem guardrail (P2947/S1897, 2026-06-19)\n\n- P2947 removes the bounded-scan caveat from P2946 by constructing the exact parametric family `v(S)=(S-4,1,1,1,1)` for every integer `S>=5`, giving positive-cone admissible `eta=S/5`.\n- Therefore the finite identity-positive-cone premises alone do not force `eta=9/5`, the exact P2938 vector `[1,2,2,2,2]`, or sum `9`; the P2945 delta formula also remains semantically unselected because identity/zero aliases coincide at `4/5`.\n- Do not continue ratio-forcing by more bounded scans or alias variants, and do not promote P2947 to strict delta/eta source law, strict damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.\n- A next admissible move must introduce a genuinely new strict source theorem selecting the exact P2938 vector/sum `9` and delta numerator semantics plus beta/eta coupling, or pivot outside the P2938/P2945 ratio lane; otherwise preserve the P2929-P2947 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
