#!/usr/bin/env python3
"""P2946/S1896: delta/eta ratio strict-forcing obstruction matrix.

P2945 produced exact arithmetic candidates delta=4/5 and eta=9/5 from the
P2938/P2944 identity-positive-cone inventory.  P2946 asks the next stricter
question: are those exact ratio formulas forced by the finite identity-positive
cone premises alone?

The audit is deliberately narrow and computable.  It preserves the P2944 finite
identity/zero counts and the P2943 positive-cone shape, then enumerates bounded
positive prime-coordinate variants.  If the premises force eta, every admissible
positive vector should have prime_vector_sum/prime_count = 9/5.  They do not:
most admissible vectors give different eta values.  Likewise the finite counts
only make several denominator-equivalent delta formulas coincide because
identity_count == zero_count == 1; they do not select the P2945 formula as a
strict theorem.

Therefore P2946 exports an obstruction/witness table for strict ratio forcing,
not a damping packet or L_total closure.
"""
from __future__ import annotations

import hashlib
import json
from fractions import Fraction
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2940_s1890_p2938_carrier_aut_orbit_selector_burden import PRIME_VECTOR
from p2945_s1895_identity_positive_cone_delta_eta_source_candidate import OUT as P2945

OUT = GEN / "p2946_s1896_delta_eta_ratio_strict_forcing_obstruction_matrix.json"
MD = GEN / "p2946_s1896_delta_eta_ratio_strict_forcing_obstruction_matrix.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

TARGET_DELTA = Fraction(4, 5)
TARGET_ETA = Fraction(9, 5)
BOUNDED_COORDINATES = [1, 2, 3]


def frac_payload(value: Fraction) -> dict[str, Any]:
    return {
        "numerator": value.numerator,
        "denominator": value.denominator,
        "as_string": f"{value.numerator}/{value.denominator}",
        "as_float": float(value),
    }


def eta_variant_rows(prime_count: int) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for vector in product(BOUNDED_COORDINATES, repeat=prime_count):
        vector_sum = sum(vector)
        eta = Fraction(vector_sum, prime_count)
        rows.append({
            "prime_vector": list(vector),
            "prime_vector_sum": vector_sum,
            "eta": frac_payload(eta),
            "matches_p2945_eta_9_5": eta == TARGET_ETA,
            "is_p2938_vector": list(vector) == PRIME_VECTOR,
        })
    return rows


def delta_formula_alias_rows(prime_count: int, identity_count: int, zero_count: int) -> list[dict[str, Any]]:
    formulas = [
        ("(prime_count-identity_count)/prime_count", prime_count - identity_count),
        ("(prime_count-zero_count)/prime_count", prime_count - zero_count),
        ("(prime_count-min(identity_count,zero_count))/prime_count", prime_count - min(identity_count, zero_count)),
        ("(prime_count-max(identity_count,zero_count))/prime_count", prime_count - max(identity_count, zero_count)),
        ("(prime_count-(identity_count+zero_count)/2)/prime_count", Fraction(prime_count, 1) - Fraction(identity_count + zero_count, 2)),
    ]
    rows = []
    for formula, numerator in formulas:
        value = Fraction(numerator, prime_count)
        rows.append({
            "formula": formula,
            "value": frac_payload(value),
            "matches_p2945_delta_4_5": value == TARGET_DELTA,
            "selected_by_strict_theorem": False,
        })
    return rows


def forcing_obligation_rows(eta_rows: list[dict[str, Any]], delta_aliases: list[dict[str, Any]]) -> list[dict[str, Any]]:
    eta_match_count = sum(row["matches_p2945_eta_9_5"] for row in eta_rows)
    eta_distinct = sorted({row["eta"]["as_string"] for row in eta_rows}, key=lambda text: Fraction(text))
    return [
        {
            "obligation": "finite_identity_positive_cone_model_class_enumerated",
            "satisfied": True,
            "evidence": f"{len(eta_rows)} positive prime-coordinate vectors over {BOUNDED_COORDINATES}",
        },
        {
            "obligation": "eta_9_5_forced_by_positive_cone_premises",
            "satisfied": eta_match_count == len(eta_rows),
            "evidence": f"{eta_match_count} of {len(eta_rows)} bounded positive vectors have eta=9/5; distinct eta values: {eta_distinct}",
        },
        {
            "obligation": "p2938_vector_uniquely_forced_by_identity_positive_cone_premises",
            "satisfied": False,
            "evidence": "the premises allow many positive vectors; P2938 [1,2,2,2,2] is one admissible vector, not uniquely forced here",
        },
        {
            "obligation": "delta_formula_uniquely_selected_by_strict_theorem",
            "satisfied": False,
            "evidence": f"{sum(row['matches_p2945_delta_4_5'] for row in delta_aliases)} denominator-equivalent formulas match 4/5 because identity_count=zero_count=1",
        },
        {
            "obligation": "strict_beta_eta_coupling_theorem_exported",
            "satisfied": False,
            "evidence": "no theorem couples either the candidate eta or beta to the strict damping/compression term",
        },
    ]


def build_payload(p2945: dict[str, Any]) -> dict[str, Any]:
    inventory = p2945["constructed_theoretical_objects"]["source_inventory"]
    prime_count = inventory["prime_count"]
    eta_rows = eta_variant_rows(prime_count)
    delta_aliases = delta_formula_alias_rows(prime_count, inventory["identity_count"], inventory["zero_count"])
    obligations = forcing_obligation_rows(eta_rows, delta_aliases)
    eta_match_count = sum(row["matches_p2945_eta_9_5"] for row in eta_rows)
    accepted = all(row["satisfied"] for row in obligations)
    return {
        "status": "P2946_DELTA_ETA_RATIO_STRICT_FORCING_OBSTRUCTION_MATRIX_NO_STRICT_EXPORT",
        "input_hashes": {"P2945": hashlib.sha256(P2945.read_bytes()).hexdigest() if P2945.exists() else None},
        "constructed_theoretical_objects": {
            "candidate_object": "Identity_Positive_Cone_Ratio_Forcing_Theorem_Obstruction_Matrix",
            "bounded_positive_coordinate_domain": BOUNDED_COORDINATES,
            "eta_variant_summary": {
                "variant_count": len(eta_rows),
                "eta_9_5_match_count": eta_match_count,
                "eta_9_5_nonmatch_count": len(eta_rows) - eta_match_count,
                "p2938_vector_present": any(row["is_p2938_vector"] for row in eta_rows),
            },
            "eta_variant_rows": eta_rows,
            "delta_formula_alias_rows": delta_aliases,
            "forcing_obligation_rows": obligations,
        },
        "strict_forcing_certificate": {
            "finite_model_class_enumerated": True,
            "eta_9_5_forced_by_positive_cone_premises": obligations[1]["satisfied"],
            "p2938_vector_uniquely_forced": False,
            "delta_formula_uniquely_selected": False,
            "strict_ratio_source_theorem_exported": False,
            "strict_beta_eta_coupling_theorem_exported": False,
            "accepted_strict_delta_eta_source_law": accepted,
        },
        "decision": {
            "positive_witnesses": {
                "p2945_exact_ratios_reused": True,
                "bounded_positive_model_class_computed": True,
                "p2938_vector_is_in_model_class": any(row["is_p2938_vector"] for row in eta_rows),
            },
            "negative_export_flags": {
                "strict_ratio_source_theorem_exported": False,
                "strict_prime_log_value_source_exported": False,
                "strict_delta_eta_source_law_exported": False,
                "strict_beta_eta_coupling_theorem_exported": False,
                "strict_damping_beta_eta_source_packet_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_hamiltonian_exported": False,
                "bridge_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2946 tests whether the P2945 ratio formulas are forced by finite identity-positive-cone premises alone.  They are not: eta=9/5 holds only on a proper subset of bounded positive prime-coordinate variants, P2938 is admissible but not uniquely forced by these premises, and the delta formula has denominator-equivalent aliases because identity_count equals zero_count.",
            "next_honest_step": "Do not add more ratio aliases or bounded positive-vector scans.  A next proof-grade move must supply an independent strict theorem selecting the exact P2938 prime-coordinate vector and the P2945 ratio formulas, together with beta/eta damping coupling; otherwise pivot to a genuinely new typed object outside the P2938/P2945 ratio-forcing lane or preserve the P2929-P2946 no-strict-export boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["strict_forcing_certificate"]
    summary = payload["constructed_theoretical_objects"]["eta_variant_summary"]
    lines = [
        "# P2946/S1896 delta/eta ratio strict-forcing obstruction matrix",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Strict-forcing certificate",
        f"- finite model class enumerated: `{cert['finite_model_class_enumerated']}`",
        f"- eta=9/5 forced by positive-cone premises: `{cert['eta_9_5_forced_by_positive_cone_premises']}`",
        f"- P2938 vector uniquely forced: `{cert['p2938_vector_uniquely_forced']}`",
        f"- delta formula uniquely selected: `{cert['delta_formula_uniquely_selected']}`",
        f"- strict ratio source theorem exported: `{cert['strict_ratio_source_theorem_exported']}`",
        f"- strict beta/eta coupling theorem exported: `{cert['strict_beta_eta_coupling_theorem_exported']}`",
        f"- accepted strict delta/eta source law: `{cert['accepted_strict_delta_eta_source_law']}`",
        "",
        "## Finite model-class summary",
        f"- positive coordinate variant count: `{summary['variant_count']}`",
        f"- eta=9/5 match count: `{summary['eta_9_5_match_count']}`",
        f"- eta=9/5 nonmatch count: `{summary['eta_9_5_nonmatch_count']}`",
        f"- P2938 vector present: `{summary['p2938_vector_present']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2945))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2946/S1896 delta/eta ratio strict-forcing obstruction", "## P2946/S1896 delta/eta ratio strict-forcing obstruction\n\n`P2946/S1896` tests whether the P2945 formulas are forced by the finite identity-positive-cone premises.  In the bounded positive prime-coordinate model class, `eta=9/5` holds only on a proper subset of admissible vectors, the exact P2938 vector is admissible but not uniquely forced, and the `delta=4/5` formula has denominator-equivalent aliases because `identity_count=zero_count=1`.  Thus P2946 exports a ratio-forcing obstruction, not a strict delta/eta source law, beta/eta coupling theorem, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2946/S1896 ratio-forcing obstruction `L_total` guard", "## P2946/S1896 ratio-forcing obstruction `L_total` guard\n\n`P2946/S1896` shows that the P2945 ratio formulas are not forced by the finite identity-positive-cone premises alone.  Since eta, the exact P2938 vector, and the delta formula selection still lack an independent strict source theorem and beta/eta coupling, this obstruction cannot enter nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current delta/eta ratio strict-forcing obstruction guardrail (P2946/S1896, 2026-06-19)", "## Current delta/eta ratio strict-forcing obstruction guardrail (P2946/S1896, 2026-06-19)\n\n- P2946 tests whether the P2945 exact formulas are forced by the finite identity-positive-cone premises alone.\n- The bounded positive-coordinate model class leaves eta nonunique: `eta=9/5` occurs only on a proper subset of admissible vectors, and the exact P2938 vector is admissible but not uniquely forced by these premises.\n- The `delta=4/5` formula also remains theorem-unselected because denominator-equivalent aliases coincide under `identity_count=zero_count=1`; no beta/eta coupling theorem is exported.\n- Do not continue ratio-alias or bounded positive-vector scans as primary strategy, and do not promote P2946 to strict delta/eta source law, strict damping packet, nonproxy `L_total`, bridge closure, role transfer, or ToE.  A next admissible move must supply an independent strict theorem selecting the exact P2938 vector and P2945 ratio formulas plus beta/eta coupling, or pivot outside this ratio-forcing lane; otherwise preserve the P2929-P2946 boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
