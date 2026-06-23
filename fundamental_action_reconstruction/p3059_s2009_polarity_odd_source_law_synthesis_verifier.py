#!/usr/bin/env python3
"""P3059/S2009: polarity-odd source-law synthesis verifier.

P3058 identified the exact missing atom more sharply: a strict polarity-odd
source-law boundary condition coupled to G_selector.  P3059 constructs the
smallest honest synthesis problem for that object instead of adding another
polarity-even compatibility check.  It combines the currently available odd
sign clues as a finite coefficient module and asks whether any current rule
selects one non-premise sign for the boundary condition.
"""
from __future__ import annotations

import hashlib, itertools, json, subprocess
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3058_s2008_unique_polarity_coupling_constraint_verifier import OUT as P3058

OUT = GEN / "p3059_s2009_polarity_odd_source_law_synthesis_verifier.json"
MD = GEN / "p3059_s2009_polarity_odd_source_law_synthesis_verifier.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "polarity_odd_source_law_boundary_condition": r"polarity-odd source-law|polarity_odd_source_law|strict_polarity_odd_source_law_boundary_condition|polarity-odd.*boundary",
    "signed_clue_module": r"receiver.*winding|chiral.*bispectrum|boundary.*cocycle|inversion-odd.*sign",
    "coefficient_sign_pair_obstruction": r"sign-pair|global.*sign|lambda.*unfixed|opposite polarity|unique.*polarity",
}

UNITS = [1, 5, 7, 11]
REVERSING_UNITS = [7, 11]
COEFF_RANGE = range(-2, 3)

ODD_CLUE_BASIS = [
    {"basis": "boundary_cocycle_orientation", "base_value": 1, "character": {1: 1, 5: 1, 7: -1, 11: -1}},
    {"basis": "chiral_bispectrum_sign", "base_value": 1, "character": {1: 1, 5: 1, 7: -1, 11: -1}},
    {"basis": "receiver_winding_sign", "base_value": 1, "character": {1: 1, 5: 1, 7: -1, 11: -1}},
]


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def sign(value: int) -> int:
    return (value > 0) - (value < 0)


def candidate_value(coefficients: tuple[int, ...], unit: int = 1) -> int:
    total = 0
    for coeff, basis in zip(coefficients, ODD_CLUE_BASIS):
        total += coeff * basis["base_value"] * basis["character"][unit]
    return total


def is_inversion_odd(coefficients: tuple[int, ...]) -> bool:
    base = candidate_value(coefficients, 1)
    return all(candidate_value(coefficients, u) == (-base if u in REVERSING_UNITS else base) for u in UNITS)


def enumerate_synthesis_module() -> list[dict[str, Any]]:
    rows = []
    for coeffs in itertools.product(COEFF_RANGE, repeat=len(ODD_CLUE_BASIS)):
        if all(c == 0 for c in coeffs):
            continue
        value = candidate_value(coeffs, 1)
        rows.append({
            "coefficients": dict(zip([b["basis"] for b in ODD_CLUE_BASIS], coeffs)),
            "coefficient_tuple": list(coeffs),
            "base_value": value,
            "base_sign": sign(value),
            "nonzero_signed_value": value != 0,
            "inversion_odd": is_inversion_odd(coeffs),
            "paired_tuple": [-c for c in coeffs],
            "paired_base_value": -value,
            "same_support_pair_present": True,
            "accepted_as_unique_nonpremise_boundary_condition": False,
            "blocker": "global coefficient sign remains paired unless a strict sign-normalization/source law is added",
        })
    return rows


def quotient_summary(rows: list[dict[str, Any]]) -> dict[str, Any]:
    nonzero = [r for r in rows if r["nonzero_signed_value"]]
    positive = [r for r in nonzero if r["base_sign"] > 0]
    negative = [r for r in nonzero if r["base_sign"] < 0]
    zero = [r for r in rows if not r["nonzero_signed_value"]]
    canonical_pairs = set()
    for row in nonzero:
        tup = tuple(row["coefficient_tuple"])
        neg = tuple(-c for c in tup)
        canonical_pairs.add(min(tup, neg))
    return {
        "total_nonzero_coefficient_vectors": len(rows),
        "nonzero_signed_candidates": len(nonzero),
        "zero_value_candidates": len(zero),
        "positive_polarity_candidates": len(positive),
        "negative_polarity_candidates": len(negative),
        "sign_pair_orbits": len(canonical_pairs),
        "accepted_unique_boundary_conditions": 0,
        "reason": "every nonzero synthesized odd law has the equally admissible negative coefficient vector on the same support",
    }


def construct_missing_object() -> dict[str, Any]:
    return {
        "object": "StrictPolarityOddSourceLawBoundaryConditionSynthesisModule",
        "carrier_symbol": "G_selector",
        "target_atom": "strict_polarity_odd_source_law_boundary_condition",
        "basis_clues": [b["basis"] for b in ODD_CLUE_BASIS],
        "coefficient_domain": "integer coefficients in [-2,2]^3 excluding the zero vector",
        "acceptance_rule": "accept only a nonzero inversion-odd synthesized law whose global coefficient sign is fixed by a strict non-premise rule rather than by convention",
    }


def build_payload() -> dict[str, Any]:
    p3058 = read_json(P3058)
    grep_rows = content_grep()
    candidates = enumerate_synthesis_module()
    summary = quotient_summary(candidates)
    obligations = [
        {"obligation": "content_first_grep_before_synthesis", "satisfied": True, "detail": "three content patterns searched for polarity-odd boundary-condition and sign-pair obstruction material"},
        {"obligation": "construct_polarity_odd_source_law_synthesis_module", "satisfied": True, "detail": "current odd clues are represented as a finite coefficient module on G_selector"},
        {"obligation": "exhaust_bounded_integer_coefficient_domain", "satisfied": True, "detail": "all nonzero coefficient vectors in [-2,2]^3 are enumerated"},
        {"obligation": "export_nonpremise_global_coefficient_sign_rule", "satisfied": False, "detail": "every nonzero law is paired with its negative and no current artifact selects one member"},
        {"obligation": "discharge_selector_or_ltotal", "satisfied": False, "detail": "no QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure follows"},
    ]
    return {
        "status": "P3059_POLARITY_ODD_SOURCE_LAW_SYNTHESIS_SIGN_PAIR_OBSTRUCTION_NO_EXPORT",
        "input_hashes": {"P3058": hashlib.sha256(P3058.read_bytes()).hexdigest() if P3058.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "synthesis_module": construct_missing_object(),
            "odd_clue_basis": ODD_CLUE_BASIS,
            "candidate_rows": candidates,
            "quotient_summary": summary,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(row["hit_count"] for row in grep_rows),
            "odd_basis_size": len(ODD_CLUE_BASIS),
            "coefficient_values_per_basis": len(list(COEFF_RANGE)),
            "total_nonzero_coefficient_vectors": summary["total_nonzero_coefficient_vectors"],
            "nonzero_signed_candidates": summary["nonzero_signed_candidates"],
            "zero_value_candidates": summary["zero_value_candidates"],
            "positive_polarity_candidates": summary["positive_polarity_candidates"],
            "negative_polarity_candidates": summary["negative_polarity_candidates"],
            "sign_pair_orbits": summary["sign_pair_orbits"],
            "accepted_unique_boundary_conditions": summary["accepted_unique_boundary_conditions"],
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "p3058_status_seen": p3058.get("status"),
        },
        "proof_obligations": obligations,
        "decision": {
            "bounded_no_go": "P3059 constructs the P3058-requested polarity-odd source-law boundary condition as a finite synthesis module over three current inversion-odd clues.  The bounded integer search over [-2,2]^3 has 124 nonzero coefficient vectors: 106 give nonzero signed odd candidates, split evenly into 53 positive and 53 negative polarities, and 53 sign-pair orbits.  Each nonzero candidate is paired with its negative on the same support, so the current clues supply oddness and signed values but not a non-premise global coefficient-sign rule.  Thus no unique boundary condition, G_selector, QW-2191 discharge, selector closure, L_total, bridge, role transfer, or ToE closure is exported.",
            "negative_export_flags": {k: False for k in ["strict_polarity_odd_source_law_boundary_condition_exported", "unique_polarity_coupling_exported", "actual_g_selector_formula_exported", "selector_gluing_object_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "The next proof-grade move should not merely add another inversion-odd clue.  It must construct a strict non-premise global coefficient-sign normalization/source rule for the P3059 synthesis module, or prove that such a rule is impossible in a larger named coefficient/source class.  If that cannot be done, pivot to a different P3057 atom while preserving no selector/readout/L_total/bridge/ToE export.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3059/S2009 polarity-odd source-law synthesis verifier", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- odd basis size: `{c['odd_basis_size']}`",
        f"- coefficient values per basis: `{c['coefficient_values_per_basis']}`",
        f"- total nonzero coefficient vectors: `{c['total_nonzero_coefficient_vectors']}`",
        f"- nonzero signed candidates: `{c['nonzero_signed_candidates']}`",
        f"- zero-value candidates: `{c['zero_value_candidates']}`",
        f"- positive polarity candidates: `{c['positive_polarity_candidates']}`",
        f"- negative polarity candidates: `{c['negative_polarity_candidates']}`",
        f"- sign-pair orbits: `{c['sign_pair_orbits']}`",
        f"- accepted unique boundary conditions: `{c['accepted_unique_boundary_conditions']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3059/S2009 polarity-odd source-law synthesis verifier", "## P3059/S2009 polarity-odd source-law synthesis verifier\n\n`P3059/S2009` responds to the P3058 recommendation by constructing a finite synthesis module for a strict polarity-odd source-law boundary condition on `G_selector`.  It combines three current inversion-odd clues — boundary-cocycle orientation, chiral-bispectrum sign, and receiver-winding sign — over integer coefficients in `[-2,2]^3` and enumerates all `124` nonzero coefficient vectors.  `106` vectors give nonzero signed odd candidates, split into `53` positive and `53` negative polarities across `53` sign-pair orbits; `0` candidates export a non-premise global coefficient-sign rule.  No unique boundary condition, actual `G_selector`, `QW-2191` discharge, selector closure, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3059/S2009 polarity-odd source-law synthesis `L_total` guard", "## P3059/S2009 polarity-odd source-law synthesis `L_total` guard\n\n`P3059/S2009` adds no physical `L_total` term.  Its finite synthesis module produces paired polarity-odd signed candidates but no strict non-premise global coefficient-sign source rule, no unit-bearing action/EOM carrier, and no nonproxy variational chain rule.\n")
    append_once(AGENTS, "Current polarity-odd source-law synthesis guardrail (P3059/S2009, 2026-06-23)", "## Current polarity-odd source-law synthesis guardrail (P3059/S2009, 2026-06-23)\n\n- P3059 constructs the P3058-requested strict polarity-odd source-law boundary-condition synthesis module on `G_selector` from three current inversion-odd clues: boundary-cocycle orientation, chiral-bispectrum sign, and receiver-winding sign.\n- The bounded coefficient search over `[-2,2]^3` enumerates all `124` nonzero coefficient vectors; `106` are nonzero signed odd candidates, split into `53` positive and `53` negative polarities across `53` sign-pair orbits, with `0` accepted unique boundary conditions.\n- Do not promote another inversion-odd clue, signed receiver/chiral observable, coefficient convention, or sign-paired source law to `QW-2191` discharge, selector closure, observed physics, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move in this lane must construct a strict non-premise global coefficient-sign normalization/source rule for this synthesis module, or prove impossibility in a larger named coefficient/source class; otherwise attack a different P3057 atom while preserving the P3048-P3059 bounded no-export boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
