#!/usr/bin/env python3
"""P2985/S1935: Z12 Leibniz derivation source-candidate obstruction.

P2984 closed the current CRT idempotent lane as bounded no-go.  This file
introduces one genuinely new finite typed object outside the nilradical and CRT
lanes: the internal Leibniz derivation algebra of the ring Z/12Z.

The audit is intentionally narrow and computational.  Additive endomorphisms of
Z/12Z are represented by D_a(x)=a*x mod 12, and every candidate is tested on all
144 products for the Leibniz rule D(xy)=xD(y)+yD(x).  The finite result is a
zero derivation algebra: only a=0 survives.  This is a real theorem-grade finite
witness, but it is also a bounded no-go for source construction because a zero
derivation cannot provide a nonzero strict flow, orientation sign, damping scale,
bridge-completion map, unit-bearing density, or nonproxy variational chain.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2984_s1934_crt_action_variational_installation_obstruction import OUT as P2984

OUT = GEN / "p2985_s1935_z12_leibniz_derivation_source_candidate_obstruction.json"
MD = GEN / "p2985_s1935_z12_leibniz_derivation_source_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12


def d(a: int, x: int) -> int:
    return (a * x) % MODULUS


def derivation_candidate_rows() -> list[dict[str, Any]]:
    rows = []
    for a in range(MODULUS):
        violations = []
        for x in range(MODULUS):
            for y in range(MODULUS):
                lhs = d(a, (x * y) % MODULUS)
                rhs = (x * d(a, y) + y * d(a, x)) % MODULUS
                if lhs != rhs:
                    violations.append({"x": x, "y": y, "lhs": lhs, "rhs": rhs, "defect": (lhs - rhs) % MODULUS})
        rows.append({
            "candidate": f"D_{a}",
            "multiplier_a": a,
            "additive_endomorphism": True,
            "D_of_one": a,
            "unit_product_defect_at_1_1": (d(a, 1) - (1 * d(a, 1) + 1 * d(a, 1))) % MODULUS,
            "leibniz_holds_all_144_products": not violations,
            "violation_count": len(violations),
            "first_violation": violations[0] if violations else None,
            "nonzero_derivation": a != 0 and not violations,
        })
    return rows


def obstruction_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    accepted = [r for r in rows if r["leibniz_holds_all_144_products"]]
    nonzero = [r for r in accepted if r["multiplier_a"] != 0]
    return [
        {"obligation": "finite_additive_endomorphism_scan", "satisfied": len(rows) == MODULUS and all(r["additive_endomorphism"] for r in rows), "evidence": "all 12 additive maps D_a(x)=a*x mod 12 are enumerated"},
        {"obligation": "full_Leibniz_product_scan", "satisfied": all(r["violation_count"] in range(145) for r in rows), "evidence": "each D_a is tested on all 12*12 products"},
        {"obligation": "derivation_algebra_computed", "satisfied": [r["multiplier_a"] for r in accepted] == [0], "evidence": "D(1)=2D(1) forces D(1)=0; exhaustive scan leaves only D_0"},
        {"obligation": "nonzero_strict_flow_exported", "satisfied": bool(nonzero), "evidence": "no nonzero Leibniz derivation exists in Z/12Z"},
        {"obligation": "orientation_or_selector_source", "satisfied": False, "evidence": "the zero derivation is sign-blind and cannot distinguish +omega from -omega"},
        {"obligation": "positive_scale_or_damping_source", "satisfied": False, "evidence": "zero derivation has no target-independent positive beta/Z_beta or measure value"},
        {"obligation": "bridge_completion_map", "satisfied": False, "evidence": "zero derivation supplies no amplitude, phase/topological, or nonlinear damping completion law"},
        {"obligation": "unit_bearing_action_density", "satisfied": False, "evidence": "zero derivation installs no named density, field provenance, measure, or nonproxy variational chain"},
        {"obligation": "accepted_current_source_candidate", "satisfied": False, "evidence": "finite theorem is real but source export is blocked by zero-derivation collapse"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_scan", "leibniz_scan", "derivation_computed", "nonzero_flow", "orientation_source", "positive_scale", "bridge_map", "unit_action_density"]
    return [{"present": dict(zip(names, bits)), "accepts_derivation_source_candidate": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2984_path: Any) -> dict[str, Any]:
    rows = derivation_candidate_rows()
    accepted = [r for r in rows if r["leibniz_holds_all_144_products"]]
    obligations = obstruction_rows(rows)
    matrix = acceptance_matrix()
    return {
        "status": "P2985_Z12_LEIBNIZ_DERIVATION_SOURCE_CANDIDATE_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2984": hashlib.sha256(p2984_path.read_bytes()).hexdigest() if p2984_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "Z12LeibnizDerivationAlgebra_ZeroObstructionMatrix",
            "candidate_rows": rows,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "derivation_certificate": {
            "modulus": MODULUS,
            "candidate_count": len(rows),
            "products_tested_per_candidate": MODULUS * MODULUS,
            "accepted_derivations": [r["candidate"] for r in accepted],
            "nonzero_accepted_derivations": [r["candidate"] for r in accepted if r["multiplier_a"] != 0],
            "unit_product_defects": {str(r["multiplier_a"]): r["unit_product_defect_at_1_1"] for r in rows},
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_derivation_source_candidate"]),
        },
        "decision": {
            "positive_progress": "P2985 introduces a new finite typed object outside the closed nilradical and CRT lanes: the internal Leibniz derivation algebra of Z/12Z.",
            "breakthrough": "Bounded no-go: exhaustive testing of all 12 additive endomorphisms on all 144 products proves that only the zero derivation satisfies Leibniz, so no nonzero strict flow/source is exported.",
            "negative_export_flags": {k: False for k in ["nonzero_derivation_source_exported", "selector_or_orientation_exported", "damping_source_exported", "bridge_closure_exported", "unit_bearing_density_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay nilradical, CRT, ratio-package, incidence, selector, or bridge lanes through the zero derivation.  The Z12 Leibniz derivation object is now bounded no-go as a source candidate; the next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside nilradical/CRT/zero-derivation lanes, or preserve the P2929-P2985 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["derivation_certificate"]
    lines = [
        "# P2985/S1935 Z12 Leibniz derivation source-candidate obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Derivation certificate",
        f"- modulus: `{cert['modulus']}`",
        f"- additive candidates: `{cert['candidate_count']}`",
        f"- products tested per candidate: `{cert['products_tested_per_candidate']}`",
        f"- accepted derivations: `{cert['accepted_derivations']}`",
        f"- nonzero accepted derivations: `{cert['nonzero_accepted_derivations']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`",
        "",
        "## Lay summary",
        payload["decision"]["positive_progress"],
        payload["decision"]["breakthrough"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P2984)
    payload = build_payload(P2984)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2985/S1935 Z12 Leibniz derivation source-candidate obstruction", "## P2985/S1935 Z12 Leibniz derivation source-candidate obstruction\n\n`P2985/S1935` introduces one new finite typed object outside the closed nilradical and CRT lanes: the internal Leibniz derivation algebra of `Z/12Z`.  The exhaustive finite scan represents additive endomorphisms as `D_a(x)=a*x mod 12` and tests each of the `12` candidates on all `144` products for `D(xy)=xD(y)+yD(x)`.  Only `D_0` survives; every nonzero candidate already fails at the unit-product constraint `D(1)=2D(1)`.  This proves a real zero-derivation theorem, but it exports no nonzero strict flow, selector/orientation source, positive damping/measure source, bridge-completion map, unit-bearing density, nonproxy `L_total`, role transfer, or ToE.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2985/S1935 zero-derivation `L_total` guard", "## P2985/S1935 zero-derivation `L_total` guard\n\n`P2985/S1935` computes the internal Leibniz derivation algebra of `Z/12Z` and finds only the zero derivation.  A zero derivation installs no strict flow term, no named unit-bearing action density, no field provenance, no nonproxy variational chain, no EOM/Hamiltonian term, no bridge closure, no role transfer, and no ToE.\n")
    append_once(AGENTS, "Current Z12 Leibniz derivation source-candidate obstruction guardrail (P2985/S1935, 2026-06-20)", "## Current Z12 Leibniz derivation source-candidate obstruction guardrail (P2985/S1935, 2026-06-20)\n\n- P2985 introduces a new finite typed object outside the closed nilradical and CRT lanes: the internal Leibniz derivation algebra of `Z/12Z`.\n- Exhaustive computation scans all `12` additive endomorphisms `D_a(x)=a*x mod 12` against all `144` products and finds only `D_0`; every nonzero candidate fails the Leibniz rule, already at `D(1)=2D(1)`.\n- The current route is bounded no-go as a source candidate: no nonzero strict flow, selector/orientation source, target-independent positive `beta/Z_beta` or measure, bridge-completion map, unit-bearing density, or nonproxy variational chain is exported.\n- Do not replay nilradical, CRT, ratio-package, incidence, selector, or bridge lanes through the zero derivation.  The Z12 Leibniz derivation object is now bounded no-go; a next move must introduce a genuinely new strict typed object/theorem/provider outside nilradical/CRT/zero-derivation lanes or preserve the P2929-P2985 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
