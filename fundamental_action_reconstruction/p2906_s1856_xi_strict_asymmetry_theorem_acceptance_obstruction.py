#!/usr/bin/env python3
"""P2906/S1856: Xi strict asymmetry theorem acceptance obstruction.

P2905 found no current provenance export for Xi_{0,+}.  P2906 constructs the
next missing object as an acceptance theorem schema: what would a strict
asymmetry/chiral/defect-generation theorem have to output in order to select
Xi_{0,+} over the 23 translated/sign-flipped alternatives?

The finite gate separates sign-breaking from origin-breaking.  A chiral sign
source can reduce 24 Xi alternatives to 12, but without a nonimported origin law
there is still no distinguished basepoint.  A translation-equivariant theorem
from translation-neutral strict data to the pointed Xi target has no fixed point,
so it cannot export Xi_{0,+}.  The remaining missing object is sharpened to a
joint origin-and-sign strict source theorem, not another Xi candidate mention.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
import hashlib
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN

P2905 = GEN / "p2905_s1855_xi_dirac_source_provenance_alternative_audit.json"
OUT = GEN / "p2906_s1856_xi_strict_asymmetry_theorem_acceptance_obstruction.json"
MD = GEN / "p2906_s1856_xi_strict_asymmetry_theorem_acceptance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
N = 12


def xi_family() -> list[dict[str, Any]]:
    return [
        {
            "xi": f"Xi_{{{b},{'+' if s > 0 else '-'}}}",
            "basepoint": b,
            "sign": s,
            "defect_edge": [b, (b + s * 5) % N],
            "coupled_axiom": f"A({b},{'+' if s > 0 else '-'})",
        }
        for b in range(N)
        for s in (-1, 1)
    ]


def translate(row: dict[str, Any], shift: int) -> tuple[int, int]:
    return ((row["basepoint"] + shift) % N, row["sign"])


def orbit(seed: dict[str, Any]) -> list[tuple[int, int]]:
    return sorted({translate(seed, shift) for shift in range(N)})


def fixed_points(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [row for row in rows if all(translate(row, shift) == (row["basepoint"], row["sign"]) for shift in range(N))]


def sign_slice(rows: list[dict[str, Any]], sign: int) -> list[dict[str, Any]]:
    return [row for row in rows if row["sign"] == sign]


def build_payload(p2905: dict[str, Any]) -> dict[str, Any]:
    rows = xi_family()
    plus = sign_slice(rows, 1)
    minus = sign_slice(rows, -1)
    fps = fixed_points(rows)
    theorem_schema = {
        "target": "Xi = Z12 x {-,+}",
        "required_outputs": [
            "origin/basepoint b0 in Z12",
            "polarity sigma0 in {-,+}",
            "coupling Xi_{b0,sigma0} -> A(b0,sigma0) -> D=(b0,b0+5*sigma0)",
            "strict nadsoliton provenance for both b0 and sigma0",
        ],
        "acceptance_rule": "A theorem sourced only by translation-neutral strict data must be translation-equivariant; a selected point must therefore be a translation fixed point of Xi.",
    }
    obstruction_rows = [
        {
            "candidate_source_kind": "translation-neutral strict invariant",
            "computed_remaining_choices": len(rows),
            "fixed_points": len(fps),
            "can_select_xi_0_plus": False,
            "reason": "Xi has no translation fixed point.",
        },
        {
            "candidate_source_kind": "chiral/sign source only, no origin source",
            "computed_remaining_choices": len(plus),
            "fixed_points_within_positive_sign_slice": len(fixed_points(plus)),
            "can_select_xi_0_plus": False,
            "reason": "The positive sign slice is one free translation orbit of 12 basepoints.",
        },
        {
            "candidate_source_kind": "origin source only, no sign source",
            "computed_remaining_choices": 2,
            "can_select_xi_0_plus": False,
            "reason": "A basepoint without polarity leaves Xi_{b,+} versus Xi_{b,-}.",
        },
        {
            "candidate_source_kind": "joint origin-and-sign strict theorem",
            "computed_remaining_choices": 1,
            "can_select_xi_0_plus": "only_if_provenance_and_coupling_are_exported",
            "reason": "This is the missing object; it is not present in P2905.",
        },
    ]
    return {
        "status": "P2906_XI_STRICT_ASYMMETRY_THEOREM_ACCEPTANCE_OBSTRUCTION_NO_EXPORT",
        "input_hashes": {"P2905": hashlib.sha256(P2905.read_bytes()).hexdigest() if P2905.exists() else None},
        "constructed_theoretical_objects": {
            "xi_target_family": rows,
            "strict_asymmetry_theorem_schema": theorem_schema,
            "positive_sign_slice": plus,
            "negative_sign_slice": minus,
            "translation_orbits": [orbit(plus[0]), orbit(minus[0])],
            "acceptance_obstruction_rows": obstruction_rows,
        },
        "acceptance_matrix": {
            "p2905_rechecked_no_positive_provenance": p2905.get("acceptance_matrix", {}).get("positive_provenance_hit_count") == 0,
            "xi_target_count": len(rows),
            "translation_orbit_count": 2,
            "translation_orbit_sizes": [len(orbit(plus[0])), len(orbit(minus[0]))],
            "translation_fixed_point_count": len(fps),
            "chiral_sign_only_remaining_basepoints": len(plus),
            "origin_only_remaining_polarities": 2,
            "joint_origin_sign_theorem_required": True,
            "joint_origin_sign_theorem_exported": False,
            "accepted_as_strict_source_theorem": False,
        },
        "decision": {
            "positive_witnesses": {
                "acceptance_theorem_schema_constructed": True,
                "sign_vs_origin_obligation_separated": True,
                "finite_fixed_point_obstruction_recomputed": True,
            },
            "negative_export_flags": {
                "strict_asymmetry_chiral_theorem_exported": False,
                "strict_origin_source_exported": False,
                "strict_sign_source_exported": False,
                "xi_0_plus_strict_source_exported": False,
                "unit_bearing_strict_density_exported": False,
                "nonproxy_ltotal_exported": False,
                "eom_closure_exported": False,
                "hamiltonian_closure_exported": False,
                "role_transfer_exported": False,
                "toe_closure_exported": False,
            },
            "reason": "P2906 constructs the acceptance theorem required after P2905 and computes that sign-breaking alone leaves 12 translated Xi candidates, origin-breaking alone leaves two polarities, and translation-neutral strict data has 0 fixed points in the Xi target.  Thus no strict Xi_{0,+} source theorem is exported; the exact missing object is a joint origin-and-sign strict provenance theorem with coupling.",
            "next_honest_step": "Either provide one concrete strict nadsoliton-derived origin-and-sign theorem that computes b=0 and sigma=+ with coupling to A(0,+), or pivot outside the Xi/defect-placement lane.  More sign-only, origin-only, translation-neutral, or inventory variants should be treated as repetition-gated; without the joint theorem, preserve no-new-live-frontier.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    acc = payload["acceptance_matrix"]
    lines = [
        "# P2906/S1856 Xi strict asymmetry theorem acceptance obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite gate",
        f"- Xi target count: `{acc['xi_target_count']}`",
        f"- translation orbit sizes: `{acc['translation_orbit_sizes']}`",
        f"- translation fixed points: `{acc['translation_fixed_point_count']}`",
        f"- chiral/sign-only remaining basepoints: `{acc['chiral_sign_only_remaining_basepoints']}`",
        f"- origin-only remaining polarities: `{acc['origin_only_remaining_polarities']}`",
        f"- joint origin-sign theorem exported: `{acc['joint_origin_sign_theorem_exported']}`",
        "",
        "## Boundary",
        payload["decision"]["reason"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(read_json(P2905))
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2906/S1856 Xi strict asymmetry theorem acceptance obstruction", "## P2906/S1856 Xi strict asymmetry theorem acceptance obstruction\n\n`P2906/S1856` constructs the missing strict asymmetry/chiral theorem as an acceptance schema over the `Xi = Z12 x {-,+}` target.  The finite gate separates the obligations: a sign/chiral source alone leaves `12` translated basepoints, an origin source alone leaves `2` polarities, and translation-neutral strict data has `0` fixed points in the `Xi` target.  Therefore the exact missing object is a joint origin-and-sign strict provenance theorem with coupling to `A(0,+)`, not another sign-only, origin-only, translation-neutral, or inventory variant.  No strict source, unit-bearing density, nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2906/S1856 Xi theorem `L_total` guard", "## P2906/S1856 Xi theorem `L_total` guard\n\n`P2906/S1856` adds no action density.  It proves at the finite target level that sign-only provenance still leaves `12` possible basepoints and origin-only provenance still leaves `2` polarities; only a joint origin-and-sign strict theorem could source `Xi_{0,+}`.  Since no such theorem or unit-bearing coupling is exported, no nonproxy `L_total`, EOM, Hamiltonian, role transfer, or ToE closure follows.\n")
    append_once(AGENTS, "Current Xi strict asymmetry theorem acceptance obstruction guardrail (P2906/S1856, 2026-06-19)", "## Current Xi strict asymmetry theorem acceptance obstruction guardrail (P2906/S1856, 2026-06-19)\n\n- P2906 constructs the exact acceptance schema for a strict asymmetry/chiral theorem selecting `Xi_{0,+}` over the `23` alternatives.\n- The finite target gate finds `0` translation fixed points; sign/chiral sourcehood alone leaves `12` translated basepoints, while origin sourcehood alone leaves `2` polarities.\n- Do not promote sign-only, origin-only, translation-neutral, candidate, inventory, symbolic `rho_9/5`, or `U_9_5` data to strict sourcehood, nonproxy `L_total`, EOM, Hamiltonian, bridge closure, role transfer, or ToE closure.\n- A next admissible proof-grade move must provide a joint origin-and-sign strict provenance theorem computing `b=0` and `sigma=+` with coupling to `A(0,+)`, pivot outside the Xi/defect-placement lane, or preserve no-new-live-frontier.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
