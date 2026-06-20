#!/usr/bin/env python3
"""P2979/S1929: nilradical action/variational installation obstruction.

P2978 left the final nilradical-lane route: install the Z12 nilpotent anchor in
a genuinely unit-bearing named action density with a nonproxy variational chain.
This audit attacks that route only.  It does not replay nilradical provenance,
named source-atom coupling, ratio/Gamma/incidence lanes, selector closure,
generic bridge completion, role transfer, or formal L_total promotion.

The finite calculation constructs the strongest formal toy receivers for the
nilpotent 6: a linear scalar weight in a quadratic density and a nilpotent-square
weight.  The first has a nonzero formal Euler coefficient, and the second is
identically zero because 6^2=0 mod 12.  Neither is a strict action theorem:
there is no unit-bearing measure, strict field provenance, boundary/integration
map, nonproxy continuum lift, or named variational chain.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from math import gcd
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2978_s1928_nilradical_named_source_atom_coupling_obstruction import OUT as P2978

OUT = GEN / "p2979_s1929_nilradical_action_variational_installation_obstruction.json"
MD = GEN / "p2979_s1929_nilradical_action_variational_installation_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MODULUS = 12
NILPOTENT = 6
UNITS = [x for x in range(MODULUS) if gcd(x, MODULUS) == 1]
FIELDS = ["psi"]


def formal_variational_receivers() -> dict[str, Any]:
    unit_images = {str(u): (u * NILPOTENT) % MODULUS for u in UNITS}
    linear_weight = NILPOTENT
    square_weight = (NILPOTENT * NILPOTENT) % MODULUS
    receivers = [
        {
            "receiver": "formal_nilpotent_linear_quadratic_density",
            "density_template": "L_nil = mu_nil * 6 * psi^2 / 2",
            "weight_mod_12": linear_weight,
            "formal_euler_coefficient": linear_weight,
            "formal_hessian_coefficient": linear_weight,
            "nonzero_formal_variation": linear_weight != 0,
            "unit_bearing_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_variational_installation": False,
        },
        {
            "receiver": "nilpotent_square_quadratic_density",
            "density_template": "L_nil2 = mu_nil * (6^2 mod 12) * psi^2 / 2",
            "weight_mod_12": square_weight,
            "formal_euler_coefficient": square_weight,
            "formal_hessian_coefficient": square_weight,
            "nonzero_formal_variation": square_weight != 0,
            "unit_bearing_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_variational_installation": False,
        },
        {
            "receiver": "symbolic_L_total_nilradical_slot",
            "density_template": "L_total += lambda_nil * N_nil(6)",
            "weight_mod_12": linear_weight,
            "formal_euler_coefficient": 0,
            "formal_hessian_coefficient": 0,
            "nonzero_formal_variation": False,
            "unit_bearing_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_variational_installation": False,
        },
        {
            "receiver": "completed_strict_nilradical_action_variational_schema",
            "density_template": "named unit-bearing density plus nonproxy variational chain",
            "weight_mod_12": linear_weight,
            "formal_euler_coefficient": linear_weight,
            "formal_hessian_coefficient": linear_weight,
            "nonzero_formal_variation": True,
            "unit_bearing_measure": True,
            "strict_field_provenance": True,
            "boundary_integration_theorem": True,
            "nonproxy_continuum_lift": True,
            "accepted_action_variational_installation": False,
        },
    ]
    return {
        "modulus": MODULUS,
        "nilpotent": NILPOTENT,
        "nilpotent_square_mod_12": square_weight,
        "units": UNITS,
        "unit_images": unit_images,
        "unit_fixed": all(v == NILPOTENT for v in unit_images.values()),
        "field_symbols": FIELDS,
        "receiver_rows": receivers,
    }


def obligation_rows(receivers: list[dict[str, Any]]) -> list[dict[str, Any]]:
    current = [r for r in receivers if r["receiver"] != "completed_strict_nilradical_action_variational_schema"]
    return [
        {"obligation": "finite_nilradical_weight", "satisfied": all(r["weight_mod_12"] in [0, 6] for r in current), "evidence": "all rows use the computed nilpotent 6 or its square 0 mod 12"},
        {"obligation": "nonzero_formal_variation", "satisfied": any(r["nonzero_formal_variation"] for r in current), "evidence": "the linear quadratic toy receiver has a nonzero formal Euler coefficient 6"},
        {"obligation": "unit_bearing_measure", "satisfied": any(r["unit_bearing_measure"] for r in current), "evidence": "mu_nil is symbolic and no strict measure/unit theorem is exported"},
        {"obligation": "strict_field_provenance", "satisfied": any(r["strict_field_provenance"] for r in current), "evidence": "psi is a formal placeholder, not a strict field sourced by the nilradical object"},
        {"obligation": "boundary_integration_theorem", "satisfied": any(r["boundary_integration_theorem"] for r in current), "evidence": "no boundary, integration, or variational domain theorem is attached"},
        {"obligation": "nonproxy_continuum_lift", "satisfied": any(r["nonproxy_continuum_lift"] for r in current), "evidence": "finite Z12 weights are not lifted to a tensor-resolved nonproxy continuum action"},
        {"obligation": "accepted_action_variational_installation", "satisfied": any(r["accepted_action_variational_installation"] for r in current), "evidence": "formal density receivers fail strict installation premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_weight", "nonzero_variation", "unit_measure", "field_provenance", "boundary_integration", "nonproxy_lift", "named_density"]
    return [{"present": dict(zip(names, bits)), "accepts_action_variational_installation": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2978_path: Any) -> dict[str, Any]:
    witness = formal_variational_receivers()
    receivers = witness["receiver_rows"]
    obligations = obligation_rows(receivers)
    matrix = acceptance_matrix()
    return {
        "status": "P2979_NILRADICAL_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2978": hashlib.sha256(p2978_path.read_bytes()).hexdigest() if p2978_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "NilradicalActionVariationalInstallation_ObstructionMatrix",
            "formal_variational_witness": witness,
            "receiver_rows": receivers,
            "proof_obligation_rows": obligations,
            "finite_acceptance_matrix": matrix,
        },
        "installation_certificate": {
            "receiver_count": len(receivers),
            "unit_fixed": witness["unit_fixed"],
            "nilpotent_square_mod_12": witness["nilpotent_square_mod_12"],
            "nonzero_formal_variation_rows": [r["receiver"] for r in receivers if r["nonzero_formal_variation"] and r["receiver"] != "completed_strict_nilradical_action_variational_schema"],
            "accepted_current_installations": [r["receiver"] for r in receivers if r["accepted_action_variational_installation"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_action_variational_installation"]),
        },
        "decision": {
            "positive_progress": "P2979 attacks the final nilradical-lane route: action/variational installation of the unit-fixed nilpotent anchor 6.",
            "breakthrough": "Bounded no-go: a formal quadratic receiver can produce a nonzero toy Euler coefficient 6, while the nilpotent-square receiver vanishes, but no current row exports a unit-bearing measure, strict field provenance, boundary/integration theorem, nonproxy continuum lift, or named density theorem.",
            "negative_export_flags": {k: False for k in ["nilradical_action_installation_exported", "unit_bearing_density_exported", "nonproxy_variational_chain_exported", "strict_nilradical_source_exported", "damping_source_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay nilradical-lane provenance, named source-atom coupling, formal quadratic receivers, nilpotent-square zero rows, symbolic L_total slots, ratio/Gamma/incidence lanes, selector closure, or bridge maps.  The nilradical lane is now bounded no-go on current artifacts; the next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside this lane, or preserve the P2929-P2979 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["installation_certificate"]
    lines = [
        "# P2979/S1929 nilradical action/variational installation obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Installation certificate",
        f"- receiver rows: `{cert['receiver_count']}`",
        f"- unit-fixed nilpotent: `{cert['unit_fixed']}`",
        f"- nilpotent square mod 12: `{cert['nilpotent_square_mod_12']}`",
        f"- nonzero formal variation rows: `{cert['nonzero_formal_variation_rows']}`",
        f"- accepted current installations: `{cert['accepted_current_installations']}`",
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
    read_json(P2978)
    payload = build_payload(P2978)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2979/S1929 nilradical action/variational installation obstruction", "## P2979/S1929 nilradical action/variational installation obstruction\n\n`P2979/S1929` attacks the final nilradical-lane route left by P2978: action/variational installation of the unit-fixed nilpotent anchor `6`.  The finite matrix builds formal receivers `L_nil = mu_nil*6*psi^2/2`, `L_nil2 = mu_nil*(6^2 mod 12)*psi^2/2`, and a symbolic `L_total` slot.  The strongest toy receiver has nonzero formal Euler/Hessian coefficient `6`, while the nilpotent-square receiver vanishes because `6^2=0`; however, no row exports a unit-bearing measure, strict field provenance, boundary/integration theorem, nonproxy continuum lift, or named density theorem.  No nilradical action installation, nonproxy variational chain, damping source, bridge closure, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2979/S1929 nilradical action/variational installation `L_total` guard", "## P2979/S1929 nilradical action/variational installation `L_total` guard\n\n`P2979/S1929` explicitly tests formal nilradical action receivers but adds no sourced term to `L_total`.  The nonzero toy Euler coefficient `6` is formal only, the nilpotent-square receiver is zero, and no unit-bearing measure, strict field provenance, boundary/integration theorem, nonproxy lift, EOM, Hamiltonian, role transfer, bridge closure, or ToE is exported.\n")
    append_once(AGENTS, "Current nilradical action/variational installation obstruction guardrail (P2979/S1929, 2026-06-20)", "## Current nilradical action/variational installation obstruction guardrail (P2979/S1929, 2026-06-20)\n\n- P2979 attacks the final nilradical-lane route: action/variational installation of the unit-fixed nilpotent anchor `6`.\n- Formal receiver positives are limited: `L_nil = mu_nil*6*psi^2/2` has toy Euler/Hessian coefficient `6`, while the nilpotent-square receiver vanishes because `6^2=0`; the acceptance matrix has `128` profiles with only the full profile accepting.\n- The current route is bounded no-go: no unit-bearing measure, strict field provenance, boundary/integration theorem, nonproxy continuum lift, or named density theorem is exported.\n- Do not promote nilradical formal receivers, toy Euler coefficients, nilpotent-square zero rows, symbolic `L_total` slots, provenance/coupling replay, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  The nilradical lane is now bounded no-go on current artifacts; a next move must introduce a genuinely new strict typed object/theorem/provider outside this lane or preserve the P2929-P2979 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
