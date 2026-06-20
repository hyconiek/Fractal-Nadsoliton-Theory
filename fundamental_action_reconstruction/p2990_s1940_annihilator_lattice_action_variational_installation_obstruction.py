#!/usr/bin/env python3
"""P2990/S1940: annihilator-lattice action/variational installation obstruction.

P2989 left the final annihilator-lattice route: action installation with a
genuinely unit-bearing measure, named density theorem, boundary/integration map,
and nonproxy continuum lift.  This audit attacks that route only.  It does not
replay strict provenance, nonpremise source-localizer, named source-atom
coupling, nilradical/CRT/zero-derivation lanes, selector closure, bridge
completion, role transfer, or L_total promotion.

The finite calculation constructs the strongest formal toy receivers available
from the six Z/12Z ideal-annihilator rows: one quadratic receiver per row using
annihilator size as formal weight, a split-sum density, and a zero-product
constraint receiver using the ideal-annihilator product table.  These receivers
produce exact nonzero toy Euler/Hessian coefficients and a finite zero-product
constraint count.  They still are not strict action terms: no current artifact
exports a unit-bearing measure, strict field provenance, boundary/integration
theorem, named density theorem, or nonproxy continuum lift.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2986_s1936_z12_annihilator_lattice_source_candidate_obstruction import annihilator_lattice_witness
from p2989_s1939_annihilator_lattice_named_source_atom_coupling_obstruction import OUT as P2989

OUT = GEN / "p2990_s1940_annihilator_lattice_action_variational_installation_obstruction.json"
MD = GEN / "p2990_s1940_annihilator_lattice_action_variational_installation_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"


def row_label(index: int, row: dict[str, Any]) -> str:
    return f"I{index}_size{row['size']}_ann{row['annihilator_size']}"


def formal_action_receivers() -> dict[str, Any]:
    lattice = annihilator_lattice_witness()
    ideal_rows = lattice["annihilator_rows"]
    row_receivers = []
    hessian_weights = []
    zero_product_total = 0
    for index, row in enumerate(ideal_rows):
        label = row_label(index, row)
        field = f"psi_{label}"
        weight = row["annihilator_size"]
        zero_pairs = row["size"] * row["annihilator_size"]
        zero_product_total += zero_pairs
        hessian_weights.append(weight)
        row_receivers.append({
            "receiver": f"formal_annihilator_row_quadratic_density_{label}",
            "density_template": f"L_{label} = mu_{label} * {weight} * {field}^2 / 2",
            "ideal": row["ideal"],
            "annihilator": row["annihilator"],
            "field_symbol": field,
            "weight": weight,
            "zero_product_pairs": zero_pairs,
            "formal_euler_coefficient": weight,
            "formal_hessian_coefficient": weight,
            "nonzero_formal_variation": weight != 0,
            "unit_bearing_measure": False,
            "strict_field_provenance": False,
            "boundary_integration_theorem": False,
            "named_density_theorem": False,
            "nonproxy_continuum_lift": False,
            "accepted_action_installation": False,
        })
    split_sum = {
        "receiver": "formal_annihilator_lattice_split_sum_density",
        "density_template": "L_ann = sum_i mu_i*|Ann(I_i)|*psi_i^2/2",
        "field_symbols": [r["field_symbol"] for r in row_receivers],
        "formal_hessian_diagonal": hessian_weights,
        "formal_hessian_rank_over_Q": sum(1 for w in hessian_weights if w != 0),
        "nonzero_formal_variation": any(w != 0 for w in hessian_weights),
        "unit_bearing_measure": False,
        "strict_field_provenance": False,
        "boundary_integration_theorem": False,
        "named_density_theorem": False,
        "nonproxy_continuum_lift": False,
        "accepted_action_installation": False,
    }
    zero_constraint = {
        "receiver": "formal_ideal_annihilator_zero_product_constraint_density",
        "density_template": "C_ann = sum_{x in I, y in Ann(I)} lambda_I * ((x*y) mod 12)^2",
        "zero_product_pair_count": zero_product_total,
        "all_products_zero": all(row["annihilator_product_zero"] for row in ideal_rows),
        "formal_constraint_residual": 0,
        "nonzero_formal_variation": False,
        "unit_bearing_measure": False,
        "strict_field_provenance": False,
        "boundary_integration_theorem": False,
        "named_density_theorem": False,
        "nonproxy_continuum_lift": False,
        "accepted_action_installation": False,
    }
    completed_schema = {
        "receiver": "completed_strict_annihilator_action_installation_schema",
        "density_template": "named unit-bearing annihilator-lattice density plus boundary/integration and nonproxy lift",
        "required_positive_inputs": ["unit_bearing_measure", "strict_field_provenance", "boundary_integration_theorem", "named_density_theorem", "nonproxy_continuum_lift"],
        "nonzero_formal_variation": True,
        "unit_bearing_measure": True,
        "strict_field_provenance": True,
        "boundary_integration_theorem": True,
        "named_density_theorem": True,
        "nonproxy_continuum_lift": True,
        "accepted_action_installation": False,
    }
    receivers = row_receivers + [split_sum, zero_constraint, completed_schema]
    return {
        "row_count": len(ideal_rows),
        "row_weights": hessian_weights,
        "zero_product_total": zero_product_total,
        "formal_hessian_rank_over_Q": split_sum["formal_hessian_rank_over_Q"],
        "receiver_rows": receivers,
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    current = [r for r in witness["receiver_rows"] if r["receiver"] != "completed_strict_annihilator_action_installation_schema"]
    return [
        {"obligation": "finite_annihilator_action_receivers", "satisfied": witness["row_count"] == 6, "evidence": "one formal quadratic receiver is built for each of the six ideal-annihilator rows"},
        {"obligation": "nonzero_formal_variation", "satisfied": any(r["nonzero_formal_variation"] for r in current), "evidence": f"row weights {witness['row_weights']} give nonzero toy Euler/Hessian coefficients"},
        {"obligation": "zero_product_constraint_witness", "satisfied": any(r.get("all_products_zero") and r.get("formal_constraint_residual") == 0 for r in current), "evidence": f"all ideal-annihilator products vanish across {witness['zero_product_total']} finite pairs"},
        {"obligation": "unit_bearing_measure", "satisfied": any(r["unit_bearing_measure"] for r in current), "evidence": "mu_i and lambda_I are symbolic; no strict unit/measure theorem is exported"},
        {"obligation": "strict_field_provenance", "satisfied": any(r["strict_field_provenance"] for r in current), "evidence": "psi_I fields are formal placeholders, not strict fields sourced by annihilator rows"},
        {"obligation": "boundary_integration_theorem", "satisfied": any(r["boundary_integration_theorem"] for r in current), "evidence": "finite row sums are not a boundary/integration theorem"},
        {"obligation": "named_density_theorem", "satisfied": any(r["named_density_theorem"] for r in current), "evidence": "formal receiver names are not exported named action densities"},
        {"obligation": "nonproxy_continuum_lift", "satisfied": any(r["nonproxy_continuum_lift"] for r in current), "evidence": "finite annihilator weights are not lifted to tensor-resolved nonproxy continuum action"},
        {"obligation": "accepted_action_installation", "satisfied": any(r["accepted_action_installation"] for r in current), "evidence": "formal annihilator receivers fail strict installation premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_receivers", "nonzero_variation", "zero_product_constraint", "unit_measure", "field_provenance", "boundary_integration", "named_density", "nonproxy_lift"]
    return [{"present": dict(zip(names, bits)), "accepts_annihilator_action_installation": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2989_path: Any) -> dict[str, Any]:
    witness = formal_action_receivers()
    matrix = acceptance_matrix()
    return {
        "status": "P2990_ANNIHILATOR_LATTICE_ACTION_VARIATIONAL_INSTALLATION_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2989": hashlib.sha256(p2989_path.read_bytes()).hexdigest() if p2989_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "AnnihilatorLatticeActionVariationalInstallation_ObstructionMatrix",
            "formal_action_witness": witness,
            "receiver_rows": witness["receiver_rows"],
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "installation_certificate": {
            "receiver_count": len(witness["receiver_rows"]),
            "row_count": witness["row_count"],
            "row_weights": witness["row_weights"],
            "zero_product_total": witness["zero_product_total"],
            "formal_hessian_rank_over_Q": witness["formal_hessian_rank_over_Q"],
            "accepted_current_installations": [r["receiver"] for r in witness["receiver_rows"] if r["accepted_action_installation"]],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_annihilator_action_installation"]),
        },
        "decision": {
            "positive_progress": "P2990 attacks the final annihilator-lattice route after P2989: action/variational installation with unit-bearing measure and nonproxy lift obligations.",
            "breakthrough": "Bounded no-go: six row quadratic receivers and a zero-product constraint receiver give exact finite toy Euler/Hessian and constraint data, but no current row exports a unit-bearing measure, strict field provenance, boundary/integration theorem, named density theorem, or nonproxy continuum lift.",
            "negative_export_flags": {k: False for k in ["annihilator_action_installation_exported", "unit_bearing_density_exported", "nonproxy_variational_chain_exported", "strict_field_provenance_exported", "named_density_theorem_exported", "bridge_closure_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay annihilator formal receivers, zero-product constraints, named-atom rows, source-localizer/provenance rows, nilradical/CRT/zero-derivation lanes, selector replay, bridge maps, role transfer, or L_total placeholders.  The annihilator-lattice lane is now bounded no-go on current artifacts; the next proof-grade move must introduce a genuinely new strict typed object/theorem/provider outside this lane, or preserve the P2929-P2990 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["installation_certificate"]
    lines = [
        "# P2990/S1940 annihilator-lattice action/variational installation obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Installation certificate",
        f"- receiver rows: `{cert['receiver_count']}`",
        f"- row count: `{cert['row_count']}`",
        f"- row weights: `{cert['row_weights']}`",
        f"- zero-product finite pairs: `{cert['zero_product_total']}`",
        f"- formal Hessian rank over Q: `{cert['formal_hessian_rank_over_Q']}`",
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
    read_json(P2989)
    payload = build_payload(P2989)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2990/S1940 annihilator-lattice action/variational installation obstruction", "## P2990/S1940 annihilator-lattice action/variational installation obstruction\n\n`P2990/S1940` attacks the final annihilator-lattice route left by P2989: action/variational installation with a genuinely unit-bearing measure, named density theorem, boundary/integration map, and nonproxy continuum lift.  The finite matrix builds one formal quadratic receiver per ideal-annihilator row using annihilator size as weight, a split-sum density, and a zero-product constraint receiver.  The toy Hessian has exact row weights and nonzero rank, and the zero-product constraint has zero residual across the finite ideal-annihilator product table; however, no row exports a unit-bearing measure, strict field provenance, boundary/integration theorem, named density theorem, or nonproxy continuum lift.  No annihilator action installation, nonproxy variational chain, bridge closure, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2990/S1940 annihilator action/variational installation `L_total` guard", "## P2990/S1940 annihilator action/variational installation `L_total` guard\n\n`P2990/S1940` explicitly tests formal annihilator-lattice action receivers but adds no sourced annihilator term to `L_total`.  The nonzero toy Euler/Hessian coefficients and zero-product constraint are formal finite algebra only; no unit-bearing measure, strict field provenance, boundary/integration theorem, named density theorem, nonproxy lift, EOM, Hamiltonian, role transfer, bridge closure, or ToE is exported.\n")
    append_once(AGENTS, "Current annihilator-lattice action/variational installation obstruction guardrail (P2990/S1940, 2026-06-20)", "## Current annihilator-lattice action/variational installation obstruction guardrail (P2990/S1940, 2026-06-20)\n\n- P2990 attacks the final annihilator-lattice route left by P2989: action/variational installation with unit-bearing measure and nonproxy lift obligations.\n- Formal receiver positives are limited: one quadratic receiver per ideal-annihilator row, a split-sum density, and a zero-product constraint receiver give exact finite toy Euler/Hessian and constraint data; the acceptance matrix has `256` profiles with only the full profile accepting.\n- The current route is bounded no-go: no unit-bearing measure, strict field provenance, boundary/integration theorem, named density theorem, or nonproxy continuum lift is exported.\n- Do not promote annihilator formal receivers, zero-product constraints, named-atom rows, provenance/localizer rows, symbolic `L_total` slots, selector replay, bridge maps, role transfer, nonproxy `L_total`, or ToE.  The annihilator-lattice lane is now bounded no-go on current artifacts; a next move must introduce a genuinely new strict typed object/theorem/provider outside this lane or preserve the P2929-P2990 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
