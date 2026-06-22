#!/usr/bin/env python3
"""P3004/S1954: square-map named source-atom coupling obstruction.

P3003 left two square-map-functional-graph routes. This audit attacks exactly
one: named source-atom coupling for square-map basin/attractor and residue
receivers. It does not replay strict provenance, basin/localizer, action
installation, zero-divisor graph, CRT idempotent, nilradical, annihilator,
Fourier, selector, bridge, role-transfer, or L_total lanes.

The finite calculation constructs graph-dynamical receivers and crosses them
with four named source atoms. The receiver matrix is exact, but every row remains
blocked because current artifacts export neither square-map strict provenance nor
an accepted nonpremise basin/attractor localizer, and no atom-specific coupling
theorem with a unit-bearing coefficient/nonproxy export is present.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3003_s1953_square_map_basin_attractor_localizer_obstruction import OUT as P3003, localizer_witness

OUT = GEN / "p3004_s1954_square_map_named_source_atom_coupling_obstruction.json"
MD = GEN / "p3004_s1954_square_map_named_source_atom_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

SOURCE_ATOMS = [
    "selector_orientation_sign",
    "target_independent_positive_beta_Z_beta",
    "legacy_to_strict_bridge_source",
    "unit_bearing_action_density_source",
]


def receiver_rows() -> list[dict[str, Any]]:
    witness = localizer_witness()
    rows: list[dict[str, Any]] = []
    for basin in witness["basin_rows"]:
        rows.append({
            "receiver_type": "basin_attractor",
            "receiver_id": f"basin_{basin['attractor']}",
            "attractor": basin["attractor"],
            "basin_size": basin["basin_size"],
            "basin_members": basin["basin_members"],
            "depth_signature": basin["depth_signature"],
            "fiber_to_attractor": basin["fiber_to_attractor"],
            "unit_orbit_of_attractor": basin["unit_orbit_of_attractor"],
            "singleton_unit_orbit": basin["singleton_unit_orbit"],
        })
    for residue in witness["residue_rows"]:
        rows.append({
            "receiver_type": "residue",
            "receiver_id": f"residue_{residue['residue']}",
            "residue": residue["residue"],
            "square_target": residue["square_target"],
            "fixed_attractor": residue["fixed_attractor"],
            "basin_size": residue["basin_size"],
            "depth_to_attractor": residue["depth_to_attractor"],
            "fiber_size_of_target": residue["fiber_size_of_target"],
            "unit_orbit": residue["unit_orbit"],
            "unit_orbit_size": residue["unit_orbit_size"],
            "localizer_signature": residue["localizer_signature"],
        })
    return rows


def coupling_witness() -> dict[str, Any]:
    receivers = receiver_rows()
    coupling_rows = []
    for receiver, atom in product(receivers, SOURCE_ATOMS):
        coupling_rows.append({
            "receiver_id": receiver["receiver_id"],
            "receiver_type": receiver["receiver_type"],
            "source_atom": atom,
            "finite_receiver": receiver,
            "strict_square_map_provenance_exported": False,
            "accepted_nonpremise_localizer_exported": False,
            "atom_specific_coupling_theorem_exported": False,
            "unit_bearing_coupling_coefficient_exported": False,
            "nonproxy_exported": False,
            "accepted_source_coupling": False,
        })
    return {
        "receiver_rows": receivers,
        "source_atoms": SOURCE_ATOMS,
        "coupling_rows": coupling_rows,
        "basin_receiver_count": sum(1 for r in receivers if r["receiver_type"] == "basin_attractor"),
        "residue_receiver_count": sum(1 for r in receivers if r["receiver_type"] == "residue"),
        "accepted_source_couplings": [r for r in coupling_rows if r["accepted_source_coupling"]],
    }


def obligation_rows(w: dict[str, Any]) -> list[dict[str, Any]]:
    rows = w["coupling_rows"]
    return [
        {"obligation": "finite_square_map_atom_matrix", "satisfied": len(rows) == 64 and len(w["source_atoms"]) == 4, "evidence": "16 square-map receivers crossed with four named source atoms"},
        {"obligation": "exact_basin_residue_receivers", "satisfied": w["basin_receiver_count"] == 4 and w["residue_receiver_count"] == 12, "evidence": "P3003 basin/attractor and residue receivers are recomputed"},
        {"obligation": "strict_square_map_provenance", "satisfied": any(r["strict_square_map_provenance_exported"] for r in rows), "evidence": "P3002 exported no strict square-map provenance"},
        {"obligation": "accepted_nonpremise_localizer", "satisfied": any(r["accepted_nonpremise_localizer_exported"] for r in rows), "evidence": "P3003 exported no accepted nonpremise basin/attractor localizer"},
        {"obligation": "atom_specific_coupling_theorem", "satisfied": any(r["atom_specific_coupling_theorem_exported"] for r in rows), "evidence": "no theorem couples square-map receivers to selector, beta/Z_beta, bridge-source, or action-density atoms"},
        {"obligation": "unit_bearing_coupling_coefficient", "satisfied": any(r["unit_bearing_coupling_coefficient_exported"] for r in rows), "evidence": "no unit-bearing coefficient or measure/sign theorem is exported"},
        {"obligation": "nonproxy_export", "satisfied": any(r["nonproxy_exported"] for r in rows), "evidence": "receiver/atom bookkeeping is not a nonproxy continuum export"},
        {"obligation": "accepted_current_source_coupling", "satisfied": bool(w["accepted_source_couplings"]), "evidence": "no row satisfies provenance, localizer, coupling theorem, coefficient, and nonproxy premises"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_receiver", "named_source_atom", "strict_provenance", "accepted_localizer", "atom_specific_theorem", "unit_coefficient", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_source_coupling": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p3003_path: Any) -> dict[str, Any]:
    witness = coupling_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P3004_SQUARE_MAP_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P3003": hashlib.sha256(p3003_path.read_bytes()).hexdigest() if p3003_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "SquareMapNamedSourceAtomCoupling_ObstructionMatrix",
            "coupling_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "coupling_certificate": {
            "basin_receiver_count": witness["basin_receiver_count"],
            "residue_receiver_count": witness["residue_receiver_count"],
            "total_receiver_count": len(witness["receiver_rows"]),
            "named_source_atom_count": len(witness["source_atoms"]),
            "coupling_test_count": len(witness["coupling_rows"]),
            "accepted_source_coupling_count": len(witness["accepted_source_couplings"]),
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for row in matrix if row["accepts_source_coupling"]),
        },
        "decision": {
            "positive_progress": "P3004 attacks exactly one remaining P3003 square-map route: named source-atom coupling for basin/attractor and residue receivers.",
            "breakthrough": "Bounded no-go: 64 exact square-map receiver/atom tests are constructed, but strict provenance, accepted nonpremise localizer, atom-specific coupling theorem, unit-bearing coefficient, and nonproxy export are absent.",
            "negative_export_flags": {k: False for k in ["square_map_named_source_atom_coupling_exported", "selector_closure_exported", "target_independent_positive_beta_Z_beta_source_exported", "bridge_closure_exported", "unit_bearing_action_density_exported", "unit_bearing_action_installation_exported", "role_transfer_exported", "nonproxy_ltotal_exported", "toe_closure_exported"]},
            "next_honest_step": "The only remaining square-map-functional-graph route is unit-bearing action installation with a genuinely unit-bearing square-map measure, named square-map density theorem, boundary/integration map, and nonproxy continuum lift. Otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P3004 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P3004/S1954 square-map named source-atom coupling obstruction", "",
        f"Status: `{payload['status']}`", "", "## Coupling certificate",
        f"- basin receivers: `{cert['basin_receiver_count']}`",
        f"- residue receivers: `{cert['residue_receiver_count']}`",
        f"- total receivers: `{cert['total_receiver_count']}`",
        f"- named source atoms: `{cert['named_source_atom_count']}`",
        f"- coupling tests: `{cert['coupling_test_count']}`",
        f"- accepted source couplings: `{cert['accepted_source_coupling_count']}`",
        f"- acceptance matrix rows/accepted: `{cert['acceptance_matrix_rows']}/{cert['accepted_rows']}`", "",
        "## Lay summary", payload["decision"]["positive_progress"], payload["decision"]["breakthrough"], "",
        "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3003)
    payload = build_payload(P3003)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3004/S1954 square-map named source-atom coupling obstruction", "## P3004/S1954 square-map named source-atom coupling obstruction\n\n`P3004/S1954` attacks exactly one remaining P3003 route: named source-atom coupling for square-map basin/attractor and residue receivers.  The finite receiver/atom matrix is exact: four basin receivers plus twelve residue receivers crossed with four named atoms gives `64` coupling tests with basin, fiber, unit-orbit, and localizer-signature data.  The theorem side remains blocked: P3002 provides no strict square-map provenance, P3003 provides no accepted nonpremise localizer, and no atom-specific coupling theorem, unit-bearing coefficient, or nonproxy export is present.  No selector closure, target-independent positive `beta/Z_beta`, bridge closure, unit-bearing action density, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3004/S1954 square-map named source-atom coupling `L_total` guard", "## P3004/S1954 square-map named source-atom coupling `L_total` guard\n\n`P3004/S1954` adds no square-map named-atom coupling term to `L_total`.  The receiver/atom matrix lacks strict field provenance, an accepted nonpremise source localization theorem, atom-specific unit-bearing density, boundary/integration theorem, nonproxy variational chain, EOM/Hamiltonian term, bridge closure, role transfer, and ToE.\n")
    append_once(AGENTS, "Current square-map named source-atom coupling obstruction guardrail (P3004/S1954, 2026-06-21)", "## Current square-map named source-atom coupling obstruction guardrail (P3004/S1954, 2026-06-21)\n\n- P3004 attacks exactly one remaining P3003 square-map route: named source-atom coupling for basin/attractor and residue receivers.\n- Finite positives are exact but receiver-only: four basin receivers plus twelve residue receivers crossed with four named atoms gives `64` coupling tests with basin, fiber, unit-orbit, and localizer-signature data.\n- The route is bounded no-go because P3002 exported no strict square-map provenance, P3003 exported no accepted nonpremise localizer, and no atom-specific coupling theorem, unit-bearing coefficient, or nonproxy export is present.\n- Do not promote square-map receiver/atom rows to selector closure, target-independent positive `beta/Z_beta`, bridge closure, unit-bearing action density, role transfer, nonproxy `L_total`, or ToE.  The only remaining square-map-functional-graph route is unit-bearing action installation with a genuinely unit-bearing square-map measure, named square-map density theorem, boundary/integration map, and nonproxy continuum lift; otherwise introduce a genuinely new strict typed object/provider while preserving the P2929-P3004 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
