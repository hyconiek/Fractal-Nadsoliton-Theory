#!/usr/bin/env python3
"""P2999/S1949: zero-divisor graph named source-atom coupling obstruction.

P2998 left two zero-divisor-graph routes.  This audit attacks exactly one:
a named source-atom coupling theorem for the Z/12Z zero-divisor graph.  It does
not replay strict provenance, nonpremise localizer, action installation,
Fourier, annihilator, nilradical, CRT, zero-derivation, selector, bridge,
role-transfer, or L_total lanes.

The finite calculation is exact: seven vertex receivers plus eight edge
receivers are crossed with four named source atoms, giving 60 receiver tests.
Each test carries graph-internal degree/orbit or zero-product data.  The
obstruction is theorem-side: no current artifact exports strict provenance,
nonpremise vertex/edge localizer, atom-specific coupling law, unit-bearing
coupling coefficient, or nonproxy export required to make these receivers a
strict source-coupling theorem.
"""
from __future__ import annotations

import hashlib, json
from itertools import product
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p2998_s1948_zero_divisor_graph_nonpremise_vertex_edge_localizer_obstruction import OUT as P2998, localizer_witness

OUT = GEN / "p2999_s1949_zero_divisor_graph_named_source_atom_coupling_obstruction.json"
MD = GEN / "p2999_s1949_zero_divisor_graph_named_source_atom_coupling_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

NAMED_SOURCE_ATOMS = [
    "selector_orientation_sign",
    "target_independent_positive_beta_Z_beta",
    "legacy_to_strict_bridge_source",
    "unit_bearing_action_density_source",
]


def graph_receiver_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for row in witness["vertex_rows"]:
        rows.append({
            "receiver_type": "vertex",
            "receiver_id": row["vertex"],
            "finite_receiver": {
                "degree": row["degree"],
                "neighbors": row["neighbors"],
                "neighbor_degree_signature": row["neighbor_degree_signature"],
                "automorphism_orbit": row["automorphism_orbit"],
                "orbit_size": row["orbit_size"],
                "singleton_graph_signature": row["singleton_graph_signature"],
            },
        })
    for row in witness["edge_rows"]:
        rows.append({
            "receiver_type": "edge",
            "receiver_id": row["edge"],
            "finite_receiver": {
                "zero_product_edge": True,
                "endpoint_degrees": row["endpoint_degrees"],
                "shared_neighbor_count": row["shared_neighbor_count"],
                "shared_neighbors": row["shared_neighbors"],
                "automorphism_orbit": row["automorphism_orbit"],
                "orbit_size": row["orbit_size"],
                "singleton_graph_signature": row["singleton_graph_signature"],
            },
        })
    return rows


def coupling_row(receiver: dict[str, Any], atom: str) -> dict[str, Any]:
    return {
        "receiver_type": receiver["receiver_type"],
        "receiver_id": receiver["receiver_id"],
        "named_source_atom": atom,
        "finite_graph_receiver": receiver["finite_receiver"],
        "zero_product_incidence_receiver_available": receiver["receiver_type"] == "edge",
        "orbit_signature_receiver_available": True,
        "strict_graph_provenance_available": False,
        "nonpremise_vertex_edge_localizer_available": False,
        "atom_specific_coupling_theorem": False,
        "unit_bearing_coupling_coefficient": False,
        "nonproxy_export_available": False,
        "accepted_source_coupling": False,
    }


def coupling_witness() -> dict[str, Any]:
    witness = localizer_witness()
    receivers = graph_receiver_rows(witness)
    rows = [coupling_row(receiver, atom) for receiver in receivers for atom in NAMED_SOURCE_ATOMS]
    return {
        "vertex_receiver_count": witness["vertex_count"],
        "edge_receiver_count": witness["edge_count"],
        "receiver_count": len(receivers),
        "named_source_atoms": NAMED_SOURCE_ATOMS,
        "coupling_test_count": len(rows),
        "coupling_rows": rows,
        "all_rows_have_exact_graph_receivers": all(r["orbit_signature_receiver_available"] for r in rows),
        "accepted_source_couplings": [
            {"receiver_type": r["receiver_type"], "receiver_id": r["receiver_id"], "atom": r["named_source_atom"]}
            for r in rows if r["accepted_source_coupling"]
        ],
    }


def obligation_rows(witness: dict[str, Any]) -> list[dict[str, Any]]:
    rows = witness["coupling_rows"]
    return [
        {"obligation": "finite_graph_atom_matrix_present", "satisfied": witness["coupling_test_count"] == (witness["vertex_receiver_count"] + witness["edge_receiver_count"]) * len(NAMED_SOURCE_ATOMS), "evidence": "7 vertex receivers plus 8 edge receivers crossed with 4 named source atoms gives 60 tests"},
        {"obligation": "exact_graph_receivers_present", "satisfied": witness["all_rows_have_exact_graph_receivers"], "evidence": "each test retains finite zero-product, degree, and automorphism-orbit receiver data"},
        {"obligation": "strict_graph_provenance_available", "satisfied": any(r["strict_graph_provenance_available"] for r in rows), "evidence": "P2997 found no strict graph provenance"},
        {"obligation": "nonpremise_vertex_edge_localizer_available", "satisfied": any(r["nonpremise_vertex_edge_localizer_available"] for r in rows), "evidence": "P2998 found no accepted nonpremise vertex/edge localizer"},
        {"obligation": "atom_specific_coupling_theorem", "satisfied": any(r["atom_specific_coupling_theorem"] for r in rows), "evidence": "no theorem couples a zero-divisor graph receiver to selector sign, beta/Z_beta, bridge-source, or action-density source atoms"},
        {"obligation": "unit_bearing_coupling_coefficient", "satisfied": any(r["unit_bearing_coupling_coefficient"] for r in rows), "evidence": "finite graph receivers carry no unit-bearing coefficient"},
        {"obligation": "nonproxy_export_available", "satisfied": any(r["nonproxy_export_available"] for r in rows), "evidence": "no nonproxy export to source law, action, EOM, or continuum lift is present"},
        {"obligation": "accepted_current_source_coupling", "satisfied": bool(witness["accepted_source_couplings"]), "evidence": "no current graph/atom pair satisfies the full source-coupling profile"},
    ]


def acceptance_matrix() -> list[dict[str, Any]]:
    names = ["finite_matrix", "exact_graph_receivers", "strict_provenance", "nonpremise_localizer", "atom_coupling_theorem", "unit_coefficient", "nonproxy_export"]
    return [{"present": dict(zip(names, bits)), "accepts_source_coupling_theorem": all(bits)} for bits in product([False, True], repeat=len(names))]


def build_payload(p2998_path: Any) -> dict[str, Any]:
    witness = coupling_witness()
    matrix = acceptance_matrix()
    return {
        "status": "P2999_ZERO_DIVISOR_GRAPH_NAMED_SOURCE_ATOM_COUPLING_OBSTRUCTION_BOUNDED_NO_GO",
        "input_hashes": {"P2998": hashlib.sha256(p2998_path.read_bytes()).hexdigest() if p2998_path.exists() else None},
        "constructed_theoretical_objects": {
            "object": "ZeroDivisorGraphNamedSourceAtomCoupling_ObstructionMatrix",
            "coupling_witness": witness,
            "proof_obligation_rows": obligation_rows(witness),
            "finite_acceptance_matrix": matrix,
        },
        "coupling_certificate": {
            "vertex_receiver_count": witness["vertex_receiver_count"],
            "edge_receiver_count": witness["edge_receiver_count"],
            "receiver_count": witness["receiver_count"],
            "named_source_atom_count": len(witness["named_source_atoms"]),
            "coupling_test_count": witness["coupling_test_count"],
            "all_rows_have_exact_graph_receivers": witness["all_rows_have_exact_graph_receivers"],
            "accepted_source_couplings": witness["accepted_source_couplings"],
            "acceptance_matrix_rows": len(matrix),
            "accepted_rows": sum(1 for r in matrix if r["accepts_source_coupling_theorem"]),
        },
        "decision": {
            "positive_progress": "P2999 attacks exactly one remaining zero-divisor-graph route after P2998: named source-atom coupling for vertex/edge graph receivers.",
            "breakthrough": "Bounded no-go: all 60 graph/atom tests have exact finite graph receivers, but no strict graph provenance, accepted nonpremise localizer, atom-specific coupling theorem, unit-bearing coefficient, or nonproxy export is available.",
            "negative_export_flags": {k: False for k in ["source_atom_coupling_exported", "strict_graph_provenance_exported", "nonpremise_vertex_edge_localizer_exported", "selector_closure_exported", "positive_beta_source_exported", "bridge_closure_exported", "unit_bearing_action_density_exported", "nonproxy_ltotal_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "The only remaining zero-divisor-graph route is unit-bearing action installation.  Attack it only as a bounded audit requiring a genuinely unit-bearing graph measure, named graph density theorem, boundary/integration map, and nonproxy continuum lift; otherwise introduce a genuinely new strict typed object/provider or preserve the P2929-P2999 no-strict-export certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    cert = payload["coupling_certificate"]
    lines = [
        "# P2999/S1949 zero-divisor graph named source-atom coupling obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Coupling certificate",
        f"- vertex/edge receivers: `{cert['vertex_receiver_count']}/{cert['edge_receiver_count']}`",
        f"- total receivers: `{cert['receiver_count']}`",
        f"- named source atom count: `{cert['named_source_atom_count']}`",
        f"- coupling tests: `{cert['coupling_test_count']}`",
        f"- all rows have exact graph receivers: `{cert['all_rows_have_exact_graph_receivers']}`",
        f"- accepted source couplings: `{cert['accepted_source_couplings']}`",
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
    read_json(P2998)
    payload = build_payload(P2998)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2999/S1949 zero-divisor graph named source-atom coupling obstruction", "## P2999/S1949 zero-divisor graph named source-atom coupling obstruction\n\n`P2999/S1949` attacks exactly one remaining zero-divisor-graph route after P2998: named source-atom coupling for vertex/edge graph receivers.  The finite receiver matrix is exact: seven vertex receivers plus eight edge receivers crossed with four named atoms (`selector_orientation_sign`, `target_independent_positive_beta_Z_beta`, `legacy_to_strict_bridge_source`, `unit_bearing_action_density_source`) gives `60` graph/atom tests with zero-product, degree, and automorphism-orbit receiver data.  The theorem side remains blocked: P2997 provides no strict graph provenance, P2998 provides no accepted nonpremise vertex/edge localizer, and no atom-specific coupling theorem, unit-bearing coefficient, or nonproxy export is present.  No selector closure, positive beta source, bridge closure, unit-bearing action density, nonproxy `L_total`, role transfer, or ToE follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2999/S1949 zero-divisor graph source-coupling `L_total` guard", "## P2999/S1949 zero-divisor graph source-coupling `L_total` guard\n\n`P2999/S1949` adds no zero-divisor graph source-coupling term to `L_total`.  The `60` finite graph/atom receivers are algebraic incidence bookkeeping only; they do not supply strict field provenance, unit-bearing coefficient, named density theorem, variational chain, EOM/Hamiltonian term, bridge closure, role transfer, or ToE.\n")
    append_once(AGENTS, "Current zero-divisor graph named source-atom coupling obstruction guardrail (P2999/S1949, 2026-06-20)", "## Current zero-divisor graph named source-atom coupling obstruction guardrail (P2999/S1949, 2026-06-20)\n\n- P2999 attacks exactly one remaining zero-divisor-graph route after P2998: named source-atom coupling for vertex/edge graph receivers.\n- Finite positives are exact but receiver-only: seven vertex receivers plus eight edge receivers crossed with four named source atoms gives `60` graph/atom tests with zero-product, degree, and automorphism-orbit data.\n- The route is bounded no-go because P2997 exported no strict graph provenance, P2998 exported no accepted nonpremise vertex/edge localizer, and no atom-specific coupling theorem, unit-bearing coefficient, or nonproxy export is present.\n- Do not promote graph/atom receivers to selector closure, target-independent positive `beta/Z_beta`, bridge closure, unit-bearing action density, role transfer, nonproxy `L_total`, or ToE.\n- The only remaining zero-divisor-graph route is unit-bearing action installation with a genuinely unit-bearing graph measure, named graph density theorem, boundary/integration map, and nonproxy continuum lift; otherwise introduce a genuinely new strict typed object/provider while preserving the P2929-P2999 no-strict-export certificate.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
