#!/usr/bin/env python3
"""P3023/S1973: kernel-dissipation time-order candidate obstruction.

After P3022, introduce one genuinely new typed object outside the internal T_K
unit-normalization lane: a directed time-order candidate from the monotone
dissipation profile of K_strict_gate(d) on d=1..12.

The finite object is real: the sampled strict kernel strictly decreases along
1 -> 2 -> ... -> 12, so it defines an acyclic directed chain with one cyclic
reset obstruction.  The bounded obstruction is also finite: the chain depends on
the integer label chart, is not U(12)-equivariant, and has no strict physical
unit theorem for the step d -> d+1.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3022_s1972_unit_atom_reconciliation_no_new_live_frontier import OUT as P3022

OUT = GEN / "p3023_s1973_kernel_dissipation_time_order_candidate_obstruction.json"
MD = GEN / "p3023_s1973_kernel_dissipation_time_order_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = [1, 5, 7, 11]
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8


def k_strict(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA * (float(d) ** ETA))


def unit_image_edge(edge: tuple[int, int], unit: int) -> tuple[int, int]:
    a, b = edge
    return (((unit * a - 1) % N) + 1, ((unit * b - 1) % N) + 1)


def build_time_order_matrix() -> dict[str, Any]:
    values = {d: k_strict(d) for d in range(1, N + 1)}
    chain_edges = [(d, d + 1) for d in range(1, N)]
    cyclic_edge = (N, 1)
    descent_rows = []
    for a, b in chain_edges + [cyclic_edge]:
        delta = values[b] - values[a]
        descent_rows.append({
            "edge": [a, b],
            "K_source": round(values[a], 15),
            "K_target": round(values[b], 15),
            "delta_K": round(delta, 15),
            "strict_descent": delta < 0.0,
            "cyclic_reset": (a, b) == cyclic_edge,
        })
    chain_edge_set = set(chain_edges)
    equivariance_rows = []
    for unit in UNITS:
        images = [unit_image_edge(edge, unit) for edge in chain_edges]
        preserved = sum(1 for edge in images if edge in chain_edge_set)
        equivariance_rows.append({
            "unit": unit,
            "preserved_chain_edges": preserved,
            "total_chain_edges": len(chain_edges),
            "chain_equivariant": preserved == len(chain_edges),
            "sample_images": [list(edge) for edge in images[:4]],
        })
    obligations = [
        {"obligation": "new_typed_time_order_object", "satisfied": True, "detail": "directed chain from monotone K_strict_gate dissipation, not a T_K unit normalization"},
        {"obligation": "strict_descent_chain", "satisfied": all(row["strict_descent"] for row in descent_rows if not row["cyclic_reset"]), "detail": "K(d+1)<K(d) for d=1..11"},
        {"obligation": "no_cyclic_reset_obstruction", "satisfied": descent_rows[-1]["strict_descent"], "detail": "the edge 12->1 is an ascent/reset"},
        {"obligation": "U12_equivariant_directed_successor", "satisfied": all(row["chain_equivariant"] for row in equivariance_rows), "detail": "unit relabeling does not preserve the integer-label chain"},
        {"obligation": "strict_chart_or_selector_source", "satisfied": False, "detail": "no strict theorem selects the integer label chart/order as physical time"},
        {"obligation": "physical_unit_theorem", "satisfied": False, "detail": "no physical tick/action/Hamiltonian unit is attached to d->d+1"},
    ]
    return {
        "object": "KernelDissipationTimeOrderCandidate_EquivarianceUnitObstructionMatrix",
        "kernel": "K_strict_gate(d)=cos(omega*d+phi)/(1+beta*d^eta)",
        "directed_successor_candidate": "d -> d+1 on labels 1..12, stopping before cyclic reset",
        "sampled_kernel_values": {str(d): round(values[d], 15) for d in range(1, N + 1)},
        "descent_rows": descent_rows,
        "unit_equivariance_rows": equivariance_rows,
        "proof_obligations": obligations,
        "accepted_as_strict_time_order_with_physical_unit": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3022_path: Any) -> dict[str, Any]:
    read_json(p3022_path)
    matrix = build_time_order_matrix()
    return {
        "status": "P3023_KERNEL_DISSIPATION_TIME_ORDER_CANDIDATE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3022": hashlib.sha256(p3022_path.read_bytes()).hexdigest() if p3022_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "chain_edge_count": N - 1,
            "strict_descent_chain_edges": sum(1 for row in matrix["descent_rows"] if row["strict_descent"] and not row["cyclic_reset"]),
            "cyclic_reset_strict_descent": matrix["descent_rows"][-1]["strict_descent"],
            "unit_equivariant_rows": sum(1 for row in matrix["unit_equivariance_rows"] if row["chain_equivariant"]),
            "unit_row_count": len(matrix["unit_equivariance_rows"]),
            "accepted_as_strict_time_order_with_physical_unit": matrix["accepted_as_strict_time_order_with_physical_unit"],
        },
        "decision": {
            "breakthrough": "A new typed time-order candidate was constructed from K_strict_gate dissipation: K(d) strictly decreases along the finite chain 1->2->...->12.  The obstruction is that this chain is chart-dependent, has a cyclic reset at 12->1, is not U(12)-equivariant, and has no strict physical unit theorem for the directed step.",
            "negative_export_flags": {k: False for k in ["strict_time_order_object_exported", "physical_unit_theorem_exported", "unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_arrow_exported", "observed_physics_exported", "ltotal_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote kernel-dissipation label chains to physical time.  A next proof-grade move may attack exactly one missing theorem for this object: a strict chart/selector source for the integer order or an independent physical tick theorem; otherwise pivot to a different genuinely new typed object while preserving the P3017-P3023 no-closure boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3023/S1973 kernel-dissipation time-order candidate obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- strict descent chain edges / total: `{c['strict_descent_chain_edges']}/{c['chain_edge_count']}`",
        f"- cyclic reset is strict descent: `{c['cyclic_reset_strict_descent']}`",
        f"- U(12)-equivariant rows / total: `{c['unit_equivariant_rows']}/{c['unit_row_count']}`",
        f"- accepted as strict time-order with physical unit: `{c['accepted_as_strict_time_order_with_physical_unit']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3022)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3023/S1973 kernel-dissipation time-order candidate obstruction", "## P3023/S1973 kernel-dissipation time-order candidate obstruction\n\n`P3023/S1973` introduces one new typed object after the P3022 no-new-live-frontier certificate: a kernel-dissipation directed time-order candidate from sampled `K_strict_gate(d)` on labels `1..12`.  The finite positive is real: `K(d+1)<K(d)` for all `11/11` chain edges `1->2->...->12`.  The bounded obstruction is also explicit: the cyclic edge `12->1` is a reset/ascent, only the identity `U(12)` unit preserves the full directed chain, no strict chart/selector source exports this integer order as physical time, and no physical tick/action/Hamiltonian unit theorem is attached.  No strict time-order object with physical unit, time arrow, unit-bearing EOM/Hamiltonian, observed-physics export, `L_total`, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3023/S1973 kernel-dissipation time-order `L_total` guard", "## P3023/S1973 kernel-dissipation time-order `L_total` guard\n\n`P3023/S1973` adds no physical `L_total` term.  Its monotone `K_strict_gate` descent chain is a finite directed-order scaffold only; the chain is chart-dependent, not nontrivially `U(12)`-equivariant, has a cyclic reset, and lacks a physical tick/action/Hamiltonian unit theorem.\n")
    append_once(AGENTS, "Current kernel-dissipation time-order candidate guardrail (P3023/S1973, 2026-06-22)", "## Current kernel-dissipation time-order candidate guardrail (P3023/S1973, 2026-06-22)\n\n- P3023 introduces a new typed time-order candidate outside internal `T_K` unit normalizations: the monotone `K_strict_gate` dissipation chain on labels `1..12`.\n- The finite positive is exact on the sampled chain (`11/11` strict descent edges), but the object remains blocked: `12->1` is a cyclic reset/ascent, the directed chain is not nontrivially `U(12)`-equivariant, and no strict chart/selector source or physical tick/action/Hamiltonian unit theorem is exported.\n- Do not promote kernel-dissipation label chains to time arrow, strict time-order with physical unit, unit-bearing EOM/Hamiltonian, observed-physics, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- A next move may attack exactly one missing theorem for this object (strict chart/selector source or independent physical tick theorem), or pivot to a different genuinely new typed object while preserving the P3017-P3023 no-closure boundary.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
