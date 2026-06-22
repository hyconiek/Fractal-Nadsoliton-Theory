#!/usr/bin/env python3
"""P3021/S1971: action-quantum/reference-cell source obstruction.

After P3020, attack the remaining independent unit atom for the typed action:
can the formal T_K action receive a strict action quantum/reference-cell source?

We build finite reference-cell candidates from the already typed quadratic action:
all lattice edges, nonzero-gradient support, Laplacian rank, and DFT/Parseval
mode cells.  Each yields a computable positive action-per-cell candidate.  The
obstruction is that these are partitions of the same dimensionless formal action:
under T_K -> c T_K all action quanta rescale by c^2, and no independent strict
reference cell, physical action unit, or Hamiltonian coupling theorem is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3018_s1968_lambda_action_unit_normalization_candidate_obstruction import t_obs, cyclic_gradient_energy
from p3020_s1970_clock_unit_theorem_candidate_obstruction import OUT as P3020, dft_magnitudes

OUT = GEN / "p3021_s1971_action_quantum_reference_cell_source_obstruction.json"
MD = GEN / "p3021_s1971_action_quantum_reference_cell_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
SCALE_FACTORS = [0.5, 1.0, 2.0, 3.0]


def action_cell_counts(q: list[float]) -> dict[str, int]:
    diffs = [q[(i + 1) % len(q)] - q[i] for i in range(len(q))]
    mags = dft_magnitudes(q)
    nonzero_modes = [k for k in range(1, N) if mags[k] > 1e-12]
    return {
        "all_cyclic_edge_cells": N,
        "nonzero_gradient_support_cells": sum(1 for x in diffs if abs(x) > 1e-12),
        "cyclic_laplacian_rank_cells": N - 1,
        "nonzero_dft_mode_cells": len(nonzero_modes),
    }


def build_action_quantum_matrix() -> dict[str, Any]:
    q = [t_obs(d) for d in range(1, N + 1)]
    base_energy = cyclic_gradient_energy(q)
    base_action_lambda_one = 0.5 * base_energy
    counts = action_cell_counts(q)
    cell_rows = []
    for name, count in counts.items():
        quantum = base_action_lambda_one / count
        cell_rows.append({
            "reference_cell_candidate": name,
            "cell_count": count,
            "action_per_cell_lambda_1": round(quantum, 15),
            "positive": count > 0 and quantum > 0.0 and math.isfinite(quantum),
            "independent_of_T_K_amplitude": False,
            "exports_physical_action_quantum": False,
        })
    scale_rows = []
    for c in SCALE_FACTORS:
        scaled_q = [c * x for x in q]
        scaled_action = 0.5 * cyclic_gradient_energy(scaled_q)
        scaled_counts = action_cell_counts(scaled_q)
        scale_rows.append({
            "observable_scale_c": c,
            "action_ratio_to_base": round(scaled_action / base_action_lambda_one, 12),
            "expected_c_square": round(c * c, 12),
            "cell_counts_same_as_base": scaled_counts == counts,
            "quanta_rescale_by_c_square": all(math.isclose((scaled_action / scaled_counts[name]) / (base_action_lambda_one / counts[name]), c * c, rel_tol=0.0, abs_tol=1e-12) for name in counts),
        })
    obligations = [
        {"obligation": "typed_formal_action_input", "satisfied": True, "detail": "uses P3017 quadratic action S=(1/2)sum(ΔT_K)^2 with lambda=1"},
        {"obligation": "positive_finite_reference_cell_candidates", "satisfied": all(row["positive"] for row in cell_rows), "detail": "edge, gradient-support, rank, and DFT-mode cells are positive"},
        {"obligation": "cell_counts_stable_under_observable_rescaling", "satisfied": all(row["cell_counts_same_as_base"] for row in scale_rows), "detail": "the combinatorial partitions are stable"},
        {"obligation": "action_quantum_invariant_under_observable_rescaling", "satisfied": not all(row["quanta_rescale_by_c_square"] for row in scale_rows), "detail": "all action-per-cell candidates rescale as c^2"},
        {"obligation": "strict_reference_cell_source_theorem", "satisfied": False, "detail": "no theorem selects one cell partition as a physical reference cell"},
        {"obligation": "physical_action_unit_and_hamiltonian_coupling", "satisfied": False, "detail": "no independent action quantum or Hamiltonian normalization is exported"},
    ]
    return {
        "object": "ActionQuantumReferenceCellSource_RescalingPartitionObstructionMatrix",
        "typed_action": "S_formal[T_K]=(lambda/2) sum_d (T_K(d+1)-T_K(d))^2, audited at lambda=1",
        "base_action_lambda_one": round(base_action_lambda_one, 15),
        "reference_cell_rows": cell_rows,
        "scale_orbit_rows": scale_rows,
        "proof_obligations": obligations,
        "accepted_as_strict_action_quantum_reference_cell_source": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3020_path: Any) -> dict[str, Any]:
    matrix = build_action_quantum_matrix()
    return {
        "status": "P3021_ACTION_QUANTUM_REFERENCE_CELL_SOURCE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3020": hashlib.sha256(p3020_path.read_bytes()).hexdigest() if p3020_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "reference_cell_candidate_count": len(matrix["reference_cell_rows"]),
            "positive_reference_cell_candidates": sum(1 for row in matrix["reference_cell_rows"] if row["positive"]),
            "scale_row_count": len(matrix["scale_orbit_rows"]),
            "stable_cell_count_rows": sum(1 for row in matrix["scale_orbit_rows"] if row["cell_counts_same_as_base"]),
            "quanta_rescale_rows": sum(1 for row in matrix["scale_orbit_rows"] if row["quanta_rescale_by_c_square"]),
            "accepted_as_strict_action_quantum_reference_cell_source": matrix["accepted_as_strict_action_quantum_reference_cell_source"],
        },
        "decision": {
            "breakthrough": "Four finite reference-cell partitions of the typed T_K action were constructed and give positive action-per-cell candidates.  The obstruction is that all candidates partition the same dimensionless formal action: under T_K -> c T_K the cell counts remain stable but every action quantum rescales by c^2, so no independent strict action unit or physical reference cell is sourced.",
            "negative_export_flags": {k: False for k in ["strict_action_quantum_reference_cell_source_exported", "strict_clock_unit_theorem_exported", "strict_observable_unit_source_exported", "strict_lambda_action_unit_source_exported", "unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_arrow_exported", "ltotal_exported", "observed_physics_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay internal partitions of the T_K formal action as action quanta.  The next proof-grade move should either construct a genuinely new strict time-order object carrying a directed successor plus a physical unit theorem, or perform a post-P3017-P3021 unit-atom reconciliation/no-new-live-frontier certificate before selecting a different typed object.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3021/S1971 action-quantum/reference-cell source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- positive reference-cell candidates / total: `{c['positive_reference_cell_candidates']}/{c['reference_cell_candidate_count']}`",
        f"- stable cell-count rows / total scale rows: `{c['stable_cell_count_rows']}/{c['scale_row_count']}`",
        f"- action-quanta rescale rows / total scale rows: `{c['quanta_rescale_rows']}/{c['scale_row_count']}`",
        f"- accepted as strict action-quantum/reference-cell source: `{c['accepted_as_strict_action_quantum_reference_cell_source']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3020)
    payload = build_payload(P3020)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3021/S1971 action-quantum/reference-cell source obstruction", "## P3021/S1971 action-quantum/reference-cell source obstruction\n\n`P3021/S1971` attacks the remaining P3020 unit atom: a strict action quantum/reference-cell source coupled to the typed formal action for `T_K`.  It constructs four positive finite reference-cell partitions of the lambda-one action: all cyclic edges, nonzero-gradient support, cyclic Laplacian rank, and nonzero DFT modes.  The bounded obstruction is rescaling-partition dependence: the cell counts are stable under `T_K -> c T_K`, but every action-per-cell candidate rescales by `c^2`, and no independent strict reference-cell theorem, physical action unit, or Hamiltonian normalization is exported.  No strict action quantum/reference-cell source, unit-bearing EOM/Hamiltonian, time arrow, `L_total`, observed-physics export, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3021/S1971 action-quantum/reference-cell `L_total` guard", "## P3021/S1971 action-quantum/reference-cell `L_total` guard\n\n`P3021/S1971` adds no physical `L_total` term.  Its edge, gradient-support, rank, and DFT-mode reference cells partition the already dimensionless formal `T_K` action; since their action-per-cell candidates rescale as `c^2` under observable rescaling and lack an independent physical action unit/Hamiltonian coupling, they cannot install unit-bearing EOM or Hamiltonian terms.\n")
    append_once(AGENTS, "Current action-quantum/reference-cell source guardrail (P3021/S1971, 2026-06-22)", "## Current action-quantum/reference-cell source guardrail (P3021/S1971, 2026-06-22)\n\n- P3021 attacks the remaining P3020 unit atom: a strict action quantum/reference-cell source coupled to the typed formal action for `T_K`.\n- Four finite reference-cell partitions are positive (all cyclic edges, nonzero-gradient support, cyclic Laplacian rank, nonzero DFT modes), but they partition a dimensionless formal action: under `T_K -> c T_K` every action-per-cell candidate rescales by `c^2`.\n- Do not promote internal action partitions to strict action quantum/reference-cell source, unit-bearing EOM/Hamiltonian, time arrow, `L_total`, observed-physics, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should either construct a genuinely new strict time-order object carrying a directed successor plus a physical unit theorem, or perform a post-P3017-P3021 unit-atom reconciliation/no-new-live-frontier certificate before selecting a different typed object.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
