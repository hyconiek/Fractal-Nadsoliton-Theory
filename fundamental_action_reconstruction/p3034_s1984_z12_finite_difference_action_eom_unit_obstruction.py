#!/usr/bin/env python3
"""P3034/S1984: Z12 finite-difference action/EOM unit obstruction.

Attack exactly one P3028 foundation atom: unit-bearing action/EOM/Hamiltonian.
This is not a replay of prior graph/source-specific action-installation lanes;
it is a narrow P3028 transition-foundation receiver for sampled K_strict_gate on
Z12.  The formal cyclic Dirichlet+mass action, EOM residual, and Hamiltonian
proxy are computable and cyclic-shift covariant, but no physical action unit,
field provenance, clock unit, or nonproxy continuum lift is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3028_s1978_nadsoliton_information_to_classical_transition_foundation_lattice import OUT as P3028
from p3033_s1983_z12_block_coarse_graining_classical_limit_obstruction import OUT as P3033

OUT = GEN / "p3034_s1984_z12_finite_difference_action_eom_unit_obstruction.json"
MD = GEN / "p3034_s1984_z12_finite_difference_action_eom_unit_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

MASS_SQUARED = 1.0
RESCALE_C = 3.0


def kernel_vector(scale: float = 1.0, shift: int = 0) -> list[float]:
    base = [scale * k_strict(label) for label in range(1, N + 1)]
    return [base[(j - shift) % N] for j in range(N)]


def cyclic_gradient(values: list[float]) -> list[float]:
    return [values[(i + 1) % N] - values[i] for i in range(N)]


def formal_action(values: list[float], mass_squared: float = MASS_SQUARED) -> float:
    grad = cyclic_gradient(values)
    dirichlet = 0.5 * sum(g * g for g in grad)
    mass_term = 0.5 * mass_squared * sum(v * v for v in values)
    return dirichlet + mass_term


def eom_residual(values: list[float], mass_squared: float = MASS_SQUARED) -> list[float]:
    # Euler row for S = 1/2 sum (x_{i+1}-x_i)^2 + m^2/2 sum x_i^2.
    return [2 * values[i] - values[(i - 1) % N] - values[(i + 1) % N] + mass_squared * values[i] for i in range(N)]


def l2_norm(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values))


def scaling_exponent(base: float, scaled: float, factor: float = RESCALE_C) -> float | None:
    if base == 0 or scaled == 0 or factor <= 0:
        return None
    return round(math.log(abs(scaled / base), factor), 12)


def build_matrix() -> dict[str, Any]:
    values = kernel_vector()
    shifted_values = kernel_vector(shift=1)
    scaled_values = kernel_vector(scale=RESCALE_C)
    action = formal_action(values)
    shifted_action = formal_action(shifted_values)
    scaled_action = formal_action(scaled_values)
    residual = eom_residual(values)
    scaled_residual = eom_residual(scaled_values)
    hamiltonian_proxy = action / N
    scaled_hamiltonian_proxy = scaled_action / N
    receiver_rows = [
        {
            "receiver": "cyclic_dirichlet_plus_mass_action",
            "value": round(action, 12),
            "shift_covariant": abs(action - shifted_action) <= 1e-12,
            "rescaling_exponent": scaling_exponent(action, scaled_action),
            "unit_bearing_exported": False,
            "failure": "formal action scales as K^2 and has no physical action quantum/unit source",
        },
        {
            "receiver": "cyclic_euler_lagrange_residual",
            "l2_norm": round(l2_norm(residual), 12),
            "rescaling_exponent": scaling_exponent(l2_norm(residual), l2_norm(scaled_residual)),
            "unit_bearing_exported": False,
            "failure": "finite residual is computable but is not a nonproxy field EOM with units or boundary/continuum lift",
        },
        {
            "receiver": "hamiltonian_proxy_action_per_label_tick",
            "value": round(hamiltonian_proxy, 12),
            "rescaling_exponent": scaling_exponent(hamiltonian_proxy, scaled_hamiltonian_proxy),
            "unit_bearing_exported": False,
            "failure": "division by the dimensionless label count N is not a physical clock/frequency theorem",
        },
    ]
    obligations = [
        {"obligation": "attacks_single_P3028_foundation_atom", "satisfied": True, "detail": "only unit-bearing action/EOM/Hamiltonian is tested"},
        {"obligation": "explicit_formal_action_receiver", "satisfied": True, "detail": "cyclic Dirichlet+mass action on sampled K_strict_gate is constructed"},
        {"obligation": "finite_eom_residual_computable", "satisfied": True, "detail": "Euler residual vector is computed on Z12"},
        {"obligation": "cyclic_shift_covariant_receiver", "satisfied": all(row.get("shift_covariant", True) for row in receiver_rows), "detail": "formal action is invariant under cyclic relabeling"},
        {"obligation": "physical_action_unit_source", "satisfied": False, "detail": "no action quantum/reference-cell/source theorem labels the action value"},
        {"obligation": "field_provenance_and_boundary_map", "satisfied": False, "detail": "K samples are not exported as a physical field with boundary/integration map"},
        {"obligation": "physical_clock_or_hamiltonian_unit", "satisfied": False, "detail": "action/N is a label-tick proxy, not a physical Hamiltonian/frequency theorem"},
        {"obligation": "nonproxy_continuum_eom_lift", "satisfied": False, "detail": "finite Z12 residual is not lifted to tensor-resolved continuum EOM"},
    ]
    return {
        "object": "Z12FiniteDifferenceActionEOMHamiltonian_UnitObstructionMatrix",
        "tested_foundation_atom": "unit_bearing_action_eom_hamiltonian",
        "formal_action": "S[x]=1/2 sum_i (x_{i+1}-x_i)^2 + m^2/2 sum_i x_i^2 on cyclic Z12, m^2=1",
        "receiver_rows": receiver_rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "receiver_rows": len(receiver_rows),
            "unit_bearing_receiver_rows": sum(1 for row in receiver_rows if row["unit_bearing_exported"]),
            "action_rescaling_exponent": receiver_rows[0]["rescaling_exponent"],
            "residual_rescaling_exponent": receiver_rows[1]["rescaling_exponent"],
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "unit_bearing_action_eom_hamiltonian_exported": all(row["satisfied"] for row in obligations) and all(row["unit_bearing_exported"] for row in receiver_rows),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3028)
    read_json(P3033)
    matrix = build_matrix()
    return {
        "status": "P3034_Z12_FINITE_DIFFERENCE_ACTION_EOM_UNIT_OBSTRUCTION_NO_LTOTAL_EXPORT",
        "input_hashes": {
            "P3028": hashlib.sha256(P3028.read_bytes()).hexdigest() if P3028.exists() else None,
            "P3033": hashlib.sha256(P3033.read_bytes()).hexdigest() if P3033.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "A formal cyclic Dirichlet+mass action, Euler residual, and Hamiltonian proxy are computable for sampled K_strict_gate on Z12.  The action scales as K^2, the residual scales as K, and the Hamiltonian proxy is action per dimensionless label tick.  Without a physical action unit, field provenance/boundary map, clock unit, or nonproxy continuum EOM lift, this receiver does not export unit-bearing action/EOM/Hamiltonian or L_total.",
            "negative_export_flags": {k: False for k in ["unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "hamiltonian_exported", "classical_transition_exported", "observed_physics_exported", "matter_sector_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay finite Z12 quadratic actions, action/N Hamiltonian proxies, or internal action normalizations as unit-bearing physics.  A next move must supply an actual action quantum/reference-cell theorem, field provenance plus boundary/integration map, or pivot to another single P3028 atom such as strict selector/branch source.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3034/S1984 Z12 finite-difference action/EOM unit obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- receiver rows: `{c['receiver_rows']}`",
        f"- unit-bearing receiver rows: `{c['unit_bearing_receiver_rows']}`",
        f"- action rescaling exponent: `{c['action_rescaling_exponent']}`",
        f"- residual rescaling exponent: `{c['residual_rescaling_exponent']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- unit-bearing action/EOM/Hamiltonian exported: `{c['unit_bearing_action_eom_hamiltonian_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3034/S1984 Z12 finite-difference action/EOM unit obstruction", "## P3034/S1984 Z12 finite-difference action/EOM unit obstruction\n\n`P3034/S1984` attacks exactly one P3028 foundation atom: unit-bearing action/EOM/Hamiltonian.  It constructs a formal cyclic Dirichlet+mass action on sampled `K_strict_gate` over `Z12`, computes the Euler residual, and forms an action-per-label Hamiltonian proxy.  These receivers are finite and computable, but `0/3` receiver rows are unit-bearing: the action scales as `K^2`, the residual scales as `K`, the Hamiltonian proxy divides by a dimensionless label count, and no physical action unit, field provenance/boundary map, clock unit, or nonproxy continuum EOM lift is exported.  No `L_total`, observed physics, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3034/S1984 finite-difference action/EOM `L_total` guard", "## P3034/S1984 finite-difference action/EOM `L_total` guard\n\n`P3034/S1984` adds no physical `L_total` term.  The cyclic Dirichlet+mass action and Euler residual are formal finite receivers for sampled strict-kernel data, but they lack action quantum/reference-cell provenance, field/boundary integration, physical clock/Hamiltonian units, and nonproxy continuum EOM lift.\n")
    append_once(AGENTS, "Current Z12 finite-difference action/EOM unit guardrail (P3034/S1984, 2026-06-23)", "## Current Z12 finite-difference action/EOM unit guardrail (P3034/S1984, 2026-06-23)\n\n- P3034 attacks exactly one P3028 foundation atom: unit-bearing action/EOM/Hamiltonian.\n- A cyclic Dirichlet+mass action, Euler residual, and action-per-label Hamiltonian proxy are computable for sampled `K_strict_gate`, but none is unit-bearing.\n- Do not promote finite Z12 quadratic actions, Euler residuals, action-per-label proxies, or internal action normalizations to physical `L_total`, EOM, Hamiltonian, observed physics, selector, bridge/role-transfer, or ToE closure.\n- A next move must supply an actual action quantum/reference-cell theorem, field provenance plus boundary/integration map, or pivot to another single P3028 atom such as strict selector/branch source.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
