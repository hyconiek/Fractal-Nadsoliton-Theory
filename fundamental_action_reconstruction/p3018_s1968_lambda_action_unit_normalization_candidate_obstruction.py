#!/usr/bin/env python3
"""P3018/S1968: lambda/action-unit normalization candidate obstruction.

This follows P3017 by attacking exactly one named unit-source atom: can the
formal quadratic action for the typed time observable receive a strict
lambda/action-unit normalization?

A concrete normalization candidate is built: choose lambda so the formal action
of the audited T_K vector equals one action quantum in dimensionless units,

    lambda_* = 2 / sum_d (T_K(d+1)-T_K(d))^2.

This is computationally real and positive.  The bounded obstruction is that it
is a target-setting convention, not a strict unit source: rescaling the observable
changes lambda_* by c^-2, and no repository artifact supplies the independent
action quantum, observable unit, clock unit, or Hamiltonian normalization needed
to interpret lambda_* physically.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3017_s1967_time_observable_action_eom_source_obstruction import OUT as P3017

OUT = GEN / "p3018_s1968_lambda_action_unit_normalization_candidate_obstruction.json"
MD = GEN / "p3018_s1968_lambda_action_unit_normalization_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
OMEGA = 0.18575
PHI = 0.16250
BETA = 1.0
ETA = 1.8
SCALE_FACTORS = [0.5, 1.0, 2.0, 3.0]
TARGET_ACTIONS = [0.5, 1.0, 2.0]


def k_strict(d: int) -> float:
    return math.cos(OMEGA * d + PHI) / (1.0 + BETA * (float(d) ** ETA))


def t_obs(d: int) -> float:
    return k_strict(d + 1) - k_strict(d)


def cyclic_gradient_energy(q: list[float]) -> float:
    return sum((q[(i + 1) % len(q)] - q[i]) ** 2 for i in range(len(q)))


def build_lambda_normalization_matrix() -> dict[str, Any]:
    q = [t_obs(d) for d in range(1, N + 1)]
    energy = cyclic_gradient_energy(q)
    lambda_star = 2.0 / energy
    normalized_action = 0.5 * lambda_star * energy
    scale_rows = []
    for c in SCALE_FACTORS:
        scaled_energy = cyclic_gradient_energy([c * x for x in q])
        scaled_lambda = 2.0 / scaled_energy
        scale_rows.append({
            "observable_scale_c": c,
            "gradient_energy": round(scaled_energy, 12),
            "lambda_star_for_action_1": round(scaled_lambda, 12),
            "lambda_ratio_to_unscaled": round(scaled_lambda / lambda_star, 12),
            "expected_c_inverse_square": round(c ** -2, 12),
            "same_lambda_as_unscaled": math.isclose(scaled_lambda, lambda_star, rel_tol=0.0, abs_tol=1e-12),
        })
    target_rows = []
    for target in TARGET_ACTIONS:
        target_lambda = 2.0 * target / energy
        target_rows.append({
            "target_action_quantum": target,
            "lambda_for_target": round(target_lambda, 12),
            "ratio_to_lambda_star": round(target_lambda / lambda_star, 12),
            "requires_external_action_quantum": target != 1.0,
        })
    obligations = [
        {"obligation": "positive_lambda_candidate", "satisfied": lambda_star > 0.0 and math.isfinite(lambda_star), "detail": f"lambda_*={lambda_star:.12g}"},
        {"obligation": "normalizes_formal_action_to_one", "satisfied": math.isclose(normalized_action, 1.0, rel_tol=0.0, abs_tol=1e-12), "detail": f"S_formal[T_K]={normalized_action:.12g}"},
        {"obligation": "observable_scale_invariant_lambda", "satisfied": all(row["same_lambda_as_unscaled"] for row in scale_rows), "detail": "lambda_* changes under q -> c q by c^-2"},
        {"obligation": "strict_action_quantum_source", "satisfied": False, "detail": "setting S=1 is a convention unless a strict action quantum is exported"},
        {"obligation": "strict_observable_unit_source", "satisfied": False, "detail": "T_K remains dimensionless kernel-difference without a physical readout unit"},
        {"obligation": "clock_and_hamiltonian_unit_source", "satisfied": False, "detail": "no strict clock unit maps action normalization to Hamiltonian/energy"},
    ]
    return {
        "object": "LambdaActionUnitNormalization_CandidateObstructionMatrix",
        "typed_observable": "T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)",
        "normalization_rule": "lambda_* = 2 / sum_d (T_K(d+1)-T_K(d))^2 so that S_formal[T_K]=1",
        "gradient_energy": round(energy, 15),
        "lambda_star": round(lambda_star, 15),
        "normalized_action": round(normalized_action, 15),
        "scale_rows": scale_rows,
        "target_action_rows": target_rows,
        "proof_obligations": obligations,
        "accepted_as_strict_lambda_action_unit_source": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3017_path: Any) -> dict[str, Any]:
    matrix = build_lambda_normalization_matrix()
    return {
        "status": "P3018_LAMBDA_ACTION_UNIT_NORMALIZATION_CANDIDATE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3017": hashlib.sha256(p3017_path.read_bytes()).hexdigest() if p3017_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "gradient_energy_positive": matrix["gradient_energy"] > 0.0,
            "lambda_star_positive": matrix["lambda_star"] > 0.0,
            "normalized_action": matrix["normalized_action"],
            "scale_row_count": len(matrix["scale_rows"]),
            "scale_invariant_lambda_rows": sum(1 for row in matrix["scale_rows"] if row["same_lambda_as_unscaled"]),
            "target_action_row_count": len(matrix["target_action_rows"]),
            "accepted_as_strict_lambda_action_unit_source": matrix["accepted_as_strict_lambda_action_unit_source"],
        },
        "decision": {
            "breakthrough": "A concrete positive lambda normalization exists for the typed observable: choose lambda_* so the formal quadratic action of T_K equals one.  The obstruction is that this is target/action-quantum setting, not a strict unit source: lambda_* changes under observable rescaling and no independent action, observable, clock, or Hamiltonian unit theorem is exported.",
            "negative_export_flags": {k: False for k in ["strict_lambda_action_unit_source_exported", "unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_arrow_exported", "ltotal_exported", "observed_physics_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not treat action-to-one normalization as a strict physical unit theorem.  The next proof-grade move should attack exactly one independent unit atom: a strict action quantum/reference-cell source, a strict observable-unit readout source, or a strict clock-unit theorem; otherwise preserve the P3017-P3018 no-closure boundary.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    obj = payload["constructed_theoretical_objects"]
    lines = [
        "# P3018/S1968 lambda/action-unit normalization candidate obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- gradient energy positive: `{c['gradient_energy_positive']}`",
        f"- lambda_star positive: `{c['lambda_star_positive']}`",
        f"- lambda_star: `{obj['lambda_star']}`",
        f"- normalized action: `{c['normalized_action']}`",
        f"- scale-invariant lambda rows / total: `{c['scale_invariant_lambda_rows']}/{c['scale_row_count']}`",
        f"- accepted as strict lambda/action-unit source: `{c['accepted_as_strict_lambda_action_unit_source']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3017)
    payload = build_payload(P3017)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3018/S1968 lambda/action-unit normalization candidate obstruction", "## P3018/S1968 lambda/action-unit normalization candidate obstruction\n\n`P3018/S1968` attacks exactly one P3017 unit source atom: `lambda`/action-unit normalization for the typed observable `T_K`.  It constructs a positive computational candidate `lambda_* = 2 / sum_d (T_K(d+1)-T_K(d))^2`, so the formal quadratic action of `T_K` is normalized to one.  This is real progress as a finite normalization witness, but the acceptance matrix is bounded no-go: `lambda_*` changes under observable rescaling, and setting action to one is a target convention unless an independent strict action quantum/reference-cell, observable-unit readout, clock-unit, and Hamiltonian normalization theorem is exported.  No strict lambda/action-unit source, unit-bearing EOM/Hamiltonian, `L_total`, observed-physics export, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3018/S1968 lambda/action-unit normalization `L_total` guard", "## P3018/S1968 lambda/action-unit normalization `L_total` guard\n\n`P3018/S1968` adds no physical `L_total` term.  The positive `lambda_*` that sets the formal action of `T_K` to one is a normalization convention unless a strict action quantum/reference-cell and observable/clock/Hamiltonian unit theorem is supplied; it changes under observable rescaling and cannot by itself install a unit-bearing action or EOM.\n")
    append_once(AGENTS, "Current lambda/action-unit normalization guardrail (P3018/S1968, 2026-06-22)", "## Current lambda/action-unit normalization guardrail (P3018/S1968, 2026-06-22)\n\n- P3018 constructs a positive `lambda_*` normalization for the typed time observable `T_K` by setting its formal quadratic action to one.\n- The finite obstruction is normalization-gauge dependence: `lambda_*` rescales as `c^-2` under `T_K -> c T_K`, and no independent strict action quantum/reference-cell, observable-unit, clock-unit, or Hamiltonian normalization theorem is exported.\n- Do not promote action-to-one normalization to strict lambda/action-unit source, unit-bearing EOM/Hamiltonian, `L_total`, observed-physics, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should attack exactly one independent unit atom: strict action quantum/reference-cell, strict observable-unit readout, or strict clock-unit theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
