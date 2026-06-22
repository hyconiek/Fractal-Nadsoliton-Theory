#!/usr/bin/env python3
"""P3025/S1975: physical tick/action/Hamiltonian unit obstruction.

After P3024 blocks a non-premise chart source for the P3023 dissipation chain,
attack the remaining single theorem atom for that same object: whether an
independent physical tick/action/Hamiltonian unit theorem can be constructed.

The finite matrix builds four explicit tick candidates from the chain and asks
whether they are (i) positive, (ii) independent of the observable amplitude
rescaling K -> cK, and (iii) coupled to an independent action/Hamiltonian unit.
The candidates are real internal calibrations, but none provides an absolute
physical unit or Hamiltonian source.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3024_s1974_dissipation_chart_selector_source_obstruction import OUT as P3024

OUT = GEN / "p3025_s1975_dissipation_physical_tick_hamiltonian_unit_obstruction.json"
MD = GEN / "p3025_s1975_dissipation_physical_tick_hamiltonian_unit_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

RESCALE = 3.0


def chain_values(scale: float = 1.0) -> list[float]:
    return [scale * k_strict(d) for d in range(1, N + 1)]


def chain_drops(values: list[float]) -> list[float]:
    return [values[d - 1] - values[d] for d in range(1, N)]


def quadratic_chain_action(values: list[float]) -> float:
    return 0.5 * sum((values[d] - values[d - 1]) ** 2 for d in range(1, N))


def tick_candidates(values: list[float]) -> dict[str, float]:
    drops = chain_drops(values)
    total_drop = values[0] - values[-1]
    rms_drop = math.sqrt(sum(drop * drop for drop in drops) / len(drops))
    return {
        "label_tick": 1.0,
        "mean_positive_drop_tick": total_drop / len(drops),
        "rms_drop_tick": rms_drop,
        "inverse_mean_slope_tick": len(drops) / total_drop,
    }


def scaling_exponent(base: float, scaled: float) -> str:
    if base == 0.0:
        return "undefined"
    ratio = scaled / base
    if math.isclose(ratio, 1.0, rel_tol=1e-12, abs_tol=1e-12):
        return "c^0"
    if math.isclose(ratio, RESCALE, rel_tol=1e-12, abs_tol=1e-12):
        return "c^1"
    if math.isclose(ratio, 1.0 / RESCALE, rel_tol=1e-12, abs_tol=1e-12):
        return "c^-1"
    if math.isclose(ratio, RESCALE * RESCALE, rel_tol=1e-12, abs_tol=1e-12):
        return "c^2"
    return f"ratio={ratio:.12g}"


def build_unit_matrix() -> dict[str, Any]:
    values = chain_values()
    scaled_values = chain_values(RESCALE)
    base_ticks = tick_candidates(values)
    scaled_ticks = tick_candidates(scaled_values)
    base_action = quadratic_chain_action(values)
    scaled_action = quadratic_chain_action(scaled_values)

    rows = []
    for name, tick in base_ticks.items():
        scaled_tick = scaled_ticks[name]
        hamiltonian = base_action / tick
        scaled_hamiltonian = scaled_action / scaled_tick
        rows.append({
            "candidate_tick": name,
            "tick_value": round(tick, 15),
            "positive": tick > 0.0,
            "tick_rescaling": scaling_exponent(tick, scaled_tick),
            "formal_action_value": round(base_action, 15),
            "formal_action_rescaling": scaling_exponent(base_action, scaled_action),
            "formal_hamiltonian_value": round(hamiltonian, 15),
            "formal_hamiltonian_rescaling": scaling_exponent(hamiltonian, scaled_hamiltonian),
            "independent_physical_tick_source": False,
            "independent_action_quantum_source": False,
            "hamiltonian_unit_theorem": False,
        })

    obligations = [
        {"obligation": "attacks_single_P3024_remaining_theorem", "satisfied": True, "detail": "only the independent physical tick/action/Hamiltonian unit theorem is tested"},
        {"obligation": "positive_tick_candidates_constructed", "satisfied": all(row["positive"] for row in rows), "detail": "label, mean-drop, RMS-drop, and inverse-slope ticks are positive"},
        {"obligation": "observable_rescaling_independent_tick", "satisfied": all(row["tick_rescaling"] == "c^0" for row in rows), "detail": "drop-based ticks scale as c or c^-1, while label tick is only dimensionless"},
        {"obligation": "independent_action_quantum_source", "satisfied": False, "detail": "formal chain action rescales as c^2 and no strict action quantum is supplied"},
        {"obligation": "physical_Hamiltonian_unit_coupling", "satisfied": False, "detail": "S/tau is an internal ratio; no energy/frequency/unit theorem maps it to a physical Hamiltonian"},
        {"obligation": "strict_physical_tick_theorem_exported", "satisfied": False, "detail": "no candidate supplies an absolute physical tick independent of chart and amplitude conventions"},
    ]
    return {
        "object": "DissipationChainPhysicalTickHamiltonianUnit_RescalingCouplingObstructionMatrix",
        "tested_theorem_atom": "independent physical tick/action/Hamiltonian unit theorem for the P3023/P3024 dissipation chain",
        "rescaling_test": f"K -> {RESCALE} K",
        "tick_rows": rows,
        "proof_obligations": obligations,
        "accepted_as_physical_tick_hamiltonian_unit_theorem": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3024_path: Any) -> dict[str, Any]:
    read_json(p3024_path)
    matrix = build_unit_matrix()
    return {
        "status": "P3025_DISSIPATION_PHYSICAL_TICK_HAMILTONIAN_UNIT_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3024": hashlib.sha256(p3024_path.read_bytes()).hexdigest() if p3024_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "tick_candidate_count": len(matrix["tick_rows"]),
            "positive_tick_candidates": sum(1 for row in matrix["tick_rows"] if row["positive"]),
            "rescaling_independent_tick_rows": sum(1 for row in matrix["tick_rows"] if row["tick_rescaling"] == "c^0"),
            "hamiltonian_unit_rows": sum(1 for row in matrix["tick_rows"] if row["hamiltonian_unit_theorem"]),
            "accepted_as_physical_tick_hamiltonian_unit_theorem": matrix["accepted_as_physical_tick_hamiltonian_unit_theorem"],
        },
        "decision": {
            "breakthrough": "Four explicit tick/action/Hamiltonian couplings were constructed for the P3023 dissipation chain.  They are positive internal calibrations, but the label tick is dimensionless, the drop-derived ticks depend on K-amplitude rescaling, the formal action rescales as c^2, and S/tau remains an internal ratio without an exported physical energy/frequency theorem.",
            "negative_export_flags": {k: False for k in ["physical_tick_theorem_exported", "action_quantum_source_exported", "hamiltonian_unit_theorem_exported", "strict_time_order_object_exported", "unit_bearing_action_eom_source_exported", "time_arrow_exported", "observed_physics_exported", "ltotal_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay internal tick/action/Hamiltonian ratios for the P3023 chain.  With chart-source and physical-unit atoms both bounded no-go, the next proof-grade move should emit a P3023-P3025 time-order no-new-live-frontier reconciliation, or pivot to a genuinely new strict typed object with an external physical unit source.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3025/S1975 dissipation physical tick/Hamiltonian unit obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- tick candidates: `{c['tick_candidate_count']}`",
        f"- positive tick candidates: `{c['positive_tick_candidates']}`",
        f"- rescaling-independent tick rows: `{c['rescaling_independent_tick_rows']}`",
        f"- Hamiltonian unit theorem rows: `{c['hamiltonian_unit_rows']}`",
        f"- accepted as physical tick/Hamiltonian unit theorem: `{c['accepted_as_physical_tick_hamiltonian_unit_theorem']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload(P3024)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3025/S1975 dissipation physical tick/Hamiltonian unit obstruction", "## P3025/S1975 dissipation physical tick/Hamiltonian unit obstruction\n\n`P3025/S1975` attacks the remaining non-replay theorem atom for the P3023/P3024 dissipation chain: an independent physical tick/action/Hamiltonian unit theorem.  It constructs four explicit tick couplings (label tick, mean positive drop, RMS drop, inverse mean slope) and the formal ratio `H=S/tau`.  The bounded obstruction is finite: all ticks are positive internal calibrations, but the label tick is dimensionless, drop-derived ticks scale under `K -> cK`, the chain action rescales as `c^2`, and no strict action quantum or physical energy/frequency coupling theorem is exported.  No strict physical tick theorem, Hamiltonian, time arrow, unit-bearing EOM, observed-physics export, `L_total`, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3025/S1975 physical tick/Hamiltonian `L_total` guard", "## P3025/S1975 physical tick/Hamiltonian `L_total` guard\n\n`P3025/S1975` adds no physical `L_total` term.  The label, mean-drop, RMS-drop, and inverse-slope ticks for the P3023 chain are internal calibrations only; the formal Hamiltonian ratio `S/tau` has no strict physical unit, action quantum, or energy/frequency coupling theorem.\n")
    append_once(AGENTS, "Current dissipation physical tick/Hamiltonian unit guardrail (P3025/S1975, 2026-06-22)", "## Current dissipation physical tick/Hamiltonian unit guardrail (P3025/S1975, 2026-06-22)\n\n- P3025 attacks the remaining non-replay P3023/P3024 theorem atom: an independent physical tick/action/Hamiltonian unit theorem for the kernel-dissipation chain.\n- Four tick couplings are positive internal calibrations, but the label tick is dimensionless, drop-derived ticks scale under `K -> cK`, the formal action scales as `c^2`, and `H=S/tau` has no strict energy/frequency coupling theorem.\n- Do not promote these internal tick/action/Hamiltonian ratios to physical time arrow, unit-bearing EOM/Hamiltonian, observed-physics, `L_total`, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should emit a P3023-P3025 time-order no-new-live-frontier reconciliation, or pivot to a genuinely new strict typed object with an external physical unit source.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
