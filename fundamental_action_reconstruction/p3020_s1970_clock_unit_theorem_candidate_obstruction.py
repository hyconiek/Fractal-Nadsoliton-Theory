#!/usr/bin/env python3
"""P3020/S1970: clock-unit theorem candidate obstruction.

After P3019, attack exactly one remaining independent unit atom: can the typed
time observable T_K(d) receive a strict clock-unit theorem?

We construct three finite clock candidates: the label tick Δd=1, the full cyclic
period 12, and the dominant discrete Fourier period of the T_K vector.  These
are real dimensionless clock scaffolds, but the audit keeps them from becoming
strict physical clock units: the directed tick is not U(12)-invariant, the full
cycle has no local directed successor/unit scale, and the spectral period is a
feature of a dimensionless sampled vector without a physical tick, action, or
Hamiltonian coupling theorem.
"""
from __future__ import annotations

import cmath, hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3018_s1968_lambda_action_unit_normalization_candidate_obstruction import t_obs
from p3019_s1969_observable_unit_readout_source_obstruction import OUT as P3019

OUT = GEN / "p3020_s1970_clock_unit_theorem_candidate_obstruction.json"
MD = GEN / "p3020_s1970_clock_unit_theorem_candidate_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = [1, 5, 7, 11]


def dft_magnitudes(q: list[float]) -> list[float]:
    mags = []
    for k in range(N):
        z = sum(q[n] * cmath.exp(-2j * math.pi * k * n / N) for n in range(N))
        mags.append(abs(z))
    return mags


def orbit_of_step(step: int) -> list[int]:
    return sorted({(u * step) % N for u in UNITS})


def build_clock_unit_matrix() -> dict[str, Any]:
    q = [t_obs(d) for d in range(1, N + 1)]
    mags = dft_magnitudes(q)
    nonzero_modes = list(range(1, N))
    dominant_mode = max(nonzero_modes, key=lambda k: mags[k])
    dominant_modes = [k for k in nonzero_modes if math.isclose(mags[k], mags[dominant_mode], rel_tol=0.0, abs_tol=1e-12)]
    dominant_periods = sorted({N // math.gcd(N, k) for k in dominant_modes})
    tick_images = [{"unit": u, "image_of_plus_one_tick": u % N, "preserves_plus_one_tick": (u % N) == 1} for u in UNITS]
    clock_rows = [
        {
            "candidate_clock": "label_tick_delta_d_equals_1",
            "finite_value": 1,
            "positive": True,
            "unit_action_images": tick_images,
            "unit_invariant_directed_tick": all(row["preserves_plus_one_tick"] for row in tick_images),
            "exports_physical_clock_unit": False,
        },
        {
            "candidate_clock": "full_z12_cycle_period",
            "finite_value": N,
            "positive": True,
            "unit_invariant_period": all((u * N) % N == 0 for u in UNITS),
            "has_local_directed_tick": False,
            "exports_physical_clock_unit": False,
        },
        {
            "candidate_clock": "dominant_dft_period_of_T_K",
            "dominant_modes": dominant_modes,
            "dominant_periods_in_label_steps": dominant_periods,
            "dominant_magnitude": round(mags[dominant_mode], 15),
            "positive": bool(dominant_periods),
            "mode_orbits_under_units": {str(k): orbit_of_step(k) for k in dominant_modes},
            "has_physical_frequency_unit": False,
            "exports_physical_clock_unit": False,
        },
    ]
    obligations = [
        {"obligation": "typed_time_observable_input", "satisfied": True, "detail": "uses the P3014/P3017/P3019 T_K(d) vector"},
        {"obligation": "positive_finite_clock_candidates", "satisfied": all(row["positive"] for row in clock_rows), "detail": "label tick, cycle period, and DFT period exist"},
        {"obligation": "unit_invariant_directed_tick", "satisfied": clock_rows[0]["unit_invariant_directed_tick"], "detail": "U(12) sends +1 tick to 1,5,7,11; orientation is not preserved"},
        {"obligation": "local_successor_from_cycle_period", "satisfied": clock_rows[1]["has_local_directed_tick"], "detail": "the period 12 is invariant but does not choose a local directed tick"},
        {"obligation": "physical_frequency_or_clock_unit_source", "satisfied": False, "detail": "DFT period is in label steps only; no strict physical tick/frequency theorem is exported"},
        {"obligation": "action_hamiltonian_coupling", "satisfied": False, "detail": "no action quantum or Hamiltonian normalization converts the clock candidate into physical time"},
    ]
    return {
        "object": "ClockUnitTheoremCandidate_UnitActionSpectralObstructionMatrix",
        "typed_observable": "T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)",
        "clock_candidate_rows": clock_rows,
        "proof_obligations": obligations,
        "accepted_as_strict_clock_unit_theorem": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3019_path: Any) -> dict[str, Any]:
    matrix = build_clock_unit_matrix()
    return {
        "status": "P3020_CLOCK_UNIT_THEOREM_CANDIDATE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3019": hashlib.sha256(p3019_path.read_bytes()).hexdigest() if p3019_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "clock_candidate_count": len(matrix["clock_candidate_rows"]),
            "positive_clock_candidates": sum(1 for row in matrix["clock_candidate_rows"] if row["positive"]),
            "plus_one_tick_preserving_unit_count": sum(1 for row in matrix["clock_candidate_rows"][0]["unit_action_images"] if row["preserves_plus_one_tick"]),
            "unit_count": len(UNITS),
            "dominant_dft_modes": matrix["clock_candidate_rows"][2]["dominant_modes"],
            "accepted_as_strict_clock_unit_theorem": matrix["accepted_as_strict_clock_unit_theorem"],
        },
        "decision": {
            "breakthrough": "Three finite clock-unit candidates were constructed for the typed T_K observable: label tick, full Z12 period, and dominant DFT period.  The obstruction is that these are dimensionless scaffolds: the directed tick is not U(12)-invariant, the invariant cycle period does not choose a local tick, and the spectral period has no strict physical frequency/action/Hamiltonian coupling.",
            "negative_export_flags": {k: False for k in ["strict_clock_unit_theorem_exported", "strict_observable_unit_source_exported", "strict_lambda_action_unit_source_exported", "unit_bearing_action_eom_source_exported", "hamiltonian_exported", "time_arrow_exported", "ltotal_exported", "observed_physics_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay label-tick, cycle-period, or DFT-period clock-unit candidates.  The next proof-grade move should attack the remaining independent action quantum/reference-cell source coupled to the typed action, or introduce a genuinely new strict time-order object carrying both a directed successor and its own physical unit theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3020/S1970 clock-unit theorem candidate obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- positive clock candidates / total: `{c['positive_clock_candidates']}/{c['clock_candidate_count']}`",
        f"- +1 tick preserving units / total units: `{c['plus_one_tick_preserving_unit_count']}/{c['unit_count']}`",
        f"- dominant DFT modes: `{c['dominant_dft_modes']}`",
        f"- accepted as strict clock-unit theorem: `{c['accepted_as_strict_clock_unit_theorem']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3019)
    payload = build_payload(P3019)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3020/S1970 clock-unit theorem candidate obstruction", "## P3020/S1970 clock-unit theorem candidate obstruction\n\n`P3020/S1970` attacks exactly one P3019 remaining unit atom: a strict clock-unit theorem for the typed time observable `T_K`.  It constructs three finite dimensionless clock candidates: the label tick `Δd=1`, the full `Z12` cycle period `12`, and the dominant DFT period of the `T_K` vector.  The bounded obstruction is finite and explicit: the directed `+1` tick is preserved by only one of four `U(12)` units, the invariant cycle period does not choose a local directed tick, and the DFT period remains a label-step feature without a strict physical frequency, action quantum, or Hamiltonian coupling theorem.  No strict clock-unit theorem, time arrow, unit-bearing EOM/Hamiltonian, `L_total`, observed-physics export, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3020/S1970 clock-unit theorem `L_total` guard", "## P3020/S1970 clock-unit theorem `L_total` guard\n\n`P3020/S1970` adds no physical `L_total` term.  The label tick, full-cycle period, and dominant DFT period are dimensionless clock scaffolds only; without a unit-invariant directed tick plus strict physical frequency/action/Hamiltonian coupling, they cannot install a unit-bearing EOM or Hamiltonian.\n")
    append_once(AGENTS, "Current clock-unit theorem candidate guardrail (P3020/S1970, 2026-06-22)", "## Current clock-unit theorem candidate guardrail (P3020/S1970, 2026-06-22)\n\n- P3020 attacks one P3019 remaining unit atom: a strict clock-unit theorem for the typed time observable `T_K`.\n- Three finite candidates are constructed (label tick `Δd=1`, full `Z12` cycle period `12`, and dominant DFT period), but they remain dimensionless scaffolds: the `+1` tick is not `U(12)`-invariant, the invariant cycle period does not choose a local directed tick, and the DFT period has no physical frequency/action/Hamiltonian coupling.\n- Do not promote label-tick, cycle-period, or DFT-period candidates to strict clock-unit theorem, time arrow, unit-bearing EOM/Hamiltonian, `L_total`, observed-physics, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should attack the remaining independent action quantum/reference-cell source coupled to the typed action, or introduce a genuinely new strict time-order object carrying both a directed successor and its own physical unit theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
