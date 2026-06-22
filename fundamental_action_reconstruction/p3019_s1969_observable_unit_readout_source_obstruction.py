#!/usr/bin/env python3
"""P3019/S1969: observable-unit readout source obstruction.

After P3018, attack exactly one independent unit atom: can the typed time
observable T_K(d) receive a strict observable-unit readout source?

We construct several finite normalization candidates (RMS, L1, L∞, and total
variation units) and audit their behavior under observable rescaling.  The
candidates are computationally real, positive, and observer-independent as
functions of the strict-kernel vector, but they are internal calibration choices:
under T_K -> c T_K the chosen observable unit rescales by c, so the normalized
readout is unchanged while no absolute physical readout unit, detector map, or
Hamiltonian/clock coupling is exported.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3018_s1968_lambda_action_unit_normalization_candidate_obstruction import OUT as P3018, t_obs

OUT = GEN / "p3019_s1969_observable_unit_readout_source_obstruction.json"
MD = GEN / "p3019_s1969_observable_unit_readout_source_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
SCALE_FACTORS = [0.5, 1.0, 2.0, 3.0]


def unit_candidates(q: list[float]) -> dict[str, float]:
    diffs = [q[(i + 1) % len(q)] - q[i] for i in range(len(q))]
    return {
        "rms_unit": math.sqrt(sum(x * x for x in q) / len(q)),
        "l1_mean_abs_unit": sum(abs(x) for x in q) / len(q),
        "linf_unit": max(abs(x) for x in q),
        "total_variation_unit": sum(abs(x) for x in diffs),
    }


def build_observable_unit_matrix() -> dict[str, Any]:
    q = [t_obs(d) for d in range(1, N + 1)]
    base_units = unit_candidates(q)
    unit_rows = []
    for name, value in base_units.items():
        normalized = [x / value for x in q]
        unit_rows.append({
            "candidate_unit": name,
            "unit_value": round(value, 15),
            "positive": value > 0.0 and math.isfinite(value),
            "normalized_rms": round(math.sqrt(sum(x * x for x in normalized) / len(normalized)), 12),
            "normalized_linf": round(max(abs(x) for x in normalized), 12),
            "exports_physical_observable_unit": False,
        })
    scale_rows = []
    for c in SCALE_FACTORS:
        scaled_units = unit_candidates([c * x for x in q])
        scale_rows.append({
            "observable_scale_c": c,
            "unit_ratios_to_base": {name: round(scaled_units[name] / base_units[name], 12) for name in base_units},
            "all_units_rescale_by_c": all(math.isclose(scaled_units[name] / base_units[name], c, rel_tol=0.0, abs_tol=1e-12) for name in base_units),
            "absolute_unit_fixed": all(math.isclose(scaled_units[name], base_units[name], rel_tol=0.0, abs_tol=1e-12) for name in base_units),
        })
    obligations = [
        {"obligation": "typed_time_observable_input", "satisfied": True, "detail": "uses P3014/P3017/P3018 T_K(d) vector"},
        {"obligation": "positive_finite_candidate_units", "satisfied": all(row["positive"] for row in unit_rows), "detail": "RMS, L1, L∞, and total-variation units are positive"},
        {"obligation": "observer_independent_formula", "satisfied": True, "detail": "candidate units are computed from the strict-kernel observable vector only"},
        {"obligation": "observable_rescaling_invariant_absolute_unit", "satisfied": all(row["absolute_unit_fixed"] for row in scale_rows), "detail": "all candidates rescale with T_K -> c T_K"},
        {"obligation": "strict_physical_readout_unit_source", "satisfied": False, "detail": "no independent detector/readout calibration theorem exports an absolute unit"},
        {"obligation": "clock_hamiltonian_coupling", "satisfied": False, "detail": "no clock/Hamiltonian theorem maps this observable unit to energy or time units"},
    ]
    return {
        "object": "ObservableUnitReadoutSource_RescalingOrbitObstructionMatrix",
        "typed_observable": "T_K(d)=K_strict_gate(d+1)-K_strict_gate(d)",
        "candidate_units": unit_rows,
        "scale_orbit_rows": scale_rows,
        "proof_obligations": obligations,
        "accepted_as_strict_observable_unit_source": all(row["satisfied"] for row in obligations),
    }


def build_payload(p3018_path: Any) -> dict[str, Any]:
    matrix = build_observable_unit_matrix()
    return {
        "status": "P3019_OBSERVABLE_UNIT_READOUT_SOURCE_OBSTRUCTION_NO_CLOSURE",
        "input_hashes": {"P3018": hashlib.sha256(p3018_path.read_bytes()).hexdigest() if p3018_path.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": {
            "candidate_unit_count": len(matrix["candidate_units"]),
            "positive_candidate_units": sum(1 for row in matrix["candidate_units"] if row["positive"]),
            "scale_row_count": len(matrix["scale_orbit_rows"]),
            "all_units_rescale_rows": sum(1 for row in matrix["scale_orbit_rows"] if row["all_units_rescale_by_c"]),
            "absolute_unit_fixed_rows": sum(1 for row in matrix["scale_orbit_rows"] if row["absolute_unit_fixed"]),
            "accepted_as_strict_observable_unit_source": matrix["accepted_as_strict_observable_unit_source"],
        },
        "decision": {
            "breakthrough": "Four observer-independent observable-unit candidates for T_K were constructed and are positive.  The obstruction is that every candidate is homogeneous in the observable: rescaling T_K rescales the unit, leaving only an internal calibration convention rather than an absolute physical readout unit.",
            "negative_export_flags": {k: False for k in ["strict_observable_unit_source_exported", "strict_lambda_action_unit_source_exported", "unit_bearing_action_eom_source_exported", "clock_unit_exported", "hamiltonian_exported", "time_arrow_exported", "ltotal_exported", "observed_physics_exported", "selector_source_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay internal observable-unit normalizations.  The next proof-grade move should attack exactly one remaining independent unit atom: a strict clock-unit theorem for the typed time observable, or a strict action quantum/reference-cell source coupled to the already typed action; otherwise introduce a genuinely new strict time-order object with its own unit theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3019/S1969 observable-unit readout source obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- positive candidate units / total: `{c['positive_candidate_units']}/{c['candidate_unit_count']}`",
        f"- all-units-rescale rows / total: `{c['all_units_rescale_rows']}/{c['scale_row_count']}`",
        f"- absolute-unit-fixed rows / total: `{c['absolute_unit_fixed_rows']}/{c['scale_row_count']}`",
        f"- accepted as strict observable-unit source: `{c['accepted_as_strict_observable_unit_source']}`", "",
        "## Decision", payload["decision"]["breakthrough"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    read_json(P3018)
    payload = build_payload(P3018)
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3019/S1969 observable-unit readout source obstruction", "## P3019/S1969 observable-unit readout source obstruction\n\n`P3019/S1969` attacks exactly one P3018 independent unit atom: strict observable-unit readout for the typed time observable `T_K`.  It constructs four positive observer-independent unit candidates from the finite `T_K` vector: RMS, mean absolute value, `L∞`, and total variation.  The finite obstruction is rescaling-orbit dependence: under `T_K -> c T_K` every candidate unit rescales by `c`, so the normalized readout is only an internal calibration convention.  No independent strict detector/readout unit, clock/Hamiltonian coupling, unit-bearing EOM, `L_total`, observed-physics export, selector, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3019/S1969 observable-unit readout `L_total` guard", "## P3019/S1969 observable-unit readout `L_total` guard\n\n`P3019/S1969` adds no physical `L_total` term.  Its RMS, mean-absolute, `L∞`, and total-variation units normalize the typed observable `T_K` only internally; because each unit rescales with `T_K -> c T_K`, no strict observable-unit theorem, clock/Hamiltonian coupling, action quantum, or unit-bearing EOM is installed.\n")
    append_once(AGENTS, "Current observable-unit readout source guardrail (P3019/S1969, 2026-06-22)", "## Current observable-unit readout source guardrail (P3019/S1969, 2026-06-22)\n\n- P3019 attacks one P3018 independent unit atom: strict observable-unit readout for the typed time observable `T_K`.\n- Four finite observer-independent candidate units are positive (RMS, mean absolute value, `L∞`, total variation), but all rescale with `T_K -> c T_K`; they are internal calibration choices, not an absolute strict physical readout unit.\n- Do not promote internal observable-unit normalizations to strict observable-unit source, unit-bearing EOM/Hamiltonian, time arrow, `L_total`, observed-physics, selector, bridge/role-transfer, or ToE closure.\n- The next honest move should attack exactly one remaining independent unit atom: a strict clock-unit theorem, a strict action quantum/reference-cell source coupled to the typed action, or a genuinely new strict time-order object with its own unit theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
