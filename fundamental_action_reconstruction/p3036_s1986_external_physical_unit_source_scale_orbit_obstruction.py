#!/usr/bin/env python3
"""P3036/S1986: external physical unit-source scale-orbit obstruction.

Attack exactly one P3028 foundation atom: external physical unit source.  The
new theoretical object is a two-axis scale torsor for sampled K_strict_gate on
Z12: amplitude scaling K -> cK and label-length scaling ell -> a ell.  Finite
candidate unit normalizations can choose dimensionless representatives, but none
exports a physical length/action/clock unit independent of the scale orbit.
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
from p3035_s1985_z12_directional_branch_selector_source_obstruction import OUT as P3035

OUT = GEN / "p3036_s1986_external_physical_unit_source_scale_orbit_obstruction.json"
MD = GEN / "p3036_s1986_external_physical_unit_source_scale_orbit_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

AMP_SCALE = 3.0
LENGTH_SCALE = 5.0
D_F = 9.0 / 5.0


def kernel_vector(scale: float = 1.0) -> list[float]:
    return [scale * k_strict(label) for label in range(1, N + 1)]


def l2_norm(values: list[float]) -> float:
    return math.sqrt(sum(v * v for v in values))


def dirichlet_action(values: list[float], label_length: float = 1.0) -> float:
    # If label spacing is assigned length a, a continuum-like Dirichlet density
    # scales as (Delta K / a)^2 * a = (Delta K)^2 / a on the finite cycle.
    grad = [values[(i + 1) % N] - values[i] for i in range(N)]
    return 0.5 * sum(g * g for g in grad) / label_length


def entropy_scale_candidate(bits: int = 1, h0: float = 0.0) -> float:
    return math.exp((bits * math.log(2.0) - h0) / D_F)


def build_matrix() -> dict[str, Any]:
    base = kernel_vector()
    amp = kernel_vector(AMP_SCALE)
    rows = [
        {
            "candidate": "label_step_d_equals_one_length_unit",
            "base_value": 1.0,
            "after_label_length_rescale": LENGTH_SCALE,
            "dimensionless_representative_available": True,
            "scale_orbit_invariant": False,
            "accepted_as_external_physical_unit_source": False,
            "failure": "choosing d=1 fixes a chart label step, not an external physical length unit",
        },
        {
            "candidate": "kernel_l2_norm_set_to_one_amplitude_unit",
            "base_value": round(l2_norm(base), 12),
            "after_amplitude_rescale": round(l2_norm(amp), 12),
            "rescaling_exponent": round(math.log(l2_norm(amp) / l2_norm(base), AMP_SCALE), 12),
            "dimensionless_representative_available": True,
            "scale_orbit_invariant": False,
            "accepted_as_external_physical_unit_source": False,
            "failure": "normalizing ||K||_2 to one is an amplitude gauge choice; it rescales with K -> cK",
        },
        {
            "candidate": "dirichlet_action_set_to_one_action_unit",
            "base_value": round(dirichlet_action(base), 12),
            "after_amplitude_rescale": round(dirichlet_action(amp), 12),
            "after_label_length_rescale": round(dirichlet_action(base, LENGTH_SCALE), 12),
            "amplitude_rescaling_exponent": round(math.log(dirichlet_action(amp) / dirichlet_action(base), AMP_SCALE), 12),
            "length_rescaling_exponent": round(math.log(dirichlet_action(base, LENGTH_SCALE) / dirichlet_action(base), LENGTH_SCALE), 12),
            "dimensionless_representative_available": True,
            "scale_orbit_invariant": False,
            "accepted_as_external_physical_unit_source": False,
            "failure": "action-to-one imports both amplitude and label-length gauges; no external action quantum is sourced",
        },
        {
            "candidate": "entropy_one_bit_scale_reference_cell",
            "base_value": round(entropy_scale_candidate(), 12),
            "after_reference_cell_shift_h0_plus_log2": round(entropy_scale_candidate(h0=math.log(2.0)), 12),
            "dimensionless_representative_available": True,
            "scale_orbit_invariant": False,
            "accepted_as_external_physical_unit_source": False,
            "failure": "one-bit entropy scale is conditional on an unsourced reference cell and bit-to-length/action map",
        },
    ]
    obligations = [
        {"obligation": "attacks_single_P3028_foundation_atom", "satisfied": True, "detail": "only external physical unit source is tested"},
        {"obligation": "explicit_scale_torsor", "satisfied": True, "detail": "amplitude K->cK and label-length ell->a ell axes are represented"},
        {"obligation": "finite_candidate_unit_receivers_computable", "satisfied": True, "detail": "label, L2, action, and entropy-reference candidates are computed"},
        {"obligation": "dimensionful_external_unit_export", "satisfied": False, "detail": "no candidate carries physical dimension beyond a normalized representative"},
        {"obligation": "amplitude_scale_orbit_quotient", "satisfied": False, "detail": "L2/action candidates move under K->cK"},
        {"obligation": "label_length_scale_orbit_quotient", "satisfied": False, "detail": "d=1/action candidates move under ell->a ell"},
        {"obligation": "reference_cell_and_bit_to_unit_map", "satisfied": False, "detail": "entropy candidate remains conditional on reference-cell and bit-to-unit premises"},
        {"obligation": "coupling_to_unit_bearing_action_or_readout", "satisfied": False, "detail": "no theorem couples a unit source to action/EOM/Hamiltonian or classical readout rows"},
    ]
    return {
        "object": "ExternalPhysicalUnitSource_ScaleOrbitObstructionMatrix",
        "tested_foundation_atom": "external_physical_unit_source",
        "scale_torsor": {"amplitude_axis": "K -> cK", "label_length_axis": "ell -> a ell", "tested_c": AMP_SCALE, "tested_a": LENGTH_SCALE},
        "candidate_rows": rows,
        "proof_obligations": obligations,
        "finite_certificate": {
            "candidate_rows": len(rows),
            "dimensionless_representative_rows": sum(1 for row in rows if row["dimensionless_representative_available"]),
            "scale_orbit_invariant_rows": sum(1 for row in rows if row["scale_orbit_invariant"]),
            "accepted_external_unit_source_rows": sum(1 for row in rows if row["accepted_as_external_physical_unit_source"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "external_physical_unit_source_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_external_physical_unit_source"] for row in rows),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3028)
    read_json(P3035)
    matrix = build_matrix()
    return {
        "status": "P3036_EXTERNAL_PHYSICAL_UNIT_SOURCE_SCALE_ORBIT_OBSTRUCTION_NO_UNIT_EXPORT",
        "input_hashes": {
            "P3028": hashlib.sha256(P3028.read_bytes()).hexdigest() if P3028.exists() else None,
            "P3035": hashlib.sha256(P3035.read_bytes()).hexdigest() if P3035.exists() else None,
        },
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "A two-axis scale torsor for external unit sourcing is explicit: amplitude K->cK and label-length ell->a ell. Four finite unit candidates are computable, but each only fixes a dimensionless representative or conditional reference cell. None exports a physical length/action/clock unit or couples that unit to action/EOM/Hamiltonian/readout rows.",
            "negative_export_flags": {k: False for k in ["external_physical_unit_source_exported", "unit_bearing_action_eom_hamiltonian_exported", "classical_transition_exported", "observed_physics_exported", "strict_selector_branch_source_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay d=1 labels, norm-to-one gauges, action-to-one gauges, or conditional entropy reference cells as external physical unit sources. A next move should reconcile the P3028/P3029-P3036 foundation atoms in a no-new-live-frontier certificate, unless a genuinely new unit theorem or readout-coupling source is supplied.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3036/S1986 external physical unit-source scale-orbit obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- candidate rows: `{c['candidate_rows']}`",
        f"- dimensionless representative rows: `{c['dimensionless_representative_rows']}`",
        f"- scale-orbit invariant rows: `{c['scale_orbit_invariant_rows']}`",
        f"- accepted external unit-source rows: `{c['accepted_external_unit_source_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- external physical unit source exported: `{c['external_physical_unit_source_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3036/S1986 external physical unit-source scale-orbit obstruction", "## P3036/S1986 external physical unit-source scale-orbit obstruction\n\n`P3036/S1986` attacks exactly one P3028 foundation atom: external physical unit source.  It constructs a two-axis scale torsor for sampled `K_strict_gate` on `Z12`: amplitude scaling `K -> cK` and label-length scaling `ell -> a ell`.  Four finite candidates are computable (`d=1`, `||K||_2=1`, Dirichlet action-to-one, and an entropy one-bit reference-cell scale), but all remain dimensionless representatives or conditional gauges.  No scale-orbit-invariant physical length/action/clock unit, unit-bearing action/EOM/Hamiltonian coupling, classical transition, `L_total`, bridge/role transfer, or ToE closure follows.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3036/S1986 external unit-source `L_total` guard", "## P3036/S1986 external unit-source `L_total` guard\n\n`P3036/S1986` adds no physical `L_total` term.  Label-step, norm-to-one, action-to-one, and entropy reference-cell candidates fix internal representatives or conditional gauges, not an external physical unit coupled to a variational/action/EOM theorem.\n")
    append_once(AGENTS, "Current external physical unit-source scale-orbit guardrail (P3036/S1986, 2026-06-23)", "## Current external physical unit-source scale-orbit guardrail (P3036/S1986, 2026-06-23)\n\n- P3036 attacks exactly one P3028 foundation atom: external physical unit source.\n- The two-axis scale torsor `K -> cK` and `ell -> a ell` is explicit; `d=1`, `||K||_2=1`, action-to-one, and entropy one-bit reference-cell candidates are computable but remain internal representatives or conditional gauges.\n- Do not promote label steps, norm-to-one gauges, action-to-one gauges, entropy reference-cell choices, or internal scale representatives to physical length/action/clock units, unit-bearing action/EOM, classical transition, observed physics, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move should reconcile the P3028/P3029-P3036 foundation atoms in a no-new-live-frontier certificate unless a genuinely new unit theorem or readout-coupling source is supplied.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
