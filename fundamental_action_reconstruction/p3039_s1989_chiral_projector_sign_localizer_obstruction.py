#!/usr/bin/env python3
"""P3039/S1989: chiral projector sign/localizer obstruction.

Follow-up to P3038.  The integrated selector candidate used
chi_i = sin(2*pi*i/12) as an inversion-odd chiral projection.  This script
attacks exactly that missing premise: can the finite Z12 data source a
non-premise phase origin and polarity for this projector?

The result is deliberately bounded: the sine torsor is computable and has
inversion-odd sign, but translations move the phase origin and Aut(Z12)
inversion flips polarity.  K_strict_gate-weighted and memory-weighted phase
receivers pick chart representatives, not Aut-compatible nonpremise localizers.
"""
from __future__ import annotations

import hashlib, json, math
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N, k_strict
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import OUT as P3038, memory_viscosity_trace

OUT = GEN / "p3039_s1989_chiral_projector_sign_localizer_obstruction.json"
MD = GEN / "p3039_s1989_chiral_projector_sign_localizer_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-12
UNITS = [1, 5, 7, 11]


def labels() -> list[int]:
    return list(range(N))


def kernel_vector() -> list[float]:
    return [k_strict(i + 1) for i in labels()]


def chi(phase: int = 0, polarity: int = 1) -> list[float]:
    return [polarity * math.sin(2.0 * math.pi * (i - phase) / N) for i in labels()]


def transform_chi(values: list[float], unit: int, shift: int = 0) -> list[float]:
    # Pull back along i -> unit*i + shift on Z12.
    return [values[(unit * i + shift) % N] for i in labels()]


def max_abs_diff(a: list[float], b: list[float]) -> float:
    return max(abs(x - y) for x, y in zip(a, b))


def phase_score(weights: list[float], phase: int, polarity: int = 1) -> float:
    return sum(w * c for w, c in zip(weights, chi(phase, polarity)))


def best_phase_rows(weights: list[float], name: str) -> dict[str, Any]:
    rows = []
    for phase in labels():
        plus = phase_score(weights, phase, +1)
        minus = phase_score(weights, phase, -1)
        rows.append({"phase": phase, "plus_score": round(plus, 15), "minus_score": round(minus, 15), "abs_score": round(abs(plus), 15)})
    max_abs = max(row["abs_score"] for row in rows)
    winners = [row["phase"] for row in rows if abs(row["abs_score"] - max_abs) <= TOL]
    return {
        "receiver": name,
        "phase_rows": rows,
        "max_abs_score": round(max_abs, 15),
        "winner_phases": winners,
        "unique_chart_winner": len(winners) == 1,
        "accepted_as_nonpremise_localizer": False,
        "failure": "phase scoring can choose representatives in the current chart, but translations move the phase labels and polarity remains externally chosen",
    }


def build_matrix() -> dict[str, Any]:
    base = chi(0, +1)
    inversion = transform_chi(base, 11)
    inversion_flips = max_abs_diff(inversion, [-x for x in base]) < 1e-12
    translation_orbit_size = len({tuple(round(x, 12) for x in transform_chi(base, 1, s)) for s in labels()})
    polarity_torsor_size = len({tuple(round(x, 12) for x in chi(p, sig)) for p in labels() for sig in (+1, -1)})
    k = kernel_vector()
    mem = memory_viscosity_trace(k)
    abs_k = [abs(x) for x in k]
    receivers = [
        {
            "receiver": "bare_chiral_projector_torsor",
            "formula": "chi_{b,sigma}(i)=sigma*sin(2*pi*(i-b)/12)",
            "finite_nonzero": any(abs(x) > TOL for x in base),
            "inversion_odd": inversion_flips,
            "accepted_as_nonpremise_localizer": False,
            "failure": "the formula defines a full phase/polarity torsor; no row fixes b or sigma",
        },
        best_phase_rows(k, "K_strict_gate_weighted_phase_receiver"),
        best_phase_rows(mem, "memory_viscosity_weighted_phase_receiver"),
        best_phase_rows(abs_k, "absolute_kernel_weighted_phase_receiver"),
    ]
    obligations = [
        {"obligation": "exact_p3038_missing_premise_targeted", "satisfied": True, "detail": "only the chi_i phase-origin/sign-localizer premise is audited"},
        {"obligation": "finite_chiral_projector_torsor_constructed", "satisfied": True, "detail": f"{polarity_torsor_size} distinct phase/polarity projectors are enumerated"},
        {"obligation": "inversion_odd_sign_verified", "satisfied": inversion_flips, "detail": "unit 11 sends the base projector to its negative"},
        {"obligation": "current_hint_receivers_tested", "satisfied": True, "detail": "K, memory/viscosity, and absolute-K weighted phase receivers are scanned"},
        {"obligation": "nonpremise_phase_origin_localizer", "satisfied": False, "detail": "translations generate the full phase orbit; any selected phase is chart-relative"},
        {"obligation": "nonpremise_polarity_sign_source", "satisfied": False, "detail": "Aut inversion exchanges the two polarities and no strict sign law fixes one"},
        {"obligation": "aut_z12_compatible_selector_source", "satisfied": False, "detail": "candidate receivers are not invariant localizer theorems on the phase/polarity torsor"},
        {"obligation": "coupling_back_to_p3038_as_source", "satisfied": False, "detail": "without sourced phase and polarity, P3038 remains a branch-separating candidate, not a strict selector"},
    ]
    return {
        "object": "ChiralProjectorSignLocalizer_ObstructionMatrix",
        "target_formula": "chi_i = sin(2*pi*i/12) from P3038",
        "torsor": {"translation_phase_orbit_size": translation_orbit_size, "phase_polarity_projector_count": polarity_torsor_size, "aut_units": UNITS, "inversion_unit_11_flips_polarity": inversion_flips},
        "receiver_rows": receivers,
        "proof_obligations": obligations,
        "finite_certificate": {
            "receiver_rows": len(receivers),
            "finite_nonzero_rows": sum(1 for row in receivers if row.get("finite_nonzero", row.get("max_abs_score", 0.0) > TOL)),
            "inversion_odd_rows": sum(1 for row in receivers if row.get("inversion_odd", False)),
            "accepted_nonpremise_localizer_rows": sum(1 for row in receivers if row["accepted_as_nonpremise_localizer"]),
            "translation_phase_orbit_size": translation_orbit_size,
            "phase_polarity_projector_count": polarity_torsor_size,
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for row in obligations if row["satisfied"]),
            "nonpremise_chiral_sign_localizer_exported": all(row["satisfied"] for row in obligations) and all(row["accepted_as_nonpremise_localizer"] for row in receivers),
        },
    }


def build_payload() -> dict[str, Any]:
    read_json(P3038)
    matrix = build_matrix()
    return {
        "status": "P3039_CHIRAL_PROJECTOR_SIGN_LOCALIZER_OBSTRUCTION_NO_SOURCE_EXPORT",
        "input_hashes": {"P3038": hashlib.sha256(P3038.read_bytes()).hexdigest() if P3038.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "The P3038 chiral projector has a real finite inversion-odd torsor, and unit 11 flips its sign.  However translations move the phase origin, inversion exchanges polarity, and K/memory/absolute-K weighted receivers choose only chart-relative representatives.  Therefore no nonpremise chiral sign/localizer theorem is exported.",
            "negative_export_flags": {k: False for k in ["nonpremise_chiral_sign_localizer_exported", "p3038_promoted_to_selector_source", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not replay sine/chiral phase receivers as a selector source.  The next proof-grade move should attack one different P3038 missing source premise: either a sourced retardation path-anisotropy theorem for the c-retarded split, or a physical unit/readout coupling theorem for the integrated density.  Continue only with a genuinely new strict source law; otherwise preserve P3038-P3039 as a branch-separating-but-unsourced certificate.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3039/S1989 chiral projector sign/localizer obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- receiver rows: `{c['receiver_rows']}`",
        f"- finite nonzero rows: `{c['finite_nonzero_rows']}`",
        f"- inversion-odd rows: `{c['inversion_odd_rows']}`",
        f"- accepted nonpremise localizer rows: `{c['accepted_nonpremise_localizer_rows']}`",
        f"- translation phase orbit size: `{c['translation_phase_orbit_size']}`",
        f"- phase/polarity projector count: `{c['phase_polarity_projector_count']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- nonpremise chiral sign/localizer exported: `{c['nonpremise_chiral_sign_localizer_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3039/S1989 chiral projector sign/localizer obstruction", "## P3039/S1989 chiral projector sign/localizer obstruction\n\n`P3039/S1989` attacks exactly one missing P3038 premise: a nonpremise phase-origin/sign-localizer theorem for `chi_i = sin(2*pi*i/12)`.  The finite torsor is real and inversion-odd: the `U(12)` inversion unit sends the base projector to its negative.  But translations move the phase origin, inversion exchanges polarity, and `K_strict_gate`/memory/absolute-`K` weighted receivers choose only chart-relative representatives.  Thus P3039 does not export a nonpremise chiral sign/localizer, does not promote P3038 to selector closure, and does not discharge `QW-2191`, observed physics, `L_total`, bridge/role transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3039/S1989 chiral projector localizer `L_total` guard", "## P3039/S1989 chiral projector localizer `L_total` guard\n\n`P3039/S1989` adds no physical `L_total` term.  It verifies a finite inversion-odd chiral projector torsor, but fails to source its phase origin or polarity nonpremise; consequently the P3038 integrated density remains a candidate receiver rather than a unit-bearing variational term.\n")
    append_once(AGENTS, "Current chiral projector sign-localizer guardrail (P3039/S1989, 2026-06-23)", "## Current chiral projector sign-localizer guardrail (P3039/S1989, 2026-06-23)\n\n- P3039 attacks exactly one P3038 missing source premise: a nonpremise phase-origin/sign-localizer theorem for `chi_i = sin(2*pi*i/12)`.\n- The finite chiral projector torsor is real and inversion-odd, but translations move phase origins, Aut inversion exchanges polarity, and `K`/memory/absolute-`K` weighted phase receivers are chart-relative rather than strict source laws.\n- Do not promote sine/chiral projectors, phase-score winners, memory-weighted phase receivers, or inversion-oddness alone to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may attack a different P3038 missing premise: sourced retardation path-anisotropy theorem or physical unit/readout coupling theorem, but only with a genuinely new strict source law.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
