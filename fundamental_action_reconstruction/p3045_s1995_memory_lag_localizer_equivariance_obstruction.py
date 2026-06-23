#!/usr/bin/env python3
"""P3045/S1995: memory-lag localizer equivariance obstruction.

P3044 produced a concrete memory-lag commutator
    C_s(i)=K_i*M_{i+s}-M_i*K_{i+s}
and left exactly one admissible next premise: a chart-independent lag/localizer
for that new object.  P3045 attacks only that premise.

The finite test is deliberately source-law conservative.  It expands P3044's
positive lags to the oriented lag torsor s in Z12 minus {0}, verifies cyclic-origin
invariance of integrated commutator scores, and tests Aut(Z12) inversion.  The
same branch score that selects an oriented positive lag is inverted by unit 11:
S_{-s}=-S_s.  Thus receiver scores can choose a lag in the current chart, but
they do not export an Aut-compatible nonpremise lag/localizer theorem.
"""
from __future__ import annotations

import hashlib, json
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3023_s1973_kernel_dissipation_time_order_candidate_obstruction import N
from p3044_s1994_memory_lag_commutator_source_candidate import OUT as P3044, commutator, kernel_vector
from p3038_s1988_viscous_retarded_chiral_selector_candidate_obstruction import memory_viscosity_trace

OUT = GEN / "p3045_s1995_memory_lag_localizer_equivariance_obstruction.json"
MD = GEN / "p3045_s1995_memory_lag_localizer_equivariance_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
TOL = 1e-12
UNITS = [1, 5, 7, 11]
ORIENTED_LAGS = [s for s in range(-N // 2 + 1, N // 2 + 1) if s != 0]


def signed_score(lag: int, shift: int = 0) -> float:
    k = kernel_vector()
    m = memory_viscosity_trace(k)
    if shift:
        k = k[shift:] + k[:shift]
        m = m[shift:] + m[:shift]
    return sum(commutator(k, m, lag % N))


def abs_lag(lag: int) -> int:
    r = lag % N
    return min(r, (-r) % N)


def aut_lag(unit: int, lag: int) -> int:
    r = (unit * lag) % N
    if r > N // 2:
        r -= N
    return r


def oriented_lag_rows() -> list[dict[str, Any]]:
    rows = []
    for lag in ORIENTED_LAGS:
        score = signed_score(lag)
        inv_lag = aut_lag(11, lag)
        inv_score = signed_score(inv_lag)
        shift_scores = [round(signed_score(lag, shift), 15) for shift in range(N)]
        rows.append({
            "oriented_lag": lag,
            "absolute_lag": abs_lag(lag),
            "signed_score": round(score, 15),
            "inversion_lag": inv_lag,
            "inversion_signed_score": round(inv_score, 15),
            "inversion_flips_score": abs(score + inv_score) < TOL,
            "translation_score_orbit_size": len(set(shift_scores)),
            "translation_origin_invariant": max(abs(x - shift_scores[0]) for x in shift_scores) < TOL,
        })
    return rows


def localizer_candidate_rows(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    positive = [r for r in rows if r["oriented_lag"] > 0]
    candidates = []
    specs = [
        ("max_signed_score", max(positive, key=lambda r: r["signed_score"])),
        ("max_abs_signed_score", max(rows, key=lambda r: abs(r["signed_score"]))),
        ("max_positive_absolute_lag_score", max(positive, key=lambda r: (r["signed_score"], -r["absolute_lag"]))),
        ("translation_integrated_receiver", max(positive, key=lambda r: r["signed_score"])),
    ]
    for name, winner in specs:
        inv_lag = winner["inversion_lag"]
        inv_score = next(r["signed_score"] for r in rows if r["oriented_lag"] == inv_lag)
        candidates.append({
            "candidate": name,
            "winner_oriented_lag": winner["oriented_lag"],
            "winner_absolute_lag": winner["absolute_lag"],
            "winner_signed_score": winner["signed_score"],
            "inversion_image_lag": inv_lag,
            "inversion_image_signed_score": inv_score,
            "translation_invariant_score": winner["translation_origin_invariant"],
            "aut_inversion_compatible": winner["oriented_lag"] == inv_lag and abs(winner["signed_score"] - inv_score) < TOL,
            "accepted_chart_independent_lag_localizer": False,
            "failure": "receiver selects an oriented lag only after a chart polarity convention; Aut inversion sends the signed row to the opposite lag/sign",
        })
    return candidates


def build_matrix() -> dict[str, Any]:
    read_json(P3044)
    rows = oriented_lag_rows()
    candidates = localizer_candidate_rows(rows)
    obligations = [
        {"obligation": "p3044_commutator_used_without_replay", "satisfied": True, "detail": "P3045 attacks only the lag/localizer premise left open by P3044"},
        {"obligation": "oriented_lag_torsor_constructed", "satisfied": True, "detail": "finite oriented lag torsor Z12\\{0} is explicit"},
        {"obligation": "translation_origin_invariance_checked", "satisfied": all(r["translation_origin_invariant"] for r in rows), "detail": "integrated scores are invariant under cyclic source-origin shifts"},
        {"obligation": "inversion_action_checked", "satisfied": all(r["inversion_flips_score"] for r in rows), "detail": "Aut unit 11 sends S_s to S_-s=-S_s"},
        {"obligation": "aut_compatible_nonpremise_lag_localizer", "satisfied": False, "detail": "all receiver winners require a polarity/chart convention under inversion"},
        {"obligation": "strict_source_law_for_lag", "satisfied": False, "detail": "no strict nadsoliton law exports a lag or lag polarity"},
        {"obligation": "selector_or_readout_coupling", "satisfied": False, "detail": "no theorem couples the lag torsor to QW-2191, selector torsor, or physical readout"},
    ]
    return {
        "object": "MemoryLagLocalizer_EquivarianceObstructionMatrix",
        "tested_premise": "chart-independent lag/localizer theorem for C_s(i)=K_i*M_{i+s}-M_i*K_{i+s}",
        "aut_units": UNITS,
        "oriented_lag_rows": rows,
        "localizer_candidate_rows": candidates,
        "proof_obligations": obligations,
        "finite_certificate": {
            "oriented_lag_rows": len(rows),
            "translation_invariant_rows": sum(1 for r in rows if r["translation_origin_invariant"]),
            "inversion_flip_rows": sum(1 for r in rows if r["inversion_flips_score"]),
            "localizer_candidate_rows": len(candidates),
            "accepted_chart_independent_lag_localizer_rows": sum(1 for r in candidates if r["accepted_chart_independent_lag_localizer"]),
            "aut_compatible_candidate_rows": sum(1 for r in candidates if r["aut_inversion_compatible"]),
            "proof_obligations": len(obligations),
            "satisfied_proof_obligations": sum(1 for r in obligations if r["satisfied"]),
            "chart_independent_lag_localizer_exported": False,
        },
    }


def build_payload() -> dict[str, Any]:
    matrix = build_matrix()
    return {
        "status": "P3045_MEMORY_LAG_LOCALIZER_EQUIVARIANCE_OBSTRUCTION_BOUNDED_NO_EXPORT",
        "input_hashes": {"P3044": hashlib.sha256(P3044.read_bytes()).hexdigest() if P3044.exists() else None},
        "constructed_theoretical_objects": matrix,
        "finite_certificate": matrix["finite_certificate"],
        "decision": {
            "bounded_no_go": "P3045 constructs the oriented lag torsor for the P3044 memory-lag commutator and verifies that integrated scores are cyclic-origin invariant.  The same computation shows the obstruction: Aut inversion sends each oriented lag score to the opposite lag with opposite sign, so receiver winners are chart-polarity choices rather than a nonpremise lag/localizer theorem.",
            "negative_export_flags": {k: False for k in ["chart_independent_lag_localizer_exported", "strict_memory_lag_source_law_exported", "selector_or_readout_coupling_exported", "qw2191_discharged", "strict_selector_closure_exported", "observed_physics_exported", "unit_bearing_action_eom_hamiltonian_exported", "ltotal_exported", "bridge_closure_exported", "role_transfer_exported", "toe_closure_exported"]},
            "next_honest_step": "Do not promote translation-invariant lag scores or lag-2 winners to selector closure.  A next proof-grade move may attack exactly one remaining P3044 premise: either a strict nadsoliton source law for the memory-lag commutator/lag polarity, or an explicit coupling theorem from this signed lag torsor to a selector torsor/readout row.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3045/S1995 memory-lag localizer equivariance obstruction", "",
        f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- oriented lag rows: `{c['oriented_lag_rows']}`",
        f"- translation-invariant rows: `{c['translation_invariant_rows']}`",
        f"- inversion-flip rows: `{c['inversion_flip_rows']}`",
        f"- localizer candidate rows: `{c['localizer_candidate_rows']}`",
        f"- accepted chart-independent lag-localizer rows: `{c['accepted_chart_independent_lag_localizer_rows']}`",
        f"- Aut-compatible candidate rows: `{c['aut_compatible_candidate_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        f"- chart-independent lag localizer exported: `{c['chart_independent_lag_localizer_exported']}`", "",
        "## Decision", payload["decision"]["bounded_no_go"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3045/S1995 memory-lag localizer equivariance obstruction", "## P3045/S1995 memory-lag localizer equivariance obstruction\n\n`P3045/S1995` attacks exactly one P3044 missing premise: a chart-independent lag/localizer theorem for the memory-lag commutator `C_s(i)=K_i*M_{i+s}-M_i*K_{i+s}`.  It constructs the oriented lag torsor, verifies cyclic-origin invariance of integrated scores, and shows the obstruction: Aut inversion sends each signed lag row to the opposite lag/sign.  Therefore lag-score winners, including the finite lag-2 receiver winner, are chart-polarity choices rather than strict localizer exports.  No strict source law, selector/readout coupling, `QW-2191` discharge, unit-bearing action/EOM, `L_total`, bridge/role transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3045/S1995 memory-lag localizer `L_total` guard", "## P3045/S1995 memory-lag localizer `L_total` guard\n\n`P3045/S1995` adds no physical `L_total` term.  Translation-invariant integrated lag scores are finite diagnostics, but Aut inversion flips the oriented lag sign and no strict lag-polarity source, selector/readout coupling, or unit-bearing variational/action/EOM installation is exported.\n")
    append_once(AGENTS, "Current memory-lag localizer equivariance guardrail (P3045/S1995, 2026-06-23)", "## Current memory-lag localizer equivariance guardrail (P3045/S1995, 2026-06-23)\n\n- P3045 attacks exactly one P3044 missing premise: chart-independent lag/localizer for the memory-lag commutator.\n- The finite oriented-lag torsor has cyclic-origin-invariant integrated scores, but Aut inversion maps each signed lag score to the opposite lag/sign; lag-score winners remain chart-polarity receiver choices.\n- Do not promote translation-invariant lag scores, lag-2 winners, oriented lag polarity, or memory-lag localizer receivers to `QW-2191` discharge, selector closure, observed physics, unit-bearing action/EOM, `L_total`, bridge/role-transfer, or ToE closure.\n- A next move may attack exactly one remaining P3044 premise: strict nadsoliton source law for the commutator/lag polarity, or explicit selector/readout coupling theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
