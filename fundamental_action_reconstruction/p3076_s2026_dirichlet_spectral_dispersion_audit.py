#!/usr/bin/env python3
"""P3076/S2026: Dirichlet/Laplacian spectral-dispersion audit.

P3075 exported a scoped internal local Dirichlet gradient source for the accepted
cycle-Laplacian scale-flow rows.  P3076 constructs the next missing object: a
finite spectral-dispersion interface for that source on Z12.  The audit is
honest and bounded: it diagonalizes the Z12 cycle Laplacian by exact Fourier-mode
labels, tests the first-order gradient-flow amplification factors at the P3074
fractional steps, and separates diffusive/internal smoothing behavior from a
hypothetical lightlike/wave-compatible branch.  It exports no observed light,
spacetime EOM, gauge photon sector, unit system, or empirical physics.
"""
from __future__ import annotations

import hashlib, json, math, subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3075_s2025_variational_source_obstruction_table import OUT as P3075

OUT = GEN / "p3076_s2026_dirichlet_spectral_dispersion_audit.json"
MD = GEN / "p3076_s2026_dirichlet_spectral_dispersion_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12_SIZE = 12
STEPS = (Fraction(1, 12), Fraction(1, 6), Fraction(1, 4))
CONTENT_PATTERNS = {
    "dirichlet_laplacian_source": r"local Dirichlet|cycle-Laplacian|Laplacian flow|gradient source",
    "spectral_dispersion_obligation": r"spectral-dispersion|continuum-limit|mode spectrum|small-k dispersion|wave-compatible|diffusive",
    "no_observed_physics_promotion": r"observed light|gauge photon|spacetime EOM|unit-bearing|empirical physics|L_total|ToE|selector closure",
}


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def eigenvalue_label(j: int) -> str:
    # lambda_j = 2 - 2 cos(2*pi*j/12); exact labels by Z12 symmetry.
    exact = {
        0: "0", 1: "2-sqrt(3)", 2: "1", 3: "2", 4: "3", 5: "2+sqrt(3)", 6: "4",
        7: "2+sqrt(3)", 8: "3", 9: "2", 10: "1", 11: "2-sqrt(3)",
    }
    return exact[j]


def eigenvalue_float(j: int) -> float:
    return 2.0 - 2.0 * math.cos(2.0 * math.pi * j / Z12_SIZE)


def mode_rows() -> list[dict[str, Any]]:
    rows = []
    for j in range(Z12_SIZE):
        lam = eigenvalue_float(j)
        q = min(j, Z12_SIZE - j)
        rows.append({
            "mode_j": j,
            "undirected_wavenumber_q": q,
            "multiplicity_role": "constant_kernel" if j == 0 else ("nyquist_singlet" if j == 6 else "paired_cos_sin_doublet"),
            "laplacian_eigenvalue_exact": eigenvalue_label(j),
            "laplacian_eigenvalue_decimal": round(lam, 12),
            "gradient_flow_generator_eigenvalue": f"-{eigenvalue_label(j)}",
            "small_k_quadratic_proxy": "lambda_j ~ (2*pi*j/12)^2 for j << 12" if j in (0, 1, 11) else "finite Z12 mode; not a continuum small-k row",
            "wave_branch_if_second_order_extra_structure_added": f"omega_j^2 = {eigenvalue_label(j)}",
            "accepted_as_lightlike_branch": False,
            "blocked_by": "P3075 source is first-order dissipative gradient flow, not a second-order/unit-bearing wave equation with time coordinate and Lorentzian metric",
        })
    return rows


def amplification_rows() -> list[dict[str, Any]]:
    rows = []
    for j in range(Z12_SIZE):
        lam = eigenvalue_float(j)
        for eps in STEPS:
            amp = 1.0 - float(eps) * lam
            rows.append({
                "mode_j": j,
                "step_fraction": f"{eps.numerator}/{eps.denominator}",
                "eigenvalue_exact": eigenvalue_label(j),
                "explicit_gradient_amplification_decimal": round(amp, 12),
                "nonexpansive_on_this_step": abs(amp) <= 1.0 + 1e-12,
                "strictly_contracting_nonconstant_mode": j != 0 and abs(amp) < 1.0,
                "oscillatory_phase_generated": False,
            })
    return rows


def build_payload() -> dict[str, Any]:
    p3075 = read_json(P3075)
    grep_rows = content_grep()
    modes = mode_rows()
    amps = amplification_rows()
    nonconstant_modes = [r for r in modes if r["mode_j"] != 0]
    proof_obligations = [
        {"obligation": "content_first_grep_before_spectral_audit", "satisfied": True, "detail": "searched Dirichlet/Laplacian source, spectral-dispersion obligation, and no-physics-promotion lanes"},
        {"obligation": "construct_z12_dirichlet_spectral_interface", "satisfied": True, "detail": "12 exact Fourier-mode eigenvalue labels for the cycle Laplacian are tabulated"},
        {"obligation": "test_fractional_step_gradient_amplifications", "satisfied": True, "detail": "12 modes x 3 P3074-style fractional steps = 36 amplification rows"},
        {"obligation": "export_internal_diffusive_smoothing_branch", "satisfied": all(r["strictly_contracting_nonconstant_mode"] for r in amps if r["mode_j"] != 0), "detail": "all nonconstant modes are contracted for steps 1/12, 1/6, and 1/4"},
        {"obligation": "export_lightlike_or_wave_branch", "satisfied": False, "detail": "the audited object is first-order gradient flow; a wave branch would require a new second-order/unit/time/metric source theorem"},
        {"obligation": "export_observed_light_or_gauge_photon_sector", "satisfied": False, "detail": "no spacetime embedding, gauge field, units, or empirical map is constructed"},
    ]
    return {
        "status": "P3076_INTERNAL_DIRICHLET_DIFFUSIVE_SPECTRAL_BRANCH_WAVE_LIGHTLIKE_OBSTRUCTION",
        "input_hashes": {"P3075": hashlib.sha256(P3075.read_bytes()).hexdigest() if P3075.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "spectral_dispersion_interface": {
                "object": "Z12DirichletSpectralDispersionAudit",
                "source_reused": "P3075 local Dirichlet negative-gradient source for accepted cycle-Laplacian rows",
                "diagonal_basis": "Z12 Fourier modes j=0..11 with lambda_j = 2 - 2 cos(2*pi*j/12)",
                "acceptance_predicate_for_lightlike_branch": "requires oscillatory/second-order or unit-bearing time-coordinate source, not merely a first-order dissipative gradient eigenvalue",
            },
            "mode_spectrum_rows": modes,
            "fractional_step_amplification_rows": amps,
            "representative_nonconstant_rows": nonconstant_modes[:6],
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "p3075_local_dirichlet_accepted_rows": p3075.get("finite_certificate", {}).get("local_dirichlet_accepted_rows"),
            "z12_modes": len(modes),
            "nonconstant_modes": len(nonconstant_modes),
            "fractional_steps": len(STEPS),
            "amplification_rows": len(amps),
            "nonexpansive_amplification_rows": sum(1 for r in amps if r["nonexpansive_on_this_step"]),
            "strictly_contracting_nonconstant_amplification_rows": sum(1 for r in amps if r["strictly_contracting_nonconstant_mode"]),
            "oscillatory_phase_rows": sum(1 for r in amps if r["oscillatory_phase_generated"]),
            "accepted_lightlike_branch_rows": sum(1 for r in modes if r["accepted_as_lightlike_branch"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3076 diagonalizes the P3075 local Dirichlet/Laplacian source on Z12 and finds an internal diffusive smoothing branch: the constant mode is neutral and all 11 nonconstant modes contract under the three audited fractional steps.  The spectrum has a formal lambda_j ~ k^2 small-k proxy, but the exported dynamics is first-order dissipative gradient flow, not a second-order, unit-bearing, Lorentzian, wave/lightlike equation.",
            "negative_export_flags": {k: False for k in ["lightlike_branch_exported", "observed_light_exported", "gauge_photon_sector_exported", "strict_spacetime_embedding_exported", "unit_bearing_time_coordinate_exported", "second_order_wave_eom_exported", "empirical_physics_exported", "hamiltonian_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"internal_dirichlet_spectrum_diagonalized": True, "internal_diffusive_smoothing_branch_exported": True, "wave_lightlike_obstruction_table_executed": True},
            "next_honest_step": "Construct one bounded second-order lift obstruction table for the same Z12 Dirichlet source: add the minimal candidate phase-space variables (rho, pi), symplectic form, and Hamiltonian H = 1/2*pi^2 + E_D, then test exactly which premises are new and whether they are internally sourced or merely imported.  Keep it scoped as a missing-source audit; do not promote it to observed light, gauge photons, spacetime EOM, units, empirical physics, L_total, bridge/role-transfer, selector closure, or ToE.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3076/S2026 Dirichlet spectral-dispersion audit", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3075 local Dirichlet accepted rows: `{c['p3075_local_dirichlet_accepted_rows']}`",
        f"- Z12 modes: `{c['z12_modes']}`",
        f"- nonconstant modes: `{c['nonconstant_modes']}`",
        f"- fractional steps: `{c['fractional_steps']}`",
        f"- amplification rows: `{c['amplification_rows']}`",
        f"- nonexpansive amplification rows: `{c['nonexpansive_amplification_rows']}`",
        f"- strictly contracting nonconstant amplification rows: `{c['strictly_contracting_nonconstant_amplification_rows']}`",
        f"- oscillatory phase rows: `{c['oscillatory_phase_rows']}`",
        f"- accepted lightlike branch rows: `{c['accepted_lightlike_branch_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3076/S2026 Dirichlet spectral-dispersion audit", "## P3076/S2026 Dirichlet spectral-dispersion audit\n\n`P3076/S2026` diagonalizes the `P3075` local Dirichlet/Laplacian source on `Z12` by Fourier-mode eigenvalues `lambda_j = 2 - 2 cos(2*pi*j/12)`.  The finite table has `12` spectral modes and `36` fractional-step amplification rows.  The constant mode is neutral and all nonconstant rows contract for steps `1/12`, `1/6`, and `1/4`, yielding an internal diffusive smoothing branch.  The formal small-`k` `lambda ~ k^2` proxy is not a lightlike branch: no second-order wave equation, unit-bearing time coordinate, Lorentzian metric, observed-light/gauge sector, empirical physics, `L_total`, selector closure, bridge/role-transfer, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3076/S2026 first-order spectrum is not wave EOM", "## P3076/S2026 first-order spectrum is not wave EOM\n\n`P3076/S2026` gives an exact finite spectral audit of the internal Dirichlet gradient source.  Its accepted dynamics is first-order dissipative smoothing: Fourier modes have generator eigenvalues `-lambda_j` and no oscillatory phase rows.  A wave/Hamiltonian reading would require new phase-space, symplectic, second-order/time-coordinate, and unit-source premises that are not exported here.\n")
    append_once(AGENTS, "Current Dirichlet spectral-dispersion guardrail (P3076/S2026, 2026-06-24)", "## Current Dirichlet spectral-dispersion guardrail (P3076/S2026, 2026-06-24)\n\n- P3076 follows the P3075 recommendation and constructs a bounded spectral-dispersion audit for the internal local Dirichlet/Laplacian source on `Z12`.\n- The audit diagonalizes `12` Fourier modes and checks `36` fractional-step amplification rows.  The constant mode is neutral; all nonconstant modes contract for steps `1/12`, `1/6`, and `1/4`, so the exported branch is internal diffusive smoothing.\n- Do not promote the formal small-`k` `lambda ~ k^2` proxy to observed light, gauge photons, spacetime EOM, units, Hamiltonian dynamics, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one bounded second-order lift obstruction table for the same source, explicitly auditing whether phase-space/symplectic/Hamiltonian/time-coordinate premises are internally sourced or merely imported.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
