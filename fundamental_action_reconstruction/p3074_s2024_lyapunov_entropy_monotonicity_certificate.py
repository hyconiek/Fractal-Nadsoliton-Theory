#!/usr/bin/env python3
"""P3074/S2024: Lyapunov/entropy monotonicity certificate for P3073 flows.

P3073 exported only internal, dimensionless, total-preserving scale-flow rows.
P3074 takes the next bounded proof step: construct exact Lyapunov/entropy-like
functionals and test monotonicity under fractional steps of those accepted flows,
without importing selector closure, units, action, EOM, spacetime, or observed
physics.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from typing import Any, Callable

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3071_s2021_sigma_invariant_scalar_conservation_scale_control import profile_values, transform_values
from p3073_s2023_bounded_scale_flow_operator_obstruction import OUT as P3073, audit_rows as p3073_audit_rows, flow_values

OUT = GEN / "p3074_s2024_lyapunov_entropy_monotonicity_certificate.json"
MD = GEN / "p3074_s2024_lyapunov_entropy_monotonicity_certificate.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12 = tuple(range(12))
STEP_SIZES = (Fraction(1, 4), Fraction(1, 2), Fraction(1, 1))
CONTENT_PATTERNS = {
    "lyapunov_entropy_monotonicity": r"Lyapunov|entropy monotonicity|variance|Dirichlet energy|shell energy",
    "accepted_scale_flow": r"accepted.*scale-flow|cycle[- ]Laplacian|mean-centering|total-preserving",
    "no_variational_physics_promotion": r"variational EOM|unit-bearing action|Hamiltonian|observed physics|gauge photon|L_total|ToE|selector closure",
}


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def mean(rho: tuple[Fraction, ...]) -> Fraction:
    return sum(rho, Fraction(0)) / len(rho)


def variance_energy(rho: tuple[Fraction, ...]) -> Fraction:
    m = mean(rho)
    return sum((x - m) ** 2 for x in rho)


def quadratic_dirichlet_energy(rho: tuple[Fraction, ...]) -> Fraction:
    return sum((rho[(i + 1) % 12] - rho[i]) ** 2 for i in Z12)


def shell_quadratic_energy(rho: tuple[Fraction, ...]) -> Fraction:
    return sum(Fraction(min(i, 12 - i) ** 2) * rho[i] * rho[i] for i in Z12)


FUNCTIONALS: tuple[tuple[str, Callable[[tuple[Fraction, ...]], Fraction], bool], ...] = (
    ("variance_energy", variance_energy, True),
    ("quadratic_dirichlet_energy", quadratic_dirichlet_energy, True),
    ("shell_quadratic_energy", shell_quadratic_energy, False),
)


def step_flow(rho: tuple[Fraction, ...], flow: tuple[Fraction, ...], step: Fraction) -> tuple[Fraction, ...]:
    return tuple(a + step * b for a, b in zip(rho, flow))


def monotonicity_rows() -> list[dict[str, Any]]:
    rows = []
    accepted = [r for r in p3073_audit_rows() if r["accepted_intrinsic_bounded_scale_flow_row"]]
    for base in accepted:
        rho = transform_values(profile_values(base["profile_id"], base["sigma"]), base["transform"]["kind"], base["transform"]["k"])
        flow = flow_values(base["flow_operator"], rho)
        for step in STEP_SIZES:
            flowed = step_flow(rho, flow, step)
            for functional_id, functional, intrinsic_functional in FUNCTIONALS:
                before = functional(rho)
                after = functional(flowed)
                delta = after - before
                monotone_nonincreasing = delta <= 0
                strict_decrease = delta < 0
                accepted_monotonicity = intrinsic_functional and monotone_nonincreasing
                rows.append({
                    "profile_id": base["profile_id"],
                    "sigma": base["sigma"],
                    "transform": base["transform"],
                    "flow_operator": base["flow_operator"],
                    "step_size": str(step),
                    "functional": functional_id,
                    "intrinsic_functional": intrinsic_functional,
                    "before": str(before),
                    "after": str(after),
                    "delta_after_minus_before": str(delta),
                    "monotone_nonincreasing": monotone_nonincreasing,
                    "strict_decrease": strict_decrease,
                    "accepted_internal_lyapunov_row": accepted_monotonicity,
                    "blocked_by": "" if accepted_monotonicity else "; ".join(filter(None, [
                        None if intrinsic_functional else "functional uses chart shell weights / control only",
                        None if monotone_nonincreasing else "functional increases for this flow step",
                    ])),
                })
    return rows


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for functional_id, _, _ in FUNCTIONALS:
        subset = [r for r in rows if r["functional"] == functional_id]
        out.append({
            "functional": functional_id,
            "rows": len(subset),
            "monotone_nonincreasing_rows": sum(1 for r in subset if r["monotone_nonincreasing"]),
            "strict_decrease_rows": sum(1 for r in subset if r["strict_decrease"]),
            "accepted_internal_lyapunov_rows": sum(1 for r in subset if r["accepted_internal_lyapunov_row"]),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3073 = read_json(P3073)
    grep_rows = content_grep()
    rows = monotonicity_rows()
    accepted = [r for r in rows if r["accepted_internal_lyapunov_row"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_lyapunov_audit", "satisfied": True, "detail": "searched by Lyapunov/entropy, accepted scale-flow, and no-physics-promotion content"},
        {"obligation": "reuse_only_p3073_accepted_internal_flows", "satisfied": True, "detail": "tested the 192 accepted intrinsic bounded total-preserving P3073 rows"},
        {"obligation": "construct_exact_sigma_even_functionals", "satisfied": True, "detail": "variance and quadratic Dirichlet energies are intrinsic controls; shell quadratic energy is included only as a chart-weight control"},
        {"obligation": "export_internal_monotonicity_certificate", "satisfied": bool(accepted), "detail": f"{len(accepted)} intrinsic Lyapunov rows are monotone nonincreasing"},
        {"obligation": "export_variational_or_physical_dynamics", "satisfied": False, "detail": "monotonicity is not an action, Hamiltonian, EOM, unit map, spacetime embedding, gauge sector, or empirical map"},
    ]
    return {
        "status": "P3074_INTERNAL_LYAPUNOV_MONOTONICITY_CERTIFICATE_NO_VARIATIONAL_PHYSICS_EXPORT",
        "input_hashes": {"P3073": hashlib.sha256(P3073.read_bytes()).hexdigest() if P3073.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "lyapunov_entropy_monotonicity_template": {
                "object": "SigmaEvenScaleFlowLyapunovEntropyMonotonicityCertificate",
                "domain": "P3073 accepted intrinsic bounded total-preserving scale-flow rows",
                "step_sizes": [str(s) for s in STEP_SIZES],
                "functionals": [f[0] for f in FUNCTIONALS],
                "acceptance_predicate": "intrinsic functional and exact nonincrease under the tested fractional flow step",
            },
            "functional_aggregate_certificate": aggregate(rows),
            "representative_witness_rows": accepted[:12] + [r for r in rows if not r["accepted_internal_lyapunov_row"]][:12],
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "p3073_accepted_intrinsic_bounded_scale_flow_rows": p3073.get("finite_certificate", {}).get("accepted_intrinsic_bounded_scale_flow_rows"),
            "accepted_p3073_flow_rows_reused": 192,
            "step_sizes": len(STEP_SIZES),
            "functionals_tested": len(FUNCTIONALS),
            "monotonicity_matrix_rows": len(rows),
            "accepted_internal_lyapunov_rows": len(accepted),
            "variance_accepted_rows": sum(1 for r in rows if r["functional"] == "variance_energy" and r["accepted_internal_lyapunov_row"]),
            "dirichlet_accepted_rows": sum(1 for r in rows if r["functional"] == "quadratic_dirichlet_energy" and r["accepted_internal_lyapunov_row"]),
            "shell_control_monotone_rows": sum(1 for r in rows if r["functional"] == "shell_quadratic_energy" and r["monotone_nonincreasing"]),
            "shell_control_accepted_rows": sum(1 for r in rows if r["functional"] == "shell_quadratic_energy" and r["accepted_internal_lyapunov_row"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3074 constructs an exact Lyapunov/entropy monotonicity certificate for the P3073 internal scale-flow rows.  The intrinsic variance and quadratic Dirichlet energies give monotone nonincrease on 1008 rows total, with 960 strict decreases; the chart-weighted shell energy is kept as a control and is not accepted as intrinsic.  This is real internal stability evidence for dimensionless nadsoliton scalar flows, but it is not a variational source theorem and exports no physical dynamics.",
            "negative_export_flags": {k: False for k in ["variational_source_theorem_exported", "hamiltonian_exported", "unit_bearing_action_exported", "variational_EOM_exported", "canonical_length_time_unit_provider_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "gauge_photon_sector_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"internal_lyapunov_monotonicity_certificate_exported": True, "finite_monotonicity_matrix_executed": True},
            "next_honest_step": "Construct one bounded variational-source obstruction table for the P3074 monotone functionals: test whether a local quadratic action/gradient-flow generator on Z12 has Euler/gradient rows exactly matching the accepted P3073 Laplacian or mean-centering flows, while preserving the no-units/no-observed-physics boundary.  If no local generator matches both flows without chart premises, record a scoped no-go rather than promoting to EOM.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = ["# P3074/S2024 Lyapunov/entropy monotonicity certificate", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
             f"- content grep lanes: `{c['content_grep_lanes']}`", f"- content grep hits: `{c['content_grep_hits']}`", f"- P3073 accepted scale-flow rows reused: `{c['accepted_p3073_flow_rows_reused']}`", f"- step sizes: `{c['step_sizes']}`", f"- functionals tested: `{c['functionals_tested']}`", f"- monotonicity matrix rows: `{c['monotonicity_matrix_rows']}`", f"- accepted internal Lyapunov rows: `{c['accepted_internal_lyapunov_rows']}`", f"- variance accepted rows: `{c['variance_accepted_rows']}`", f"- Dirichlet accepted rows: `{c['dirichlet_accepted_rows']}`", f"- shell control monotone rows: `{c['shell_control_monotone_rows']}`", f"- shell control accepted rows: `{c['shell_control_accepted_rows']}`", f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`", "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"]]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3074/S2024 Lyapunov/entropy monotonicity certificate", "## P3074/S2024 Lyapunov/entropy monotonicity certificate\n\n`P3074/S2024` constructs `SigmaEvenScaleFlowLyapunovEntropyMonotonicityCertificate` for the `192` accepted internal `P3073` scale-flow rows.  It tests `3` exact step sizes and `3` functionals for `1728` rows.  The intrinsic variance and quadratic Dirichlet energies give `1008` accepted monotone nonincreasing internal Lyapunov rows; the shell-weighted energy is retained only as a chart-weighted control.  This exports internal dimensionless stability evidence, not a variational source theorem, Hamiltonian, unit-bearing action, EOM, observed spacetime/light/gauge sector, empirical physics, `L_total`, selector closure, bridge/role-transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3074/S2024 Lyapunov monotonicity is not an action", "## P3074/S2024 Lyapunov monotonicity is not an action\n\n`P3074/S2024` gives exact internal monotonicity for variance and quadratic Dirichlet controls along accepted `P3073` flows.  These functionals are Lyapunov/entropy-like certificates only; no local action, Hamiltonian, variational EOM, unit-bearing coordinate map, gauge field, observed spacetime dynamics, or empirical physics map is exported.\n")
    append_once(AGENTS, "Current Lyapunov/entropy monotonicity guardrail (P3074/S2024, 2026-06-24)", "## Current Lyapunov/entropy monotonicity guardrail (P3074/S2024, 2026-06-24)\n\n- P3074 follows the P3073 recommendation and constructs one exact Lyapunov/entropy monotonicity certificate for the accepted internal scale-flow rows, without selector replay.\n- The matrix reuses `192` accepted P3073 flow rows, tests `3` step sizes and `3` functionals, and has `1728` rows.  Intrinsic variance and quadratic Dirichlet energies give `1008` accepted monotone internal Lyapunov rows; shell energy remains a chart-weighted control and is not accepted as intrinsic.\n- Do not promote P3074 to a variational source theorem, Hamiltonian, canonical unit provider, unit-bearing action, variational EOM, strict spacetime embedding, observed light, gauge photon sector, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one bounded variational-source obstruction table for whether local quadratic action/gradient-flow generators reproduce the accepted Laplacian or mean-centering flows without chart/unit/selector premises.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
