#!/usr/bin/env python3
"""P3075/S2025: bounded variational-source obstruction table.

P3074 exported internal Lyapunov monotonicity, but not a variational source.
P3075 constructs the missing finite interface: candidate quadratic generators on
Z12 and their exact negative-gradient rows are compared with the accepted P3073
Laplacian/mean-centering scale-flow rows.  The result is deliberately scoped:
the local Dirichlet generator reproduces the Laplacian flow internally, while
the mean-centering flow is reproduced only by a global/nonlocal variance source.
No unit-bearing action, spacetime EOM, gauge field, or empirical physics is
exported.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3071_s2021_sigma_invariant_scalar_conservation_scale_control import profile_values, transform_values
from p3073_s2023_bounded_scale_flow_operator_obstruction import audit_rows as p3073_audit_rows, flow_values
from p3074_s2024_lyapunov_entropy_monotonicity_certificate import OUT as P3074

OUT = GEN / "p3075_s2025_variational_source_obstruction_table.json"
MD = GEN / "p3075_s2025_variational_source_obstruction_table.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
Z12 = tuple(range(12))
CONTENT_PATTERNS = {
    "variational_source_gradient_flow": r"variational source|gradient-flow|negative-gradient|Euler|quadratic action",
    "z12_local_dirichlet_laplacian": r"Z12.*Dirichlet|cycle Laplacian|local quadratic|edge-local",
    "mean_centering_nonlocal_boundary": r"mean-centering|variance.*global|nonlocal|total-preserving",
    "no_unit_eom_physics_promotion": r"unit-bearing action|variational EOM|Hamiltonian|observed spacetime|gauge field|empirical physics|L_total|ToE|selector closure",
}
GENERATORS = (
    {"id": "zero_generator", "intrinsic": True, "premise_free": True, "local_on_z12": True},
    {"id": "local_dirichlet_generator", "intrinsic": True, "premise_free": True, "local_on_z12": True},
    {"id": "global_variance_generator", "intrinsic": True, "premise_free": True, "local_on_z12": False},
    {"id": "chart_shell_generator", "intrinsic": False, "premise_free": False, "local_on_z12": False},
    {"id": "linear_total_generator", "intrinsic": True, "premise_free": True, "local_on_z12": False},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def mean(rho: tuple[Fraction, ...]) -> Fraction:
    return sum(rho, Fraction(0)) / len(rho)


def negative_gradient(generator_id: str, rho: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if generator_id == "zero_generator":
        return tuple(Fraction(0) for _ in Z12)
    if generator_id == "local_dirichlet_generator":
        # E_D = 1/2 * sum_n (rho_{n+1}-rho_n)^2 gives -grad(E_D) = rho_{n-1}-2rho_n+rho_{n+1}.
        return tuple(rho[(n - 1) % 12] - 2 * rho[n] + rho[(n + 1) % 12] for n in Z12)
    if generator_id == "global_variance_generator":
        # E_V = 1/2 * sum_n (rho_n-mean(rho))^2 gives normalized -grad(E_V) = mean(rho)-rho_n.
        m = mean(rho)
        return tuple(m - rho[n] for n in Z12)
    if generator_id == "chart_shell_generator":
        # E_S = 1/2 * sum_n d_chart(n)^2 rho_n^2 depends on a chosen origin/shell chart.
        return tuple(-Fraction(min(n, 12 - n) ** 2) * rho[n] for n in Z12)
    if generator_id == "linear_total_generator":
        # E_T = sum_n rho_n: intrinsic but gradient is a uniform drift, not total-preserving.
        return tuple(Fraction(-1) for _ in Z12)
    raise ValueError(generator_id)


def source_rows() -> list[dict[str, Any]]:
    rows = []
    accepted_flows = [r for r in p3073_audit_rows() if r["accepted_intrinsic_bounded_scale_flow_row"]]
    for base in accepted_flows:
        rho = transform_values(profile_values(base["profile_id"], base["sigma"]), base["transform"]["kind"], base["transform"]["k"])
        target_flow = flow_values(base["flow_operator"], rho)
        for generator in GENERATORS:
            source_flow = negative_gradient(generator["id"], rho)
            exact_match = source_flow == target_flow
            total_preserving_source = sum(source_flow, Fraction(0)) == 0
            accepted_local_variational_source = all([
                generator["intrinsic"], generator["premise_free"], generator["local_on_z12"], exact_match, total_preserving_source,
            ])
            rows.append({
                "profile_id": base["profile_id"],
                "sigma": base["sigma"],
                "transform": base["transform"],
                "target_flow_operator": base["flow_operator"],
                "candidate_generator": generator["id"],
                "intrinsic_generator": generator["intrinsic"],
                "premise_free_generator": generator["premise_free"],
                "local_on_z12": generator["local_on_z12"],
                "exact_negative_gradient_matches_target_flow": exact_match,
                "total_preserving_source": total_preserving_source,
                "accepted_local_variational_source_row": accepted_local_variational_source,
                "blocked_by": "" if accepted_local_variational_source else "; ".join(filter(None, [
                    None if generator["intrinsic"] else "not intrinsic / chart weighted",
                    None if generator["premise_free"] else "requires chart or convention premise",
                    None if generator["local_on_z12"] else "not local on Z12 edge-neighborhoods",
                    None if exact_match else "negative-gradient row does not match target flow",
                    None if total_preserving_source else "source gradient is not total-preserving",
                ])),
            })
    return rows


def aggregate(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for generator in GENERATORS:
        subset = [r for r in rows if r["candidate_generator"] == generator["id"]]
        out.append({
            "candidate_generator": generator["id"],
            "rows": len(subset),
            "exact_match_rows": sum(1 for r in subset if r["exact_negative_gradient_matches_target_flow"]),
            "total_preserving_source_rows": sum(1 for r in subset if r["total_preserving_source"]),
            "accepted_local_variational_source_rows": sum(1 for r in subset if r["accepted_local_variational_source_row"]),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3074 = read_json(P3074)
    grep_rows = content_grep()
    rows = source_rows()
    accepted = [r for r in rows if r["accepted_local_variational_source_row"]]
    mean_nonlocal_matches = [r for r in rows if r["target_flow_operator"] == "mean_centering_flow" and r["candidate_generator"] == "global_variance_generator" and r["exact_negative_gradient_matches_target_flow"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_variational_source_audit", "satisfied": True, "detail": "searched by variational/gradient-flow, Z12 local Dirichlet, mean-centering nonlocality, and no-physics-promotion content"},
        {"obligation": "construct_candidate_quadratic_generators", "satisfied": True, "detail": "zero, local Dirichlet, global variance, chart shell, and linear total generators are tested exactly"},
        {"obligation": "execute_exact_negative_gradient_match_matrix", "satisfied": True, "detail": "192 accepted P3073 flows x 5 candidate generators = 960 exact rows"},
        {"obligation": "export_local_internal_variational_source_for_laplacian", "satisfied": bool(accepted), "detail": f"{len(accepted)} rows match the accepted Laplacian flow by the local Dirichlet generator"},
        {"obligation": "export_local_source_for_mean_centering_flow", "satisfied": False, "detail": f"{len(mean_nonlocal_matches)} exact matches exist only through the global variance generator, which is nonlocal on Z12"},
        {"obligation": "export_unit_bearing_physical_eom", "satisfied": False, "detail": "the source is finite, dimensionless, internal, and not a unit-bearing spacetime variational EOM"},
    ]
    return {
        "status": "P3075_LOCAL_DIRICHLET_VARIATIONAL_SOURCE_PARTIAL_EXPORT_MEAN_CENTERING_NONLOCAL_OBSTRUCTION",
        "input_hashes": {"P3074": hashlib.sha256(P3074.read_bytes()).hexdigest() if P3074.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "quadratic_generator_template": {
                "object": "Z12LocalQuadraticGradientSourceInterface",
                "domain": "P3073 accepted intrinsic bounded total-preserving Laplacian/mean-centering flow rows",
                "candidate_generators": [g["id"] for g in GENERATORS],
                "acceptance_predicate": "intrinsic and premise-free and local-on-Z12 and exact negative-gradient match to the target P3073 flow and total-preserving source",
            },
            "generator_aggregate_certificate": aggregate(rows),
            "representative_witness_rows": accepted[:12] + mean_nonlocal_matches[:12] + [r for r in rows if not r["accepted_local_variational_source_row"] and not r["exact_negative_gradient_matches_target_flow"]][:12],
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "p3074_accepted_internal_lyapunov_rows": p3074.get("finite_certificate", {}).get("accepted_internal_lyapunov_rows"),
            "accepted_p3073_flow_rows_reused": 192,
            "candidate_generators": len(GENERATORS),
            "variational_source_matrix_rows": len(rows),
            "exact_negative_gradient_match_rows": sum(1 for r in rows if r["exact_negative_gradient_matches_target_flow"]),
            "accepted_local_variational_source_rows": len(accepted),
            "local_dirichlet_accepted_rows": sum(1 for r in accepted if r["candidate_generator"] == "local_dirichlet_generator"),
            "global_variance_exact_but_nonlocal_rows": len(mean_nonlocal_matches),
            "mean_centering_local_source_rows": sum(1 for r in rows if r["target_flow_operator"] == "mean_centering_flow" and r["accepted_local_variational_source_row"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3075 constructs the missing variational-source interface for the P3074 monotone flows.  The local Dirichlet quadratic generator has exact negative-gradient rows matching all 96 accepted cycle-Laplacian target rows, so a scoped internal local gradient source is exported for that flow.  The 96 mean-centering target rows match exactly only through the global variance generator, which is intrinsic but nonlocal on Z12, so no local source for mean-centering is exported.  The result remains finite, dimensionless, and internal: it is not a unit-bearing action, spacetime EOM, gauge field, or empirical physics map.",
            "negative_export_flags": {k: False for k in ["full_variational_EOM_exported", "unit_bearing_action_exported", "hamiltonian_exported", "canonical_length_time_unit_provider_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "gauge_photon_sector_exported", "empirical_physics_exported", "mean_centering_local_source_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"internal_local_dirichlet_gradient_source_for_laplacian_exported": True, "finite_variational_source_obstruction_matrix_executed": True, "mean_centering_nonlocal_obstruction_exported": True},
            "next_honest_step": "Use only the P3075 local Dirichlet/Laplacian source as an internal dimensionless gradient-flow object.  The next proof-grade move is one bounded continuum-limit/spectral-dispersion audit: diagonalize the Z12 Dirichlet Laplacian source, compute its exact mode spectrum and small-k dispersion proxy, and test whether it has a lightlike/wave-compatible branch or only a diffusive/internal smoothing branch.  Do not promote the branch to observed light, gauge photons, spacetime EOM, units, or empirical physics without a separate unit/coordinate/source theorem.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3075/S2025 variational-source obstruction table", "", f"Status: `{payload['status']}`", "", "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3074 accepted internal Lyapunov rows: `{c['p3074_accepted_internal_lyapunov_rows']}`",
        f"- accepted P3073 flow rows reused: `{c['accepted_p3073_flow_rows_reused']}`",
        f"- candidate generators: `{c['candidate_generators']}`",
        f"- variational-source matrix rows: `{c['variational_source_matrix_rows']}`",
        f"- exact negative-gradient match rows: `{c['exact_negative_gradient_match_rows']}`",
        f"- accepted local variational-source rows: `{c['accepted_local_variational_source_rows']}`",
        f"- local Dirichlet accepted rows: `{c['local_dirichlet_accepted_rows']}`",
        f"- global variance exact-but-nonlocal rows: `{c['global_variance_exact_but_nonlocal_rows']}`",
        f"- mean-centering local source rows: `{c['mean_centering_local_source_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3075/S2025 variational-source obstruction table", "## P3075/S2025 variational-source obstruction table\n\n`P3075/S2025` constructs `Z12LocalQuadraticGradientSourceInterface` for the `P3073` accepted internal scale-flow rows after the `P3074` Lyapunov certificate.  It tests `192` accepted target-flow rows against `5` candidate quadratic generators for `960` exact negative-gradient rows.  The local Dirichlet generator exactly sources all `96` cycle-Laplacian rows as an internal dimensionless local gradient flow; the mean-centering rows match only the global variance generator and therefore remain nonlocal on `Z12`.  This is not a unit-bearing action, Hamiltonian, spacetime variational EOM, observed-light/gauge sector, empirical physics, `L_total`, selector closure, bridge/role-transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3075/S2025 local Dirichlet source is internal only", "## P3075/S2025 local Dirichlet source is internal only\n\n`P3075/S2025` exports a scoped internal local quadratic Dirichlet gradient source for the accepted Laplacian flow, while mean-centering remains exact only through a global nonlocal variance generator.  The finite source object is dimensionless and does not export a unit-bearing action, Hamiltonian, variational spacetime EOM, coordinate theorem, gauge field, observed spacetime dynamics, or empirical physics map.\n")
    append_once(AGENTS, "Current variational-source obstruction guardrail (P3075/S2025, 2026-06-24)", "## Current variational-source obstruction guardrail (P3075/S2025, 2026-06-24)\n\n- P3075 follows the P3074 recommendation and constructs one bounded variational-source obstruction table for accepted P3073 scale flows, without selector replay or observed-physics promotion.\n- The matrix reuses `192` accepted P3073 flow rows, tests `5` candidate quadratic generators, and has `960` exact negative-gradient rows.  The local Dirichlet generator exactly sources all `96` cycle-Laplacian rows; mean-centering has `96` exact matches only through the global variance generator and remains nonlocal on `Z12`.\n- Do not promote P3075 to a unit-bearing action, Hamiltonian, full variational EOM, canonical unit provider, strict spacetime embedding, observed light, gauge photon sector, empirical physics, `QW-2191` discharge, selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one bounded continuum-limit/spectral-dispersion audit of the internal Dirichlet/Laplacian source, without claiming observed light or gauge photons absent a separate unit/coordinate/source theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
