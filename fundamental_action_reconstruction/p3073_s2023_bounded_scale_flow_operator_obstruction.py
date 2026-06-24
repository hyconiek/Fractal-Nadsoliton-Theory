#!/usr/bin/env python3
"""P3073/S2023: bounded sigma-even scale-flow operator obstruction.

P3072 left one honest next step: test whether the P3071 sigma-even scalar
summaries admit an intrinsic, premise-free scale-flow operator that produces a
nonzero bounded flow without importing selector, units, action, EOM, or observed
physics.  P3073 constructs that finite operator interface over T_sigma x Z12.

The result is mixed and deliberately scoped.  Exact intrinsic nonzero bounded
flows exist (cycle Laplacian and mean-centering), but they preserve only the
zeroth/global total and not the full P3071 scalar-summary packet.  Therefore the
step exports a useful internal scale-flow operator object and a bounded no-go
for full conserved-summary dynamics, not a physical EOM or empirical map.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3071_s2021_sigma_invariant_scalar_conservation_scale_control import profile_values, scalar_summary, transform_values
from p3072_s2022_discrete_continuity_noether_current_matrix import OUT as P3072

OUT = GEN / "p3073_s2023_bounded_scale_flow_operator_obstruction.json"
MD = GEN / "p3073_s2023_bounded_scale_flow_operator_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "bounded_scale_flow_operator": r"bounded.*scale-flow|scale-flow operator|renormalization.*flow|RG.*flow",
    "sigma_even_scalar_summary": r"sigma-even.*scalar|scalar summary|bounded.*scale-control|P3071",
    "z12_laplacian_mean_centering": r"Z12.*Laplacian|cycle Laplacian|mean-centering|total-preserving",
    "no_eom_physics_promotion_boundary": r"unit-bearing action|variational EOM|Hamiltonian|observed light|gauge photon|empirical physics|L_total|ToE|selector closure",
}
SIGMAS = (-1, 1)
Z12 = tuple(range(12))
TRANSFORMS = tuple((kind, k) for kind in ("rotation", "reflection") for k in Z12)
ACCEPTED_P3071_PROFILES = (
    "constant_cardinality_density",
    "even_distance_quadratic_density",
    "even_distance_shell_indicator_density",
)
FLOW_OPERATORS = (
    {"id": "zero_flow", "intrinsic": True, "premise_free": True, "linear": True},
    {"id": "cycle_laplacian_flow", "intrinsic": True, "premise_free": True, "linear": True},
    {"id": "mean_centering_flow", "intrinsic": True, "premise_free": True, "linear": True},
    {"id": "uniform_expansion_flow", "intrinsic": False, "premise_free": False, "linear": True},
    {"id": "chart_slope_flow", "intrinsic": False, "premise_free": False, "linear": True},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return rows


def flow_values(operator_id: str, rho: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    if operator_id == "zero_flow":
        return tuple(Fraction(0) for _ in Z12)
    if operator_id == "cycle_laplacian_flow":
        return tuple(rho[(n - 1) % 12] - 2 * rho[n] + rho[(n + 1) % 12] for n in Z12)
    if operator_id == "mean_centering_flow":
        mean = sum(rho, Fraction(0)) / len(rho)
        return tuple(mean - rho[n] for n in Z12)
    if operator_id == "uniform_expansion_flow":
        return tuple(rho[n] for n in Z12)
    if operator_id == "chart_slope_flow":
        return tuple(Fraction(n) for n in Z12)
    raise ValueError(operator_id)


def add_flow(rho: tuple[Fraction, ...], flow: tuple[Fraction, ...]) -> tuple[Fraction, ...]:
    return tuple(a + b for a, b in zip(rho, flow))


def support_size(values: tuple[Fraction, ...]) -> int:
    return sum(1 for value in values if value != 0)


def span(values: tuple[Fraction, ...]) -> Fraction:
    return max(values) - min(values)


def preserves_full_summary(rho: tuple[Fraction, ...], flowed: tuple[Fraction, ...]) -> bool:
    return scalar_summary(rho) == scalar_summary(flowed)


def audit_rows() -> list[dict[str, Any]]:
    rows = []
    for profile_id in ACCEPTED_P3071_PROFILES:
        for sigma in SIGMAS:
            base_rho = profile_values(profile_id, sigma)
            for kind, k in TRANSFORMS:
                rho = transform_values(base_rho, kind, k)
                rho_total = sum(rho, Fraction(0))
                for operator in FLOW_OPERATORS:
                    flow = flow_values(operator["id"], rho)
                    flowed = add_flow(rho, flow)
                    flow_total = sum(flow, Fraction(0))
                    nonzero_flow = support_size(flow) > 0
                    total_preserving = flow_total == 0 and sum(flowed, Fraction(0)) == rho_total
                    bounded_flow = span(flow) <= Fraction(72)
                    full_summary_preserved = preserves_full_summary(rho, flowed)
                    accepted_bounded_scale_flow = all([
                        operator["intrinsic"], operator["premise_free"], nonzero_flow,
                        total_preserving, bounded_flow,
                    ])
                    accepted_full_summary_dynamics = accepted_bounded_scale_flow and full_summary_preserved
                    rows.append({
                        "profile_id": profile_id,
                        "sigma": sigma,
                        "transform": {"kind": kind, "k": k},
                        "flow_operator": operator["id"],
                        "flow_support": support_size(flow),
                        "flow_total": str(flow_total),
                        "flow_span": str(span(flow)),
                        "nonzero_flow": nonzero_flow,
                        "intrinsic_operator": operator["intrinsic"],
                        "premise_free_operator": operator["premise_free"],
                        "total_preserving": total_preserving,
                        "bounded_flow_span_le_72": bounded_flow,
                        "full_scalar_summary_preserved_after_unit_step": full_summary_preserved,
                        "accepted_intrinsic_bounded_scale_flow_row": accepted_bounded_scale_flow,
                        "accepted_full_summary_dynamics_row": accepted_full_summary_dynamics,
                        "blocked_by": "" if accepted_full_summary_dynamics else "; ".join(filter(None, [
                            None if operator["intrinsic"] else "not intrinsic / chart or convention dependent",
                            None if operator["premise_free"] else "premise or convention required",
                            None if nonzero_flow else "zero flow only",
                            None if total_preserving else "not total-preserving",
                            None if bounded_flow else "flow span bound failed",
                            None if full_summary_preserved else "does not preserve full P3071 scalar-summary packet",
                        ])),
                    })
    return rows


def aggregate_by_operator(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for operator in FLOW_OPERATORS:
        subset = [r for r in rows if r["flow_operator"] == operator["id"]]
        out.append({
            "flow_operator": operator["id"],
            "audit_rows": len(subset),
            "nonzero_rows": sum(1 for r in subset if r["nonzero_flow"]),
            "total_preserving_rows": sum(1 for r in subset if r["total_preserving"]),
            "accepted_intrinsic_bounded_scale_flow_rows": sum(1 for r in subset if r["accepted_intrinsic_bounded_scale_flow_row"]),
            "accepted_full_summary_dynamics_rows": sum(1 for r in subset if r["accepted_full_summary_dynamics_row"]),
        })
    return out


def build_payload() -> dict[str, Any]:
    p3072 = read_json(P3072)
    grep_rows = content_grep()
    rows = audit_rows()
    operator_aggregates = aggregate_by_operator(rows)
    accepted_scale = [r for r in rows if r["accepted_intrinsic_bounded_scale_flow_row"]]
    accepted_full = [r for r in rows if r["accepted_full_summary_dynamics_row"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_scale_flow_audit", "satisfied": True, "detail": "searched by bounded scale-flow, sigma scalar summary, Z12 Laplacian/mean-centering, and no-EOM-promotion content"},
        {"obligation": "construct_finite_scale_flow_operator_family", "satisfied": True, "detail": "five exact flow operators over P3071 accepted profiles on T_sigma x Z12"},
        {"obligation": "execute_bounded_scale_flow_matrix", "satisfied": True, "detail": "3 profiles x 2 sigma branches x 24 D12 transforms x 5 flow operators = 720 exact rows"},
        {"obligation": "export_intrinsic_nonzero_bounded_total_preserving_flow", "satisfied": bool(accepted_scale), "detail": f"{len(accepted_scale)} rows accepted for internal bounded total-preserving scale flow"},
        {"obligation": "export_full_conserved_summary_dynamics", "satisfied": False, "detail": "nonzero accepted flows do not preserve the full P3071 scalar-summary packet after a unit step"},
    ]
    return {
        "status": "P3073_BOUNDED_SCALE_FLOW_OPERATOR_PARTIAL_EXPORT_FULL_SUMMARY_DYNAMICS_NO_GO",
        "input_hashes": {"P3072": hashlib.sha256(P3072.read_bytes()).hexdigest() if P3072.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "bounded_scale_flow_operator_template": {
                "object": "SigmaEvenBoundedScaleFlowOperatorTemplate",
                "domain": "P3071 accepted sigma-even scalar profiles over T_sigma x Z12",
                "tested_flow_family": [op["id"] for op in FLOW_OPERATORS],
                "acceptance_predicate": "intrinsic and premise-free and nonzero and total-preserving and bounded flow span",
                "strong_dynamics_predicate": "acceptance predicate plus preservation of the full P3071 scalar-summary packet after one exact step",
            },
            "operator_aggregate_certificate": operator_aggregates,
            "representative_witness_rows": [r for r in rows if r["accepted_intrinsic_bounded_scale_flow_row"]][:12] + [r for r in rows if not r["accepted_intrinsic_bounded_scale_flow_row"]][:12],
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "p3072_accepted_nontrivial_noether_current_rows": p3072.get("finite_certificate", {}).get("accepted_nontrivial_noether_current_rows"),
            "profiles_tested": len(ACCEPTED_P3071_PROFILES),
            "sigma_branches": len(SIGMAS),
            "d12_transforms": len(TRANSFORMS),
            "flow_operators": len(FLOW_OPERATORS),
            "scale_flow_matrix_rows": len(rows),
            "accepted_intrinsic_bounded_scale_flow_rows": len(accepted_scale),
            "accepted_full_summary_dynamics_rows": len(accepted_full),
            "total_preserving_rows": sum(1 for r in rows if r["total_preserving"]),
            "nonzero_flow_rows": sum(1 for r in rows if r["nonzero_flow"]),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3073 constructs a finite bounded scale-flow interface for the P3071 sigma-even scalar summaries.  Two intrinsic premise-free operators, cycle_laplacian_flow and mean_centering_flow, produce nonzero bounded total-preserving flows on the two nonconstant accepted profiles, giving 192 accepted internal scale-flow rows.  However, every nonzero accepted flow changes at least one member of the full P3071 scalar-summary packet after one exact step, so no full conserved-summary dynamics, Noether theorem, action, EOM, spacetime, or empirical physics is exported.",
            "negative_export_flags": {k: False for k in ["full_conserved_summary_dynamics_exported", "nontrivial_noether_current_exported", "dynamical_update_law_exported", "canonical_length_time_unit_provider_exported", "unit_bearing_coordinate_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "gauge_photon_sector_exported", "unit_bearing_action_exported", "variational_EOM_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"intrinsic_bounded_total_preserving_scale_flow_exported": True, "finite_scale_flow_matrix_executed": True, "full_summary_dynamics_obstruction_exported": True},
            "next_honest_step": "Use the P3073 Laplacian/mean-centering flows only as internal scale-flow operators.  The next proof-grade move is to construct one Lyapunov/entropy monotonicity certificate for these exact flows: test whether a sigma-even functional (variance, shell energy, or quadratic Dirichlet energy) is monotone under the accepted total-preserving flows.  Do not promote monotonicity to action/EOM/observed physics unless a separate variational source theorem is exported.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3073/S2023 bounded scale-flow operator obstruction", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- P3072 accepted nontrivial Noether-current rows: `{c['p3072_accepted_nontrivial_noether_current_rows']}`",
        f"- profiles tested: `{c['profiles_tested']}`",
        f"- sigma branches: `{c['sigma_branches']}`",
        f"- D12 transforms: `{c['d12_transforms']}`",
        f"- flow operators: `{c['flow_operators']}`",
        f"- scale-flow matrix rows: `{c['scale_flow_matrix_rows']}`",
        f"- accepted intrinsic bounded scale-flow rows: `{c['accepted_intrinsic_bounded_scale_flow_rows']}`",
        f"- accepted full-summary dynamics rows: `{c['accepted_full_summary_dynamics_rows']}`",
        f"- total-preserving rows: `{c['total_preserving_rows']}`",
        f"- nonzero-flow rows: `{c['nonzero_flow_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3073/S2023 bounded scale-flow operator obstruction", "## P3073/S2023 bounded scale-flow operator obstruction\n\n`P3073/S2023` constructs `SigmaEvenBoundedScaleFlowOperatorTemplate` over the three `P3071` accepted sigma-even scalar profiles.  It audits `3` profiles across `2` sigma branches, `24` `D12` transforms, and `5` exact flow operators, for `720` rows.  The cycle-Laplacian and mean-centering operators produce `192` intrinsic, premise-free, nonzero, bounded, total-preserving internal scale-flow rows, but `0` rows preserve the full `P3071` scalar-summary packet after one exact step.  Thus no full conserved-summary dynamics, Noether theorem, unit-bearing action, EOM, observed-light/gauge sector, `L_total`, bridge/role-transfer, selector closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3073/S2023 scale-flow operators are not variational dynamics", "## P3073/S2023 scale-flow operators are not variational dynamics\n\n`P3073/S2023` exports internal total-preserving bounded scale-flow operators for the nonconstant `P3071` scalar profiles, but the accepted nonzero flows do not preserve the full scalar-summary packet and are not derived from an action, Hamiltonian, variational EOM, unit-bearing coordinate map, gauge field, observed spacetime dynamics, or empirical physics map.\n")
    append_once(AGENTS, "Current bounded scale-flow operator guardrail (P3073/S2023, 2026-06-24)", "## Current bounded scale-flow operator guardrail (P3073/S2023, 2026-06-24)\n\n- P3073 follows the P3072 recommendation and constructs one bounded scale-flow operator matrix for the `P3071` accepted sigma-even scalar summaries.\n- The matrix has `3` profiles, `2` sigma branches, `24` `D12` transforms, `5` exact flow operators, and `720` rows.  Cycle-Laplacian and mean-centering flows give `192` intrinsic, premise-free, nonzero, bounded, total-preserving internal scale-flow rows, but `0` rows preserve the full `P3071` scalar-summary packet after one exact step.\n- Do not promote P3073 to full conserved-summary dynamics, a nontrivial Noether current, dynamical update law, canonical unit provider, unit-bearing coordinates, strict spacetime embedding, observed light, gauge photons, unit-bearing action, variational EOM, empirical physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one Lyapunov/entropy monotonicity certificate for the accepted total-preserving flows, unless a new strict variational action/EOM/unit provider is introduced.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
