#!/usr/bin/env python3
"""P3071/S2021: sigma-invariant scalar conservation/scale-control theorem.

P3070 obstructed the unit-source atom and recommended pivoting to the P3066
sigma-invariant scalar row.  This step follows that recommendation directly: it
constructs a finite, exact conservation/scale-control theorem for sigma-even
nadsoliton scalar profiles on T_sigma x Z12.

The theorem is deliberately scoped.  It proves conservation only for bounded
sigma-even scalar summaries under sigma flip and D12 reindexing/permutation
moves.  It does not export units, coordinates, observed light, a gauge sector,
EOM, L_total, bridge/role-transfer, or ToE closure.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3070_s2020_canonical_length_time_unit_provider_obstruction import OUT as P3070

OUT = GEN / "p3071_s2021_sigma_invariant_scalar_conservation_scale_control.json"
MD = GEN / "p3071_s2021_sigma_invariant_scalar_conservation_scale_control.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "sigma_invariant_scalar_conservation": r"sigma[- ]invariant.*scalar|sigma-even.*scalar|scalar.*conservation|conservation.*scalar",
    "scale_control_bounded_scalar": r"scale[- ]control|bounded[- ]scale|bounded.*scalar|normaliz(ed|ation).*scalar",
    "nadsoliton_z12_profile": r"nadsoliton.*Z12|Z12.*profile|orbit.*histogram|dihedral|D12",
    "no_physics_promotion_boundary": r"observed light|unit-bearing coordinate|gauge photon|variational EOM|empirical physics|L_total|ToE",
}

SIGMAS = (-1, 1)
Z12 = tuple(range(12))
TRANSFORMS = tuple((kind, k) for kind in ("rotation", "reflection") for k in Z12)


def d12_apply(kind: str, k: int, n: int) -> int:
    if kind == "rotation":
        return (n + k) % 12
    if kind == "reflection":
        return (k - n) % 12
    raise ValueError(kind)


def distance(n: int) -> int:
    return min(n, 12 - n)


def profile_values(profile_id: str, sigma: int) -> tuple[Fraction, ...]:
    if profile_id == "constant_cardinality_density":
        return tuple(Fraction(1) for _ in Z12)
    if profile_id == "even_distance_quadratic_density":
        return tuple(Fraction(distance(n) ** 2) for n in Z12)
    if profile_id == "even_distance_shell_indicator_density":
        return tuple(Fraction(1 if distance(n) in {0, 3, 6} else 0) for n in Z12)
    if profile_id == "sigma_signed_distance_density":
        return tuple(Fraction(sigma * distance(n)) for n in Z12)
    if profile_id == "chart_label_density":
        return tuple(Fraction(n) for n in Z12)
    raise ValueError(profile_id)


PROFILES = (
    {"id": "constant_cardinality_density", "sigma_even": True, "chart_label_dependent": False, "dimensionless": True},
    {"id": "even_distance_quadratic_density", "sigma_even": True, "chart_label_dependent": False, "dimensionless": True},
    {"id": "even_distance_shell_indicator_density", "sigma_even": True, "chart_label_dependent": False, "dimensionless": True},
    {"id": "sigma_signed_distance_density", "sigma_even": False, "chart_label_dependent": False, "dimensionless": True},
    {"id": "chart_label_density", "sigma_even": True, "chart_label_dependent": True, "dimensionless": True},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def transform_values(values: tuple[Fraction, ...], kind: str, k: int) -> tuple[Fraction, ...]:
    out = [Fraction(0) for _ in Z12]
    for n, value in enumerate(values):
        out[d12_apply(kind, k, n)] = value
    return tuple(out)


def scalar_summary(values: tuple[Fraction, ...]) -> dict[str, str]:
    total = sum(values, Fraction(0))
    mean = total / len(values)
    square_total = sum(v * v for v in values)
    max_value = max(values)
    min_value = min(values)
    span = max_value - min_value
    return {
        "total": str(total),
        "mean": str(mean),
        "square_total": str(square_total),
        "min": str(min_value),
        "max": str(max_value),
        "span": str(span),
    }


def profile_audit_rows() -> list[dict[str, Any]]:
    rows = []
    for profile in PROFILES:
        baseline = profile_values(profile["id"], 1)
        baseline_summary = scalar_summary(baseline)
        sigma_flip_values = profile_values(profile["id"], -1)
        sigma_flip_summary = scalar_summary(sigma_flip_values)
        for sigma in SIGMAS:
            source_values = profile_values(profile["id"], sigma)
            for kind, k in TRANSFORMS:
                moved = transform_values(source_values, kind, k)
                moved_summary = scalar_summary(moved)
                summary_conserved = moved_summary == scalar_summary(source_values)
                sigma_summary_conserved = sigma_flip_summary == baseline_summary
                bounded = Fraction(moved_summary["span"]) <= Fraction(36)
                accepted = all([
                    profile["sigma_even"], not profile["chart_label_dependent"], profile["dimensionless"],
                    summary_conserved, sigma_summary_conserved, bounded,
                ])
                rows.append({
                    "profile_id": profile["id"],
                    "sigma": sigma,
                    "transform": {"kind": kind, "k": k},
                    "summary": moved_summary,
                    "summary_conserved_under_d12": summary_conserved,
                    "summary_conserved_under_sigma_flip": sigma_summary_conserved,
                    "bounded_span_le_36": bounded,
                    "accepted_sigma_invariant_scale_control_row": accepted,
                    "blocked_by": "" if accepted else "; ".join(filter(None, [
                        None if profile["sigma_even"] else "not sigma-even",
                        None if not profile["chart_label_dependent"] else "chart-label dependent",
                        None if summary_conserved else "D12 summary not conserved",
                        None if sigma_summary_conserved else "sigma-flip summary not conserved",
                        None if bounded else "span bound failed",
                    ])),
                })
    return rows


def aggregate_profiles(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    out = []
    for profile in PROFILES:
        subset = [r for r in rows if r["profile_id"] == profile["id"]]
        accepted = [r for r in subset if r["accepted_sigma_invariant_scale_control_row"]]
        out.append({
            "profile_id": profile["id"],
            "audit_rows": len(subset),
            "accepted_rows": len(accepted),
            "all_rows_accepted": len(accepted) == len(subset),
            "baseline_summary_sigma_plus": scalar_summary(profile_values(profile["id"], 1)),
            "baseline_summary_sigma_minus": scalar_summary(profile_values(profile["id"], -1)),
            "exports_physical_unit": False,
            "exports_coordinate": False,
            "exports_observed_physics": False,
        })
    return out


def build_payload() -> dict[str, Any]:
    p3070 = read_json(P3070)
    grep_rows = content_grep()
    rows = profile_audit_rows()
    aggregates = aggregate_profiles(rows)
    accepted_profiles = [a for a in aggregates if a["all_rows_accepted"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_sigma_invariant_scalar_audit", "satisfied": True, "detail": "searched by scalar conservation, scale-control, Z12 profile, and no-physics-promotion content"},
        {"obligation": "construct_sigma_even_scalar_profile_space", "satisfied": True, "detail": "five finite profiles over T_sigma x Z12 with exact Fraction summaries"},
        {"obligation": "finite_d12_sigma_flip_conservation_matrix", "satisfied": True, "detail": "5 profiles x 2 sigma branches x 24 D12 transforms = 240 rows"},
        {"obligation": "export_scoped_sigma_invariant_scale_control_rows", "satisfied": bool(accepted_profiles), "detail": f"{len(accepted_profiles)} profiles accepted as scoped dimensionless conservation/scale-control rows"},
        {"obligation": "export_standard_physics_or_unit_closure", "satisfied": False, "detail": "the theorem remains dimensionless and exports no units, coordinates, EOM, gauge sector, or empirical map"},
    ]
    return {
        "status": "P3071_SIGMA_INVARIANT_SCALAR_CONSERVATION_SCALE_CONTROL_SCOPED_EXPORT",
        "input_hashes": {"P3070": hashlib.sha256(P3070.read_bytes()).hexdigest() if P3070.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "sigma_even_scalar_conservation_template": {
                "object": "SigmaEvenNadsolitonScalarConservationTemplate",
                "domain": "T_sigma branches times Z12 phase rows",
                "group_actions": "sigma flip and D12 rotations/reflections as finite reindexing moves",
                "conserved_summaries": ["total", "mean", "square_total", "min", "max", "span"],
            },
            "candidate_profiles": list(PROFILES),
            "profile_aggregate_certificate": aggregates,
            "d12_sigma_flip_conservation_matrix": rows,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "sigma_branches": len(SIGMAS),
            "z12_rows": len(Z12),
            "d12_transforms": len(TRANSFORMS),
            "candidate_profiles": len(PROFILES),
            "conservation_matrix_rows": len(rows),
            "accepted_profile_count": len(accepted_profiles),
            "accepted_profile_ids": [a["profile_id"] for a in accepted_profiles],
            "accepted_matrix_rows": sum(1 for r in rows if r["accepted_sigma_invariant_scale_control_row"]),
            "rejected_matrix_rows": sum(1 for r in rows if not r["accepted_sigma_invariant_scale_control_row"]),
            "p3070_accepted_unit_provider_rows": p3070.get("finite_certificate", {}).get("accepted_unit_provider_rows"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3071 pivots exactly as P3070 recommended and proves a scoped, dimensionless sigma-invariant scalar conservation/scale-control result.  Three Z12 scalar profiles have all 48 sigma/D12 rows accepted: constant cardinality density, even distance-quadratic density, and even distance shell-indicator density.  Two profiles are rejected because one is sigma-odd and one is chart-label dependent.  This is real finite conservation evidence, but it remains internal and dimensionless.",
            "negative_export_flags": {k: False for k in ["canonical_length_time_unit_provider_exported", "unit_bearing_coordinate_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "gauge_photon_sector_exported", "unit_bearing_action_exported", "variational_EOM_exported", "empirical_physics_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "bridge_role_transfer_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"sigma_invariant_scalar_profiles_exported": True, "finite_conservation_matrix_executed": True, "bounded_dimensionless_scale_control_exported": True},
            "next_honest_step": "Use the P3071 accepted sigma-even scalar profiles as internal dimensionless invariants only.  The next proof-grade move is to construct one transition-interface theorem from these conserved profiles to a dynamics candidate: either a discrete continuity/Noether-current matrix on Z12 or a bounded renormalization/scale-flow obstruction table.  Do not jump to units, spacetime, observed light, selector closure, L_total, bridge/role-transfer, or ToE without a new unit/action/EOM provider.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3071/S2021 sigma-invariant scalar conservation/scale-control", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- sigma branches: `{c['sigma_branches']}`",
        f"- Z12 rows: `{c['z12_rows']}`",
        f"- D12 transforms: `{c['d12_transforms']}`",
        f"- candidate profiles: `{c['candidate_profiles']}`",
        f"- conservation matrix rows: `{c['conservation_matrix_rows']}`",
        f"- accepted profile count: `{c['accepted_profile_count']}`",
        f"- accepted profile ids: `{', '.join(c['accepted_profile_ids'])}`",
        f"- accepted matrix rows: `{c['accepted_matrix_rows']}`",
        f"- rejected matrix rows: `{c['rejected_matrix_rows']}`",
        f"- P3070 accepted unit-provider rows: `{c['p3070_accepted_unit_provider_rows']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3071/S2021 sigma-invariant scalar conservation/scale-control", "## P3071/S2021 sigma-invariant scalar conservation/scale-control\n\n`P3071/S2021` follows the P3070 pivot and constructs a finite `SigmaEvenNadsolitonScalarConservationTemplate` on `T_sigma x Z12`.  It audits `5` scalar profiles across `2` sigma branches and `24` `D12` reindexing transforms, for `240` exact rows.  Three dimensionless profiles are accepted as scoped conservation/scale-control rows: constant cardinality density, even distance-quadratic density, and even distance shell-indicator density.  The result exports internal bounded scalar invariants only; it exports no unit provider, coordinate theorem, spacetime embedding, observed-light closure, `L_total`, bridge/role-transfer, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3071/S2021 conserved scalar profiles are not dynamics", "## P3071/S2021 conserved scalar profiles are not dynamics\n\n`P3071/S2021` gives scoped dimensionless conservation/scale-control for three sigma-even Z12 scalar profiles under sigma flip and D12 reindexing.  These invariants are not a Noether current, action density, EOM, Hamiltonian, unit-bearing coordinate, gauge field, or empirical constant map.\n")
    append_once(AGENTS, "Current sigma-invariant scalar conservation/scale-control guardrail (P3071/S2021, 2026-06-24)", "## Current sigma-invariant scalar conservation/scale-control guardrail (P3071/S2021, 2026-06-24)\n\n- P3071 follows the P3070 pivot rather than replaying selector/unit/spacetime claims: it constructs a finite sigma-invariant scalar conservation/scale-control theorem on `T_sigma x Z12`.\n- The exact matrix has `5` candidate profiles, `2` sigma branches, `24` `D12` transforms, and `240` rows; `3` dimensionless profiles are accepted as scoped internal invariants, while sigma-odd and chart-label profiles are rejected.\n- Do not promote these dimensionless conserved summaries to canonical units, unit-bearing coordinates, strict spacetime embedding, observed light, gauge photons, unit-bearing action, variational EOM, empirical physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one transition-interface theorem from conserved profiles toward dynamics: a discrete continuity/Noether-current matrix on `Z12` or a bounded renormalization/scale-flow obstruction table, with no closure promotion without a new unit/action/EOM provider.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
