#!/usr/bin/env python3
"""P3069/S2019: coordinate-pair source rank obstruction.

P3068 identified the next single atom: a strict coordinate-pair source theorem
for t_sigma and x_sigma from nadsoliton data, with units tracked.  P3069 attacks
that atom directly.  It does not try to close the selector and does not promote
P3067/P3068 proxies to observed light.

The constructed object is a finite feature-rank/provider audit over the
T_sigma branches and the Z12 phase index.  It separates raw algebraic ability to
span two coordinate-like columns from the stronger strict requirement: two
nonconventional, unit-bearing, source-provenanced coordinate maps.
"""
from __future__ import annotations

import hashlib, json, subprocess
from fractions import Fraction
from itertools import combinations
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once
from p2812_s1762_two_wl_refined_collision_audit import read_json
from p2855_s1805_z12_rational_phase_lattice_source_candidate_audit import GEN
from p3068_s2018_strict_spacetime_embedding_unit_metric_obstruction import OUT as P3068

OUT = GEN / "p3069_s2019_coordinate_pair_source_rank_obstruction.json"
MD = GEN / "p3069_s2019_coordinate_pair_source_rank_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "coordinate_pair_source_theorem": r"coordinate[- ]pair|time_coordinate|space_coordinate|t_sigma|x_sigma|coordinate.*source",
    "nadsoliton_data_to_coordinates": r"nadsoliton.*coordinate|Z12.*coordinate|phase.*index|kernel distance|chart.*label",
    "unit_tracked_embedding_blocker": r"unit.*coordinate|unit-bearing.*coordinate|speed[- ]of[- ]light|unit-normalized.*metric|spacetime.*embedding",
}

DOMAIN = tuple((sigma, n) for sigma in (-1, 1) for n in range(12))

FEATURES = (
    {"id": "branch_sign_sigma", "kind": "orientation_branch", "unit_bearing": False, "nonconventional_source": False, "chart_dependent": False},
    {"id": "cyclic_phase_index_n", "kind": "phase_chart_label", "unit_bearing": False, "nonconventional_source": False, "chart_dependent": True},
    {"id": "reflection_even_distance", "kind": "kernel_distance_proxy", "unit_bearing": False, "nonconventional_source": False, "chart_dependent": False},
    {"id": "phase_parity", "kind": "discrete_parity", "unit_bearing": False, "nonconventional_source": False, "chart_dependent": True},
    {"id": "sigma_weighted_phase_index", "kind": "orientation_times_chart_label", "unit_bearing": False, "nonconventional_source": False, "chart_dependent": True},
    {"id": "formal_unit_slot", "kind": "empty_unit_placeholder", "unit_bearing": False, "nonconventional_source": False, "chart_dependent": False},
)


def content_grep() -> list[dict[str, Any]]:
    rows = []
    for lane, pattern in CONTENT_PATTERNS.items():
        proc = subprocess.run(["rg", "-n", pattern, "fundamental_action_reconstruction", "AGENTS.md", "-S"], cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
        hits = [line for line in proc.stdout.splitlines() if line.strip()]
        rows.append({"content_lane": lane, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:10]})
    return rows


def feature_value(feature_id: str, sigma: int, n: int) -> Fraction:
    if feature_id == "branch_sign_sigma":
        return Fraction(sigma)
    if feature_id == "cyclic_phase_index_n":
        return Fraction(n)
    if feature_id == "reflection_even_distance":
        return Fraction(min(n, 12 - n))
    if feature_id == "phase_parity":
        return Fraction(1 if n % 2 else 0)
    if feature_id == "sigma_weighted_phase_index":
        return Fraction(sigma * n)
    if feature_id == "formal_unit_slot":
        return Fraction(0)
    raise ValueError(feature_id)


def rational_rank(columns: list[list[Fraction]]) -> int:
    if not columns:
        return 0
    matrix = [list(row) for row in zip(*columns)]
    rows, cols = len(matrix), len(matrix[0])
    rank = 0
    pivot_col = 0
    while rank < rows and pivot_col < cols:
        pivot = next((r for r in range(rank, rows) if matrix[r][pivot_col] != 0), None)
        if pivot is None:
            pivot_col += 1
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        pv = matrix[rank][pivot_col]
        matrix[rank] = [x / pv for x in matrix[rank]]
        for r in range(rows):
            if r != rank and matrix[r][pivot_col] != 0:
                factor = matrix[r][pivot_col]
                matrix[r] = [x - factor * y for x, y in zip(matrix[r], matrix[rank])]
        rank += 1
        pivot_col += 1
    return rank


def subset_rows() -> list[dict[str, Any]]:
    rows = []
    features_by_id = {f["id"]: f for f in FEATURES}
    for size in range(1, len(FEATURES) + 1):
        for ids in combinations([f["id"] for f in FEATURES], size):
            columns = [[feature_value(fid, sigma, n) for sigma, n in DOMAIN] for fid in ids]
            rank = rational_rank(columns)
            unit_bearing = [fid for fid in ids if features_by_id[fid]["unit_bearing"]]
            strict_sources = [fid for fid in ids if features_by_id[fid]["nonconventional_source"]]
            chart_dependent = [fid for fid in ids if features_by_id[fid]["chart_dependent"]]
            accepted = rank >= 2 and len(unit_bearing) >= 2 and len(strict_sources) >= 2 and not chart_dependent
            rows.append({
                "feature_subset": list(ids),
                "subset_size": size,
                "rational_rank": rank,
                "raw_two_column_span": rank >= 2,
                "unit_bearing_features": unit_bearing,
                "strict_source_features": strict_sources,
                "chart_dependent_features": chart_dependent,
                "accepted_coordinate_pair_source": accepted,
            })
    return rows


def feature_rows() -> list[dict[str, Any]]:
    rows = []
    for f in FEATURES:
        col = [feature_value(f["id"], sigma, n) for sigma, n in DOMAIN]
        rows.append({
            **f,
            "nonzero_values": sorted({str(v) for v in col}),
            "single_feature_rank": rational_rank([col]),
            "accepted_as_coordinate_source": False,
            "blocked_by": "not a nonconventional unit-bearing coordinate source theorem",
        })
    return rows


def build_payload() -> dict[str, Any]:
    p3068 = read_json(P3068)
    grep_rows = content_grep()
    features = feature_rows()
    subsets = subset_rows()
    raw_rank2_rows = [r for r in subsets if r["raw_two_column_span"]]
    accepted_rows = [r for r in subsets if r["accepted_coordinate_pair_source"]]
    proof_obligations = [
        {"obligation": "content_first_grep_before_coordinate_pair_source_audit", "satisfied": True, "detail": "searched by coordinate-pair, nadsoliton-data-to-coordinate, and unit-tracked embedding content"},
        {"obligation": "construct_coordinate_pair_source_candidate_space", "satisfied": True, "detail": "constructed six current feature/provider channels over 24 T_sigma x Z12 domain rows"},
        {"obligation": "finite_rational_rank_audit", "satisfied": True, "detail": "enumerated all 63 nonempty feature subsets and computed exact rational ranks"},
        {"obligation": "separate_raw_rank_from_strict_coordinate_source", "satisfied": True, "detail": "raw rank-two chart spans are not accepted without unit-bearing nonconventional source maps"},
        {"obligation": "export_coordinate_pair_source_theorem", "satisfied": False, "detail": "no current feature supplies a unit-bearing nonconventional t_sigma/x_sigma coordinate source theorem"},
    ]
    return {
        "status": "P3069_COORDINATE_PAIR_SOURCE_RANK_OBSTRUCTION_NO_EXPORT",
        "input_hashes": {"P3068": hashlib.sha256(P3068.read_bytes()).hexdigest() if P3068.exists() else None},
        "constructed_theoretical_objects": {
            "content_first_repo_grep": grep_rows,
            "coordinate_pair_source_candidate_space": {
                "object": "CoordinatePairSourceCandidateSpace",
                "domain": "T_sigma branches times Z12 phase rows: 2 x 12 = 24 rows",
                "candidate_features": features,
            },
            "coordinate_pair_rank_obstruction_matrix": subsets,
        },
        "finite_certificate": {
            "content_grep_lanes": len(grep_rows),
            "content_grep_hits": sum(r["hit_count"] for r in grep_rows),
            "domain_rows": len(DOMAIN),
            "candidate_features": len(FEATURES),
            "feature_subset_rows": len(subsets),
            "raw_rank_two_or_more_rows": len(raw_rank2_rows),
            "accepted_coordinate_pair_source_rows": len(accepted_rows),
            "unit_bearing_candidate_features": sum(1 for f in FEATURES if f["unit_bearing"]),
            "nonconventional_source_features": sum(1 for f in FEATURES if f["nonconventional_source"]),
            "max_raw_rank": max(r["rational_rank"] for r in subsets),
            "p3068_current_missing_atoms": p3068.get("finite_certificate", {}).get("current_missing_atoms"),
            "p3068_current_artifact_accepted": p3068.get("finite_certificate", {}).get("current_artifact_accepted"),
            "proof_obligations": len(proof_obligations),
            "satisfied_proof_obligations": sum(1 for o in proof_obligations if o["satisfied"]),
        },
        "proof_obligations": proof_obligations,
        "decision": {
            "bounded_result": "P3069 attacks the P3068 coordinate-pair atom.  The exact 24-row T_sigma x Z12 feature audit finds raw algebraic rank-two spans among chart-dependent features, but 0 accepted coordinate-pair source rows because no current feature is both unit-bearing and nonconventionally sourced as t_sigma or x_sigma.  Thus the P3068 spacetime embedding remains blocked at the coordinate-source layer.",
            "negative_export_flags": {k: False for k in ["coordinate_pair_source_theorem_exported", "time_coordinate_exported", "space_coordinate_exported", "unit_bearing_coordinate_exported", "strict_spacetime_embedding_exported", "observed_light_exported", "variational_EOM_exported", "QW_2191_discharged", "strict_selector_closure_exported", "ltotal_exported", "toe_closure_exported"]},
            "positive_scoped_flags": {"coordinate_candidate_space_constructed": True, "rational_rank_matrix_executed": True, "raw_rank_vs_strict_source_boundary_separated": True},
            "next_honest_step": "Do not use chart-rank spans as spacetime.  The next admissible move is exactly one unit-source atom: construct a canonical length/time unit provider that converts one nonconventional nadsoliton scalar into a unit-bearing coordinate, then rerun P3069.  If no unit provider is exported, pivot to the P3066 sigma-invariant scalar conservation/scale-control theorem and keep light compatibility unclaimed.",
        },
    }


def write_markdown(payload: dict[str, Any]) -> None:
    c = payload["finite_certificate"]
    lines = [
        "# P3069/S2019 coordinate-pair source rank obstruction", "", f"Status: `{payload['status']}`", "",
        "## Finite certificate",
        f"- content grep lanes: `{c['content_grep_lanes']}`",
        f"- content grep hits: `{c['content_grep_hits']}`",
        f"- domain rows: `{c['domain_rows']}`",
        f"- candidate features: `{c['candidate_features']}`",
        f"- feature subset rows: `{c['feature_subset_rows']}`",
        f"- raw rank-two-or-more rows: `{c['raw_rank_two_or_more_rows']}`",
        f"- accepted coordinate-pair source rows: `{c['accepted_coordinate_pair_source_rows']}`",
        f"- unit-bearing candidate features: `{c['unit_bearing_candidate_features']}`",
        f"- nonconventional source features: `{c['nonconventional_source_features']}`",
        f"- max raw rank: `{c['max_raw_rank']}`",
        f"- P3068 current missing atoms: `{c['p3068_current_missing_atoms']}`",
        f"- P3068 current artifact accepted: `{c['p3068_current_artifact_accepted']}`",
        f"- satisfied proof obligations: `{c['satisfied_proof_obligations']}/{c['proof_obligations']}`",
        "", "## Decision", payload["decision"]["bounded_result"], "", "## Recommendation", payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P3069/S2019 coordinate-pair source rank obstruction", "## P3069/S2019 coordinate-pair source rank obstruction\n\n`P3069/S2019` attacks one P3068 atom: a strict coordinate-pair source theorem for `t_sigma` and `x_sigma`.  It constructs a `CoordinatePairSourceCandidateSpace` over `24` `T_sigma x Z12` rows and audits all `63` nonempty subsets of six current features by exact rational rank.  Some chart-dependent subsets have raw rank at least two, but `0` subsets are accepted as coordinate-pair sources because current artifacts provide `0` unit-bearing candidate features and `0` nonconventional coordinate-source features.  Thus no strict spacetime embedding, observed-light closure, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P3069/S2019 coordinate rank is not a coordinate/action source", "## P3069/S2019 coordinate rank is not a coordinate/action source\n\n`P3069/S2019` separates raw finite rank from physical coordinate sourcing.  Chart-dependent `T_sigma x Z12` feature spans do not add a unit-bearing coordinate map, Lorentzian metric, light-sector action density, EOM, Hamiltonian, or empirical constant map.\n")
    append_once(AGENTS, "Current coordinate-pair source rank obstruction guardrail (P3069/S2019, 2026-06-24)", "## Current coordinate-pair source rank obstruction guardrail (P3069/S2019, 2026-06-24)\n\n- P3069 attacks exactly one P3068 atom: a strict coordinate-pair source theorem for `t_sigma` and `x_sigma` from nadsoliton data with units tracked.\n- The finite `T_sigma x Z12` rank audit has `24` domain rows and `63` feature-subset rows; raw chart-dependent rank-two spans exist, but `0` rows are accepted because current features include `0` unit-bearing and `0` nonconventional coordinate-source channels.\n- Do not promote chart-rank spans to spacetime, observed light, Lorentz closure, variational EOM, empirical physics, `QW-2191` discharge, strict selector closure, `L_total`, bridge/role-transfer, or ToE closure.\n- The next honest move is one unit-source atom: construct a canonical length/time unit provider for a nonconventional nadsoliton scalar, then rerun P3069; otherwise pivot to the P3066 sigma-invariant scalar conservation/scale-control theorem.\n")
    return payload


if __name__ == "__main__":
    print(json.dumps(main(), indent=2, sort_keys=True, ensure_ascii=False))
