#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import itertools
import json
import math
import re
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.json"
MD = GEN / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.md"

P2653 = GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json"
P2654 = GEN / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.json"
P2655 = GEN / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

ETA = 9.0 / 5.0
BASE_BETA = 1.0
BASE_POINTS: dict[str, tuple[float, float]] = {
    "n0_nadsoliton_origin": (0.0, 0.0),
    "n1_light_axis": (1.0, 0.0),
    "n2_matter_axis": (0.0, 1.0),
    "n3_observer_downstream": (1.0, 1.0),
    "n4_compression_probe": (2.0, 1.0),
}
AUDITED_SCALES = [0.25, 0.5, 1.0, 2.0, 3.0]
TOL = 1e-12

NEGATIVE_EXPORT_FLAGS = [
    "typed_metric_uv_source_theorem_exported",
    "nadsoliton_dynamics_unit_selected",
    "laplacian_action_identity_exported_as_source",
    "dimensionless_operator_identity_exported",
    "scale_orbit_discharged_as_beta_selector",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
    "absolute_spectral_anchor_promoted_to_source",
    "bridge_completion_exported",
    "role_transfer_revalidated",
    "q_w_2191_discharged",
    "role_bearing_ltotal_reenabled",
    "toe_closure_claimed",
]


def sha256_file(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run(
        [
            "rg", "-n", pattern, ".",
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json", "-g", "*.lean", "-g", "*.csv", "-g", "*.tsv",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:70]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "typed_dynamics_operator_content": (
            "typed nadsoliton dynamics|nadsoliton dynamics|Laplacian|operator identity|"
            "dimensionless operator|conserved.*action|action quantum"
        ),
        "scale_source_content": (
            "scale orbit|scale quotient|UV unit|canonical unit|spectral gap|"
            "normalization convention|absolute spectral"
        ),
        "beta_source_content": (
            "beta source|beta=1|unique beta|beta-selecting|target-independent beta|"
            "strict denominator"
        ),
        "nonclosure_guard_content": (
            "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|"
            "source theorem"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for typed dynamics/operator scale-source audit, not packet-name search",
        "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()},
    }


def latest_commit_audit() -> list[dict[str, Any]]:
    proc = subprocess.run(["git", "log", "-n", "12", "--oneline", "--name-only"], cwd=REPO, text=True, stdout=subprocess.PIPE, check=True)
    rows: list[dict[str, Any]] = []
    current: dict[str, Any] | None = None
    for line in proc.stdout.splitlines():
        if not line.strip():
            continue
        if re.match(r"^[0-9a-f]{7,12} ", line):
            if current:
                rows.append(current)
            sha, subject = line.split(" ", 1)
            current = {"sha": sha, "subject": subject, "files": []}
        elif current is not None:
            current["files"].append(line)
    if current:
        rows.append(current)
    return rows


def base_distance(a: str, b: str) -> float:
    ax, ay = BASE_POINTS[a]
    bx, by = BASE_POINTS[b]
    return math.hypot(ax - bx, ay - by)


def laplacian_matrix(scale: float) -> list[list[float]]:
    names = list(BASE_POINTS)
    n = len(names)
    matrix = [[0.0 for _ in names] for _ in names]
    for i, j in itertools.combinations(range(n), 2):
        d = scale * base_distance(names[i], names[j])
        weight = 1.0 / (d * d)
        matrix[i][i] += weight
        matrix[j][j] += weight
        matrix[i][j] -= weight
        matrix[j][i] -= weight
    return matrix


def matmul(a: list[list[float]], b: list[list[float]]) -> list[list[float]]:
    n = len(a)
    return [[sum(a[i][k] * b[k][j] for k in range(n)) for j in range(n)] for i in range(n)]


def trace(matrix: list[list[float]]) -> float:
    return sum(matrix[i][i] for i in range(len(matrix)))


def max_matrix_diff(a: list[list[float]], b: list[list[float]]) -> float:
    return max(abs(a[i][j] - b[i][j]) for i in range(len(a)) for j in range(len(a)))


def operator_invariants(scale: float) -> dict[str, float]:
    l1 = laplacian_matrix(scale)
    l2 = matmul(l1, l1)
    l3 = matmul(l2, l1)
    tr1 = trace(l1)
    tr2 = trace(l2)
    tr3 = trace(l3)
    return {
        "trace_l": tr1,
        "trace_l2": tr2,
        "trace_l3": tr3,
        "dimensionless_m2": tr2 / (tr1 * tr1),
        "dimensionless_m3": tr3 / (tr1 * tr1 * tr1),
    }


def laplacian_scale_covariance_audit() -> dict[str, Any]:
    reference = laplacian_matrix(1.0)
    reference_invariants = operator_invariants(1.0)
    rows = []
    max_covariance_error = 0.0
    max_invariant_error = 0.0
    for scale in AUDITED_SCALES:
        current = laplacian_matrix(scale)
        rescaled = [[(scale * scale) * value for value in row] for row in current]
        covariance_error = max_matrix_diff(rescaled, reference)
        invariants = operator_invariants(scale)
        invariant_error = max(
            abs(invariants[key] - reference_invariants[key])
            for key in ["dimensionless_m2", "dimensionless_m3"]
        )
        max_covariance_error = max(max_covariance_error, covariance_error)
        max_invariant_error = max(max_invariant_error, invariant_error)
        rows.append({
            "scale": scale,
            "trace_l": invariants["trace_l"],
            "trace_l_times_scale_squared": invariants["trace_l"] * scale * scale,
            "dimensionless_m2": invariants["dimensionless_m2"],
            "dimensionless_m3": invariants["dimensionless_m3"],
            "covariance_error_after_a_squared_rescaling": covariance_error,
            "dimensionless_invariant_error_from_scale_one": invariant_error,
        })
    return {
        "statement": "The candidate typed dynamics operator is a weighted nadsoliton graph Laplacian with w_ij=1/d_ij^2.  Under d -> a d it transforms as L -> L/a^2; normalized trace moments are invariant, while absolute trace/spectral scale is not selected.",
        "rows": rows,
        "max_covariance_error": max_covariance_error,
        "max_dimensionless_invariant_error": max_invariant_error,
        "operator_covariance_verified": max_covariance_error < TOL,
        "dimensionless_operator_rank_on_scale_orbit": 1 if max_invariant_error < TOL else len(AUDITED_SCALES),
        "absolute_operator_scale_selected": False,
    }


def beta_covariance_from_operator_scale_audit() -> dict[str, Any]:
    rows = []
    reference_beta = BASE_BETA
    reference_scale = 1.0
    reference_length_unit = 1.0 / math.sqrt(operator_invariants(reference_scale)["trace_l"])
    for scale in AUDITED_SCALES:
        inv = operator_invariants(scale)
        induced_length_unit = 1.0 / math.sqrt(inv["trace_l"])
        beta_covariant = reference_beta / (scale**ETA)
        rows.append({
            "scale": scale,
            "operator_induced_length_unit_if_trace_is_used": induced_length_unit,
            "length_unit_ratio_to_scale_one": induced_length_unit / reference_length_unit,
            "covariant_beta_representative": beta_covariant,
            "beta_equals_one_without_external_trace_anchor": abs(beta_covariant - 1.0) < TOL,
        })
    return {
        "statement": "Using an absolute Laplacian trace or spectral gap would fix a length only after declaring an external numerical action/spectral quantum.  Without that declaration, beta remains a covariant representative on the same scale orbit.",
        "rows": rows,
        "all_nonunit_scales_keep_beta_one_false": all(
            (row["scale"] == 1.0) or (not row["beta_equals_one_without_external_trace_anchor"])
            for row in rows
        ),
        "trace_anchor_is_external_until_dynamics_quantizes_it": True,
        "beta_source_exported_by_operator_now": False,
    }


def source_candidate_matrix(operator: dict[str, Any], beta_audit: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "dimensionless_laplacian_trace_moment_identity",
            "operator_covariant": operator["operator_covariance_verified"],
            "uses_external_absolute_anchor": False,
            "breaks_scale_orbit": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "dimensionless moments are rank-one on the scale orbit, so they cannot choose an absolute UV unit",
        },
        {
            "candidate": "absolute_trace_l_equals_declared_quantum",
            "operator_covariant": True,
            "uses_external_absolute_anchor": True,
            "breaks_scale_orbit": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "would fix scale only by importing a declared spectral/action quantum not derived here",
        },
        {
            "candidate": "spectral_gap_or_heat_time_normalization",
            "operator_covariant": True,
            "uses_external_absolute_anchor": True,
            "breaks_scale_orbit": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "equivalent to choosing a clock/length normalization unless a nadsoliton dynamics quantizes it",
        },
        {
            "candidate": "operator_beta_covariance_identity",
            "operator_covariant": True,
            "uses_external_absolute_anchor": False,
            "breaks_scale_orbit": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "confirms beta transforms covariantly but does not select beta=1 as a sourced representative",
        },
        {
            "candidate": "full_nadsoliton_action_quantization_theorem",
            "operator_covariant": None,
            "uses_external_absolute_anchor": False,
            "breaks_scale_orbit": None,
            "passes_as_uv_unit_source_now": False,
            "verdict": "the next theorem target: derive the absolute quantum from nadsoliton dynamics rather than declaring it",
        },
    ]


def upstream_consistency() -> dict[str, Any]:
    p2653 = load_json(P2653)
    p2654 = load_json(P2654)
    p2655 = load_json(P2655)
    return {
        "p2653_typed_metric_uv_source_not_exported": p2653.get("closure_decision", {}).get("typed_metric_uv_source_theorem_exported_now") is False,
        "p2654_scale_invariant_selector_blocked": p2654.get("closure_decision", {}).get("scale_invariant_selector_exists_now") is False,
        "p2655_metric_skeleton_constructed": p2655.get("closure_decision", {}).get("typed_metric_skeleton_constructed") is True,
        "p2655_uv_unit_not_selected": p2655.get("closure_decision", {}).get("uv_unit_selected_now") is False,
        "p2655_next_step_requests_operator_identity": "operator" in p2655.get("closure_decision", {}).get("next_honest_step", ""),
    }


def closure_decision(matrix: list[dict[str, Any]], operator: dict[str, Any], beta_audit: dict[str, Any]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_uv_unit_source_now"]]
    return {
        "decision": "CANDIDATE_OPERATOR_IDENTITY_AUDITED__ONLY_COVARIANCE_NOT_UV_SOURCE",
        "professorial_verdict": (
            "P2656 follows P2655 by trying the next concrete proof object: a typed nadsoliton Laplacian/action identity.  "
            "The computation succeeds as an equivariant operator scaffold: L scales as a^{-2} and dimensionless trace moments are invariant.  "
            "But that success is still quotient-level; an absolute trace, gap, heat-time, or action quantum would have to be dynamically quantized by the nadsoliton rather than declared as a normalization.  "
            "Therefore P2656 exports no UV unit, no beta source, no bridge completion, and no role-bearing L_total."
        ),
        "professorial_closure_path": [
            "Keep the typed Laplacian/action operator as a candidate dynamics scaffold, not as a completed source theorem.",
            "Do not promote trace(L)=1, spectral gap=1, or heat time=1 into a source; those are absolute anchors unless a nadsoliton quantization theorem derives them.",
            "The next honest proof-grade task is to derive an intrinsic action quantum or conserved integer/phase condition from nadsoliton dynamics and then rerun this operator source audit.",
            "Only after such a theorem fixes a UV unit should beta=1 be retested through P2649; empirical holdout packets remain support only.",
        ],
        "next_honest_step": (
            "Attempt a real nadsoliton action-quantization theorem: define the dynamical evolution law on the typed metric state space, prove a conserved integer/phase/action quantum intrinsic to that law, and then test whether it fixes the Laplacian trace/gap scale without importing an external normalization."
        ),
        "passing_uv_unit_source_candidates": passing,
        "operator_covariance_verified": operator["operator_covariance_verified"],
        "dimensionless_operator_rank_one": operator["dimensionless_operator_rank_on_scale_orbit"] == 1,
        "absolute_operator_scale_selected": operator["absolute_operator_scale_selected"],
        "trace_anchor_external_until_quantized": beta_audit["trace_anchor_is_external_until_dynamics_quantizes_it"],
        "uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    operator = payload["laplacian_scale_covariance_audit"]
    beta_audit = payload["beta_covariance_from_operator_scale_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2656/S1606 typed nadsoliton Laplacian action identity scale-source audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps typed dynamics/operator, scale-source, beta-source, and nonclosure guard content before adding the operator-source audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Laplacian/action operator audit",
        "",
        operator["statement"],
        f"Operator covariance verified? `{operator['operator_covariance_verified']}`.",
        f"Max covariance error: `{operator['max_covariance_error']:.3e}`.",
        f"Dimensionless operator rank on scale orbit: `{operator['dimensionless_operator_rank_on_scale_orbit']}`.",
        f"Max dimensionless invariant error: `{operator['max_dimensionless_invariant_error']:.3e}`.",
        f"Absolute operator scale selected? `{operator['absolute_operator_scale_selected']}`.",
        "",
        "## Beta/source consequence",
        "",
        beta_audit["statement"],
        f"Trace anchor external until dynamics quantizes it? `{beta_audit['trace_anchor_is_external_until_dynamics_quantizes_it']}`.",
        f"Beta source exported by operator now? `{beta_audit['beta_source_exported_by_operator_now']}`.",
        "",
        "## Source candidate matrix",
        "",
        "| candidate | external absolute anchor? | breaks scale orbit? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['uses_external_absolute_anchor']}` | "
            f"`{row['breaks_scale_orbit']}` | `{row['passes_as_uv_unit_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Passing UV-unit source candidates: `{decision['passing_uv_unit_source_candidates']}`.",
        f"Operator covariance verified? `{decision['operator_covariance_verified']}`.",
        f"Dimensionless operator rank one? `{decision['dimensionless_operator_rank_one']}`.",
        f"Absolute operator scale selected? `{decision['absolute_operator_scale_selected']}`.",
        f"UV unit selected now? `{decision['uv_unit_selected_now']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"Role-bearing L_total now? `{decision['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Professorial closure path",
        "",
    ])
    for item in decision["professorial_closure_path"]:
        lines.append(f"- {item}")
    lines.extend([
        "",
        "## Next honest step",
        "",
        decision["next_honest_step"],
        "",
        "## Negative exports",
        "",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    operator = laplacian_scale_covariance_audit()
    beta_audit = beta_covariance_from_operator_scale_audit()
    matrix = source_candidate_matrix(operator, beta_audit)
    decision = closure_decision(matrix, operator, beta_audit)
    payload: dict[str, Any] = {
        "status": "P2656_TYPED_LAPLACIAN_ACTION_IDENTITY_SCALE_SOURCE_AUDIT_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2653_TYPED_METRIC_UV_OBLIGATION": sha256_file(P2653),
            "P2654_SCALE_INVARIANT_SELECTOR_NO_GO": sha256_file(P2654),
            "P2655_TYPED_METRIC_SCALE_QUOTIENT": sha256_file(P2655),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_operator_model": {
            "nodes": BASE_POINTS,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "operator": "complete weighted graph Laplacian with w_ij = 1/d_ij^2",
            "no_sub_nadsoliton_information_layer": True,
        },
        "laplacian_scale_covariance_audit": operator,
        "beta_covariance_from_operator_scale_audit": beta_audit,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2656/S1606 typed Laplacian action identity source guard",
        "## P2656/S1606 typed Laplacian action identity source guard\n\n"
        "`P2656/S1606` tests the next concrete proof object after P2655: a typed nadsoliton weighted-Laplacian/action operator.  The operator scales covariantly as `L -> L/a^2`, and normalized trace moments are rank-one on the scale orbit.  This is a useful dynamics scaffold, but an absolute trace, spectral gap, heat-time, or action quantum would still be an external anchor unless a nadsoliton quantization theorem derives it.  Therefore this exports no typed metric/UV source theorem, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2656/S1606 typed Laplacian action identity Ltotal guard",
        "## P2656/S1606 typed Laplacian action identity Ltotal guard\n\n"
        "`P2656/S1606` does not re-open `L_total`: the candidate Laplacian/action identity proves operator covariance but not an intrinsic UV unit.  A variational damping coefficient still requires a nadsoliton-derived action quantum or conserved dimensionless identity that fixes the absolute scale before beta-source rerun, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
