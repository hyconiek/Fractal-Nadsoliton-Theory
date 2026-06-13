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
OUT = GEN / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.json"
MD = GEN / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.md"
P2658 = GEN / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.json"
P2659 = GEN / "p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NODES = {
    "n0_nadsoliton_origin": (0.0, 0.0),
    "n1_light_axis": (1.0, 0.0),
    "n2_matter_axis": (0.0, 1.0),
    "n3_observer_downstream": (1.0, 1.0),
    "n4_compression_probe": (2.0, 1.0),
}
EDGES = [
    ("n0_nadsoliton_origin", "n1_light_axis"),
    ("n1_light_axis", "n3_observer_downstream"),
    ("n3_observer_downstream", "n2_matter_axis"),
    ("n2_matter_axis", "n0_nadsoliton_origin"),
    ("n1_light_axis", "n4_compression_probe"),
    ("n4_compression_probe", "n3_observer_downstream"),
]
TRIANGLES = [("n1_light_axis", "n3_observer_downstream", "n4_compression_probe")]
SCALES = [0.25, 0.5, 1.0, 2.0, 3.0]
EXPONENTS = [1.0, 1.8, 2.0]
QUANTA = [1, 2, 3]
TOL = 1e-11
NEGATIVE_EXPORT_FLAGS = [
    "boundary_cocycle_anomaly_coefficient_sourced",
    "dimensionful_action_unit_exported_from_topology",
    "intrinsic_clock_or_flow_time_exported",
    "uv_unit_selected_by_topological_anomaly",
    "typed_metric_uv_source_theorem_exported",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
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
            "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "*.json",
            "-g", "!fundamental_action_reconstruction/generated/**", "-g", "!.git/**",
        ],
        cwd=REPO,
        text=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        check=False,
    )
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:60]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "boundary_cocycle_content": "boundary|cocycle|cycle|Betti|Euler|topological",
        "anomaly_coefficient_content": "anomaly coefficient|lambda_anomaly|anomalous.*term|clock-source",
        "dimension_unit_content": "dimensionful|dimensionless|action unit|UV unit|canonical unit",
        "nonclosure_guard_content": "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|beta source",
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for boundary/cocycle anomaly coefficient sourcing, not packet-name search",
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


def distance(a: str, b: str, scale: float = 1.0) -> float:
    ax, ay = NODES[a]
    bx, by = NODES[b]
    return scale * math.hypot(ax - bx, ay - by)


def graph_invariants() -> dict[str, Any]:
    vertex_count = len(NODES)
    edge_count = len(EDGES)
    face_count = len(TRIANGLES)
    connected_components = 1
    beta_1_graph = edge_count - vertex_count + connected_components
    beta_1_after_filling_triangles = beta_1_graph - face_count
    euler_characteristic = vertex_count - edge_count + face_count
    cycle_rank_witness = {
        "square_cycle_boundary_mod2": 0,
        "triangle_boundary_mod2": 0,
        "rank_like_cycle_count_before_2_cell_fill": beta_1_graph,
        "rank_like_cycle_count_after_2_cell_fill": beta_1_after_filling_triangles,
    }
    return {
        "vertex_count": vertex_count,
        "edge_count": edge_count,
        "filled_triangle_count": face_count,
        "connected_components": connected_components,
        "graph_beta_1": beta_1_graph,
        "beta_1_after_triangle_fill": beta_1_after_filling_triangles,
        "euler_characteristic": euler_characteristic,
        "cycle_rank_witness": cycle_rank_witness,
        "topological_numbers_are_scale_invariant": True,
    }


def trace_action(scale: float, exponent: float) -> float:
    return sum(2.0 / (distance(a, b, scale) ** exponent) for a, b in itertools.combinations(NODES, 2))


def dimension_typing_audit(invariants: dict[str, Any]) -> dict[str, Any]:
    candidate_terms = [
        {
            "candidate": "lambda_equals_beta1_graph",
            "topological_number": invariants["graph_beta_1"],
            "dimension": "dimensionless_integer",
        },
        {
            "candidate": "lambda_equals_beta1_after_triangle_fill",
            "topological_number": invariants["beta_1_after_triangle_fill"],
            "dimension": "dimensionless_integer",
        },
        {
            "candidate": "lambda_equals_euler_characteristic",
            "topological_number": invariants["euler_characteristic"],
            "dimension": "dimensionless_integer",
        },
    ]
    rows = []
    for term in candidate_terms:
        rows.append({
            **term,
            "can_be_added_to_trace_action_without_unit_map": False,
            "reason": "Tr_p(L_a) scales as length^-p, while the cocycle/topological integer is dimensionless; adding them requires an action/unit conversion coefficient.",
        })
    return {
        "statement": "Boundary/cocycle data can supply scale-invariant integers, but those integers are dimensionless.  They do not by themselves supply the dimensionful coefficient needed to add a constant anomaly to Tr_p(L_a) or to define an absolute action quantum.",
        "rows": rows,
        "all_topological_candidates_need_unit_map": all(not row["can_be_added_to_trace_action_without_unit_map"] for row in rows),
        "dimensionful_anomaly_coefficient_derived_now": False,
    }


def topological_anomaly_phase_audit(invariants: dict[str, Any]) -> dict[str, Any]:
    rows = []
    max_phase_error = 0.0
    selected_by_declared_quantum = []
    topological_lambda = float(invariants["graph_beta_1"])
    for exponent in EXPONENTS:
        for scale in SCALES:
            action_value = trace_action(scale, exponent) + topological_lambda
            for quantum in QUANTA:
                tau = (2.0 * math.pi * quantum) / action_value
                phase = tau * action_value
                phase_error = abs(phase - 2.0 * math.pi * quantum)
                max_phase_error = max(max_phase_error, phase_error)
                rows.append({
                    "exponent": exponent,
                    "scale": scale,
                    "topological_lambda_used_as_declared_additive_constant": topological_lambda,
                    "action_value": action_value,
                    "integer_quantum_n": quantum,
                    "quantized_tau": tau,
                    "phase_error": phase_error,
                })
        declared_quantum = trace_action(1.0, exponent) + topological_lambda
        passing_scales = [
            scale for scale in SCALES
            if abs((trace_action(scale, exponent) + topological_lambda) - declared_quantum) < TOL
        ]
        selected_by_declared_quantum.append({
            "exponent": exponent,
            "declared_absolute_action_quantum_from_scale_one": declared_quantum,
            "passing_scales": passing_scales,
            "selects_scale_one": passing_scales == [1.0],
            "selection_is_external_absolute_action_anchor": True,
        })
    return {
        "statement": "Even if a topological integer is provisionally inserted as an additive anomaly, integer phase quantization remains satisfiable at every scale by choosing tau.  A unique scale appears only after declaring an absolute action quantum at one representative.",
        "rows": rows,
        "declared_quantum_rows": selected_by_declared_quantum,
        "max_phase_error": max_phase_error,
        "all_integer_phase_conditions_satisfied": max_phase_error < TOL,
        "declared_quantum_selectors_are_external": all(row["selection_is_external_absolute_action_anchor"] for row in selected_by_declared_quantum),
        "uv_unit_selected_by_topological_anomaly_now": False,
    }


def upstream_consistency() -> dict[str, Any]:
    p2658 = load_json(P2658)
    p2659 = load_json(P2659)
    return {
        "p2658_homogeneous_no_go_verified": p2658.get("closure_decision", {}).get("all_homogeneous_covariances_verified") is True,
        "p2659_no_intrinsic_anomaly_source": p2659.get("closure_decision", {}).get("uv_unit_selected_now") is False,
        "p2659_next_step_requests_boundary_cocycle_phase_law": "boundary/cocycle/phase" in p2659.get("closure_decision", {}).get("next_honest_step", ""),
    }


def source_candidate_matrix(dimension_audit: dict[str, Any], phase_audit: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "raw_boundary_cocycle_integer_as_lambda",
            "covered_by_audit": True,
            "uses_external_clock_or_scale_anchor": False,
            "dimensionally_admissible_without_unit_map": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: the cocycle integer is scale-invariant but dimensionless and cannot be added to Tr_p(L) without a unit/action map",
        },
        {
            "candidate": "topological_integer_plus_integer_phase",
            "covered_by_audit": True,
            "uses_external_clock_or_scale_anchor": False,
            "dimensionally_admissible_without_unit_map": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: phase quantization still fixes tau times the candidate action at each scale",
        },
        {
            "candidate": "declared_absolute_action_quantum_after_topology",
            "covered_by_audit": True,
            "uses_external_clock_or_scale_anchor": True,
            "dimensionally_admissible_without_unit_map": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: unique scale comes from the declared absolute action quantum, not from topology alone",
        },
        {
            "candidate": "derived_unit_map_from_nadsoliton_boundary_phase_law",
            "covered_by_audit": False,
            "uses_external_clock_or_scale_anchor": False,
            "dimensionally_admissible_without_unit_map": None,
            "passes_as_uv_unit_source_now": False,
            "verdict": "open theorem target: derive a dimensionful action/unit map internally and rerun this audit",
        },
    ]


def closure_decision(dimension_audit: dict[str, Any], phase_audit: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_uv_unit_source_now"]]
    return {
        "decision": "BOUNDARY_COCYCLE_ANOMALY_COEFFICIENT_DIMENSION_AUDIT__NO_INTRINSIC_UV_UNIT",
        "professorial_verdict": (
            "P2660 addresses the precise theorem target left by P2659: can boundary/cocycle topology source the anomaly coefficient? "
            "The finite audit finds scale-invariant topological integers, but they are dimensionless.  They do not supply the unit/action conversion needed to add a dimensionful anomaly to a Laplacian trace or to define an absolute action quantum. "
            "If inserted provisionally, the integer phase condition still leaves the clock-scale freedom.  Therefore topology is a promising carrier of a discrete datum, but not yet a sourced UV unit or beta source."
        ),
        "next_honest_step": (
            "Build an explicit nadsoliton boundary-phase unit-map theorem: a chain-level law must turn the cocycle integer into a dimensionful action/clock coefficient without importing an external unit.  If no such map can be derived, formalize a no-go separating scale-invariant topological integers from dimensionful UV-unit sourcing."
        ),
        "passing_uv_unit_source_candidates": passing,
        "all_topological_candidates_need_unit_map": dimension_audit["all_topological_candidates_need_unit_map"],
        "all_integer_phase_conditions_satisfied": phase_audit["all_integer_phase_conditions_satisfied"],
        "uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    dimension_audit = payload["dimension_typing_audit"]
    phase_audit = payload["topological_anomaly_phase_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2660/S1610 boundary/cocycle anomaly coefficient dimension audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first audit",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Finite boundary/cocycle witness",
        f"Graph beta_1: `{payload['boundary_cocycle_invariants']['graph_beta_1']}`.",
        f"Beta_1 after triangle fill: `{payload['boundary_cocycle_invariants']['beta_1_after_triangle_fill']}`.",
        f"Euler characteristic: `{payload['boundary_cocycle_invariants']['euler_characteristic']}`.",
        "Topological numbers are scale-invariant, but scale-invariance is not yet a dimensionful action coefficient.",
        "",
        "## Dimension typing audit",
        dimension_audit["statement"],
        f"All topological candidates need a unit map? `{dimension_audit['all_topological_candidates_need_unit_map']}`.",
        f"Dimensionful anomaly coefficient derived now? `{dimension_audit['dimensionful_anomaly_coefficient_derived_now']}`.",
        "",
        "## Phase audit",
        phase_audit["statement"],
        f"All integer phase conditions satisfied by choosing tau? `{phase_audit['all_integer_phase_conditions_satisfied']}`.",
        f"Declared quantum selectors are external? `{phase_audit['declared_quantum_selectors_are_external']}`.",
        f"UV unit selected by topological anomaly now? `{phase_audit['uv_unit_selected_by_topological_anomaly_now']}`.",
        "",
        "## Source candidate matrix",
        "| candidate | covered? | external anchor? | dimensionally admissible without unit map? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['covered_by_audit']}` | `{row['uses_external_clock_or_scale_anchor']}` | "
            f"`{row['dimensionally_admissible_without_unit_map']}` | `{row['passes_as_uv_unit_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Passing UV-unit source candidates: `{decision['passing_uv_unit_source_candidates']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"Role-bearing L_total now? `{decision['role_bearing_ltotal_now']}`.",
        f"ToE closure now? `{decision['toe_closure_now']}`.",
        "",
        "## Next honest step",
        decision["next_honest_step"],
        "",
        "## Negative exports",
    ])
    for key, value in payload["negative_export_flags"].items():
        lines.append(f"- `{key}`: `{value}`")
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    invariants = graph_invariants()
    dimension_audit = dimension_typing_audit(invariants)
    phase_audit = topological_anomaly_phase_audit(invariants)
    matrix = source_candidate_matrix(dimension_audit, phase_audit)
    decision = closure_decision(dimension_audit, phase_audit, matrix)
    payload: dict[str, Any] = {
        "status": "P2660_BOUNDARY_COCYCLE_ANOMALY_COEFFICIENT_DIMENSION_AUDIT_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2658_HOMOGENEOUS_NO_GO": sha256_file(P2658),
            "P2659_NONHOMOGENEOUS_ANOMALY_AUDIT": sha256_file(P2659),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_boundary_model": {
            "nodes": NODES,
            "edges": EDGES,
            "filled_triangles": TRIANGLES,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "no_sub_nadsoliton_information_layer": True,
        },
        "boundary_cocycle_invariants": invariants,
        "dimension_typing_audit": dimension_audit,
        "topological_anomaly_phase_audit": phase_audit,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2660/S1610 boundary/cocycle anomaly coefficient dimension guard",
        "## P2660/S1610 boundary/cocycle anomaly coefficient dimension guard\n\n"
        "`P2660/S1610` audits whether boundary/cocycle topology can source the P2659 anomaly coefficient.  The finite complex supplies scale-invariant integers such as `beta_1` and Euler characteristic, but those are dimensionless; they do not by themselves provide the dimensionful action/unit map needed to add an anomaly to `Tr_p(L_a)` or define an absolute action quantum.  Provisionally inserting the integer still leaves integer phase quantization satisfiable by clock choice.  This exports no intrinsic UV unit, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2660/S1610 boundary/cocycle anomaly Ltotal guard",
        "## P2660/S1610 boundary/cocycle anomaly Ltotal guard\n\n"
        "`P2660/S1610` does not re-open `L_total`: boundary/cocycle integers can carry discrete topology, but a variational anomaly term still requires an internally derived dimensionful nadsoliton boundary-phase unit map.  Until that map exists, beta-source rerun, role-transfer rerun, and selector/source discharge remain blocked.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
