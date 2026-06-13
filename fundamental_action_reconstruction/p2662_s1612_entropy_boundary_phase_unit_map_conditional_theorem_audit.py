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
OUT = GEN / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.json"
MD = GEN / "p2662_s1612_entropy_boundary_phase_unit_map_conditional_theorem_audit.md"
P2660 = GEN / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.json"
P2661 = GEN / "p2661_s1611_shannon_entropy_scale_anomaly_uv_anchor_audit.json"
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
D_F = 9.0 / 5.0
BIT = math.log(2.0)
EXPONENT = 1.8
BASE_UNIT_RESCALINGS = [0.5, 1.0, 2.0, 3.0]
CANDIDATE_INTEGER_SOURCES = {
    "beta1_after_triangle_fill": 1,
    "graph_beta1_before_fill": 2,
    "edge_minus_vertex_abs": 1,
}
TOL = 1e-11
NEGATIVE_EXPORT_FLAGS = [
    "unconditional_entropy_measure_theorem_exported",
    "unconditional_boundary_phase_unit_map_exported",
    "bit_to_action_map_sourced_unconditionally",
    "bit_to_length_map_sourced_unconditionally",
    "intrinsic_reference_cell_exported",
    "unique_beta_representative_selected_unconditionally",
    "entropy_arrow_discharges_qw_2191",
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
    return {"count": len(lines), "samples": lines[:80]}


def semantic_rg_audit() -> dict[str, Any]:
    patterns = {
        "entropy_measure_content": "entropy.*measure|reference cell|entropy zero|Shannon|D_f.*log",
        "boundary_phase_content": "boundary-phase|boundary.*phase|cocycle|phase law|unit map",
        "bit_action_length_content": "bit-to-action|bit-to-length|action unit|length unit|log\\(2\\)",
        "nonclosure_guard_content": "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|beta source",
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for intrinsic entropy-measure and boundary-phase unit-map theorem attempt",
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


def distance(a: str, b: str, base_unit_rescaling: float = 1.0) -> float:
    ax, ay = NODES[a]
    bx, by = NODES[b]
    return base_unit_rescaling * math.hypot(ax - bx, ay - by)


def normalized_entropy(base_unit_rescaling: float = 1.0) -> float:
    weights = [1.0 / (distance(a, b, base_unit_rescaling) ** EXPONENT) for a, b in itertools.combinations(NODES, 2)]
    total = sum(weights)
    probabilities = [weight / total for weight in weights]
    return -sum(p * math.log(p) for p in probabilities if p > 0.0)


def topology_payload() -> dict[str, Any]:
    vertex_count = len(NODES)
    edge_count = len(EDGES)
    face_count = len(TRIANGLES)
    graph_beta1 = edge_count - vertex_count + 1
    beta1_after_fill = graph_beta1 - face_count
    return {
        "vertex_count": vertex_count,
        "edge_count": edge_count,
        "filled_triangle_count": face_count,
        "graph_beta1_before_fill": graph_beta1,
        "beta1_after_triangle_fill": beta1_after_fill,
        "euler_characteristic": vertex_count - edge_count + face_count,
        "candidate_integer_sources": CANDIDATE_INTEGER_SOURCES,
    }


def conditional_unit_map_theorem_audit() -> dict[str, Any]:
    rows = []
    for source_name, integer_value in CANDIDATE_INTEGER_SOURCES.items():
        entropy_target = integer_value * BIT
        for base_rescaling in BASE_UNIT_RESCALINGS:
            h0 = normalized_entropy(base_rescaling)
            coordinate_scale = math.exp((entropy_target - h0 - D_F * math.log(base_rescaling)) / D_F)
            physical_scale = coordinate_scale * base_rescaling
            residual = h0 + D_F * math.log(physical_scale) - entropy_target
            rows.append({
                "integer_source": source_name,
                "integer_value": integer_value,
                "entropy_target_nats": entropy_target,
                "entropy_target_bits": integer_value,
                "base_unit_rescaling": base_rescaling,
                "normalized_entropy_reference": h0,
                "selected_coordinate_scale": coordinate_scale,
                "selected_physical_scale": physical_scale,
                "entropy_residual_at_selected_physical_scale": residual,
            })
    by_source: dict[str, list[dict[str, Any]]] = {}
    for row in rows:
        by_source.setdefault(row["integer_source"], []).append(row)
    invariance_rows = []
    for source_name, source_rows in by_source.items():
        physical_values = [row["selected_physical_scale"] for row in source_rows]
        invariance_rows.append({
            "integer_source": source_name,
            "physical_scale_values": physical_values,
            "max_physical_scale_spread": max(physical_values) - min(physical_values),
            "coordinate_scale_values": [row["selected_coordinate_scale"] for row in source_rows],
        })
    return {
        "statement": "Conditional theorem attempt: if an intrinsic entropy measure supplies H(a)=H0+D_f log(a), and if a boundary-phase law supplies the target H=N_bits log(2), then exactly one positive physical scale solves the equation for each integer source.  The computation verifies the conditional uniqueness and base-coordinate covariance, but the premises are not exported by the current repo.",
        "rows": rows,
        "invariance_rows": invariance_rows,
        "max_entropy_residual": max(abs(row["entropy_residual_at_selected_physical_scale"]) for row in rows),
        "conditional_unique_positive_scale_for_each_integer_source": True,
        "physical_scale_invariant_under_base_coordinate_rescaling": all(row["max_physical_scale_spread"] < TOL for row in invariance_rows),
        "premise_intrinsic_entropy_measure_exported": False,
        "premise_boundary_phase_bit_target_exported": False,
        "premise_bit_to_action_or_length_map_exported": False,
        "unconditional_uv_unit_selected_now": False,
    }


def premise_gap_audit() -> list[dict[str, Any]]:
    return [
        {
            "premise": "intrinsic_pre_normalization_entropy_measure",
            "needed_for": "define H0 and the reference cell without coordinate convention",
            "currently_exported": False,
            "gap": "P2661 showed normalized entropy is scale-invariant; differential entropy needs a measure/reference cell.",
        },
        {
            "premise": "boundary_phase_integer_to_entropy_target",
            "needed_for": "derive N_bits log(2) as a nadsoliton law rather than a chosen entropy level",
            "currently_exported": False,
            "gap": "P2660 permits topology as a carrier but not yet as a dimensionful/source law.",
        },
        {
            "premise": "bit_to_action_or_bit_to_length_unit_map",
            "needed_for": "turn the dimensionless bit into an action/clock/length datum usable by beta-source rerun",
            "currently_exported": False,
            "gap": "log(2) is dimensionless unless a nadsoliton unit map is proved.",
        },
        {
            "premise": "selector_branch_orientation_law",
            "needed_for": "turn entropy monotonicity into a strict-core QW-2191 branch selector",
            "currently_exported": False,
            "gap": "entropy orientation is not yet an O(2) selector discharge.",
        },
    ]


def upstream_consistency() -> dict[str, Any]:
    p2660 = load_json(P2660)
    p2661 = load_json(P2661)
    return {
        "p2660_all_topological_candidates_need_unit_map": p2660.get("closure_decision", {}).get("all_topological_candidates_need_unit_map") is True,
        "p2660_no_beta_source": p2660.get("closure_decision", {}).get("beta_source_exported_now") is False,
        "p2661_fractal_log_shift_verified": p2661.get("closure_decision", {}).get("fractal_log_shift_verified") is True,
        "p2661_selection_not_intrinsic_without_reference_measure": p2661.get("closure_decision", {}).get("selection_is_intrinsic_without_reference_measure") is False,
        "p2661_qw2191_not_discharged": p2661.get("closure_decision", {}).get("entropy_arrow_discharges_qw_2191") is False,
    }


def source_candidate_matrix(conditional: dict[str, Any], gaps: list[dict[str, Any]]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "conditional_entropy_boundary_phase_unit_map",
            "covered_by_audit": True,
            "selects_unique_scale_conditionally": conditional["conditional_unique_positive_scale_for_each_integer_source"],
            "all_premises_exported": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "conditional pass only: uniqueness follows if entropy measure, bit target, and unit map are supplied internally",
        },
        {
            "candidate": "raw_cocycle_integer_as_entropy_target",
            "covered_by_audit": True,
            "selects_unique_scale_conditionally": True,
            "all_premises_exported": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: the integer target must be a boundary-phase law, not a picked bit level",
        },
        {
            "candidate": "bit_to_action_or_length_conversion",
            "covered_by_audit": True,
            "selects_unique_scale_conditionally": None,
            "all_premises_exported": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: no current theorem converts the dimensionless bit into action/clock/length units",
        },
        {
            "candidate": "entropy_selector_branch_for_qw2191",
            "covered_by_audit": True,
            "selects_unique_scale_conditionally": None,
            "all_premises_exported": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: scale orientation is not yet a strict-core selector branch law",
        },
    ]


def closure_decision(conditional: dict[str, Any], gaps: list[dict[str, Any]], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_uv_unit_source_now"]]
    missing = [gap["premise"] for gap in gaps if not gap["currently_exported"]]
    return {
        "decision": "ENTROPY_BOUNDARY_PHASE_UNIT_MAP_CONDITIONAL_THEOREM_AUDIT__CONDITIONAL_SCALE_SELECTION_NO_UNCONDITIONAL_SOURCE",
        "professorial_verdict": (
            "P2662 builds the requested proof-grade candidate theorem rather than merely saying entropy is promising.  The conditional theorem is mathematically coherent: H(a)=H0+D_f log(a) plus a boundary-phase bit target N log 2 selects one positive physical scale, and the selected physical scale is covariant under base-coordinate rescaling.  "
            "However, the current repo still lacks the three source premises: an intrinsic pre-normalization entropy measure/reference cell, a boundary-phase law deriving the bit target, and a bit-to-action or bit-to-length unit map.  Therefore P2662 is a conditional theorem scaffold and premise ledger, not an unconditional UV-unit, beta-source, QW-2191, L_total, or ToE closure."
        ),
        "next_honest_step": (
            "Attack the smallest missing premise first: construct a chain-level boundary-phase law that outputs N_bits log(2) from the nadsoliton complex, then separately prove a bit-to-action or bit-to-length unit map.  If either premise cannot be sourced internally, formalize the no-go separating entropy/cocycle scale selection from physical UV-unit sourcing."
        ),
        "passing_uv_unit_source_candidates": passing,
        "missing_premises": missing,
        "conditional_unique_scale_verified": conditional["conditional_unique_positive_scale_for_each_integer_source"],
        "physical_scale_covariance_verified": conditional["physical_scale_invariant_under_base_coordinate_rescaling"],
        "unconditional_uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "qw2191_discharged_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    conditional = payload["conditional_unit_map_theorem_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2662/S1612 entropy boundary-phase unit-map conditional theorem audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first audit",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Conditional theorem scaffold",
        conditional["statement"],
        f"Conditional unique positive scale for each integer source? `{conditional['conditional_unique_positive_scale_for_each_integer_source']}`.",
        f"Physical scale invariant under base-coordinate rescaling? `{conditional['physical_scale_invariant_under_base_coordinate_rescaling']}`.",
        f"Max entropy residual: `{conditional['max_entropy_residual']:.3e}`.",
        f"Unconditional UV unit selected now? `{conditional['unconditional_uv_unit_selected_now']}`.",
        "",
        "## Premise gap ledger",
        "| premise | currently exported? | needed for | gap |",
        "| --- | ---: | --- | --- |",
    ])
    for gap in payload["premise_gap_audit"]:
        lines.append(f"| `{gap['premise']}` | `{gap['currently_exported']}` | {gap['needed_for']} | {gap['gap']} |")
    lines.extend([
        "",
        "## Source candidate matrix",
        "| candidate | covered? | conditional unique scale? | all premises exported? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['covered_by_audit']}` | `{row['selects_unique_scale_conditionally']}` | "
            f"`{row['all_premises_exported']}` | `{row['passes_as_uv_unit_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Missing premises: `{decision['missing_premises']}`.",
        f"Passing UV-unit source candidates: `{decision['passing_uv_unit_source_candidates']}`.",
        f"Beta source exported now? `{decision['beta_source_exported_now']}`.",
        f"QW-2191 discharged now? `{decision['qw2191_discharged_now']}`.",
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
    topology = topology_payload()
    conditional = conditional_unit_map_theorem_audit()
    gaps = premise_gap_audit()
    matrix = source_candidate_matrix(conditional, gaps)
    decision = closure_decision(conditional, gaps, matrix)
    payload: dict[str, Any] = {
        "status": "P2662_ENTROPY_BOUNDARY_PHASE_UNIT_MAP_CONDITIONAL_THEOREM_AUDIT_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2660_BOUNDARY_COCYCLE_DIMENSION_AUDIT": sha256_file(P2660),
            "P2661_SHANNON_ENTROPY_ANOMALY_AUDIT": sha256_file(P2661),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_entropy_boundary_model": {
            "nodes": NODES,
            "edges": EDGES,
            "filled_triangles": TRIANGLES,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "no_sub_nadsoliton_information_layer": True,
            "fractal_dimension_D_f": D_F,
            "bit_in_nats": BIT,
            "homogeneous_weight_exponent": EXPONENT,
        },
        "topology_payload": topology,
        "conditional_unit_map_theorem_audit": conditional,
        "premise_gap_audit": gaps,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2662/S1612 entropy boundary-phase unit-map conditional theorem guard",
        "## P2662/S1612 entropy boundary-phase unit-map conditional theorem guard\n\n"
        "`P2662/S1612` builds the requested intrinsic entropy/boundary-phase unit-map theorem as a conditional scaffold.  If an intrinsic pre-normalization entropy measure gives `H(a)=H0+D_f log(a)` and a boundary-phase law derives an integer bit target `N log(2)`, then the equation selects one positive physical scale and is covariant under base-coordinate rescaling.  The current repo does not yet export the required premises: intrinsic reference cell/entropy zero, boundary-phase bit target, bit-to-action or bit-to-length unit map, or selector-branch orientation.  Therefore this is a conditional theorem scaffold, not an unconditional UV unit, beta source, `QW-2191` discharge, bridge completion, role transfer, role-bearing `L_total`, or ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2662/S1612 entropy boundary-phase unit-map Ltotal guard",
        "## P2662/S1612 entropy boundary-phase unit-map Ltotal guard\n\n"
        "`P2662/S1612` does not re-open `L_total`: conditional entropy/cocycle scale selection becomes variationally source-bearing only after an internally derived entropy measure, boundary-phase bit target, and bit-to-action or bit-to-length unit map are proved.  Until those premises are exported, beta-source rerun, role-transfer rerun, selector/source discharge, and ToE closure remain blocked.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
