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
OUT = GEN / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.json"
MD = GEN / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.md"

P2653 = GEN / "p2653_s1603_typed_nadsoliton_metric_uv_source_obligation_rank_audit.json"
P2654 = GEN / "p2654_s1604_scale_invariant_beta_selector_no_go_and_equivariant_functional_rank_certificate.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

ETA = 9.0 / 5.0
BASE_BETA = 1.0
BASE_POINTS: dict[str, tuple[float, float]] = {
    "n0_origin": (0.0, 0.0),
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
    "scale_orbit_quotient_discharged_as_beta_selector",
    "dimensionless_operator_identity_exported",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
    "raw_metric_normalization_promoted_to_source",
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
        "typed_state_space_metric_content": (
            "typed nadsoliton state space|state space|metric axioms|positive metric|distance functional|"
            "typed nadsoliton metric"
        ),
        "scale_quotient_content": (
            "scale orbit|scale quotient|normalization orbit|gauge-fixed|gauge fixed|raw coordinate|"
            "unit selector|canonical unit"
        ),
        "operator_identity_content": (
            "operator identity|conservation identity|dimensionless operator|dimensionless conservation|"
            "unique beta|beta-selecting"
        ),
        "nonclosure_guard_content": (
            "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|"
            "source theorem|beta source"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for typed metric state-space scale-quotient pretheorem, not packet-name search",
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


def scaled_distance(scale: float, a: str, b: str) -> float:
    return scale * base_distance(a, b)


def distance_matrix(scale: float) -> dict[str, dict[str, float]]:
    return {
        a: {b: scaled_distance(scale, a, b) for b in BASE_POINTS}
        for a in BASE_POINTS
    }


def normalized_distance_vector(scale: float) -> list[float]:
    raw = []
    for a, b in itertools.combinations(BASE_POINTS, 2):
        raw.append(scaled_distance(scale, a, b))
    diameter = max(raw)
    return [value / diameter for value in raw]


def metric_axiom_audit() -> dict[str, Any]:
    scale_rows = []
    max_triangle_violation = 0.0
    max_normalized_difference = 0.0
    reference_normalized = normalized_distance_vector(1.0)
    for scale in AUDITED_SCALES:
        nonnegative = True
        identity = True
        symmetric = True
        triangle = True
        local_triangle_violation = 0.0
        for a in BASE_POINTS:
            for b in BASE_POINTS:
                dab = scaled_distance(scale, a, b)
                dba = scaled_distance(scale, b, a)
                nonnegative = nonnegative and dab >= -TOL
                identity = identity and ((a == b and abs(dab) <= TOL) or (a != b and dab > TOL))
                symmetric = symmetric and abs(dab - dba) <= TOL
                for c in BASE_POINTS:
                    violation = dab - (scaled_distance(scale, a, c) + scaled_distance(scale, c, b))
                    local_triangle_violation = max(local_triangle_violation, violation)
                    triangle = triangle and violation <= TOL
        normalized = normalized_distance_vector(scale)
        max_norm_diff = max(abs(a - b) for a, b in zip(normalized, reference_normalized))
        max_triangle_violation = max(max_triangle_violation, local_triangle_violation)
        max_normalized_difference = max(max_normalized_difference, max_norm_diff)
        scale_rows.append({
            "scale": scale,
            "metric_axioms_pass": nonnegative and identity and symmetric and triangle,
            "nonnegative": nonnegative,
            "identity_of_indiscernibles": identity,
            "symmetric": symmetric,
            "triangle": triangle,
            "max_triangle_violation": local_triangle_violation,
            "raw_diameter": max(scaled_distance(scale, a, b) for a, b in itertools.combinations(BASE_POINTS, 2)),
            "normalized_distance_vector": normalized,
            "max_normalized_difference_from_scale_one": max_norm_diff,
        })
    return {
        "statement": "A typed finite nadsoliton metric skeleton can satisfy metric axioms for every positive global scale, while its normalized quotient geometry is unchanged; metric axioms alone do not select the UV unit.",
        "typed_order_guard": "Nodes are typed in the preferred order nadsoliton -> light -> matter -> emergent observer; no layer is placed underneath the nadsoliton.",
        "scale_rows": scale_rows,
        "all_scales_metric_axioms_pass": all(row["metric_axioms_pass"] for row in scale_rows),
        "normalized_geometry_rank_on_scale_orbit": 1 if max_normalized_difference < TOL else len(AUDITED_SCALES),
        "max_triangle_violation": max_triangle_violation,
        "max_normalized_difference": max_normalized_difference,
        "unit_selected_by_metric_axioms_alone": False,
    }


def envelope(beta: float, distance: float) -> float:
    return 1.0 / (1.0 + beta * distance**ETA)


def damping_covariance_audit() -> dict[str, Any]:
    pairs = [("n0_origin", "n1_light_axis"), ("n0_origin", "n3_observer_downstream"), ("n1_light_axis", "n4_compression_probe")]
    reference = [envelope(BASE_BETA, scaled_distance(1.0, a, b)) for a, b in pairs]
    rows = []
    for scale in AUDITED_SCALES:
        beta_for_scale = BASE_BETA / (scale**ETA)
        values = [envelope(beta_for_scale, scaled_distance(scale, a, b)) for a, b in pairs]
        rows.append({
            "scale": scale,
            "covariant_beta": beta_for_scale,
            "max_envelope_difference_from_scale_one": max(abs(a - b) for a, b in zip(values, reference)),
            "pair_values": values,
        })
    return {
        "statement": "Under d -> a d and beta -> beta/a^eta, the strict denominator values are unchanged on the typed metric skeleton; this is a computational scale-orbit witness, not a beta source.",
        "eta": ETA,
        "rows": rows,
        "max_covariance_error": max(row["max_envelope_difference_from_scale_one"] for row in rows),
        "covariance_verified": all(row["max_envelope_difference_from_scale_one"] < TOL for row in rows),
        "beta_one_selected_by_covariance": False,
    }


def source_candidate_matrix(metric: dict[str, Any], covariance: dict[str, Any]) -> list[dict[str, Any]]:
    rows = [
        {
            "candidate": "typed_state_space_plus_metric_axioms_only",
            "metric_axioms_pass": metric["all_scales_metric_axioms_pass"],
            "uses_external_unit_anchor": False,
            "breaks_scale_orbit": False,
            "passes_as_typed_metric_uv_source_now": False,
            "verdict": "honest typed metric skeleton, but all positive global scales remain admissible",
        },
        {
            "candidate": "diameter_equals_one_normalization",
            "metric_axioms_pass": True,
            "uses_external_unit_anchor": True,
            "breaks_scale_orbit": True,
            "passes_as_typed_metric_uv_source_now": False,
            "verdict": "chooses a representative by convention; not derived from nadsoliton dynamics",
        },
        {
            "candidate": "nearest_neighbor_equals_one_lattice_unit",
            "metric_axioms_pass": True,
            "uses_external_unit_anchor": True,
            "breaks_scale_orbit": True,
            "passes_as_typed_metric_uv_source_now": False,
            "verdict": "imports a lattice/unit premise and therefore cannot be the missing source theorem",
        },
        {
            "candidate": "strict_denominator_covariance_identity",
            "metric_axioms_pass": True,
            "uses_external_unit_anchor": False,
            "breaks_scale_orbit": False,
            "passes_as_typed_metric_uv_source_now": False,
            "verdict": "covariance is verified, but invariance across beta-distance representatives blocks unique beta=1",
        },
        {
            "candidate": "typed_nadsoliton_dynamics_plus_conserved_uv_action_quantum",
            "metric_axioms_pass": None,
            "uses_external_unit_anchor": False,
            "breaks_scale_orbit": None,
            "passes_as_typed_metric_uv_source_now": False,
            "verdict": "only acceptable theorem target; currently hypothetical because the dynamics and conserved dimensionless identity are not exported",
        },
    ]
    return rows


def upstream_consistency() -> dict[str, Any]:
    p2653 = load_json(P2653)
    p2654 = load_json(P2654)
    return {
        "p2653_typed_metric_uv_source_not_exported": p2653.get("closure_decision", {}).get("typed_metric_uv_source_theorem_exported_now") is False,
        "p2653_missing_atoms_include_uv_unit": "uv_unit_selected_by_nadsoliton_dynamics" in p2653.get("obligation_rank_audit", {}).get("currently_missing_atoms_union", []),
        "p2654_scale_invariant_selector_blocked": p2654.get("closure_decision", {}).get("scale_invariant_selector_exists_now") is False,
        "p2654_raw_anchor_not_promoted": p2654.get("closure_decision", {}).get("raw_anchor_promoted_to_source_now") is False,
    }


def closure_decision(matrix: list[dict[str, Any]], metric: dict[str, Any], covariance: dict[str, Any]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_typed_metric_uv_source_now"]]
    return {
        "decision": "TYPED_METRIC_SKELETON_CONSTRUCTED__SCALE_QUOTIENT_STILL_BLOCKS_UV_UNIT_SOURCE",
        "professorial_verdict": (
            "P2655 makes the proof route more concrete: a typed nadsoliton state-space metric skeleton can be written down and mechanically checked, and the strict damping denominator is covariant under the expected scale action.  "
            "But those successful checks are exactly quotient-level facts: they do not choose the global UV unit, do not turn beta=1 into a sourced number, and do not re-enable role-bearing dynamics.  "
            "The missing theorem must therefore add a genuine nadsoliton-dynamical unit source or conserved dimensionless operator identity, not another raw normalization convention."
        ),
        "professorial_closure_path": [
            "Keep the typed state-space and metric axioms as a useful theorem scaffold, but mark their global scale as unselected.",
            "Reject diameter=1, nearest-neighbor=1, or beta=1-by-coordinate choices unless an independent nadsoliton dynamics derives the unit.",
            "Next attempt should specify a concrete dynamics or conserved operator whose spectrum/action quantum fixes the UV unit before beta=1 is tested.",
            "Only after that source is supplied should P2649 be rerun as a beta-source route; P2652/P2647/P2648 remain empirical support only.",
        ],
        "next_honest_step": (
            "Build a candidate typed nadsoliton dynamics with an explicit conserved dimensionless action/operator identity that fixes a UV unit before any strict-kernel target is consulted; then run a source-matrix audit to see whether that identity breaks the scale orbit uniquely or collapses back to a normalization convention."
        ),
        "passing_typed_metric_uv_source_candidates": passing,
        "typed_metric_skeleton_constructed": metric["all_scales_metric_axioms_pass"],
        "scale_quotient_rank_one": metric["normalized_geometry_rank_on_scale_orbit"] == 1,
        "strict_denominator_covariance_verified": covariance["covariance_verified"],
        "uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    metric = payload["metric_axiom_audit"]
    covariance = payload["damping_covariance_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2655/S1605 typed nadsoliton metric state-space scale-quotient pretheorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps typed state-space/metric, scale quotient, operator identity, and nonclosure guard content before adding the pretheorem audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Metric skeleton audit",
        "",
        metric["statement"],
        metric["typed_order_guard"],
        f"All audited scales satisfy metric axioms? `{metric['all_scales_metric_axioms_pass']}`.",
        f"Normalized geometry rank on scale orbit: `{metric['normalized_geometry_rank_on_scale_orbit']}`.",
        f"Max normalized distance difference: `{metric['max_normalized_difference']:.3e}`.",
        f"Unit selected by metric axioms alone? `{metric['unit_selected_by_metric_axioms_alone']}`.",
        "",
        "## Damping covariance audit",
        "",
        covariance["statement"],
        f"Max covariance error: `{covariance['max_covariance_error']:.3e}`.",
        f"Covariance verified? `{covariance['covariance_verified']}`.",
        f"Beta=1 selected by covariance? `{covariance['beta_one_selected_by_covariance']}`.",
        "",
        "## Source candidate matrix",
        "",
        "| candidate | external unit anchor? | breaks scale orbit? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['uses_external_unit_anchor']}` | "
            f"`{row['breaks_scale_orbit']}` | `{row['passes_as_typed_metric_uv_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Passing typed metric/UV source candidates: `{decision['passing_typed_metric_uv_source_candidates']}`.",
        f"Typed metric skeleton constructed? `{decision['typed_metric_skeleton_constructed']}`.",
        f"Scale quotient rank one? `{decision['scale_quotient_rank_one']}`.",
        f"Strict denominator covariance verified? `{decision['strict_denominator_covariance_verified']}`.",
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
    metric = metric_axiom_audit()
    covariance = damping_covariance_audit()
    matrix = source_candidate_matrix(metric, covariance)
    decision = closure_decision(matrix, metric, covariance)
    payload: dict[str, Any] = {
        "status": "P2655_TYPED_METRIC_SCALE_QUOTIENT_PRETHEOREM_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2653_TYPED_METRIC_UV_OBLIGATION": sha256_file(P2653),
            "P2654_SCALE_INVARIANT_SELECTOR_NO_GO": sha256_file(P2654),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_state_space": {
            "nodes": BASE_POINTS,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "no_sub_nadsoliton_information_layer": True,
        },
        "metric_axiom_audit": metric,
        "damping_covariance_audit": covariance,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2655/S1605 typed metric scale-quotient pretheorem guard",
        "## P2655/S1605 typed metric scale-quotient pretheorem guard\n\n"
        "`P2655/S1605` constructs and checks a finite typed nadsoliton metric skeleton, verifies metric axioms across positive global scales, and confirms strict-denominator covariance under `d -> a d`, `beta -> beta/a^eta`.  These are useful quotient-level theorem scaffolds, but they do not select the UV unit: diameter/nearest-neighbor normalizations remain external conventions unless sourced by nadsoliton dynamics.  Therefore this exports no typed metric/UV source theorem, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2655/S1605 typed metric scale-quotient Ltotal guard",
        "## P2655/S1605 typed metric scale-quotient Ltotal guard\n\n"
        "`P2655/S1605` does not re-open `L_total`: a metric skeleton and damping covariance only show that the typed route is internally consistent modulo a scale quotient.  A variational damping coefficient still requires a real nadsoliton-dynamical UV unit plus an independent conserved dimensionless operator/action identity, followed by beta-source rerun, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
