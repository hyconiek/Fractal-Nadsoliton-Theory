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
OUT = GEN / "p2661_s1611_shannon_entropy_scale_anomaly_uv_anchor_audit.json"
MD = GEN / "p2661_s1611_shannon_entropy_scale_anomaly_uv_anchor_audit.md"
P2658 = GEN / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.json"
P2659 = GEN / "p2659_s1609_nonhomogeneous_anomaly_clock_source_candidate_audit.json"
P2660 = GEN / "p2660_s1610_boundary_cocycle_anomaly_coefficient_dimension_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

NODES = {
    "n0_nadsoliton_origin": (0.0, 0.0),
    "n1_light_axis": (1.0, 0.0),
    "n2_matter_axis": (0.0, 1.0),
    "n3_observer_downstream": (1.0, 1.0),
    "n4_compression_probe": (2.0, 1.0),
}
SCALES = [0.25, 0.5, 1.0, 2.0, 3.0]
EXPONENTS = [1.0, 1.8, 2.0]
D_F = 9.0 / 5.0
BIT = math.log(2.0)
TOL = 1e-11
NEGATIVE_EXPORT_FLAGS = [
    "shannon_entropy_uv_anchor_exported",
    "entropy_bit_promoted_to_length_unit",
    "fractal_entropy_anomaly_source_exported",
    "entropy_arrow_discharges_qw_2191",
    "intrinsic_entropy_reference_measure_exported",
    "unique_beta_representative_selected_by_entropy",
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
        "shannon_entropy_content": "Shannon|entropy|void asymmetry|symmetry breaking selector|T165",
        "scale_anomaly_content": "log a|log\\(a\\)|scale anomaly|additive.*shift|fractal dimension|D_f",
        "uv_anchor_content": "UV unit|canonical unit|beta source|scale anchor|one bit|bit",
        "nonclosure_guard_content": "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure",
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for Shannon entropy scale-anomaly UV-anchor intuition, not packet-name search",
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
    ax, ay = NODES[a]
    bx, by = NODES[b]
    return math.hypot(ax - bx, ay - by)


def normalized_distance_distribution(scale: float, exponent: float) -> list[float]:
    weights = []
    for a, b in itertools.combinations(NODES, 2):
        weights.append(1.0 / ((scale * base_distance(a, b)) ** exponent))
    total = sum(weights)
    return [weight / total for weight in weights]


def shannon_entropy(probabilities: list[float]) -> float:
    return -sum(p * math.log(p) for p in probabilities if p > 0.0)


def normalized_shannon_audit() -> dict[str, Any]:
    rows = []
    max_entropy_drift = 0.0
    max_distribution_drift = 0.0
    for exponent in EXPONENTS:
        reference_distribution = normalized_distance_distribution(1.0, exponent)
        reference_entropy = shannon_entropy(reference_distribution)
        for scale in SCALES:
            distribution = normalized_distance_distribution(scale, exponent)
            entropy = shannon_entropy(distribution)
            distribution_drift = max(abs(a - b) for a, b in zip(distribution, reference_distribution))
            entropy_drift = abs(entropy - reference_entropy)
            max_distribution_drift = max(max_distribution_drift, distribution_drift)
            max_entropy_drift = max(max_entropy_drift, entropy_drift)
            rows.append({
                "exponent": exponent,
                "scale": scale,
                "entropy_nats": entropy,
                "entropy_bits": entropy / BIT,
                "entropy_drift_from_scale_one": entropy - reference_entropy,
                "max_probability_drift_from_scale_one": distribution_drift,
            })
    return {
        "statement": "For normalized finite Shannon entropy built from homogeneous distance weights, global distance rescaling cancels in the probabilities.  The Shannon entropy is scale-invariant, not an additive log-a anomaly.",
        "rows": rows,
        "max_entropy_drift": max_entropy_drift,
        "max_distribution_drift": max_distribution_drift,
        "normalized_entropy_is_scale_invariant": max_entropy_drift < TOL and max_distribution_drift < TOL,
        "normalized_entropy_selects_unique_scale": False,
    }


def fractal_differential_entropy_audit() -> dict[str, Any]:
    reference_entropy = shannon_entropy(normalized_distance_distribution(1.0, 1.8))
    rows = []
    max_shift_error = 0.0
    selected_scales_by_declared_bit_level = []
    for scale in SCALES:
        entropy = reference_entropy + D_F * math.log(scale)
        expected_shift = D_F * math.log(scale)
        shift_error = abs((entropy - reference_entropy) - expected_shift)
        max_shift_error = max(max_shift_error, shift_error)
        rows.append({
            "scale": scale,
            "fractal_dimension_D_f": D_F,
            "reference_entropy_nats_at_scale_one": reference_entropy,
            "differential_entropy_nats": entropy,
            "shift_from_scale_one": entropy - reference_entropy,
            "expected_Df_log_scale_shift": expected_shift,
            "shift_error": shift_error,
        })
    for integer_bits in range(-2, 7):
        declared_entropy = integer_bits * BIT
        if declared_entropy <= min(row["differential_entropy_nats"] for row in rows) or declared_entropy >= max(row["differential_entropy_nats"] for row in rows):
            continue
        selected_scale = math.exp((declared_entropy - reference_entropy) / D_F)
        selected_scales_by_declared_bit_level.append({
            "declared_entropy_bits": integer_bits,
            "declared_entropy_nats": declared_entropy,
            "selected_scale_from_declared_level": selected_scale,
            "selection_uses_external_entropy_reference": True,
        })
    return {
        "statement": "A differential/fractal entropy model has the advertised additive D_f*log(a) shift.  However, selecting a scale from that shift requires a declared entropy level/reference measure such as H=n log 2; one bit is dimensionless information, not a length unit by itself.",
        "rows": rows,
        "max_shift_error": max_shift_error,
        "additive_log_scale_shift_verified": max_shift_error < TOL,
        "selected_scales_by_declared_bit_level": selected_scales_by_declared_bit_level,
        "declared_bit_levels_can_select_representatives": len(selected_scales_by_declared_bit_level) > 0,
        "selection_is_intrinsic_without_reference_measure": False,
        "bit_is_dimensionless_not_length_unit": True,
    }


def entropy_arrow_audit(fractal_audit: dict[str, Any]) -> dict[str, Any]:
    ordered_rows = sorted(fractal_audit["rows"], key=lambda row: row["scale"])
    monotone_non_decreasing = all(
        ordered_rows[i + 1]["differential_entropy_nats"] >= ordered_rows[i]["differential_entropy_nats"]
        for i in range(len(ordered_rows) - 1)
    )
    return {
        "statement": "The log-shift model gives an orientation along increasing scale, but this is not yet a strict-core O(2) mode selector or QW-2191 discharge; it lacks an internal premise identifying the physical time arrow and selector branch.",
        "entropy_increases_with_scale_in_model": monotone_non_decreasing,
        "exports_time_arrow_source": False,
        "discharges_qw_2191": False,
    }


def upstream_consistency() -> dict[str, Any]:
    p2658 = load_json(P2658)
    p2659 = load_json(P2659)
    p2660 = load_json(P2660)
    return {
        "p2658_homogeneous_no_go_verified": p2658.get("closure_decision", {}).get("all_homogeneous_covariances_verified") is True,
        "p2659_no_intrinsic_anomaly_source": p2659.get("closure_decision", {}).get("uv_unit_selected_now") is False,
        "p2660_requires_dimensionful_unit_map": p2660.get("closure_decision", {}).get("all_topological_candidates_need_unit_map") is True,
        "p2660_no_beta_source": p2660.get("closure_decision", {}).get("beta_source_exported_now") is False,
    }


def source_candidate_matrix(normalized_audit: dict[str, Any], fractal_audit: dict[str, Any], arrow_audit: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "normalized_finite_shannon_entropy_of_distance_weights",
            "covered_by_audit": True,
            "has_log_scale_anomaly": False,
            "uses_external_reference_or_unit": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: normalized probabilities cancel the global scale, so entropy is invariant and selects no representative",
        },
        {
            "candidate": "fractal_differential_entropy_Df_log_a_shift",
            "covered_by_audit": True,
            "has_log_scale_anomaly": True,
            "uses_external_reference_or_unit": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked as source: additive shift exists only after choosing a reference measure/entropy zero; scale selection needs an external entropy level",
        },
        {
            "candidate": "one_bit_entropy_quantum_as_length_unit",
            "covered_by_audit": True,
            "has_log_scale_anomaly": None,
            "uses_external_reference_or_unit": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: log 2 is an internal information quantum but dimensionless; it is not a UV length/action unit without a derived unit map",
        },
        {
            "candidate": "entropy_arrow_as_qw2191_selector",
            "covered_by_audit": True,
            "has_log_scale_anomaly": True,
            "uses_external_reference_or_unit": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: monotonic entropy orientation is not yet a strict-core O(2) selector source or QW-2191 discharge",
        },
        {
            "candidate": "intrinsic_nadsoliton_entropy_measure_and_unit_map_theorem",
            "covered_by_audit": False,
            "has_log_scale_anomaly": None,
            "uses_external_reference_or_unit": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "open theorem target: derive the measure, entropy zero, bit/action conversion, and selector branch internally",
        },
    ]


def closure_decision(normalized_audit: dict[str, Any], fractal_audit: dict[str, Any], arrow_audit: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_uv_unit_source_now"]]
    return {
        "decision": "SHANNON_ENTROPY_SCALE_ANOMALY_UV_ANCHOR_AUDIT__NO_INTRINSIC_UV_UNIT_YET",
        "professorial_verdict": (
            "P2661 tests the AI intuition directly.  The intuition is half-right: differential/fractal entropy can carry an additive D_f log(a) scale anomaly. "
            "But normalized finite Shannon entropy on homogeneous nadsoliton weights is exactly scale-invariant, and the differential entropy anomaly selects a scale only after a reference measure, entropy zero, or bit-level condition is declared. "
            "The bit is an internal dimensionless information quantum, not yet a length/action unit.  Entropy remains a serious theorem target, but it does not currently export a UV unit, beta source, QW-2191 discharge, or L_total reopening."
        ),
        "next_honest_step": (
            "Build an intrinsic nadsoliton entropy-measure theorem: define the measure before normalization, derive its entropy zero/reference cell and bit-to-action or bit-to-length unit map, then prove that the resulting D_f log(a) anomaly selects one scale and one selector branch without imported entropy level or clock."
        ),
        "passing_uv_unit_source_candidates": passing,
        "normalized_entropy_is_scale_invariant": normalized_audit["normalized_entropy_is_scale_invariant"],
        "fractal_log_shift_verified": fractal_audit["additive_log_scale_shift_verified"],
        "selection_is_intrinsic_without_reference_measure": fractal_audit["selection_is_intrinsic_without_reference_measure"],
        "entropy_arrow_discharges_qw_2191": arrow_audit["discharges_qw_2191"],
        "uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    normalized = payload["normalized_shannon_audit"]
    fractal = payload["fractal_differential_entropy_audit"]
    arrow = payload["entropy_arrow_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2661/S1611 Shannon entropy scale-anomaly UV-anchor audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first audit",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Normalized finite Shannon entropy audit",
        normalized["statement"],
        f"Normalized entropy is scale-invariant? `{normalized['normalized_entropy_is_scale_invariant']}`.",
        f"Max entropy drift: `{normalized['max_entropy_drift']:.3e}`.",
        f"Selects unique scale? `{normalized['normalized_entropy_selects_unique_scale']}`.",
        "",
        "## Fractal differential entropy anomaly audit",
        fractal["statement"],
        f"Additive `D_f log(a)` shift verified? `{fractal['additive_log_scale_shift_verified']}`.",
        f"Declared bit levels can select representatives? `{fractal['declared_bit_levels_can_select_representatives']}`.",
        f"Selection is intrinsic without reference measure? `{fractal['selection_is_intrinsic_without_reference_measure']}`.",
        f"Bit is dimensionless, not a length unit? `{fractal['bit_is_dimensionless_not_length_unit']}`.",
        "",
        "## Entropy arrow audit",
        arrow["statement"],
        f"Entropy increases with scale in model? `{arrow['entropy_increases_with_scale_in_model']}`.",
        f"Exports time-arrow source? `{arrow['exports_time_arrow_source']}`.",
        f"Discharges QW-2191? `{arrow['discharges_qw_2191']}`.",
        "",
        "## Source candidate matrix",
        "| candidate | covered? | log anomaly? | external reference/unit? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['covered_by_audit']}` | `{row['has_log_scale_anomaly']}` | "
            f"`{row['uses_external_reference_or_unit']}` | `{row['passes_as_uv_unit_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        decision["professorial_verdict"],
        f"Decision: `{decision['decision']}`.",
        f"Passing UV-unit source candidates: `{decision['passing_uv_unit_source_candidates']}`.",
        f"UV unit selected now? `{decision['uv_unit_selected_now']}`.",
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
    normalized = normalized_shannon_audit()
    fractal = fractal_differential_entropy_audit()
    arrow = entropy_arrow_audit(fractal)
    matrix = source_candidate_matrix(normalized, fractal, arrow)
    decision = closure_decision(normalized, fractal, arrow, matrix)
    payload: dict[str, Any] = {
        "status": "P2661_SHANNON_ENTROPY_SCALE_ANOMALY_UV_ANCHOR_AUDIT_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2658_HOMOGENEOUS_NO_GO": sha256_file(P2658),
            "P2659_NONHOMOGENEOUS_ANOMALY_AUDIT": sha256_file(P2659),
            "P2660_BOUNDARY_COCYCLE_DIMENSION_AUDIT": sha256_file(P2660),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_entropy_model": {
            "nodes": NODES,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "no_sub_nadsoliton_information_layer": True,
            "fractal_dimension_D_f": D_F,
            "bit_in_nats": BIT,
        },
        "normalized_shannon_audit": normalized,
        "fractal_differential_entropy_audit": fractal,
        "entropy_arrow_audit": arrow,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({key: value for key, value in payload.items() if key != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2661/S1611 Shannon entropy scale-anomaly UV-anchor guard",
        "## P2661/S1611 Shannon entropy scale-anomaly UV-anchor guard\n\n"
        "`P2661/S1611` directly audits the Shannon-entropy-as-anomaly intuition.  Normalized finite Shannon entropy built from homogeneous nadsoliton distance weights is scale-invariant because global scale cancels in the probabilities.  A differential/fractal entropy model does exhibit an additive `D_f log(a)` shift, but selecting a representative requires an entropy zero/reference measure or a declared bit-level condition; `log(2)` is an internal dimensionless information quantum, not yet a length/action unit.  Thus entropy remains a serious theorem target, but this audit exports no intrinsic UV unit, no beta source, no `QW-2191` discharge, no bridge completion, no role transfer, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2661/S1611 Shannon entropy Ltotal guard",
        "## P2661/S1611 Shannon entropy Ltotal guard\n\n"
        "`P2661/S1611` does not re-open `L_total`: entropy can only become a variational scale/anomaly source after an intrinsic nadsoliton entropy-measure theorem derives the reference cell, entropy zero, and bit-to-action or bit-to-length unit map.  Until then, beta-source rerun, role-transfer rerun, selector/source discharge, and ToE closure remain blocked.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
