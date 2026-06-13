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
OUT = GEN / "p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.json"
MD = GEN / "p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.md"

P2655 = GEN / "p2655_s1605_typed_nadsoliton_metric_state_space_scale_quotient_pretheorem.json"
P2656 = GEN / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.json"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"

BASE_POINTS: dict[str, tuple[float, float]] = {
    "n0_nadsoliton_origin": (0.0, 0.0),
    "n1_light_axis": (1.0, 0.0),
    "n2_matter_axis": (0.0, 1.0),
    "n3_observer_downstream": (1.0, 1.0),
    "n4_compression_probe": (2.0, 1.0),
}
AUDITED_SCALES = [0.25, 0.5, 1.0, 2.0, 3.0]
ACTION_QUANTA = [1, 2, 3, 5, 8]
TOL = 1e-12

NEGATIVE_EXPORT_FLAGS = [
    "action_quantization_theorem_exported",
    "intrinsic_clock_or_flow_time_exported",
    "uv_unit_selected_by_quantization",
    "typed_metric_uv_source_theorem_exported",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
    "laplacian_trace_anchor_promoted_to_source",
    "spectral_gap_anchor_promoted_to_source",
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
        "action_quantization_content": (
            "action quantization|action quantum|quantization theorem|integer.*action|phase condition|"
            "conserved integer|conserved.*phase"
        ),
        "clock_scale_anchor_content": (
            "clock|flow time|heat time|spectral gap|trace\\(L\\)|absolute trace|"
            "scale anchor|normalization convention"
        ),
        "uv_beta_source_content": (
            "UV unit|canonical unit|beta source|beta=1|target-independent beta|"
            "typed metric/UV"
        ),
        "nonclosure_guard_content": (
            "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|"
            "source theorem"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for action-quantization scale-anchor obstruction, not packet-name search",
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


def laplacian_trace(scale: float) -> float:
    total = 0.0
    for a, b in itertools.combinations(BASE_POINTS, 2):
        d = scale * base_distance(a, b)
        weight = 1.0 / (d * d)
        total += 2.0 * weight
    return total


def quantized_time(scale: float, quantum: int) -> float:
    return (2.0 * math.pi * quantum) / laplacian_trace(scale)


def action_phase(scale: float, tau: float) -> float:
    return tau * laplacian_trace(scale)


def action_quantization_family_audit() -> dict[str, Any]:
    reference_trace = laplacian_trace(1.0)
    rows = []
    max_trace_covariance_error = 0.0
    max_time_ratio_error = 0.0
    for scale in AUDITED_SCALES:
        trace_l = laplacian_trace(scale)
        trace_covariance_error = abs(trace_l * scale * scale - reference_trace)
        max_trace_covariance_error = max(max_trace_covariance_error, trace_covariance_error)
        for quantum in ACTION_QUANTA:
            tau = quantized_time(scale, quantum)
            reference_tau = quantized_time(1.0, quantum)
            phase = action_phase(scale, tau)
            phase_error = abs(phase - 2.0 * math.pi * quantum)
            time_ratio_error = abs((tau / reference_tau) - (scale * scale))
            max_time_ratio_error = max(max_time_ratio_error, time_ratio_error)
            rows.append({
                "scale": scale,
                "integer_quantum_n": quantum,
                "trace_l": trace_l,
                "trace_l_times_scale_squared": trace_l * scale * scale,
                "quantized_tau_n": tau,
                "tau_ratio_to_scale_one_same_n": tau / reference_tau,
                "expected_tau_ratio_a_squared": scale * scale,
                "phase_action_tau_trace_l": phase,
                "phase_error_from_2pi_n": phase_error,
            })
    return {
        "statement": "For every audited positive scale and integer n, tau_n(a)=2*pi*n/Tr(L_a) satisfies tau_n(a)*Tr(L_a)=2*pi*n.  Since Tr(L_a)=Tr(L_1)/a^2, the quantized clock rescales as tau_n(a)=a^2*tau_n(1); the integer phase condition fixes only a product, not the absolute UV length.",
        "rows": rows,
        "max_trace_covariance_error": max_trace_covariance_error,
        "max_tau_ratio_error_from_a_squared": max_time_ratio_error,
        "integer_phase_condition_satisfied_all_scales": all(row["phase_error_from_2pi_n"] < TOL for row in rows),
        "quantized_family_count": len(rows),
        "unique_scale_selected_by_integer_phase_alone": False,
    }


def fixed_clock_false_selector_audit() -> dict[str, Any]:
    tau_fixed = quantized_time(1.0, 1)
    rows = []
    for scale in AUDITED_SCALES:
        phase_over_2pi = action_phase(scale, tau_fixed) / (2.0 * math.pi)
        rows.append({
            "scale": scale,
            "fixed_tau_from_scale_one_n1": tau_fixed,
            "phase_over_2pi_with_fixed_tau": phase_over_2pi,
            "looks_integer_n1": abs(phase_over_2pi - 1.0) < TOL,
            "failure_reason_if_promoted": "fixed tau was imported from the scale-one representative, so this is a clock anchor rather than an intrinsic selector",
        })
    return {
        "statement": "A fixed clock tau chosen at scale=1 makes only that representative satisfy the n=1 phase condition, but the clock choice already smuggles in the scale-one unit.",
        "rows": rows,
        "scale_one_only_passes_with_imported_clock": [row["scale"] for row in rows if row["looks_integer_n1"]] == [1.0],
        "fixed_clock_is_external_anchor": True,
        "fixed_clock_selector_admissible_now": False,
    }


def source_candidate_matrix(family: dict[str, Any], fixed_clock: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "integer_phase_condition_tau_trace_l_equals_2pi_n",
            "satisfies_quantization": family["integer_phase_condition_satisfied_all_scales"],
            "uses_external_clock_or_scale_anchor": False,
            "selects_unique_scale": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "integer phase quantizes the product tau*Tr(L), leaving a continuous scale-clock family",
        },
        {
            "candidate": "fixed_tau_from_scale_one_selector",
            "satisfies_quantization": fixed_clock["scale_one_only_passes_with_imported_clock"],
            "uses_external_clock_or_scale_anchor": True,
            "selects_unique_scale": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "selects scale=1 only after importing the scale-one clock; this is a hidden normalization anchor",
        },
        {
            "candidate": "declared_trace_or_gap_action_quantum",
            "satisfies_quantization": True,
            "uses_external_clock_or_scale_anchor": True,
            "selects_unique_scale": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "would be a source only if the action quantum is derived from nadsoliton dynamics rather than declared",
        },
        {
            "candidate": "dimensionless_integer_phase_ratio_identity",
            "satisfies_quantization": True,
            "uses_external_clock_or_scale_anchor": False,
            "selects_unique_scale": False,
            "passes_as_uv_unit_source_now": False,
            "verdict": "dimensionless ratios are invariant across the scale-clock orbit and cannot select the UV unit",
        },
        {
            "candidate": "full_intrinsic_nadsoliton_quantization_theorem",
            "satisfies_quantization": None,
            "uses_external_clock_or_scale_anchor": False,
            "selects_unique_scale": None,
            "passes_as_uv_unit_source_now": False,
            "verdict": "still the required theorem target: derive both the clock/action quantum and the scale from the nadsoliton law",
        },
    ]


def upstream_consistency() -> dict[str, Any]:
    p2655 = load_json(P2655)
    p2656 = load_json(P2656)
    return {
        "p2655_uv_unit_not_selected": p2655.get("closure_decision", {}).get("uv_unit_selected_now") is False,
        "p2655_scale_quotient_rank_one": p2655.get("closure_decision", {}).get("scale_quotient_rank_one") is True,
        "p2656_operator_covariance_verified": p2656.get("closure_decision", {}).get("operator_covariance_verified") is True,
        "p2656_absolute_operator_scale_not_selected": p2656.get("closure_decision", {}).get("absolute_operator_scale_selected") is False,
        "p2656_next_step_requests_action_quantization": "action-quantization" in p2656.get("closure_decision", {}).get("next_honest_step", ""),
    }


def closure_decision(matrix: list[dict[str, Any]], family: dict[str, Any], fixed_clock: dict[str, Any]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_uv_unit_source_now"]]
    return {
        "decision": "ACTION_QUANTIZATION_FAMILY_EXISTS__CLOCK_SCALE_ANCHOR_STILL_NOT_SOURCED",
        "professorial_verdict": (
            "P2657 attempts the action-quantization theorem requested by P2656 in the most conservative finite setting.  "
            "The integer phase condition tau*Tr(L)=2*pi*n is computationally satisfiable for every audited scale by rescaling the clock as tau -> a^2 tau.  "
            "Thus the condition quantizes a scale-clock product, not the UV length itself.  A unique scale appears only if a fixed clock, trace, gap, or action quantum is imported from outside the typed nadsoliton dynamics.  "
            "Therefore no intrinsic UV unit, beta source, role-bearing L_total, or ToE closure is exported."
        ),
        "professorial_closure_path": [
            "Keep the integer phase/action condition as a useful obstruction certificate: it shows exactly where clock-scale degeneracy enters.",
            "Do not use fixed tau, trace(L)=constant, spectral gap=constant, or heat-time=constant as a selector unless the nadsoliton evolution law derives that absolute clock/action datum.",
            "The next proof-grade task must add an intrinsic clock/source premise or prove a no-go for all local Laplacian-type action quantizations under the current scale action.",
            "Only after an intrinsic clock/action datum is derived should P2649 be rerun for beta=1; empirical holdout packets remain support only.",
        ],
        "next_honest_step": (
            "Either derive an intrinsic nadsoliton clock/action quantum from a typed evolution law that is not covariant under tau -> a^2 tau, or formalize a broader no-go theorem for local Laplacian/action quantizations showing that every integer phase condition leaves the same scale-clock orbit unbroken."
        ),
        "passing_uv_unit_source_candidates": passing,
        "integer_phase_condition_satisfied_all_scales": family["integer_phase_condition_satisfied_all_scales"],
        "unique_scale_selected_by_integer_phase_alone": family["unique_scale_selected_by_integer_phase_alone"],
        "fixed_clock_selector_admissible_now": fixed_clock["fixed_clock_selector_admissible_now"],
        "uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    family = payload["action_quantization_family_audit"]
    fixed_clock = payload["fixed_clock_false_selector_audit"]
    decision = payload["closure_decision"]
    lines = [
        "# P2657/S1607 nadsoliton action quantization scale-anchor obstruction audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps action-quantization, clock/scale-anchor, UV/beta-source, and nonclosure guard content before adding the obstruction audit.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Integer phase/action family audit",
        "",
        family["statement"],
        f"Integer phase condition satisfied on all audited scales? `{family['integer_phase_condition_satisfied_all_scales']}`.",
        f"Max trace covariance error: `{family['max_trace_covariance_error']:.3e}`.",
        f"Max tau-ratio error from a^2: `{family['max_tau_ratio_error_from_a_squared']:.3e}`.",
        f"Unique scale selected by integer phase alone? `{family['unique_scale_selected_by_integer_phase_alone']}`.",
        "",
        "## Fixed-clock false selector audit",
        "",
        fixed_clock["statement"],
        f"Scale=1 only passes with imported fixed clock? `{fixed_clock['scale_one_only_passes_with_imported_clock']}`.",
        f"Fixed clock is external anchor? `{fixed_clock['fixed_clock_is_external_anchor']}`.",
        f"Fixed-clock selector admissible now? `{fixed_clock['fixed_clock_selector_admissible_now']}`.",
        "",
        "## Source candidate matrix",
        "",
        "| candidate | external clock/scale anchor? | selects unique scale? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['uses_external_clock_or_scale_anchor']}` | "
            f"`{row['selects_unique_scale']}` | `{row['passes_as_uv_unit_source_now']}` | {row['verdict']} |"
        )
    lines.extend([
        "",
        "## Verdict",
        "",
        decision["professorial_verdict"],
        "",
        f"Decision: `{decision['decision']}`.",
        f"Passing UV-unit source candidates: `{decision['passing_uv_unit_source_candidates']}`.",
        f"Integer phase condition satisfied all scales? `{decision['integer_phase_condition_satisfied_all_scales']}`.",
        f"Unique scale selected by integer phase alone? `{decision['unique_scale_selected_by_integer_phase_alone']}`.",
        f"Fixed-clock selector admissible now? `{decision['fixed_clock_selector_admissible_now']}`.",
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
    family = action_quantization_family_audit()
    fixed_clock = fixed_clock_false_selector_audit()
    matrix = source_candidate_matrix(family, fixed_clock)
    decision = closure_decision(matrix, family, fixed_clock)
    payload: dict[str, Any] = {
        "status": "P2657_ACTION_QUANTIZATION_SCALE_ANCHOR_OBSTRUCTION_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2655_TYPED_METRIC_SCALE_QUOTIENT": sha256_file(P2655),
            "P2656_LAPLACIAN_ACTION_SOURCE_AUDIT": sha256_file(P2656),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_action_model": {
            "nodes": BASE_POINTS,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "phase_condition": "tau * Tr(L) = 2*pi*n",
            "no_sub_nadsoliton_information_layer": True,
        },
        "action_quantization_family_audit": family,
        "fixed_clock_false_selector_audit": fixed_clock,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2657/S1607 action quantization scale-anchor obstruction guard",
        "## P2657/S1607 action quantization scale-anchor obstruction guard\n\n"
        "`P2657/S1607` audits the next theorem target after P2656: an integer action/phase condition `tau*Tr(L)=2*pi*n`.  The condition is satisfiable for every audited scale by rescaling the clock as `tau -> a^2 tau`, so it quantizes a scale-clock product rather than the UV length.  A unique representative appears only if a fixed clock, trace, spectral gap, heat-time, or action quantum is imported as an external anchor.  Therefore this exports no intrinsic action-quantization theorem, no typed metric/UV source theorem, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2657/S1607 action quantization Ltotal guard",
        "## P2657/S1607 action quantization Ltotal guard\n\n"
        "`P2657/S1607` does not re-open `L_total`: the integer phase/action condition leaves a scale-clock degeneracy unless the nadsoliton evolution law derives an intrinsic clock or action quantum.  A variational damping coefficient still requires that intrinsic scale source before beta-source rerun, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
