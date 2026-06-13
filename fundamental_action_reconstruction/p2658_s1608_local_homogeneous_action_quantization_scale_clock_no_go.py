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
OUT = GEN / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.json"
MD = GEN / "p2658_s1608_local_homogeneous_action_quantization_scale_clock_no_go.md"

P2656 = GEN / "p2656_s1606_typed_nadsoliton_laplacian_action_identity_scale_source_audit.json"
P2657 = GEN / "p2657_s1607_nadsoliton_action_quantization_scale_anchor_obstruction_audit.json"
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
HOMOGENEITY_EXPONENTS = [0.5, 1.0, 1.8, 2.0, 3.0]
ACTION_POWERS = [1.0, 2.0]
ACTION_QUANTA = [1, 2, 3, 5]
TOL = 1e-11

NEGATIVE_EXPORT_FLAGS = [
    "local_homogeneous_action_quantization_no_go_exported_as_full_global_theorem",
    "intrinsic_clock_or_flow_time_exported",
    "uv_unit_selected_by_local_quantization",
    "typed_metric_uv_source_theorem_exported",
    "target_independent_beta_source_exported",
    "canonical_unit_exported",
    "fixed_clock_anchor_promoted_to_source",
    "homogeneous_trace_anchor_promoted_to_source",
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
        "homogeneous_local_operator_content": (
            "homogeneous|local.*operator|Laplacian|trace moment|weighted.*graph|operator covariance"
        ),
        "scale_clock_no_go_content": (
            "scale-clock|clock-scale|flow time|fixed clock|tau -> a\\^|scale orbit|normalization anchor"
        ),
        "action_quantization_content": (
            "action quantization|integer phase|phase condition|action quantum|tau\\*Tr|conserved.*integer"
        ),
        "nonclosure_guard_content": (
            "role-bearing L_total|QW-2191|role-transfer|bridge completion|ToE closure|beta source|UV unit"
        ),
    }
    return {
        "tool": "rg",
        "mode": "content-first semantic audit for local homogeneous action-quantization no-go, not packet-name search",
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


def homogeneous_trace(scale: float, exponent: float) -> float:
    total = 0.0
    for a, b in itertools.combinations(BASE_POINTS, 2):
        d = scale * base_distance(a, b)
        total += 2.0 / (d**exponent)
    return total


def action_functional(scale: float, exponent: float, power: float) -> float:
    return homogeneous_trace(scale, exponent) ** power


def homogeneous_quantization_no_go() -> dict[str, Any]:
    rows = []
    max_covariance_error = 0.0
    max_tau_ratio_error = 0.0
    max_phase_error = 0.0
    for exponent in HOMOGENEITY_EXPONENTS:
        for power in ACTION_POWERS:
            degree = exponent * power
            reference_action = action_functional(1.0, exponent, power)
            for scale in AUDITED_SCALES:
                current_action = action_functional(scale, exponent, power)
                covariance_error = abs(current_action * (scale**degree) - reference_action)
                max_covariance_error = max(max_covariance_error, covariance_error)
                for quantum in ACTION_QUANTA:
                    tau = (2.0 * math.pi * quantum) / current_action
                    reference_tau = (2.0 * math.pi * quantum) / reference_action
                    phase = tau * current_action
                    tau_ratio_error = abs((tau / reference_tau) - (scale**degree))
                    phase_error = abs(phase - 2.0 * math.pi * quantum)
                    max_tau_ratio_error = max(max_tau_ratio_error, tau_ratio_error)
                    max_phase_error = max(max_phase_error, phase_error)
                    rows.append({
                        "homogeneity_exponent_p": exponent,
                        "action_power_q": power,
                        "total_degree_pq": degree,
                        "scale": scale,
                        "integer_quantum_n": quantum,
                        "action_value": current_action,
                        "action_times_scale_degree": current_action * (scale**degree),
                        "quantized_tau": tau,
                        "tau_ratio_to_scale_one_same_n": tau / reference_tau,
                        "expected_tau_ratio_a_degree": scale**degree,
                        "phase_error_from_2pi_n": phase_error,
                    })
    return {
        "statement": "For audited local homogeneous action functionals A_pq(a)=Tr_p(L_a)^q with total degree m=p*q, A_pq(a)=A_pq(1)/a^m.  The integer condition tau*A_pq(a)=2*pi*n is therefore satisfiable at every scale by tau(a)=a^m*tau(1), so it cannot select a unique UV length without an independent clock/action anchor.",
        "rows": rows,
        "max_covariance_error": max_covariance_error,
        "max_tau_ratio_error_from_expected_degree": max_tau_ratio_error,
        "max_phase_error": max_phase_error,
        "all_integer_phase_conditions_satisfied": max_phase_error < TOL,
        "all_homogeneous_covariances_verified": max_covariance_error < TOL,
        "audited_functional_count": len(HOMOGENEITY_EXPONENTS) * len(ACTION_POWERS),
        "unique_scale_selected_by_local_homogeneous_quantization": False,
    }


def fixed_clock_selector_generalization() -> dict[str, Any]:
    rows = []
    false_selector_count = 0
    for exponent in HOMOGENEITY_EXPONENTS:
        for power in ACTION_POWERS:
            reference_action = action_functional(1.0, exponent, power)
            tau_fixed = 2.0 * math.pi / reference_action
            passing_scales = []
            for scale in AUDITED_SCALES:
                phase_over_2pi = tau_fixed * action_functional(scale, exponent, power) / (2.0 * math.pi)
                looks_n1 = abs(phase_over_2pi - 1.0) < TOL
                if looks_n1:
                    passing_scales.append(scale)
                rows.append({
                    "homogeneity_exponent_p": exponent,
                    "action_power_q": power,
                    "fixed_tau_from_scale_one_n1": tau_fixed,
                    "scale": scale,
                    "phase_over_2pi_with_fixed_tau": phase_over_2pi,
                    "looks_integer_n1": looks_n1,
                })
            if passing_scales == [1.0]:
                false_selector_count += 1
    return {
        "statement": "For every audited homogeneous functional, a fixed clock chosen at the scale-one representative makes scale=1 look unique, but the fixed clock is exactly the imported scale anchor.",
        "rows": rows,
        "false_selector_count": false_selector_count,
        "all_audited_fixed_clock_selectors_are_external": false_selector_count == len(HOMOGENEITY_EXPONENTS) * len(ACTION_POWERS),
        "fixed_clock_selector_admissible_now": False,
    }


def source_candidate_matrix(no_go: dict[str, Any], fixed_clock: dict[str, Any]) -> list[dict[str, Any]]:
    return [
        {
            "candidate": "local_homogeneous_integer_phase_family",
            "covered_by_audit": True,
            "uses_external_clock_or_scale_anchor": False,
            "selects_unique_scale": no_go["unique_scale_selected_by_local_homogeneous_quantization"],
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: integer phase conditions are solved by compensating clock rescaling for every audited homogeneous degree",
        },
        {
            "candidate": "fixed_clock_for_homogeneous_functional",
            "covered_by_audit": True,
            "uses_external_clock_or_scale_anchor": True,
            "selects_unique_scale": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked: uniqueness comes from importing the scale-one clock",
        },
        {
            "candidate": "declared_absolute_homogeneous_trace_quantum",
            "covered_by_audit": True,
            "uses_external_clock_or_scale_anchor": True,
            "selects_unique_scale": True,
            "passes_as_uv_unit_source_now": False,
            "verdict": "blocked unless the absolute trace/action quantum is derived internally rather than declared",
        },
        {
            "candidate": "nonhomogeneous_or_boundary_anomaly_source",
            "covered_by_audit": False,
            "uses_external_clock_or_scale_anchor": False,
            "selects_unique_scale": None,
            "passes_as_uv_unit_source_now": False,
            "verdict": "not ruled out here; it is the next admissible target if genuinely derived from nadsoliton dynamics",
        },
        {
            "candidate": "intrinsic_nadsoliton_clock_source_theorem",
            "covered_by_audit": False,
            "uses_external_clock_or_scale_anchor": False,
            "selects_unique_scale": None,
            "passes_as_uv_unit_source_now": False,
            "verdict": "still required for a positive source theorem; not exported by this audit",
        },
    ]


def upstream_consistency() -> dict[str, Any]:
    p2656 = load_json(P2656)
    p2657 = load_json(P2657)
    return {
        "p2656_operator_covariance_verified": p2656.get("closure_decision", {}).get("operator_covariance_verified") is True,
        "p2656_absolute_operator_scale_not_selected": p2656.get("closure_decision", {}).get("absolute_operator_scale_selected") is False,
        "p2657_integer_phase_condition_satisfied_all_scales": p2657.get("closure_decision", {}).get("integer_phase_condition_satisfied_all_scales") is True,
        "p2657_unique_scale_not_selected": p2657.get("closure_decision", {}).get("unique_scale_selected_by_integer_phase_alone") is False,
        "p2657_next_step_allows_broader_no_go": "broader no-go" in p2657.get("closure_decision", {}).get("next_honest_step", ""),
    }


def closure_decision(no_go: dict[str, Any], fixed_clock: dict[str, Any], matrix: list[dict[str, Any]]) -> dict[str, Any]:
    passing = [row["candidate"] for row in matrix if row["passes_as_uv_unit_source_now"]]
    return {
        "decision": "LOCAL_HOMOGENEOUS_ACTION_QUANTIZATION_SCALE_CLOCK_NO_GO__NO_INTRINSIC_UV_UNIT",
        "professorial_verdict": (
            "P2658 generalizes the P2657 obstruction from the single Laplacian trace to an audited class of local homogeneous action functionals.  "
            "Every audited integer phase condition remains satisfiable across the scale orbit by compensating the clock with the matching homogeneity degree.  "
            "Thus local homogeneous action quantization does not source the UV unit; any apparent unique representative comes from a fixed clock, absolute trace, or declared action quantum.  "
            "This is a strong finite no-go certificate, not a full global theorem over all possible nonhomogeneous/anomalous nadsoliton dynamics."
        ),
        "professorial_closure_path": [
            "Treat local homogeneous Laplacian/action quantizations as blocked for UV-unit sourcing unless a nonhomogeneous anomaly or intrinsic clock source is explicitly derived.",
            "Do not promote fixed clocks, absolute trace quanta, spectral gaps, or heat times to source status; P2658 classifies them as external anchors.",
            "The next honest proof-grade move is either a derived nonhomogeneous/anomalous term that breaks clock-scale covariance, or a formal theorem proving that the admissible local dynamics class is necessarily homogeneous and therefore source-blocked.",
            "Only after such a source theorem should beta=1 be rerun through P2649; empirical holdout packets remain support only.",
        ],
        "next_honest_step": (
            "Attempt a nonhomogeneous/anomalous nadsoliton clock-source theorem: specify an internally derived boundary, anomaly, or discrete phase term whose scaling degree is not cancelled by tau -> a^m tau; if none can be derived, formalize the stronger theorem that all admissible local typed action functionals are homogeneous and hence cannot source the UV unit."
        ),
        "passing_uv_unit_source_candidates": passing,
        "audited_functional_count": no_go["audited_functional_count"],
        "all_homogeneous_covariances_verified": no_go["all_homogeneous_covariances_verified"],
        "all_integer_phase_conditions_satisfied": no_go["all_integer_phase_conditions_satisfied"],
        "fixed_clock_selector_admissible_now": fixed_clock["fixed_clock_selector_admissible_now"],
        "uv_unit_selected_now": False,
        "beta_source_exported_now": False,
        "role_bearing_ltotal_now": False,
        "toe_closure_now": False,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    no_go = payload["homogeneous_quantization_no_go"]
    fixed_clock = payload["fixed_clock_selector_generalization"]
    decision = payload["closure_decision"]
    lines = [
        "# P2658/S1608 local homogeneous action quantization scale-clock no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Content-first anti-duplication audit",
        "",
        "This packet greps homogeneous local operator, scale-clock no-go, action-quantization, and nonclosure guard content before adding the finite class no-go.",
        "",
    ]
    for name, data in payload["semantic_rg_antiduplication_audit"]["patterns"].items():
        lines.append(f"- `{name}`: {data['count']} hits")
    lines.extend([
        "",
        "## Homogeneous quantization no-go",
        "",
        no_go["statement"],
        f"Audited homogeneous functional count: `{no_go['audited_functional_count']}`.",
        f"All homogeneous covariances verified? `{no_go['all_homogeneous_covariances_verified']}`.",
        f"All integer phase conditions satisfied? `{no_go['all_integer_phase_conditions_satisfied']}`.",
        f"Max covariance error: `{no_go['max_covariance_error']:.3e}`.",
        f"Max tau-ratio error: `{no_go['max_tau_ratio_error_from_expected_degree']:.3e}`.",
        f"Unique scale selected by local homogeneous quantization? `{no_go['unique_scale_selected_by_local_homogeneous_quantization']}`.",
        "",
        "## Fixed-clock generalization",
        "",
        fixed_clock["statement"],
        f"All audited fixed-clock selectors are external? `{fixed_clock['all_audited_fixed_clock_selectors_are_external']}`.",
        f"Fixed-clock selector admissible now? `{fixed_clock['fixed_clock_selector_admissible_now']}`.",
        "",
        "## Source candidate matrix",
        "",
        "| candidate | covered? | external anchor? | selects unique scale? | source now? | verdict |",
        "| --- | ---: | ---: | ---: | ---: | --- |",
    ])
    for row in payload["source_candidate_matrix"]:
        lines.append(
            f"| `{row['candidate']}` | `{row['covered_by_audit']}` | `{row['uses_external_clock_or_scale_anchor']}` | "
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
    no_go = homogeneous_quantization_no_go()
    fixed_clock = fixed_clock_selector_generalization()
    matrix = source_candidate_matrix(no_go, fixed_clock)
    decision = closure_decision(no_go, fixed_clock, matrix)
    payload: dict[str, Any] = {
        "status": "P2658_LOCAL_HOMOGENEOUS_ACTION_QUANTIZATION_SCALE_CLOCK_NO_GO_NO_FALSE_PASS",
        "latest_commit_audit": latest_commit_audit(),
        "semantic_rg_antiduplication_audit": semantic_rg_audit(),
        "source_hashes": {
            "P2656_LAPLACIAN_ACTION_SOURCE_AUDIT": sha256_file(P2656),
            "P2657_ACTION_QUANTIZATION_OBSTRUCTION": sha256_file(P2657),
            "STRICT_EQUATION_SHEET": sha256_file(STRICT_EQUATION_SHEET),
            "STRICT_LAGRANGIAN_DRAFT": sha256_file(STRICT_LAGRANGIAN_DRAFT),
        },
        "upstream_consistency": upstream_consistency(),
        "typed_local_action_model": {
            "nodes": BASE_POINTS,
            "ontology_order": "nadsoliton -> light -> matter -> emergent observer",
            "audited_functional_form": "A_pq(a)=Tr_p(L_a)^q with local weights d^-p and integer phase tau*A=2*pi*n",
            "no_sub_nadsoliton_information_layer": True,
        },
        "homogeneous_quantization_no_go": no_go,
        "fixed_clock_selector_generalization": fixed_clock,
        "source_candidate_matrix": matrix,
        "closure_decision": decision,
        "negative_export_flags": {key: False for key in NEGATIVE_EXPORT_FLAGS},
    }
    payload["payload_sha256"] = sha256_json({k: v for k, v in payload.items() if k != "payload_sha256"})
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(
        STRICT_EQUATION_SHEET,
        "P2658/S1608 local homogeneous action quantization no-go guard",
        "## P2658/S1608 local homogeneous action quantization no-go guard\n\n"
        "`P2658/S1608` generalizes the P2657 scale-clock obstruction from `Tr(L)` to an audited class of local homogeneous action functionals `A_pq(a)=Tr_p(L_a)^q`.  Each functional scales as `A_pq(a)=A_pq(1)/a^(p*q)`, so the integer phase condition `tau*A_pq=2*pi*n` is satisfiable at every scale by `tau -> a^(p*q) tau`.  Fixed-clock uniqueness is therefore an imported clock/scale anchor.  This is a finite class no-go, not a full global theorem over possible nonhomogeneous/anomalous nadsoliton dynamics, and it exports no intrinsic UV unit, no canonical unit, no beta source, no bridge completion, no role transfer, no `QW-2191` discharge, no role-bearing `L_total`, and no ToE closure.\n",
    )
    append_once(
        STRICT_LAGRANGIAN_DRAFT,
        "P2658/S1608 local homogeneous action quantization Ltotal guard",
        "## P2658/S1608 local homogeneous action quantization Ltotal guard\n\n"
        "`P2658/S1608` does not re-open `L_total`: local homogeneous action quantizations leave a scale-clock orbit unbroken, and fixed clocks or absolute trace quanta remain external anchors.  A variational damping coefficient still requires a genuinely derived nonhomogeneous/anomalous clock-source term or an intrinsic nadsoliton clock/action theorem before beta-source rerun, role-transfer rerun, and selector/source discharge.\n",
    )
    return payload


if __name__ == "__main__":
    result = main()
    print(rel(OUT))
    print(rel(MD))
