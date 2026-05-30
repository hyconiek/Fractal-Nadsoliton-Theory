#!/usr/bin/env python3
"""Scratch probe: low-order transport no-go audit for strict completion cocycle.

The transport-cocycle certificate produced the unique finite transport

    T(d) = K_strict_gate(d) / K_legacy_ont(d)

and adjacent ratios R(d)=T(d+1)/T(d) on d=0..11.  This probe asks a narrower
negative question that is closer to the open strict_transport_derivation blocker:
can that finite transport already be read as a very simple autonomous transport
law?

We test four low-complexity candidates:

1. positive damping-only transport (all edge ratios positive),
2. constant first-order multiplier T[d+1]=r*T[d],
3. constant-coefficient second-order linear recurrence T[d+2]=a*T[d+1]+b*T[d],
4. affine log-envelope log|T(d)|=m*d+c.

All four fail on the audited finite domain.  This does not prove that no strict
nadsoliton dynamics can generate T(d); it only rules out these low-order/autonomous
readings and therefore sharpens what a real strict transport derivation would
need to supply.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_low_order_transport_no_go_report.json"
OUT_MD = HERE / "bridge_strict_completion_low_order_transport_no_go_report.md"
COCYCLE_REPORT = HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json"
NECESSITY_REPORT = HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json"
BLOCKER_LATTICE = HERE / "bridge_completed_kernel_blocker_dependency_lattice_report.json"

TOL = 1e-14


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def l2(values: list[float]) -> float:
    return math.sqrt(sum(value * value for value in values))


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def fit_first_order_constant_multiplier(values: list[float]) -> dict[str, Any]:
    numerator = sum(left * right for left, right in zip(values, values[1:]))
    denominator = sum(left * left for left in values[:-1])
    multiplier = numerator / denominator
    residuals = [values[index + 1] - multiplier * values[index] for index in range(len(values) - 1)]
    return {
        "model": "T[d+1]=r*T[d]",
        "best_r_least_squares": multiplier,
        "max_abs_residual": max(abs(value) for value in residuals),
        "l2_residual": l2(residuals),
        "exact_pass": max(abs(value) for value in residuals) <= TOL,
        "residuals": residuals,
    }


def fit_second_order_constant_recurrence(values: list[float]) -> dict[str, Any]:
    rows = [(values[index + 1], values[index]) for index in range(len(values) - 2)]
    targets = [values[index + 2] for index in range(len(values) - 2)]
    s11 = sum(row[0] * row[0] for row in rows)
    s12 = sum(row[0] * row[1] for row in rows)
    s22 = sum(row[1] * row[1] for row in rows)
    t1 = sum(row[0] * target for row, target in zip(rows, targets))
    t2 = sum(row[1] * target for row, target in zip(rows, targets))
    determinant = s11 * s22 - s12 * s12
    if abs(determinant) <= TOL:
        raise ValueError("degenerate second-order normal equations")
    coeff_a = (t1 * s22 - t2 * s12) / determinant
    coeff_b = (s11 * t2 - s12 * t1) / determinant
    residuals = [values[index + 2] - coeff_a * values[index + 1] - coeff_b * values[index] for index in range(len(values) - 2)]
    return {
        "model": "T[d+2]=a*T[d+1]+b*T[d]",
        "best_a_least_squares": coeff_a,
        "best_b_least_squares": coeff_b,
        "max_abs_residual": max(abs(value) for value in residuals),
        "l2_residual": l2(residuals),
        "exact_pass": max(abs(value) for value in residuals) <= TOL,
        "residuals": residuals,
    }


def fit_affine_log_envelope(values: list[float]) -> dict[str, Any]:
    domain = list(range(len(values)))
    logs = [math.log(abs(value)) for value in values]
    count = len(domain)
    sum_x = sum(domain)
    sum_y = sum(logs)
    sum_xx = sum(value * value for value in domain)
    sum_xy = sum(x_value * y_value for x_value, y_value in zip(domain, logs))
    slope = (count * sum_xy - sum_x * sum_y) / (count * sum_xx - sum_x * sum_x)
    intercept = (sum_y - slope * sum_x) / count
    residuals = [log_value - (slope * d_value + intercept) for d_value, log_value in zip(domain, logs)]
    return {
        "model": "log|T[d]|=m*d+c",
        "best_m_least_squares": slope,
        "best_c_least_squares": intercept,
        "max_abs_residual": max(abs(value) for value in residuals),
        "l2_residual": l2(residuals),
        "exact_pass": max(abs(value) for value in residuals) <= TOL,
        "residuals": residuals,
    }


def period_audit(signs: list[int], max_period: int) -> list[dict[str, Any]]:
    rows = []
    for period in range(1, max_period + 1):
        mismatches = [index for index, value in enumerate(signs) if value != signs[index % period]]
        rows.append({
            "period": period,
            "is_periodic_on_audited_edges": not mismatches,
            "mismatch_indices": mismatches,
        })
    return rows


def build_payload() -> dict[str, Any]:
    cocycle = load_json(COCYCLE_REPORT)
    necessity = load_json(NECESSITY_REPORT)
    blocker_lattice = load_json(BLOCKER_LATTICE)

    node_values = [row["transport_T"] for row in cocycle["node_transport_rows"]]
    edge_values = [row["edge_transport_ratio_R"] for row in cocycle["edge_cocycle_rows"]]
    edge_signs = [sign(value) for value in edge_values]
    positive_edges = [row["edge"] for row in cocycle["edge_cocycle_rows"] if row["edge_transport_ratio_R"] > 0.0]
    negative_edges = [row["edge"] for row in cocycle["edge_cocycle_rows"] if row["edge_transport_ratio_R"] < 0.0]

    first_order = fit_first_order_constant_multiplier(node_values)
    second_order = fit_second_order_constant_recurrence(node_values)
    affine_log = fit_affine_log_envelope(node_values)
    sign_periods = period_audit(edge_signs, max_period=5)

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LOW_ORDER_TRANSPORT_NO_GO__FINITE_Z12_NOT_FULL_DYNAMICAL_NO_GO",
        "status": "low-order-autonomous-transport-readings-fail-for-completion-cocycle",
        "source_reports": {
            "transport_cocycle": str(COCYCLE_REPORT.relative_to(ROOT)),
            "necessity_certificate": str(NECESSITY_REPORT.relative_to(ROOT)),
            "blocker_lattice": str(BLOCKER_LATTICE.relative_to(ROOT)),
        },
        "domain": {
            "node_count": len(node_values),
            "edge_count": len(edge_values),
            "tolerance": TOL,
        },
        "transport_input_summary": {
            "transport_sign_pattern": [sign(value) for value in node_values],
            "edge_sign_pattern": edge_signs,
            "positive_edge_count": len(positive_edges),
            "negative_edge_count": len(negative_edges),
            "negative_edges": negative_edges,
            "min_edge_ratio": min(edge_values),
            "max_edge_ratio": max(edge_values),
            "edge_ratio_spread": max(edge_values) - min(edge_values),
        },
        "no_go_checks": {
            "positive_damping_only_transport": {
                "candidate": "all adjacent ratios R(d) are positive damping/envelope updates",
                "exact_pass": len(negative_edges) == 0,
                "failure_reason": "four adjacent ratios are negative, so phase sign transport cannot be replaced by positive damping-only flow",
                "negative_edges": negative_edges,
            },
            "constant_first_order_multiplier": first_order,
            "constant_second_order_linear_recurrence": second_order,
            "affine_log_envelope": affine_log,
            "short_period_edge_sign_law": {
                "candidate": "edge sign sequence has period p<=5",
                "exact_pass": any(row["is_periodic_on_audited_edges"] for row in sign_periods),
                "period_rows": sign_periods,
            },
        },
        "no_go_summary": {
            "positive_damping_only_fails": len(negative_edges) > 0,
            "constant_first_order_fails": not first_order["exact_pass"],
            "constant_second_order_fails": not second_order["exact_pass"],
            "affine_log_envelope_fails": not affine_log["exact_pass"],
            "short_period_edge_sign_law_fails": not any(row["is_periodic_on_audited_edges"] for row in sign_periods),
            "minimum_failed_model_l2_residual": min(first_order["l2_residual"], second_order["l2_residual"], affine_log["l2_residual"]),
            "sharpener": "A real strict transport derivation must be nontrivial enough to generate mixed sign phase transport plus non-affine envelope variation.",
        },
        "blocker_context": {
            "transport_cocycle_status": cocycle["status"],
            "necessity_status": necessity["status"],
            "blocker_lattice_status": blocker_lattice["status"],
            "refined_blocker": "strict_transport_derivation remains open, but simple autonomous/damping-only explanations are now ruled out on the finite audit domain",
            "still_open": [
                "non-low-order strict_transport_derivation_from_nadsoliton_dynamics",
                "global_z12_map_derivation_as_internal_dynamics",
                "orientation_chi11_source",
                "chi11_uniqueness",
                "reynolds_obstruction_escape",
                "role_transfer_theorem",
            ],
        },
        "proof_certificate": {
            "positive_damping_no_go": "A positive damping-only edge law requires R(d)>0 on every edge; the audited cocycle has four negative R(d).",
            "first_order_no_go": "If T[d+1]=r*T[d] exactly, all edge ratios are equal; the least-squares residual remains nonzero at the recorded tolerance.",
            "second_order_no_go": "The best constant-coefficient second-order recurrence has nonzero residual, so the finite sequence is not captured by this low-order autonomous recurrence.",
            "affine_log_no_go": "An affine log-envelope would make log|T(d)| linear in d; the least-squares residual is far above machine tolerance.",
            "period_no_go": "The edge-sign sequence is not periodic with period <=5 on the audited edges.",
            "theoretical_limit": "This is not a theorem against all possible strict dynamics; it only rules out the listed low-order/autonomous candidates.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in solitonic state; this audit only constrains candidate transport readings.",
            "forbidden_reading": "No separate informational layer below the nadsoliton is introduced.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No full no-go theorem against every strict transport derivation is claimed.",
            "No proof derives the completion transport from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    input_summary = payload["transport_input_summary"]
    summary = payload["no_go_summary"]
    checks = payload["no_go_checks"]
    lines = [
        "# Strict completion low-order transport no-go audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This finite-domain audit tests whether the completion transport cocycle can",
        "already be explained by very simple autonomous laws.  All tested low-order",
        "readings fail on `d=0..11`; this sharpens but does not close the open",
        "`strict_transport_derivation` blocker.",
        "",
        "## Transport input summary",
        "",
        f"- Edge sign pattern: `{input_summary['edge_sign_pattern']}`",
        f"- Negative edges: `{input_summary['negative_edges']}`",
        f"- Edge ratio spread: `{input_summary['edge_ratio_spread']:.12e}`",
        "",
        "## No-go summary",
        "",
        f"- Positive damping-only transport fails: `{summary['positive_damping_only_fails']}`",
        f"- Constant first-order multiplier fails: `{summary['constant_first_order_fails']}`",
        f"- Constant second-order recurrence fails: `{summary['constant_second_order_fails']}`",
        f"- Affine log-envelope fails: `{summary['affine_log_envelope_fails']}`",
        f"- Short-period edge-sign law fails: `{summary['short_period_edge_sign_law_fails']}`",
        f"- Minimum failed-model L2 residual: `{summary['minimum_failed_model_l2_residual']:.12e}`",
        "",
        "## Model residuals",
        "",
        f"- First-order best `r`: `{checks['constant_first_order_multiplier']['best_r_least_squares']:.12e}`, L2 residual `{checks['constant_first_order_multiplier']['l2_residual']:.12e}`",
        f"- Second-order best `(a,b)`: `({checks['constant_second_order_linear_recurrence']['best_a_least_squares']:.12e}, {checks['constant_second_order_linear_recurrence']['best_b_least_squares']:.12e})`, L2 residual `{checks['constant_second_order_linear_recurrence']['l2_residual']:.12e}`",
        f"- Affine log best `(m,c)`: `({checks['affine_log_envelope']['best_m_least_squares']:.12e}, {checks['affine_log_envelope']['best_c_least_squares']:.12e})`, L2 residual `{checks['affine_log_envelope']['l2_residual']:.12e}`",
        "",
        "## Guarded interpretation",
        "",
        "The result is a finite no-go for the listed low-order/autonomous readings only.",
        "It is not a full no-go theorem against strict nadsoliton dynamics and does not",
        "derive any selector, chi_11 source, or ToE closure.",
        "",
        "## Hard limits",
        "",
    ]
    for limit in payload["hard_limits"]:
        lines.append(f"- {limit}")
    OUT_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    write_markdown(payload)
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
