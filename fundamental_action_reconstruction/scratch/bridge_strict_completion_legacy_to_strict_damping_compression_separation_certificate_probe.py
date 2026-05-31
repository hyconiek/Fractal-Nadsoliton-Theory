#!/usr/bin/env python3
"""Scratch probe: legacy-linear vs strict-nonlinear damping separation.

The component-gap matrix says the damping/compression bridge row is certified as
finite bookkeeping but still lacks a strict-source theorem.  This probe makes the
local obstruction sharper: on the positive Z12 tail, no single legacy-style
linear torsion denominator 1+gamma*d can equal the strict nonlinear denominator
1+d^eta at two distinct positive nodes, hence the strict d^eta compression is a
real strict-side addition rather than a hidden beta_tors re-labeling.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.md"

BETA_TORS = 0.01
STRICT_BETA = 1.0
ETA = 9.0 / 5.0
ETA_MINUS_ONE = ETA - 1.0
DOMAIN = list(range(1, 12))

SOURCE_REPORTS = {
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "damping_exact": HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
    "compression_ontology_audit": HERE / "bridge_compression_ontology_audit_report.json",
    "t15_bridge_theorem_spec": ROOT / "fundamental_action_reconstruction" / "T15_LEGACY_TO_STRICT_KERNEL_BRIDGE_THEOREM_SPEC.md",
}

SEARCH_TERMS = [
    "linear torsion denominator no-go",
    "1+gamma*d vs 1+d^eta",
    "beta_tors*d to beta*d^eta",
    "legacy linear torsion damping strict nonlinear compression",
    "no single linear gamma matches strict eta damping",
]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite text: {path}")
    return path.read_text(encoding="utf-8")


def report_path(path: Path) -> str:
    return str(path.relative_to(ROOT))


def strict_denominator(d: int) -> float:
    return 1.0 + STRICT_BETA * d**ETA


def legacy_linear_denominator(d: int, gamma: float = BETA_TORS) -> float:
    return 1.0 + gamma * d


def required_linear_gamma_for_node(d: int) -> float:
    # If 1+gamma*d = 1+d^eta and d>0, then gamma=d^(eta-1).
    return STRICT_BETA * d**ETA_MINUS_ONE


def best_l2_linear_gamma() -> float:
    # Minimize sum_d (gamma*d - d^eta)^2 over DOMAIN.
    numerator = sum(d * d**ETA for d in DOMAIN)
    denominator = sum(d * d for d in DOMAIN)
    return numerator / denominator


def build_rows(gamma_l2: float) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for d in DOMAIN:
        strict_value = strict_denominator(d)
        legacy_value = legacy_linear_denominator(d)
        l2_value = legacy_linear_denominator(d, gamma_l2)
        required_gamma = required_linear_gamma_for_node(d)
        rows.append({
            "d": d,
            "strict_denominator_1_plus_d_eta": strict_value,
            "legacy_beta_tors_denominator_1_plus_beta_tors_d": legacy_value,
            "best_l2_linear_denominator_1_plus_gamma_star_d": l2_value,
            "required_gamma_for_exact_single_node_match": required_gamma,
            "required_gamma_minus_legacy_beta_tors": required_gamma - BETA_TORS,
            "strict_minus_legacy_denominator_residual": strict_value - legacy_value,
            "strict_minus_best_l2_linear_residual": strict_value - l2_value,
            "legacy_beta_tors_matches_this_node": math.isclose(required_gamma, BETA_TORS, rel_tol=0.0, abs_tol=1e-15),
        })
    return rows


def build_payload() -> dict[str, Any]:
    loaded: dict[str, Any] = {}
    source_reports: dict[str, str] = {}
    for name, path in SOURCE_REPORTS.items():
        source_reports[name] = report_path(path)
        loaded[name] = load_json(path) if path.suffix == ".json" else {"text": read_text(path)}

    gamma_l2 = best_l2_linear_gamma()
    rows = build_rows(gamma_l2)
    required_gammas = [row["required_gamma_for_exact_single_node_match"] for row in rows]
    l2_residuals = [row["strict_minus_best_l2_linear_residual"] for row in rows]
    legacy_residuals = [row["strict_minus_legacy_denominator_residual"] for row in rows]
    gamma_diffs = [b - a for a, b in zip(required_gammas, required_gammas[1:])]

    component_gap = loaded["component_gap_matrix"]["completion_gap_summary"]
    guardrail = loaded["legacy_bridge_guardrail"]["legacy_kernel_intermediate_bridge_summary"]
    necessity = loaded["necessity"]["necessity_summary"]
    damping_exact = loaded["damping_exact"]["exact_rational_damping_summary"]
    compression_audit = loaded["compression_ontology_audit"]["compression_characteristic"]
    t15_text = loaded["t15_bridge_theorem_spec"]["text"]

    separation_summary = {
        "domain": DOMAIN,
        "eta": ETA,
        "eta_minus_one": ETA_MINUS_ONE,
        "legacy_beta_tors": BETA_TORS,
        "strict_beta": STRICT_BETA,
        "required_gamma_values_strictly_increase": all(diff > 0.0 for diff in gamma_diffs),
        "no_single_linear_gamma_matches_two_distinct_positive_nodes": all(diff > 0.0 for diff in gamma_diffs),
        "legacy_beta_tors_matches_no_positive_strict_node": not any(row["legacy_beta_tors_matches_this_node"] for row in rows),
        "best_l2_linear_gamma": gamma_l2,
        "best_l2_residual_l2_norm": math.sqrt(sum(res * res for res in l2_residuals)),
        "best_l2_residual_max_abs": max(abs(res) for res in l2_residuals),
        "legacy_beta_tors_residual_min": min(legacy_residuals),
        "legacy_beta_tors_residual_max": max(legacy_residuals),
        "required_gamma_at_d1": required_gammas[0],
        "required_gamma_at_d11": required_gammas[-1],
        "required_gamma_spread_d11_minus_d1": required_gammas[-1] - required_gammas[0],
        "d11_strict_over_legacy_denominator_ratio": rows[-1]["strict_denominator_1_plus_d_eta"] / rows[-1]["legacy_beta_tors_denominator_1_plus_beta_tors_d"],
        "component_gap_records_compression_missing": component_gap["strict_compression_recorded_as_missing_from_legacy"],
        "guardrail_records_legacy_incomplete": guardrail["legacy_kernel_incomplete_for_strict_characteristics"],
        "necessity_marks_damping_shape_critical": necessity["damping_factor_positive"] and necessity["damping_factor_strictly_decreasing_after_d0"],
        "exact_damping_monotone_inherited": damping_exact["continuous_strictly_decreasing_certificate"],
        "compression_ontology_identifies_missing_characteristic": "missing explicit characteristic" in compression_audit["missing_in_legacy_reading"]["description"],
        "t15_names_damping_compression_obligation": "linear torsion damping" in t15_text and ("nonlinear" in t15_text or "non-linear" in t15_text),
        "full_bridge_theorem_exported": False,
    }

    cross_checks = {
        "source_gap_and_guardrail_loaded": separation_summary["component_gap_records_compression_missing"] and separation_summary["guardrail_records_legacy_incomplete"],
        "required_gamma_strictly_increases": separation_summary["required_gamma_values_strictly_increase"],
        "no_single_linear_gamma_matches_two_nodes": separation_summary["no_single_linear_gamma_matches_two_distinct_positive_nodes"],
        "legacy_beta_tors_matches_no_positive_node": separation_summary["legacy_beta_tors_matches_no_positive_strict_node"],
        "best_linear_fit_still_has_residual": separation_summary["best_l2_residual_l2_norm"] > 0.0 and separation_summary["best_l2_residual_max_abs"] > 0.0,
        "damping_shape_critical_and_monotone_inherited": separation_summary["necessity_marks_damping_shape_critical"] and separation_summary["exact_damping_monotone_inherited"],
        "bridge_obligation_named_but_not_closed": separation_summary["t15_names_damping_compression_obligation"] and not separation_summary["full_bridge_theorem_exported"],
    }

    payload = {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_DAMPING_COMPRESSION_SEPARATION__LINEAR_TORSION_NO_GO",
        "status": "legacy-linear-torsion-denominator-separated-from-strict-nonlinear-d-eta-compression-no-bridge-theorem",
        "source_reports": source_reports,
        "grep_disambiguation": {
            "searched_terms": SEARCH_TERMS,
            "finding": "Existing reports certify damping necessity, monotonicity, ontology, inverse branches, and the component gap matrix; this report adds the exact finite/algebraic separation: no single linear gamma in a legacy-style denominator can match strict d^eta compression on two positive Z12 nodes.",
        },
        "model_definition": {
            "legacy_linear_family": "L_gamma(d)=1+gamma*d",
            "legacy_beta_tors_member": f"L_beta_tors(d)=1+{BETA_TORS}*d",
            "strict_nonlinear_target": "S(d)=1+d^(9/5)",
            "node_match_condition": "For d>0, L_gamma(d)=S(d) iff gamma=d^(4/5).",
        },
        "node_rows": rows,
        "separation_summary": separation_summary,
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all(cross_checks.values()),
        "proof_certificate": {
            "nonduplication_step": "rg was used to separate this linear-vs-nonlinear denominator no-go from prior damping monotonicity, inverse-branch, ontology, and component-gap reports.",
            "algebraic_step": "If 1+gamma*d=1+d^(9/5) for d>0, then gamma=d^(4/5); because d^(4/5) is strictly increasing on positive d, no constant gamma can match two distinct positive nodes.",
            "legacy_beta_tors_step": "The legacy beta_tors=0.01 member matches no positive strict node on d=1..11; the required gamma already equals 1 at d=1 and increases to 11^(4/5).",
            "least_squares_step": f"The best finite L2 linear denominator fit has gamma={gamma_l2:.12f}, but its residual norm is {separation_summary['best_l2_residual_l2_norm']:.12f}, so even the best constant linear torsion replacement is not exact on the audited domain.",
            "bridge_meaning_step": "The strict d^eta denominator is therefore a genuine strict-side compression addition that must be supplied by the completion map; it is not hidden inside beta_tors*d.",
            "theoretical_limit": "This is a finite/algebraic separation certificate for the damping-compression row, not a strict-source derivation of eta, beta, or the full legacy->strict bridge.",
        },
        "hard_limits": [
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No derivation of eta=9/5, beta=1, or strict d^eta compression from nadsoliton dynamics is exported here.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }
    return payload


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-linear vs strict-nonlinear damping compression separation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Separation summary",
        "",
    ]
    for key, value in payload["separation_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", f"All cross-checks pass: `{payload['all_cross_checks_pass']}`", "", "## Proof certificate", ""])
    for key, value in payload["proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
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
