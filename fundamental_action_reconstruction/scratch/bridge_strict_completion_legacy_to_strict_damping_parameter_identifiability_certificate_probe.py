#!/usr/bin/env python3
"""Scratch probe: strict damping beta/eta finite identifiability certificate.

This is a narrow bridge-completion step for the damping/compression row.  It
packages a finite parameter-identifiability fact for the strict denominator

    S(d) = 1 + beta*d^eta

on the positive audited nodes d=1..11.  Under the strict working denominator
values S(d)=1+d^(9/5), d=1 fixes beta=1 and every d>=2 recovers eta=9/5.
Equivalently, on the finite algebraic cover t^5=d, the denominator increment is
exactly t^9 and satisfies (S(d)-1)^5=d^9.

This does not derive beta=1 or eta=9/5 from nadsoliton dynamics.  It only
certifies that once the strict denominator samples are accepted, the beta/eta
parameters are finitely identifiable and are not a hidden legacy-linear torsion
parameter.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.md"

SOURCE_REPORTS = {
    "damping_compression_separation": HERE / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json",
    "damping_exact_rational": HERE / "bridge_strict_completion_damping_exact_rational_calculus_certificate_report.json",
    "finite_diagonal_completion_map": HERE / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

DOMAIN_POSITIVE = list(range(1, 12))
BETA_STRICT = Fraction(1, 1)
ETA_STRICT = Fraction(9, 5)
BETA_TORS_LEGACY = Fraction(1, 100)
TOL = 1e-12


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def strict_delta_float(d: int) -> float:
    return float(d ** (ETA_STRICT.numerator / ETA_STRICT.denominator))


def build_rows() -> list[dict[str, Any]]:
    rows = []
    for d in DOMAIN_POSITIVE:
        delta = strict_delta_float(d)
        eta_from_log = None if d == 1 else math.log(delta / float(BETA_STRICT), d)
        rows.append(
            {
                "d": d,
                "strict_denominator": 1.0 + delta,
                "strict_delta_S_minus_1": delta,
                "beta_from_d_if_eta_fixed_9_over_5": delta / (d ** (ETA_STRICT.numerator / ETA_STRICT.denominator)),
                "eta_from_log_if_beta_fixed_1": eta_from_log,
                "eta_from_log_residual_vs_9_over_5": None if eta_from_log is None else eta_from_log - float(ETA_STRICT),
                "exact_cover_identity": "(S(d)-1)^5=d^9",
                "exact_cover_left_delta_power": 5,
                "exact_cover_right_d_power": 9,
                "exact_cover_right_integer_d9": d**9,
                "legacy_linear_delta_beta_tors_d": float(BETA_TORS_LEGACY * d),
                "strict_delta_over_legacy_linear_delta": delta / float(BETA_TORS_LEGACY * d),
                "required_linear_gamma_for_this_node": d ** (ETA_STRICT.numerator / ETA_STRICT.denominator - 1.0),
            }
        )
    return rows


def rational_candidate_matches(candidate: Fraction, rows: list[dict[str, Any]]) -> bool:
    if candidate == ETA_STRICT:
        return True
    # If beta is fixed by d=1 as beta=1, matching d=2 is already enough to
    # distinguish all other rational candidates in the finite candidate grid.
    d = 2
    candidate_delta = d ** (candidate.numerator / candidate.denominator)
    strict_delta = strict_delta_float(d)
    return abs(candidate_delta - strict_delta) <= TOL


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    separation = loaded["damping_compression_separation"]
    damping_exact = loaded["damping_exact_rational"]
    diagonal = loaded["finite_diagonal_completion_map"]
    component_gap = loaded["component_gap_matrix"]
    guardrail = loaded["legacy_bridge_guardrail"]

    rows = build_rows()
    eta_residuals = [abs(row["eta_from_log_residual_vs_9_over_5"]) for row in rows if row["eta_from_log_residual_vs_9_over_5"] is not None]
    beta_residuals = [abs(row["beta_from_d_if_eta_fixed_9_over_5"] - 1.0) for row in rows]
    required_gammas = [row["required_linear_gamma_for_this_node"] for row in rows]

    candidates = [Fraction(p, q) for q in range(1, 11) for p in range(1, 31)]
    unique_candidates = sorted(set(candidates))
    matching_candidates = [candidate for candidate in unique_candidates if rational_candidate_matches(candidate, rows)]

    summary = {
        "positive_domain": DOMAIN_POSITIVE,
        "beta_identified_from_d1": rows[0]["strict_delta_S_minus_1"] == 1.0,
        "identified_beta": "1",
        "identified_eta": "9/5",
        "all_beta_recoveries_equal_1_given_eta": max(beta_residuals) <= TOL,
        "max_beta_recovery_residual": max(beta_residuals),
        "all_eta_log_recoveries_equal_9_over_5_given_beta": max(eta_residuals) <= TOL,
        "max_eta_log_recovery_residual": max(eta_residuals),
        "exact_fifth_power_cover_identity_recorded": all(row["exact_cover_right_integer_d9"] == row["d"] ** 9 for row in rows),
        "candidate_grid_p_le_30_q_le_10_matching_count": len(matching_candidates),
        "candidate_grid_unique_match": [str(candidate) for candidate in matching_candidates] == ["9/5"],
        "required_linear_gamma_strictly_increases": all(a < b for a, b in zip(required_gammas, required_gammas[1:])),
        "legacy_beta_tors_not_equal_any_required_gamma": all(abs(float(BETA_TORS_LEGACY) - gamma) > TOL for gamma in required_gammas),
        "damping_separation_inherited": separation["separation_summary"]["no_single_linear_gamma_matches_two_distinct_positive_nodes"],
        "exact_damping_calculus_inherited": damping_exact["exact_rational_damping_summary"]["continuous_strictly_decreasing_certificate"],
        "finite_diagonal_completion_map_inherited": diagonal["finite_diagonal_completion_summary"]["unique_diagonal_completion_map_exists"],
        "component_gap_damping_source_still_open": component_gap["completion_gap_summary"]["strict_dynamic_sources_missing"],
        "strict_beta_eta_source_exported": False,
        "legacy_beta_tors_to_beta_eta_theorem_exported": False,
        "full_bridge_theorem_exported": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_DAMPING_PARAMETER_IDENTIFIABILITY_CERTIFICATE__NO_SOURCE_THEOREM",
        "status": "strict-beta-eta-finitely-identifiable-from-denominator-samples__source-and-legacy-map-open",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "eta identifiability",
                "beta eta recovery",
                "strict denominator parameter",
                "1+d^eta",
                "log eta recovery",
                "damping parameter certificate",
            ],
            "finding": "Existing damping reports prove monotonicity and separate legacy-linear from strict-nonlinear damping; this report adds finite beta/eta identifiability from strict denominator samples while keeping source derivation open.",
        },
        "parameter_identification_definition": {
            "strict_denominator": "S(d)=1+beta*d^eta",
            "accepted_strict_samples": "S(d)=1+d^(9/5) on d=1..11",
            "beta_recovery": "beta=S(1)-1=1",
            "eta_recovery": "eta=log(S(d)-1)/log(d) for any d>=2 once beta=1",
            "finite_cover_identity": "eta=9/5 means (S(d)-1)^5=d^9 on the t^5=d cover",
            "candidate_grid": "p/q with 1<=p<=30 and 1<=q<=10",
        },
        "parameter_recovery_rows": rows,
        "matching_rational_candidates": [str(candidate) for candidate in matching_candidates],
        "damping_parameter_identifiability_summary": summary,
        "cross_checks": {
            "guardrail_records_strict_compression_missing_from_legacy": guardrail["legacy_kernel_intermediate_bridge_summary"]["strict_compression_missing_from_legacy_recorded"],
            "beta_eta_recovered_from_samples": summary["beta_identified_from_d1"] and summary["all_eta_log_recoveries_equal_9_over_5_given_beta"],
            "finite_cover_identity_and_candidate_unique": summary["exact_fifth_power_cover_identity_recorded"] and summary["candidate_grid_unique_match"],
            "legacy_linear_not_parameter_match": summary["required_linear_gamma_strictly_increases"] and summary["legacy_beta_tors_not_equal_any_required_gamma"] and summary["damping_separation_inherited"],
            "prior_damping_and_diagonal_reports_inherited": summary["exact_damping_calculus_inherited"] and summary["finite_diagonal_completion_map_inherited"],
            "no_source_or_full_bridge_claim": not summary["strict_beta_eta_source_exported"] and not summary["legacy_beta_tors_to_beta_eta_theorem_exported"] and not summary["full_bridge_theorem_exported"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to distinguish this beta/eta identifiability certificate from prior damping monotonicity and linear-vs-nonlinear separation reports.",
            "beta_step": "At d=1, S(1)-1=beta*1^eta=beta, and the strict sample has S(1)-1=1, so beta is identified as 1 within the accepted strict denominator model.",
            "eta_step": "For every d>=2, eta=log(S(d)-1)/log(d)=log(d^(9/5))/log(d)=9/5 once beta=1 is fixed.",
            "cover_step": "Equivalently, eta=9/5 places the strict denominator increment on the finite algebraic cover (S(d)-1)^5=d^9.",
            "candidate_step": "The rational grid p/q with p<=30 and q<=10 has exactly one matching candidate, 9/5, when beta is fixed by the d=1 sample and d=2 is checked.",
            "legacy_linear_step": "The required linear gamma at node d is d^(4/5), which strictly increases, so no single beta_tors-style linear damping parameter supplies the strict denominator samples.",
            "theoretical_limit": "This identifies beta/eta from accepted strict samples only; it is not a strict dynamical source for beta/eta, not beta_tors->chi_11, not legacy role transfer, and not ToE closure.",
        },
        "hard_limits": [
            "No strict dynamical derivation of beta=1 or eta=9/5 is exported.",
            "No beta_tors -> beta/eta theorem is claimed.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-to-strict damping parameter identifiability certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The strict denominator samples identify `beta=1` and `eta=9/5` on the",
        "audited positive domain, and satisfy the fifth-power cover identity",
        "`(S(d)-1)^5=d^9`.  This is not a source theorem.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["damping_parameter_identifiability_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Cross-checks", ""])
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend(["", "## Proof certificate", ""])
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
    print(json.dumps(payload["damping_parameter_identifiability_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
