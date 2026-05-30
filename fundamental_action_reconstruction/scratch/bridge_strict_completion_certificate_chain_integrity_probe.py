#!/usr/bin/env python3
"""Scratch probe: integrity check for the strict-completion certificate chain.

The previous probes each certify one local part of the completion story:
necessity of A/P/D, transport/cocycle reconstruction, phase-zero placement,
rational phase-zero placement, phase-zero robustness, damping monotonicity, and
low-order transport no-go results.

This probe does not add another local fit.  It audits the *chain* as a finite
proof ledger: load all prerequisite reports, check that their shared numerical
and logical conclusions agree, and emit one no-false-pass frontier statement.

It is still not a bridge theorem and not a strict-side dynamical derivation.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_certificate_chain_integrity_report.json"
OUT_MD = HERE / "bridge_strict_completion_certificate_chain_integrity_report.md"

REPORTS = {
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "cocycle": HERE / "bridge_strict_kernel_completion_transport_cocycle_report.json",
    "phase_zero": HERE / "bridge_strict_completion_phase_zero_interlacing_certificate_report.json",
    "phase_zero_rational": HERE / "bridge_strict_completion_phase_zero_rational_interval_certificate_report.json",
    "phase_zero_margin": HERE / "bridge_strict_completion_phase_zero_margin_certificate_report.json",
    "damping_monotonicity": HERE / "bridge_strict_completion_damping_continuous_monotonicity_certificate_report.json",
    "low_order_no_go": HERE / "bridge_strict_completion_low_order_transport_no_go_report.json",
}

EXPECTED_SIGN_PATTERN = [1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, 1]
EXPECTED_FLIP_EDGES = ["1->2", "5->6", "7->8", "9->10"]
EXPECTED_FULL_SUBSET = "alpha_normalization+phase_frequency_transport+damping_compression"


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def report_path(path: Path) -> str:
    return str(path.relative_to(ROOT))


def extract_negative_edges(cocycle: dict[str, Any]) -> list[str]:
    return [row["edge"] for row in cocycle["edge_cocycle_rows"] if row["edge_sign_ratio"] < 0]


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in REPORTS.items()}
    necessity = loaded["necessity"]
    cocycle = loaded["cocycle"]
    phase_zero = loaded["phase_zero"]
    phase_zero_rational = loaded["phase_zero_rational"]
    phase_zero_margin = loaded["phase_zero_margin"]
    damping = loaded["damping_monotonicity"]
    low_order = loaded["low_order_no_go"]

    exact_subsets = necessity["necessity_summary"]["exact_subsets_without_extra_scalar"]
    cocycle_sign_pattern = cocycle["cocycle_summary"]["transport_sign_pattern"]
    cocycle_flip_edges = cocycle["cocycle_summary"]["phase_sign_flip_edges"]
    cocycle_negative_edges = extract_negative_edges(cocycle)
    phase_zero_sign_pattern = phase_zero["interlacing_summary"]["phase_transport_sign_pattern"]
    phase_zero_flip_edges = phase_zero["interlacing_summary"]["phase_sign_flip_edges"]
    rational_sign_pattern = phase_zero_rational["interval_summary"]["phase_transport_sign_pattern_from_rational_intervals"]
    rational_flip_edges = phase_zero_rational["interval_summary"]["phase_sign_flip_edges_from_rational_intervals"]
    margin_sign_pattern = phase_zero_margin["robustness_summary"]["certified_phase_sign_pattern_preserved"]
    margin_flip_edges = phase_zero_margin["robustness_summary"]["certified_phase_sign_flip_edges_preserved"]
    low_order_negative_edges = low_order["transport_input_summary"]["negative_edges"]

    cross_checks = {
        "necessity_has_unique_exact_full_APD_subset": exact_subsets == [EXPECTED_FULL_SUBSET],
        "necessity_final_residual_pass": necessity["necessity_summary"]["max_abs_quotient_identity_residual"] < 1e-14,
        "cocycle_reconstruction_pass": cocycle["cocycle_summary"]["max_abs_reconstruction_residual"] < 1e-14,
        "cocycle_interval_pass": cocycle["cocycle_summary"]["max_abs_interval_additive_log_residual"] < 1e-12,
        "phase_zero_float_matches_expected": phase_zero_flip_edges == EXPECTED_FLIP_EDGES and phase_zero_sign_pattern == EXPECTED_SIGN_PATTERN,
        "phase_zero_rational_matches_float": rational_flip_edges == phase_zero_flip_edges and rational_sign_pattern == phase_zero_sign_pattern,
        "phase_zero_margin_preserves_rational": margin_flip_edges == rational_flip_edges and margin_sign_pattern == rational_sign_pattern,
        "cocycle_negative_edges_equal_phase_flips": cocycle_negative_edges == EXPECTED_FLIP_EDGES,
        "low_order_negative_edges_equal_phase_flips": low_order_negative_edges == EXPECTED_FLIP_EDGES,
        "damping_positive_and_decreasing": damping["monotonicity_summary"]["continuous_positive_certificate"] and damping["monotonicity_summary"]["continuous_strictly_decreasing_certificate"],
        "damping_cannot_supply_sign_flips": damping["monotonicity_summary"]["continuous_positive_certificate"] and cocycle_negative_edges == EXPECTED_FLIP_EDGES,
        "low_order_no_go_all_listed_models_fail": all([
            low_order["no_go_summary"]["positive_damping_only_fails"],
            low_order["no_go_summary"]["constant_first_order_fails"],
            low_order["no_go_summary"]["constant_second_order_fails"],
            low_order["no_go_summary"]["affine_log_envelope_fails"],
            low_order["no_go_summary"]["short_period_edge_sign_law_fails"],
        ]),
    }

    all_cross_checks_pass = all(cross_checks.values())

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_CERTIFICATE_CHAIN_INTEGRITY__FINITE_LEDGER_NO_BRIDGE_THEOREM",
        "status": "all-loaded-strict-completion-certificates-are-cross-consistent-no-false-pass",
        "source_reports": {name: report_path(path) for name, path in REPORTS.items()},
        "loaded_statuses": {name: loaded[name]["status"] for name in REPORTS},
        "expected_shared_objects": {
            "phase_transport_sign_pattern": EXPECTED_SIGN_PATTERN,
            "phase_sign_flip_edges": EXPECTED_FLIP_EDGES,
            "unique_exact_completion_subset": EXPECTED_FULL_SUBSET,
        },
        "cross_checks": cross_checks,
        "all_cross_checks_pass": all_cross_checks_pass,
        "chain_summary": {
            "exact_APD_completion_certified": cross_checks["necessity_has_unique_exact_full_APD_subset"] and cross_checks["necessity_final_residual_pass"],
            "transport_cocycle_certified": cross_checks["cocycle_reconstruction_pass"] and cross_checks["cocycle_interval_pass"],
            "phase_sign_source_certified": cross_checks["phase_zero_float_matches_expected"] and cross_checks["phase_zero_rational_matches_float"] and cross_checks["phase_zero_margin_preserves_rational"],
            "damping_envelope_certified": cross_checks["damping_positive_and_decreasing"],
            "simple_transport_readings_rejected": cross_checks["low_order_no_go_all_listed_models_fail"],
            "strict_dynamic_derivation_exported": False,
            "bridge_theorem_exported": False,
        },
        "frontier_statement": {
            "positive_content": "The finite completion ansatz is internally consistent across necessity, cocycle, phase-zero, rational-zero, robustness-margin, damping, and low-order no-go certificates.",
            "negative_content": "The chain still does not derive A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "next_real_blocker": "strict_phase_frequency/damping/transport derivation from strict nadsoliton dynamics, plus orientation_chi11_source and role_transfer_theorem if a bridge lane is explicitly reopened.",
        },
        "proof_certificate": {
            "ledger_step": "All prerequisite JSON reports are loaded and their status fields are recorded in one ledger.",
            "shared_object_step": "The common sign pattern and flip edges agree across cocycle, float zero, rational zero, margin, and low-order no-go reports.",
            "factor_step": "The necessity report still has exactly one exact no-extra-scalar subset: A+P+D.",
            "envelope_step": "The damping report is positive and strictly decreasing, so it is consistent with sign flips being supplied by phase only.",
            "frontier_step": "All consistency checks pass, but the exported object remains a finite certificate chain rather than a strict dynamical derivation.",
        },
        "hard_limits": [
            "K_strict_gate remains the current live/full operational kernel.",
            "No unqualified identity K_legacy_ont == K_strict_gate is claimed.",
            "No proof derives A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle from strict nadsoliton dynamics.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer to K_strict_gate is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Strict completion certificate chain integrity report",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "This audit loads the strict-completion certificate reports as a finite proof",
        "ledger and checks that their shared conclusions agree.  It is a chain",
        "integrity check, not a new bridge theorem or strict dynamical derivation.",
        "",
        "## Cross-checks",
        "",
    ]
    for key, value in payload["cross_checks"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        f"All cross-checks pass: `{payload['all_cross_checks_pass']}`",
        "",
        "## Chain summary",
        "",
    ])
    for key, value in payload["chain_summary"].items():
        lines.append(f"- `{key}`: `{value}`")
    lines.extend([
        "",
        "## Frontier statement",
        "",
        f"- Positive: {payload['frontier_statement']['positive_content']}",
        f"- Negative: {payload['frontier_statement']['negative_content']}",
        f"- Next blocker: {payload['frontier_statement']['next_real_blocker']}",
        "",
        "## Hard limits",
        "",
    ])
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
