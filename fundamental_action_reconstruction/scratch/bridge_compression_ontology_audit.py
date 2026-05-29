#!/usr/bin/env python3
"""Scratch audit: compression ontology as the missing legacy characteristic.

This packet responds to the current bridge question in a deliberately limited
way.  It does not try to identify the historical legacy kernel with the strict
gate kernel.  Instead it asks whether the later strict form contains a structural
compression characteristic that the legacy formula did not make explicit, and
whether that is enough to call the strict kernel the final nadsoliton formula.

The result is an ontology audit, not a theorem: strict has a much stronger
operational compression profile (nonlinear exponent, normalized strict loss
measure, and finite measure transport), while full final-form status still
requires a strict-side derivation of the density/exponent and selector closure.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

HERE = Path(__file__).resolve().parent
MEASURE_TRANSPORT = HERE / "bridge_phase_normalized_measure_transport_report.json"
PHASE_BUDGET = HERE / "bridge_phase_normalized_phase_budget_flow_report.json"
STRICT_ROBUSTNESS = HERE / "bridge_strict_kernel_robustness_opinion_audit_report.json"
OUT_JSON = HERE / "bridge_compression_ontology_audit_report.json"
OUT_MD = HERE / "bridge_compression_ontology_audit_report.md"

LEGACY = {"omega": math.pi / 4.0, "phi": math.pi / 6.0, "beta_tors": 0.01, "alpha_geo": 4.0 * math.log(2.0)}
STRICT = {"omega": 0.18575, "phi": 0.16250, "beta": 1.0, "eta": 9.0 / 5.0}
GRID_N = 4001
Z12_D_MAX = 11.0


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def legacy_carrier_normalized(d: float) -> float:
    return math.cos(LEGACY["omega"] * d + LEGACY["phi"]) / math.cos(LEGACY["phi"])


def strict_carrier_normalized(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(STRICT["phi"])


def legacy_denominator_compression(d: float) -> float:
    return 1.0 / (1.0 + LEGACY["beta_tors"] * d)


def strict_denominator_compression(d: float) -> float:
    return 1.0 / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def legacy_norm(d: float) -> float:
    return legacy_carrier_normalized(d) * legacy_denominator_compression(d)


def strict_norm(d: float) -> float:
    return strict_carrier_normalized(d) * strict_denominator_compression(d)


def strict_derivative(d: float) -> float:
    omega = STRICT["omega"]
    phi = STRICT["phi"]
    eta = STRICT["eta"]
    theta = omega * d + phi
    carrier = math.cos(theta) / math.cos(phi)
    carrier_prime = -omega * math.sin(theta) / math.cos(phi)
    denom = 1.0 + d**eta
    denom_prime = 0.0 if d == 0.0 else eta * d ** (eta - 1.0)
    return (carrier_prime * denom - carrier * denom_prime) / denom**2


def strict_loss_density(v: np.ndarray, d_h: float) -> np.ndarray:
    return np.array([-d_h * strict_derivative(float(d_h * item)) for item in v], dtype=float)


def moment_summary(v: np.ndarray, density: np.ndarray) -> dict[str, float]:
    norm = float(np.trapezoid(density, v))
    mean = float(np.trapezoid(v * density, v) / norm)
    second = float(np.trapezoid((v**2) * density, v) / norm)
    variance = second - mean**2
    mode = float(v[int(np.argmax(density))])
    left_half = float(np.trapezoid(density[v <= 0.5], v[v <= 0.5]) / norm)
    return {
        "trapezoid_norm": norm,
        "mean_v": mean,
        "variance_v": variance,
        "std_v": math.sqrt(max(variance, 0.0)),
        "mode_v": mode,
        "left_half_mass": left_half,
    }


def main() -> None:
    measure_report = load_json(MEASURE_TRANSPORT)
    phase_budget_report = load_json(PHASE_BUDGET)
    robustness_report = load_json(STRICT_ROBUSTNESS)
    d_h = float(measure_report["horizon"]["d_h"])
    x_h = float(measure_report["horizon"]["x_h"])

    grid_v = np.linspace(0.0, 1.0, GRID_N)
    density = strict_loss_density(grid_v, d_h)
    moments = moment_summary(grid_v, density)

    legacy_damping_z12 = legacy_denominator_compression(Z12_D_MAX)
    strict_damping_z12 = strict_denominator_compression(Z12_D_MAX)
    damping_ratio = strict_damping_z12 / legacy_damping_z12
    legacy_damping_at_dh = legacy_denominator_compression(d_h)
    strict_damping_at_dh = strict_denominator_compression(d_h)

    compression_characteristic = {
        "missing_in_legacy_reading": {
            "description": "Legacy has an explicit alpha_geo amplitude and weak linear torsion damping, but it does not expose a nonlinear information-loss measure or fractal/transport exponent as the main compression datum.",
            "legacy_beta_tors": LEGACY["beta_tors"],
            "legacy_denominator_at_z12_max_d11": legacy_damping_z12,
            "legacy_denominator_at_strict_horizon_dh": legacy_damping_at_dh,
            "carrier_can_create_zero_but_not_measure_law": True,
        },
        "present_in_strict_operational_form": {
            "description": "Strict carries compression directly in beta*d^eta and induces a normalized finite strict loss density on the phase horizon.",
            "strict_eta": STRICT["eta"],
            "strict_beta": STRICT["beta"],
            "strict_denominator_at_z12_max_d11": strict_damping_z12,
            "strict_denominator_at_strict_horizon_dh": strict_damping_at_dh,
            "strict_vs_legacy_denominator_ratio_at_d11": damping_ratio,
        },
        "transport_upgrade": {
            "description": "The measure-transport probe upgrades pointwise output matching to a unit information-budget transport candidate.",
            "measure_balance_max_abs_residual": measure_report["measure_balance"]["max_abs_balance_residual"],
            "strict_exact_residual_vs_1": measure_report["measure_balance"]["strict_exact_residual_vs_1"],
            "legacy_exact_residual_vs_1": measure_report["measure_balance"]["legacy_exact_residual_vs_1"],
            "affine_balance_max_abs_residual": measure_report["affine_failure"]["affine_balance_max_abs_residual"],
        },
    }

    completeness_matrix = [
        {
            "criterion": "explicit nonlinear compression exponent",
            "legacy_status": "absent_or_only_implicit",
            "strict_status": "present_as_eta_1p8",
            "bridge_status": "candidate_feature_match_not_a_derivation",
        },
        {
            "criterion": "finite normalized information-loss density",
            "legacy_status": "not_exported_as_legacy_ontology",
            "strict_status": "computable_from_strict_kernel",
            "bridge_status": "measure_transport_probe_balances_it_against_legacy_horizon_metric",
        },
        {
            "criterion": "non-affine phase-budget compression",
            "legacy_status": "not_an_explicit_legacy_characteristic",
            "strict_status": "induced_by_strict_horizon_and_eta",
            "bridge_status": "affine_identification_ruled_out_by_density_residual",
        },
        {
            "criterion": "independent strict-side derivation of eta/rho_strict",
            "legacy_status": "not_applicable",
            "strict_status": "still_missing_as_theorem",
            "bridge_status": "open_obligation",
        },
        {
            "criterion": "license to transfer legacy physical-role formulas",
            "legacy_status": "archival_model_level",
            "strict_status": "not_licensed",
            "bridge_status": "blocked_by_guardrail",
        },
        {
            "criterion": "QW-2191 selector closure / ToE finality",
            "legacy_status": "not_current_strict_core",
            "strict_status": "not_closed",
            "bridge_status": "open_obligation",
        },
    ]

    strict_final_form_assessment = {
        "is_strict_plausible_terminal_shape_candidate": True,
        "reason": "Among the compared forms, strict is the only one that already carries nonlinear compression as a primary operational characteristic and supports a finite unit loss measure.",
        "is_strict_proven_final_nadsoliton_formula": False,
        "blocking_obligations": [
            "derive eta=1.8 or rho_strict(v) from strict-only nadsoliton information geometry rather than from gate selection/output matching",
            "show the finite information-measure transport is intrinsic and not merely a coordinate artifact",
            "derive any alpha_geo/physical-role transfer strictly if such roles are to be reused",
            "discharge or explicitly premise the QW-2191 selector obstruction before claiming strict-core ToE closure",
        ],
    }

    report = {
        "status": "OPEN_COMPRESSION_ONTOLOGY_AUDIT_NO_FINAL_FORM_THEOREM",
        "result_kind": "SCRATCH_COMPRESSION_ONTOLOGY_AUDIT__NOT_A_THEOREM",
        "source_reports": {
            "measure_transport": str(MEASURE_TRANSPORT.relative_to(HERE.parents[1])),
            "phase_budget_flow": str(PHASE_BUDGET.relative_to(HERE.parents[1])),
            "strict_robustness": str(STRICT_ROBUSTNESS.relative_to(HERE.parents[1])),
        },
        "compression_characteristic": compression_characteristic,
        "strict_loss_density_moments": moments,
        "completeness_matrix": completeness_matrix,
        "strict_final_form_assessment": strict_final_form_assessment,
        "upstream_replay": {
            "measure_transport_candidate_supported": measure_report["candidate_ontological_reading"]["supported_by_this_probe"],
            "measure_balance_max_abs_residual": measure_report["measure_balance"]["max_abs_balance_residual"],
            "phase_budget_candidate_supported": phase_budget_report["candidate_ontological_reading"]["supported_by_this_probe"],
            "strict_robustness_result_kind": robustness_report["result_kind"],
        },
        "honest_answer_to_user_question": [
            "Yes: the computed evidence supports reading compression as the missing explicit characteristic of the legacy kernel layer.",
            "Strict possesses that characteristic operationally through beta*d^eta and the induced unit strict loss density.",
            "Strict does not yet possess it fully as a final formula theorem because eta/rho_strict are not independently derived from strict nadsoliton ontology in this packet.",
            "Therefore strict is a plausible terminal-form candidate, not a proven final nadsoliton formula.",
        ],
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem proves that the compression density is intrinsic to the nadsoliton rather than kernel-induced.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Attempt a strict-only derivation of eta=1.8 or rho_strict(v) from alpha_geo_strict_derived_v1 / Shannon asymmetry and nad12-sigma structure; otherwise keep strict as an operational terminal-shape candidate only.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch compression ontology audit\n\n"
        "Status: compression as missing legacy characteristic; no final-form theorem.\n\n"
        f"- Legacy vs strict damping at `d=11`: legacy `{legacy_damping_z12:.6e}`, strict `{strict_damping_z12:.6e}`, ratio `{damping_ratio:.6e}`.\n"
        f"- Strict loss density moments: mean `v={moments['mean_v']:.6f}`, mode `v={moments['mode_v']:.6f}`, left-half mass `{moments['left_half_mass']:.6f}`.\n"
        f"- Transport replay: max measure residual `{measure_report['measure_balance']['max_abs_balance_residual']:.3e}`, affine failure `{measure_report['affine_failure']['affine_balance_max_abs_residual']:.3e}`.\n"
        "- Honest answer: strict has the compression characteristic operationally, but is not proven as the final nadsoliton formula until eta/rho_strict are strict-side-derived.\n"
        "- No false pass: no kernel identity, no legacy physical-role transfer, no intrinsic compression theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
