#!/usr/bin/env python3
"""Scratch probe: legacy -> strict phase/frequency affine transport certificate.

This bridge step isolates the phase/frequency row.  It packages two finite facts:

1. The continuous affine phase-coordinate map

       x(d) = (omega_S*d + phi_S - phi_L)/omega_L

   exactly transports the legacy cosine argument to the strict cosine argument.
2. That map is not a Z12 automorphism/reindexing and the phase transport factor
   P(d)=cos_S(d)/cos_L(d) is not replaceable by one scalar on the audited domain.

So this is a bridge witness for the phase/frequency passage, but not a strict
source theorem for omega/phi, not an orientation/selector theorem, and not a
legacy-role transfer theorem.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.md"

SOURCE_REPORTS = {
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "phase_sign_z2_coboundary": HERE / "bridge_strict_completion_phase_sign_z2_coboundary_certificate_report.json",
    "gf2_linear_system": HERE / "bridge_strict_completion_phase_sign_gf2_linear_system_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "finite_diagonal_completion_map": HERE / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

DOMAIN = list(range(12))
OMEGA_L = math.pi / 4.0
PHI_L = math.pi / 6.0
OMEGA_S = 0.18575
PHI_S = 0.16250
UNITS_Z12 = [1, 5, 7, 11]
TOL = 1e-14


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def sign(value: float) -> int:
    if value > 0.0:
        return 1
    if value < 0.0:
        return -1
    return 0


def bit_from_sign(value: int) -> int:
    if value == 1:
        return 0
    if value == -1:
        return 1
    raise ValueError("zero sign cannot be converted to GF(2) bit")


def best_scalar_fit(values: list[float], targets: list[float]) -> dict[str, float]:
    denominator = sum(value * value for value in values)
    scalar = sum(value * target for value, target in zip(values, targets)) / denominator
    residuals = [scalar * value - target for value, target in zip(values, targets)]
    return {
        "best_scalar": scalar,
        "max_abs_residual": max(abs(value) for value in residuals),
        "l2_residual": math.sqrt(sum(value * value for value in residuals)),
    }


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    necessity = loaded["necessity"]
    z2 = loaded["phase_sign_z2_coboundary"]
    gf2 = loaded["gf2_linear_system"]
    component_gap = loaded["component_gap_matrix"]
    diagonal = loaded["finite_diagonal_completion_map"]
    guardrail = loaded["legacy_bridge_guardrail"]

    z2_bits = {row["d"]: row["node_bit"] for row in z2["node_bit_rows"]}
    point_rows = necessity["pointwise_quotient_certificate"]
    phase_rows = []
    legacy_cos_values = []
    strict_cos_values = []
    for row in point_rows:
        d = row["d"]
        legacy_arg = OMEGA_L * d + PHI_L
        strict_arg = OMEGA_S * d + PHI_S
        x_affine = (strict_arg - PHI_L) / OMEGA_L
        legacy_cos = math.cos(legacy_arg)
        strict_cos = math.cos(strict_arg)
        legacy_cos_at_x = math.cos(OMEGA_L * x_affine + PHI_L)
        phase_factor = row["phase_frequency_transport_factor"]
        phase_sign = sign(phase_factor)
        phase_bit = bit_from_sign(phase_sign)
        legacy_cos_values.append(legacy_cos)
        strict_cos_values.append(strict_cos)
        phase_rows.append(
            {
                "d": d,
                "legacy_argument": legacy_arg,
                "strict_argument": strict_arg,
                "affine_legacy_coordinate_x_d": x_affine,
                "x_d_minus_d": x_affine - d,
                "legacy_cos_at_integer_d": legacy_cos,
                "strict_cos_at_d": strict_cos,
                "legacy_cos_at_affine_x_d": legacy_cos_at_x,
                "affine_transport_residual": legacy_cos_at_x - strict_cos,
                "phase_frequency_transport_factor_P_d": phase_factor,
                "phase_factor_sign": phase_sign,
                "phase_factor_gf2_bit": phase_bit,
                "matches_z2_node_bit": phase_bit == z2_bits[d],
            }
        )

    legacy_sign_pattern = [sign(value) for value in legacy_cos_values]
    strict_sign_pattern = [sign(value) for value in strict_cos_values]
    automorphism_rows = []
    for unit in UNITS_Z12:
        for offset in DOMAIN:
            mapped = [legacy_sign_pattern[(unit * d + offset) % 12] for d in DOMAIN]
            mismatch_positions = [d for d, got, want in zip(DOMAIN, mapped, strict_sign_pattern) if got != want]
            automorphism_rows.append(
                {
                    "unit": unit,
                    "offset": offset,
                    "mapped_legacy_sign_pattern": mapped,
                    "mismatch_positions_against_strict_sign": mismatch_positions,
                    "mismatch_count": len(mismatch_positions),
                    "matches_strict_sign_pattern": len(mismatch_positions) == 0,
                }
            )

    scalar_fit = best_scalar_fit(legacy_cos_values, strict_cos_values)
    affine_slope = OMEGA_S / OMEGA_L
    affine_intercept = (PHI_S - PHI_L) / OMEGA_L
    x_values = [row["affine_legacy_coordinate_x_d"] for row in phase_rows]
    integer_distances = [abs(x - round(x)) for x in x_values]
    phase_bits = [row["phase_factor_gf2_bit"] for row in phase_rows]

    summary = {
        "domain_size": len(DOMAIN),
        "affine_slope_omega_s_over_omega_l": affine_slope,
        "affine_intercept": affine_intercept,
        "continuous_affine_phase_transport_exact": max(abs(row["affine_transport_residual"]) for row in phase_rows) <= TOL,
        "max_abs_affine_transport_residual": max(abs(row["affine_transport_residual"]) for row in phase_rows),
        "affine_map_is_not_z12_automorphism": affine_slope not in UNITS_Z12,
        "affine_coordinates_not_all_integers": any(distance > TOL for distance in integer_distances),
        "minimum_distance_to_integer_affine_coordinate": min(integer_distances),
        "maximum_distance_to_integer_affine_coordinate": max(integer_distances),
        "z12_unit_offset_automorphism_count_checked": len(automorphism_rows),
        "no_z12_unit_offset_reindex_matches_strict_sign_pattern": not any(row["matches_strict_sign_pattern"] for row in automorphism_rows),
        "best_z12_unit_offset_mismatch_count": min(row["mismatch_count"] for row in automorphism_rows),
        "phase_factor_bits": phase_bits,
        "phase_factor_bits_match_z2_node_bits": all(row["matches_z2_node_bit"] for row in phase_rows),
        "gf2_solution_inherited_unique": gf2["gf2_linear_system_summary"]["unique_solution"],
        "scalar_phase_replacement_fails": scalar_fit["max_abs_residual"] > 1e-3,
        "scalar_phase_best_fit_max_abs_residual": scalar_fit["max_abs_residual"],
        "scalar_phase_best_fit_l2_residual": scalar_fit["l2_residual"],
        "finite_diagonal_completion_map_inherited": diagonal["finite_diagonal_completion_summary"]["unique_diagonal_completion_map_exists"],
        "component_gap_phase_row_still_source_open": component_gap["completion_gap_summary"]["selector_source_gap_remains"],
        "strict_phase_frequency_source_exported": False,
        "orientation_selector_source_exported": False,
        "raw_kernel_identity_claimed": False,
        "full_bridge_theorem_exported": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_PHASE_FREQUENCY_AFFINE_TRANSPORT_CERTIFICATE__NO_SELECTOR_SOURCE",
        "status": "continuous-affine-phase-transport-exact__not-z12-automorphism__phase-source-open",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "phase affine",
                "affine phase transport",
                "phase/frequency bridge",
                "omega_S/omega_L",
                "Z12 automorphism phase",
                "phase scalar replacement",
            ],
            "finding": "Existing reports certify pointwise phase factors and GF(2) signs; this report separately packages the continuous affine phase-coordinate transport, rejects Z12 automorphism/scalar replacements, and keeps source/selector closure open.",
        },
        "phase_affine_definition": {
            "legacy_argument": "theta_L(x)=omega_L*x+phi_L",
            "strict_argument": "theta_S(d)=omega_S*d+phi_S",
            "affine_transport": "x(d)=(theta_S(d)-phi_L)/omega_L=(omega_S/omega_L)d+(phi_S-phi_L)/omega_L",
            "omega_L": OMEGA_L,
            "phi_L": PHI_L,
            "omega_S": OMEGA_S,
            "phi_S": PHI_S,
            "z12_units_checked": UNITS_Z12,
        },
        "phase_transport_rows": phase_rows,
        "z12_unit_offset_automorphism_rows": automorphism_rows,
        "best_scalar_phase_fit": scalar_fit,
        "phase_frequency_affine_transport_summary": summary,
        "cross_checks": {
            "guardrail_restores_legacy_as_intermediate": guardrail["legacy_kernel_intermediate_bridge_summary"]["legacy_kernel_restored_as_intermediate"],
            "continuous_affine_transport_exact": summary["continuous_affine_phase_transport_exact"],
            "not_discrete_z12_automorphism": summary["affine_map_is_not_z12_automorphism"] and summary["affine_coordinates_not_all_integers"] and summary["no_z12_unit_offset_reindex_matches_strict_sign_pattern"],
            "phase_bits_match_gf2_chain": summary["phase_factor_bits_match_z2_node_bits"] and summary["gf2_solution_inherited_unique"],
            "scalar_phase_replacement_rejected": summary["scalar_phase_replacement_fails"],
            "no_source_selector_or_bridge_claim": not summary["strict_phase_frequency_source_exported"] and not summary["orientation_selector_source_exported"] and not summary["raw_kernel_identity_claimed"] and not summary["full_bridge_theorem_exported"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to separate this affine phase-frequency transport certificate from existing pointwise phase-factor and GF(2) sign reports.",
            "affine_step": "By construction theta_L(x(d))=omega_L*((omega_S*d+phi_S-phi_L)/omega_L)+phi_L=theta_S(d), so the continuous phase-coordinate transport is exact.",
            "non_automorphism_step": "The affine slope omega_S/omega_L is not a Z12 unit and the sampled x(d) are not all integers; exhaustive unit+offset checks over Aut(Z12) find no reindexing that reproduces the strict cosine sign pattern.",
            "phase_factor_step": "The pointwise factor P(d)=cos(theta_S(d))/cos(theta_L(d)) has GF(2) sign bits matching the existing Z2/GF(2) phase-sign chain.",
            "non_scalar_step": "The best scalar fit from the legacy cosine carrier to the strict cosine carrier has nonzero residual, so phase/frequency transport is not scalar normalization.",
            "theoretical_limit": "This is a finite/affine phase-frequency bridge witness, not a strict derivation of omega/phi, not an orientation selector source, not beta_tors->chi_11, and not ToE closure.",
        },
        "hard_limits": [
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No strict dynamical source for omega/phi or phase transport is exported.",
            "No orientation/selector source is exported.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-to-strict phase/frequency affine transport certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The continuous affine coordinate `x(d)` exactly transports the legacy phase",
        "argument to the strict phase argument, but it is not a `Z12` automorphism",
        "or scalar replacement.  The selector/source question remains open.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["phase_frequency_affine_transport_summary"].items():
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
    print(json.dumps(payload["phase_frequency_affine_transport_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
