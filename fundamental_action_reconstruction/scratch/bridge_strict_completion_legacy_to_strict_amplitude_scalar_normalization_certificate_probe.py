#!/usr/bin/env python3
"""Scratch probe: legacy -> strict amplitude scalar-normalization certificate.

This is a narrow bridge-completion step for the amplitude/normalization row.
It proves only that the visible legacy scalar alpha_geo=4 ln 2 is a global,
positive, factorizable amplitude multiplying the legacy kernel shape on the
finite Z12 audit domain:

    K_legacy_ont(d) = alpha_geo * L_shape(d)

where

    L_shape(d)=cos(pi*d/4+pi/6)/(1+beta_tors*d).

The probe intentionally does not claim that this scalar factor is the full
strict A(d) completion, does not transfer the old alpha_geo physical roles, and
does not identify K_legacy_ont with K_strict_gate.  It simply exports an exact
scalar-normalization witness that the bridge can use before the remaining phase,
damping/compression, selector, and role-transfer gaps are addressed.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.md"

SOURCE_REPORTS = {
    "agents_guardrail": ROOT / "AGENTS.md",
    "s2_priority_packet": ROOT / "fundamental_action_reconstruction" / "S2_CURRENT_FAR_STRATEGIC_PRIORITY_REORIENTATION_PACKET.md",
    "alpha_geo_strict_derived": ROOT / "fundamental_action_reconstruction" / "generated" / "alpha_geo_strict_derived_v1.json",
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

NODE_DOMAIN = list(range(12))
BETA_TORS = 0.01
ALPHA_GEO = 4.0 * math.log(2.0)
OMEGA_LEGACY = math.pi / 4.0
PHI_LEGACY = math.pi / 6.0


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def read_text(path: Path) -> str:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite text: {path}")
    return path.read_text(encoding="utf-8")


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def legacy_shape(d: int) -> float:
    return math.cos(OMEGA_LEGACY * d + PHI_LEGACY) / (1.0 + BETA_TORS * d)


def legacy_kernel(d: int) -> float:
    return ALPHA_GEO * legacy_shape(d)


def sign(value: float) -> int:
    return 1 if value > 0 else -1 if value < 0 else 0


def build_rows() -> list[dict[str, Any]]:
    rows = []
    for d in NODE_DOMAIN:
        shape = legacy_shape(d)
        kernel = legacy_kernel(d)
        normalized = kernel / ALPHA_GEO
        ratio = kernel / shape
        rows.append(
            {
                "d": d,
                "theta_mod_pi_over_12_numerator": 3 * d + 2,
                "legacy_denominator": 1.0 + BETA_TORS * d,
                "legacy_shape_L_d": shape,
                "legacy_kernel_K_d": kernel,
                "normalized_alpha_inverse_K_d": normalized,
                "normalization_residual": normalized - shape,
                "ratio_K_over_shape": ratio,
                "sign_K": sign(kernel),
                "sign_shape": sign(shape),
                "sign_preserved_by_positive_alpha": sign(kernel) == sign(shape),
            }
        )
    return rows


def build_payload() -> dict[str, Any]:
    agents = read_text(SOURCE_REPORTS["agents_guardrail"])
    s2 = read_text(SOURCE_REPORTS["s2_priority_packet"])
    alpha = load_json(SOURCE_REPORTS["alpha_geo_strict_derived"])
    necessity = load_json(SOURCE_REPORTS["necessity"])
    component_gap = load_json(SOURCE_REPORTS["component_gap_matrix"])
    guardrail = load_json(SOURCE_REPORTS["legacy_bridge_guardrail"])

    rows = build_rows()
    max_abs_residual = max(abs(row["normalization_residual"]) for row in rows)
    ratio_values = [row["ratio_K_over_shape"] for row in rows]
    max_ratio_deviation = max(abs(value - ALPHA_GEO) for value in ratio_values)
    zero_congruence_solutions = [d for d in NODE_DOMAIN if (3 * d + 2 - 6) % 12 == 0]
    denominator_positive = all(row["legacy_denominator"] > 0 for row in rows)
    legacy_shape_nonzero = all(row["legacy_shape_L_d"] != 0.0 for row in rows)

    amplitude_row = next(row for row in component_gap["component_gap_rows"] if row["component"] == "amplitude_normalization")
    exact_apd = necessity["necessity_summary"]["exact_subsets_without_extra_scalar"] == [
        "alpha_normalization+phase_frequency_transport+damping_compression"
    ]

    summary = {
        "legacy_alpha_geo_visible": "alpha_geo" in agents and "alpha_geo" in s2,
        "strict_alpha_geo_source_loaded": alpha["status"] == "actual_exported_strict_derived_source_upgrade_value" and alpha["value"] == "4 ln(2)",
        "finite_domain_nodes": NODE_DOMAIN,
        "legacy_shape_nonzero_on_domain": legacy_shape_nonzero,
        "legacy_denominator_positive_on_domain": denominator_positive,
        "cos_zero_congruence_has_no_domain_solution": zero_congruence_solutions == [],
        "alpha_geo_positive": ALPHA_GEO > 0,
        "alpha_inverse_normalization_residual_zero_formally": True,
        "alpha_inverse_normalization_residual_max_abs_float": max_abs_residual,
        "ratio_K_over_shape_constant_alpha_geo_formally": True,
        "ratio_K_over_shape_max_float_deviation": max_ratio_deviation,
        "sign_pattern_preserved_by_positive_alpha": all(row["sign_preserved_by_positive_alpha"] for row in rows),
        "necessity_exact_apd_inherited": exact_apd,
        "component_gap_amplitude_row_present": amplitude_row["component"] == "amplitude_normalization",
        "scalar_normalization_witness_exported": True,
        "full_strict_A_d_derivation_exported": False,
        "legacy_role_transfer_allowed": False,
        "raw_kernel_identity_claimed": False,
        "full_bridge_theorem_exported": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_AMPLITUDE_SCALAR_NORMALIZATION_CERTIFICATE__NO_ROLE_TRANSFER",
        "status": "legacy-alpha-geo-factors-out-as-positive-global-scalar__bridge-row-witness-not-role-transfer",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "amplitude normalization",
                "alpha_geo normalization",
                "alpha_geo strict derived",
                "amplitude nonabsorption",
                "legacy strict amplitude",
                "strict_amplitude_normalization",
            ],
            "finding": "Older amplitude-nonabsorption packets discuss legacy physical-role nontransfer and missing full A_abs maps; this certificate is narrower and only exports the scalar alpha_geo^{-1} normalization witness for K_legacy_ont on the finite bridge domain.",
        },
        "normalization_definition": {
            "legacy_kernel": "K_legacy_ont(d)=alpha_geo*cos(pi*d/4+pi/6)/(1+beta_tors*d)",
            "legacy_shape_after_scalar_normalization": "L_shape(d)=alpha_geo^{-1}*K_legacy_ont(d)=cos(pi*d/4+pi/6)/(1+beta_tors*d)",
            "alpha_geo": "4 ln(2)",
            "beta_tors": BETA_TORS,
            "domain": NODE_DOMAIN,
            "cos_zero_congruence": "cos((3d+2)pi/12)=0 would require 3d+2 == 6 mod 12, i.e. 3d == 4 mod 12, impossible because gcd(3,12)=3 does not divide 4",
        },
        "node_rows": rows,
        "amplitude_scalar_normalization_summary": summary,
        "cross_checks": {
            "guardrail_requires_completion_map": guardrail["legacy_kernel_intermediate_bridge_summary"]["legacy_kernel_restored_as_intermediate"] and guardrail["legacy_kernel_intermediate_bridge_summary"]["strict_kernel_treated_as_completed_legacy_continuation"],
            "component_gap_amplitude_row_open_but_finite": amplitude_row["matrix_bits"]["finite_certificate_exported"] and not amplitude_row["matrix_bits"]["role_transfer_allowed_now"],
            "alpha_source_strict_derived_value_loaded": summary["strict_alpha_geo_source_loaded"],
            "scalar_normalization_formally_exact": summary["alpha_inverse_normalization_residual_zero_formally"] and summary["ratio_K_over_shape_constant_alpha_geo_formally"],
            "finite_domain_safe_for_ratio": summary["legacy_shape_nonzero_on_domain"] and summary["legacy_denominator_positive_on_domain"] and summary["cos_zero_congruence_has_no_domain_solution"],
            "signs_preserved": summary["sign_pattern_preserved_by_positive_alpha"],
            "no_role_transfer_or_bridge_claim": not summary["legacy_role_transfer_allowed"] and not summary["raw_kernel_identity_claimed"] and not summary["full_bridge_theorem_exported"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to distinguish this scalar-normalization bridge witness from older amplitude nonabsorption/role-transfer packets.",
            "factorization_step": "For every audited d, K_legacy_ont(d)=alpha_geo*L_shape(d) by direct algebra; applying alpha_geo^{-1} gives exactly L_shape(d).",
            "nonzero_step": "On d=0..11, cos(pi*d/4+pi/6)=cos((3d+2)pi/12) cannot vanish because 3d==4 mod 12 has no solution; the denominator 1+0.01*d is positive.",
            "positivity_step": "alpha_geo=4 ln(2)>0, so scalar normalization preserves the legacy sign pattern and cannot supply the strict GF(2) phase flips.",
            "strict_source_step": "The repo already exports alpha_geo_strict_derived_v1=4 ln(2), but this certificate uses it only as scalar provenance and not as a selector or physical-role transfer theorem.",
            "bridge_meaning_step": "The amplitude row now has an explicit scalar-normalization witness; the remaining bridge still needs strict A(d) source, phase/frequency transport source, nonlinear damping/compression source, selector source, and role-transfer audit.",
            "theoretical_limit": "This is not a full A_abs absorption theorem, not a full legacy->strict bridge, and not a transfer of sin^2(theta_W), alpha_EM, beta^N, or beta_tors->chi_11 roles.",
        },
        "hard_limits": [
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No full strict A(d) dynamical derivation is exported.",
            "No legacy alpha_geo physical-role transfer is licensed.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-to-strict amplitude scalar-normalization certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The visible legacy scalar `alpha_geo=4 ln(2)` factors out of `K_legacy_ont`",
        "as a positive global amplitude scalar on the audited finite domain.  This is",
        "a scalar-normalization witness only, not a physical-role transfer theorem.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["amplitude_scalar_normalization_summary"].items():
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
    print(json.dumps(payload["amplitude_scalar_normalization_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
