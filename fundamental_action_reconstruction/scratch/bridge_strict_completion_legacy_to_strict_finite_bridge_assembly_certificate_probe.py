#!/usr/bin/env python3
"""Scratch probe: finite legacy -> strict bridge assembly certificate.

This is the next narrow bridge-completion step after the separate amplitude,
phase/frequency, damping, and diagonal certificates.  It assembles the three
explicit component maps

    A(d)=1/alpha_geo,
    P(d)=cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L),
    D(d)=(1+beta_tors*d)/(1+beta*d^eta)

into one finite kernel-comparison witness

    K_strict_gate(d) = K_legacy_ont(d) * A(d) * P(d) * D(d)

on the audited Z12 domain.  The output is intentionally finite-scope: it is a
comparison witness for the restored intermediate bridge kernel, not a strict
source theorem, not a legacy physical-role transfer theorem, and not ToE
closure.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.md"
T99_SPEC = ROOT / "fundamental_action_reconstruction" / "T99_CURRENT_LEGACY_TO_STRICT_BRIDGE_CLOSURE_WITNESS_TARGET_SPEC.md"

SOURCE_REPORTS = {
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "amplitude_scalar_normalization": HERE / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json",
    "finite_diagonal_completion_map": HERE / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.json",
    "phase_frequency_affine_transport": HERE / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.json",
    "damping_parameter_identifiability": HERE / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

TOL = 1e-12


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
    raise ValueError("zero sign cannot be represented as a GF(2) sign bit")


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    necessity = loaded["necessity"]
    amplitude = loaded["amplitude_scalar_normalization"]
    diagonal = loaded["finite_diagonal_completion_map"]
    phase = loaded["phase_frequency_affine_transport"]
    damping = loaded["damping_parameter_identifiability"]
    component_gap = loaded["component_gap_matrix"]
    guardrail = loaded["legacy_bridge_guardrail"]
    t99_text = T99_SPEC.read_text(encoding="utf-8")

    constants = necessity["constants"]
    alpha_geo = constants["alpha_geo"]
    legacy_constants = constants["legacy"]
    strict_constants = constants["strict"]
    omega_l = legacy_constants["omega"]
    phi_l = legacy_constants["phi"]
    beta_tors = legacy_constants["beta_tors"]
    omega_s = strict_constants["omega"]
    phi_s = strict_constants["phi"]
    beta = strict_constants["beta"]
    eta = strict_constants["eta"]

    diagonal_by_d = {
        row["d"]: row
        for row in diagonal["operator_rows"]
    }
    phase_bits = phase["phase_frequency_affine_transport_summary"]["phase_factor_bits"]

    assembly_rows = []
    reconstruction_residuals = []
    q_residuals = []
    necessity_residuals = []
    log_residuals = []
    diagonal_sign_bits = []
    phase_sign_bits = []

    for necessity_row in necessity["pointwise_quotient_certificate"]:
        d = necessity_row["d"]
        legacy_cos = math.cos(omega_l * d + phi_l)
        strict_cos = math.cos(omega_s * d + phi_s)
        legacy_denominator = 1.0 + beta_tors * d
        strict_denominator = 1.0 + beta * (d ** eta)
        legacy_value = alpha_geo * legacy_cos / legacy_denominator
        strict_value = strict_cos / strict_denominator
        a_factor = 1.0 / alpha_geo
        p_factor = strict_cos / legacy_cos
        d_factor = legacy_denominator / strict_denominator
        assembled_q = a_factor * p_factor * d_factor
        diagonal_q = diagonal_by_d[d]["diagonal_q_d"]
        reconstructed_strict = legacy_value * assembled_q
        reconstruction_residual = reconstructed_strict - strict_value
        q_residual = assembled_q - diagonal_q
        necessity_residual = assembled_q - necessity_row["factor_product"]
        log_residual = math.log(abs(assembled_q)) - (math.log(a_factor) + math.log(abs(p_factor)) + math.log(d_factor))
        phase_bit = bit_from_sign(sign(p_factor))
        diagonal_bit = bit_from_sign(sign(assembled_q))

        reconstruction_residuals.append(abs(reconstruction_residual))
        q_residuals.append(abs(q_residual))
        necessity_residuals.append(abs(necessity_residual))
        log_residuals.append(abs(log_residual))
        phase_sign_bits.append(phase_bit)
        diagonal_sign_bits.append(diagonal_bit)

        assembly_rows.append(
            {
                "d": d,
                "legacy_kernel_value_from_formula": legacy_value,
                "strict_kernel_value_from_formula": strict_value,
                "A_amplitude_normalization": a_factor,
                "P_phase_frequency_transport": p_factor,
                "D_damping_compression": d_factor,
                "assembled_Q_APD": assembled_q,
                "diagonal_Q_reported": diagonal_q,
                "assembled_Q_minus_diagonal_Q": q_residual,
                "assembled_Q_minus_necessity_factor_product": necessity_residual,
                "reconstructed_strict_value": reconstructed_strict,
                "reconstruction_residual": reconstruction_residual,
                "log_abs_Q_minus_log_abs_A_plus_log_abs_P_plus_log_abs_D": log_residual,
                "phase_sign_bit": phase_bit,
                "diagonal_sign_bit": diagonal_bit,
                "phase_bit_matches_prior_affine_report": phase_bit == phase_bits[d],
            }
        )

    summary = {
        "domain_size": len(assembly_rows),
        "finite_kernel_comparison_witness_exported": True,
        "t99_positive_bridge_closure_target_detected": "Omega_legacy_strict_bridge_closure_witness_target_v1" in t99_text,
        "comparison_scope_only": "not closure of strict-core selector theory" in t99_text,
        "amplitude_scalar_normalization_inherited": amplitude["amplitude_scalar_normalization_summary"]["scalar_normalization_witness_exported"],
        "phase_affine_transport_inherited": phase["phase_frequency_affine_transport_summary"]["continuous_affine_phase_transport_exact"],
        "damping_beta_eta_identifiability_inherited": damping["damping_parameter_identifiability_summary"]["candidate_grid_unique_match"],
        "finite_diagonal_completion_map_inherited": diagonal["finite_diagonal_completion_summary"]["unique_diagonal_completion_map_exists"],
        "max_abs_reconstruction_residual": max(reconstruction_residuals),
        "max_abs_assembled_q_minus_diagonal_q": max(q_residuals),
        "max_abs_assembled_q_minus_necessity_factor_product": max(necessity_residuals),
        "max_abs_log_additive_residual": max(log_residuals),
        "assembled_map_reconstructs_strict_kernel": max(reconstruction_residuals) <= TOL,
        "assembled_map_matches_finite_diagonal_certificate": max(q_residuals) <= TOL,
        "assembled_map_matches_necessity_apd_product": max(necessity_residuals) <= TOL,
        "log_abs_decomposition_additive": max(log_residuals) <= TOL,
        "phase_sign_bits_match_diagonal_sign_bits": phase_sign_bits == diagonal_sign_bits,
        "phase_sign_bits_match_prior_affine_report": all(row["phase_bit_matches_prior_affine_report"] for row in assembly_rows),
        "component_gap_still_not_full_bridge": component_gap["completion_gap_summary"]["completion_map_partial_not_full"],
        "guardrail_role_transfer_required_after_full_bridge": guardrail["legacy_kernel_intermediate_bridge_summary"]["role_transfer_audit_required_after_full_bridge"],
        "strict_dynamic_source_exported": False,
        "legacy_role_transfer_exported": False,
        "selector_source_exported": False,
        "full_bridge_theorem_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_FINITE_BRIDGE_ASSEMBLY_CERTIFICATE__COMPARISON_SCOPE_ONLY",
        "status": "finite-apd-bridge-assembly-reconstructs-strict-kernel-on-z12__source-selector-role-transfer-open",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "source_specs": {"t99_bridge_closure_witness_target_spec": rel(T99_SPEC)},
        "grep_disambiguation": {
            "searched_terms": [
                "legacy_to_strict finite bridge assembly",
                "Omega_legacy_strict_bridge_closure_witness_target_v1",
                "A(d)*P(d)*D(d)",
                "assembled bridge witness",
                "comparison scope only",
            ],
            "finding": "Prior reports exported the separate amplitude, phase, damping, and diagonal witnesses; this report assembles them into one finite kernel-comparison witness while preserving source/selector/role-transfer hard limits.",
        },
        "assembly_definition": {
            "legacy_kernel": "K_legacy_ont(d)=alpha_geo*cos(omega_L*d+phi_L)/(1+beta_tors*d)",
            "strict_kernel": "K_strict_gate(d)=cos(omega_S*d+phi_S)/(1+beta*d^eta)",
            "assembled_completion_factor": "Q_assembly(d)=A(d)*P(d)*D(d)",
            "A_factor": "1/alpha_geo",
            "P_factor": "cos(omega_S*d+phi_S)/cos(omega_L*d+phi_L)",
            "D_factor": "(1+beta_tors*d)/(1+beta*d^eta)",
            "finite_identity_checked": "K_strict_gate(d)=K_legacy_ont(d)*Q_assembly(d) for d=0..11",
            "scope": "finite kernel-comparison witness, not strict dynamical source theorem",
        },
        "assembly_rows": assembly_rows,
        "finite_bridge_assembly_summary": summary,
        "cross_checks": {
            "source_reports_present": set(loaded) == set(SOURCE_REPORTS),
            "t99_target_scope_detected": summary["t99_positive_bridge_closure_target_detected"] and summary["comparison_scope_only"],
            "component_certificates_inherited": summary["amplitude_scalar_normalization_inherited"] and summary["phase_affine_transport_inherited"] and summary["damping_beta_eta_identifiability_inherited"] and summary["finite_diagonal_completion_map_inherited"],
            "assembled_identity_exact_on_finite_domain": summary["assembled_map_reconstructs_strict_kernel"] and summary["assembled_map_matches_finite_diagonal_certificate"] and summary["assembled_map_matches_necessity_apd_product"],
            "sign_and_log_decomposition_coherent": summary["phase_sign_bits_match_diagonal_sign_bits"] and summary["phase_sign_bits_match_prior_affine_report"] and summary["log_abs_decomposition_additive"],
            "guardrail_limits_preserved": summary["component_gap_still_not_full_bridge"] and summary["guardrail_role_transfer_required_after_full_bridge"] and not summary["strict_dynamic_source_exported"] and not summary["legacy_role_transfer_exported"] and not summary["selector_source_exported"] and not summary["full_bridge_theorem_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used before adding this probe to avoid duplicating the separate amplitude, phase-affine, damping-parameter, and diagonal certificates.",
            "assembly_step": "For every audited d, Q_assembly(d)=alpha_geo^{-1} * cos(theta_S(d))/cos(theta_L(d)) * (1+beta_tors*d)/(1+beta*d^eta).",
            "identity_step": "Substituting K_legacy_ont into K_legacy_ont(d)*Q_assembly(d) cancels alpha_geo, cos(theta_L), and the legacy denominator, leaving K_strict_gate(d).",
            "diagonal_step": "The assembled Q matches the previously certified unique diagonal Q=diag(K_strict/K_legacy) on all 12 audited nodes.",
            "log_sign_step": "Since A and D are positive, sign(Q_assembly) is carried by P; log|Q|=log(A)+log|P|+log(D) numerically closes on every audited node.",
            "scope_step": "This is a finite T99-style kernel-comparison witness only; it does not export a strict dynamical source, selector source, role-transfer theorem, QW-2191 discharge, or ToE closure.",
        },
        "hard_limits": [
            "No unqualified raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No strict dynamical derivation of A/P/D, omega/phi, beta/eta, or chi_11 is exported.",
            "No beta_tors -> beta/eta or beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-to-strict finite bridge assembly certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The separate amplitude, phase/frequency, and damping/compression witnesses",
        "assemble into one finite kernel-comparison map `Q_assembly=A*P*D` on Z12.",
        "The assembled map reconstructs `K_strict_gate` from `K_legacy_ont` on the",
        "audited finite domain, but remains below source/selector/role-transfer closure.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["finite_bridge_assembly_summary"].items():
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
    print(json.dumps(payload["finite_bridge_assembly_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
