#!/usr/bin/env python3
"""Scratch probe: completed strict-kernel factorization certificate.

Current working reading requested by the author: the live kernel is the strict
kernel, while the legacy kernel should be read as the historical/nadsoliton-
characteristic carrier that becomes completed by strict compression and the
missing orientation bit.  This probe makes that statement algebraic rather than
rhetorical: it factors the strict kernel over the legacy kernel as

    K_strict(d) = C_completion(d) * K_legacy_reduced(d)

on the nonzero legacy cosine domain, with explicit factor components for
amplitude removal, phase/frequency transport, and damping/compression
renormalization.  The factorization is exact pointwise where defined; what is
*not* derived is the internal theorem producing C_completion from nadsoliton
dynamics or beta_tors -> chi_11.

No false pass: the strict kernel is the current/live full form; legacy is not
used as the live kernel.  The completion factor is a certificate/target object,
not a bridge theorem, not a role-transfer theorem, and not QW-2191/ToE closure.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
OUT_JSON = HERE / "bridge_completed_strict_kernel_factorization_certificate_report.json"
OUT_MD = HERE / "bridge_completed_strict_kernel_factorization_certificate_report.md"
REACTIVATION = HERE / "bridge_legacy_kernel_reactivation_from_diagrams_candidate_report.json"
REYNOLDS = HERE / "bridge_strict_alpha_reynolds_annihilator_chi11_matrix_certificate_report.json"

ALPHA_GEO = 4.0 * math.log(2.0)
LEGACY = {
    "omega": math.pi / 4.0,
    "phi": math.pi / 6.0,
    "beta_tors": 0.01,
}
STRICT = {
    "omega": 0.18575,
    "phi": 0.16250,
    "beta": 1.0,
    "eta": 9.0 / 5.0,
}
DOMAIN = list(range(1, 12))
TARGET_Q_POWER = "256/243"
TARGET_ETA = "9/5"
BALANCED_LEDGER = [2, 2, 2, 1, 1]


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def legacy_reduced(d: float) -> float:
    return math.cos(LEGACY["omega"] * d + LEGACY["phi"]) / (1.0 + LEGACY["beta_tors"] * d)


def legacy_full(d: float) -> float:
    return ALPHA_GEO * legacy_reduced(d)


def strict_kernel(d: float) -> float:
    return math.cos(STRICT["omega"] * d + STRICT["phi"]) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])


def completion_factors(d: float) -> dict[str, float]:
    amplitude_factor = 1.0 / ALPHA_GEO
    phase_factor = math.cos(STRICT["omega"] * d + STRICT["phi"]) / math.cos(LEGACY["omega"] * d + LEGACY["phi"])
    damping_factor = (1.0 + LEGACY["beta_tors"] * d) / (1.0 + STRICT["beta"] * d ** STRICT["eta"])
    completion_over_legacy_full = amplitude_factor * phase_factor * damping_factor
    completion_over_legacy_reduced = phase_factor * damping_factor
    return {
        "amplitude_factor_alpha_removal": amplitude_factor,
        "phase_frequency_transport_factor": phase_factor,
        "damping_compression_factor": damping_factor,
        "completion_factor_over_legacy_full": completion_over_legacy_full,
        "completion_factor_over_legacy_reduced": completion_over_legacy_reduced,
    }


def factorization_rows() -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for d in DOMAIN:
        factors = completion_factors(float(d))
        k_legacy_full = legacy_full(float(d))
        k_legacy_reduced = legacy_reduced(float(d))
        k_strict = strict_kernel(float(d))
        reconstructed_from_full = factors["completion_factor_over_legacy_full"] * k_legacy_full
        reconstructed_from_reduced = factors["completion_factor_over_legacy_reduced"] * k_legacy_reduced
        rows.append(
            {
                "d": d,
                "legacy_full": k_legacy_full,
                "legacy_reduced_no_alpha": k_legacy_reduced,
                "strict_kernel": k_strict,
                **factors,
                "reconstructed_from_legacy_full": reconstructed_from_full,
                "reconstructed_from_legacy_reduced": reconstructed_from_reduced,
                "residual_from_legacy_full": reconstructed_from_full - k_strict,
                "residual_from_legacy_reduced": reconstructed_from_reduced - k_strict,
            }
        )
    return rows


def sign_pattern(values: list[float]) -> list[int]:
    return [1 if value > 0 else -1 if value < 0 else 0 for value in values]


def monotone_decreasing(values: list[float]) -> bool:
    return all(left > right for left, right in zip(values, values[1:]))


def build_payload() -> dict[str, Any]:
    reactivation = load_json(REACTIVATION)
    reynolds = load_json(REYNOLDS)
    rows = factorization_rows()
    residuals_full = [abs(row["residual_from_legacy_full"]) for row in rows]
    residuals_reduced = [abs(row["residual_from_legacy_reduced"]) for row in rows]
    damping = [row["damping_compression_factor"] for row in rows]
    phase = [row["phase_frequency_transport_factor"] for row in rows]
    completion_full = [row["completion_factor_over_legacy_full"] for row in rows]

    return {
        "result_kind": "SCRATCH_COMPLETED_STRICT_KERNEL_FACTORIZATION_CERTIFICATE__EXACT_FACTOR_TARGET_NOT_BRIDGE_THEOREM",
        "status": "strict-kernel-current-full-form-with-explicit-legacy-completion-factorization",
        "kernel_reading": {
            "current_live_kernel": "K_strict_gate(d)=cos(omega_S*d+phi_S)/(1+beta*d^eta)",
            "legacy_reading": "K_legacy_ont is the historical/nadsoliton-characteristic carrier, not the live current kernel.",
            "completion_reading": "K_strict_gate is represented as legacy carrier plus explicit completion factors for alpha removal, phase/frequency transport, damping/compression, and the still-open orientation bit lane.",
            "guarded_identity": "Pointwise factorization is exact where cos(omega_L*d+phi_L) != 0; it is not yet an internal theorem deriving the factors from nadsoliton dynamics.",
        },
        "constants": {
            "alpha_geo": ALPHA_GEO,
            "legacy": LEGACY,
            "strict": STRICT,
            "target_q_power": TARGET_Q_POWER,
            "target_eta": TARGET_ETA,
            "balanced_ledger": BALANCED_LEDGER,
        },
        "factorization_summary": {
            "domain_d_values": DOMAIN,
            "max_abs_residual_from_legacy_full": max(residuals_full),
            "max_abs_residual_from_legacy_reduced": max(residuals_reduced),
            "residual_tolerance_pass": max(residuals_full + residuals_reduced) < 1e-15,
            "damping_compression_factor_positive": all(value > 0 for value in damping),
            "damping_compression_factor_strictly_decreasing": monotone_decreasing(damping),
            "damping_factor_d1": damping[0],
            "damping_factor_d11": damping[-1],
            "damping_factor_d1_over_d11": damping[0] / damping[-1],
            "phase_transport_sign_pattern": sign_pattern(phase),
            "completion_factor_over_legacy_full_sign_pattern": sign_pattern(completion_full),
        },
        "factorization_rows": rows,
        "upstream_bridge_context": {
            "reactivation_status": reactivation["status"],
            "one_bit_frontier": reactivation["one_bit_frontier"],
            "reynolds_full_aut_obstruction": reynolds["interpretation"]["honest_negative"],
            "reynolds_annihilator_zero": reynolds["matrix_certificate"]["reynolds_times_chi11_numerator_is_zero_matrix"],
        },
        "completion_factor_interpretation": {
            "amplitude": "1/alpha_geo removes the legacy explicit amplitude; role-safe use still requires bridge discipline.",
            "phase_frequency": "cos(strict phase)/cos(legacy phase) is the exact pointwise transport from legacy torsion/resonance phase to strict gate phase on the sampled domain.",
            "damping_compression": "((1+beta_tors*d)/(1+beta*d^eta)) is the exact pointwise compression factor upgrading legacy hyperbolic damping into strict eta damping.",
            "orientation_bit": "The factorization does not supply beta_tors -> chi_11; that remains the one-bit theorem-target from the reactivation audit.",
        },
        "exact_proof_certificate": {
            "factor_identity": "For every listed d, K_strict(d) - C_full(d)*K_legacy_full(d) = 0 up to floating arithmetic, and similarly for the reduced no-alpha carrier.",
            "completion_objects": "C_full=(1/alpha_geo)*phase_transport*damping_compression and C_reduced=phase_transport*damping_compression are explicit completion factors.",
            "computational_positive": "The current full kernel can be represented exactly as strict completion of the legacy carrier on the finite d=1..11 audit domain.",
            "theoretical_limit": "This is an exact factor certificate, not a derivation of the factors, not the beta_tors->chi_11 theorem, and not QW-2191 closure.",
        },
        "ontology_guardrail": {
            "allowed_reading": "The nadsoliton itself remains the primordial information in a solitonic state; legacy kernel data is read as a historical internal characteristic carrier completed by the strict kernel.",
            "forbidden_reading": "No separate informational layer underneath the nadsoliton is introduced.",
        },
        "hard_limits": [
            "No claim that the legacy kernel is the current live kernel; K_strict_gate is the current full form.",
            "No unqualified identity K_legacy_ont == K_strict_gate is asserted; only explicit factorization/completion is certified.",
            "No beta_tors -> chi_11 theorem is asserted.",
            "No proof derives the completion factors from strict nadsoliton dynamics yet.",
            "No legacy physical-role transfer onto K_strict_gate is used without an explicit bridge theorem.",
            "No QW-2191 discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> str:
    summary = payload["factorization_summary"]
    reading = payload["kernel_reading"]
    lines = [
        "# Completed strict-kernel factorization certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Kernel reading",
        "",
        f"- Current live kernel: {reading['current_live_kernel']}",
        f"- Legacy reading: {reading['legacy_reading']}",
        f"- Completion reading: {reading['completion_reading']}",
        f"- Guarded identity: {reading['guarded_identity']}",
        "",
        "## Factorization summary",
        "",
        f"- Domain d values: `{summary['domain_d_values']}`",
        f"- Max abs residual from legacy full: `{summary['max_abs_residual_from_legacy_full']:.3e}`",
        f"- Max abs residual from legacy reduced: `{summary['max_abs_residual_from_legacy_reduced']:.3e}`",
        f"- Residual tolerance pass: `{summary['residual_tolerance_pass']}`",
        f"- Damping compression positive/decreasing: `{summary['damping_compression_factor_positive']}` / `{summary['damping_compression_factor_strictly_decreasing']}`",
        f"- Damping factor d=1 to d=11: `{summary['damping_factor_d1']:.12f}` -> `{summary['damping_factor_d11']:.12f}`",
        f"- Damping factor ratio d1/d11: `{summary['damping_factor_d1_over_d11']:.6f}`",
        f"- Phase transport sign pattern: `{summary['phase_transport_sign_pattern']}`",
        "",
        "## Completion factor interpretation",
        "",
    ]
    for key, value in payload["completion_factor_interpretation"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Proof certificate", ""])
    for key, value in payload["exact_proof_certificate"].items():
        lines.append(f"- `{key}`: {value}")
    lines.extend(["", "## Hard limits", ""])
    lines.extend(f"- {item}" for item in payload["hard_limits"])
    lines.append("")
    return "\n".join(lines)


def main() -> None:
    payload = build_payload()
    OUT_JSON.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    OUT_MD.write_text(write_markdown(payload), encoding="utf-8")
    print(json.dumps(payload, indent=2, sort_keys=True))
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
