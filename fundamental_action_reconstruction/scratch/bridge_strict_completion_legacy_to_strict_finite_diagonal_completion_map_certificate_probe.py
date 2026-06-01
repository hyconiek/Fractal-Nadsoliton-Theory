#!/usr/bin/env python3
"""Scratch probe: finite diagonal legacy -> strict completion map certificate.

This is a deliberately narrow bridge-completion step.  On the audited finite
Z12 domain, the legacy vector K_L(d)=K_legacy_ont(d) has no zero components, so
there is a unique diagonal operator Q with

    Q * K_L = K_strict.

The entries are the pointwise quotients q_d=K_strict(d)/K_L(d), and the existing
necessity certificate already decomposes q_d as A(d)*P(d)*D(d).  This probe
packages that fact as an explicit finite diagonal completion-map witness and
records what it does *not* solve: it is not a strict dynamical derivation of the
factors, not a role-transfer theorem, not beta_tors -> chi_11, and not ToE
closure.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.md"

SOURCE_REPORTS = {
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "amplitude_scalar_normalization": HERE / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json",
    "damping_compression_separation": HERE / "bridge_strict_completion_legacy_to_strict_damping_compression_separation_certificate_report.json",
    "positive_factor_sign_separation": HERE / "bridge_strict_completion_positive_factor_sign_separation_certificate_report.json",
    "component_gap_matrix": HERE / "bridge_strict_completion_legacy_to_strict_completion_component_gap_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

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


def gf2_rank(rows: list[list[int]]) -> int:
    work = [row[:] for row in rows if any(row)]
    if not work:
        return 0
    rank = 0
    col_count = len(work[0])
    for col in range(col_count):
        pivot = next((idx for idx in range(rank, len(work)) if work[idx][col] & 1), None)
        if pivot is None:
            continue
        work[rank], work[pivot] = work[pivot], work[rank]
        for idx, row in enumerate(work):
            if idx != rank and row[col] & 1:
                work[idx] = [a ^ b for a, b in zip(row, work[rank])]
        rank += 1
        if rank == len(work):
            break
    return rank


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    necessity = loaded["necessity"]
    amplitude = loaded["amplitude_scalar_normalization"]
    damping = loaded["damping_compression_separation"]
    positive_sign = loaded["positive_factor_sign_separation"]
    component_gap = loaded["component_gap_matrix"]
    guardrail = loaded["legacy_bridge_guardrail"]

    point_rows = necessity["pointwise_quotient_certificate"]
    diagonal_entries = [row["strict_over_legacy_quotient"] for row in point_rows]
    legacy_vector = [row["legacy_full"] for row in point_rows]
    strict_vector = [row["strict_kernel"] for row in point_rows]
    reconstructed = [q * legacy for q, legacy in zip(diagonal_entries, legacy_vector)]
    residuals = [got - want for got, want in zip(reconstructed, strict_vector)]
    factor_residuals = [row["quotient_minus_factor_product"] for row in point_rows]
    diagonal_sign_bits = [0 if sign(q) > 0 else 1 for q in diagonal_entries]

    diagonal_matrix_rows = [
        [diagonal_entries[row] if row == col else 0.0 for col in range(len(diagonal_entries))]
        for row in range(len(diagonal_entries))
    ]
    support_matrix_gf2 = [
        [1 if row == col and abs(diagonal_entries[row]) > TOL else 0 for col in range(len(diagonal_entries))]
        for row in range(len(diagonal_entries))
    ]
    determinant_nonzero_witness = math.prod(diagonal_entries)
    quotient_first = diagonal_entries[0]
    scalar_only_residuals = [quotient_first * legacy - strict for legacy, strict in zip(legacy_vector, strict_vector)]

    operator_rows = []
    for row, reconstructed_value, residual, sign_bit in zip(point_rows, reconstructed, residuals, diagonal_sign_bits):
        operator_rows.append(
            {
                "d": row["d"],
                "legacy_value": row["legacy_full"],
                "strict_value": row["strict_kernel"],
                "diagonal_q_d": row["strict_over_legacy_quotient"],
                "A_factor": row["alpha_normalization_factor"],
                "P_factor": row["phase_frequency_transport_factor"],
                "D_factor": row["damping_compression_factor"],
                "APD_product": row["factor_product"],
                "q_minus_APD": row["quotient_minus_factor_product"],
                "reconstructed_strict_value": reconstructed_value,
                "reconstruction_residual": residual,
                "diagonal_sign_bit": sign_bit,
                "legacy_sign": row["legacy_sign"],
                "strict_sign": row["strict_sign"],
            }
        )

    summary = {
        "domain_size": len(point_rows),
        "legacy_vector_has_no_zero_components": all(abs(value) > TOL for value in legacy_vector),
        "strict_vector_has_no_zero_components": all(abs(value) > TOL for value in strict_vector),
        "unique_diagonal_completion_map_exists": all(abs(value) > TOL for value in legacy_vector),
        "diagonal_support_rank_over_gf2": gf2_rank(support_matrix_gf2),
        "diagonal_operator_full_rank_on_finite_domain": gf2_rank(support_matrix_gf2) == len(point_rows),
        "diagonal_determinant_nonzero_witness_float": determinant_nonzero_witness,
        "max_abs_reconstruction_residual": max(abs(value) for value in residuals),
        "max_abs_q_minus_apd_residual": max(abs(value) for value in factor_residuals),
        "apd_factorization_inherited_exact": max(abs(value) for value in factor_residuals) <= TOL,
        "scalar_only_completion_fails": max(abs(value) for value in scalar_only_residuals) > 1e-3,
        "scalar_only_max_abs_residual_using_q0": max(abs(value) for value in scalar_only_residuals),
        "diagonal_sign_bits": diagonal_sign_bits,
        "diagonal_sign_bits_match_positive_factor_phase_bits": diagonal_sign_bits == [row["phase_z2_bit"] for row in positive_sign["node_sign_separation_rows"]],
        "amplitude_scalar_normalization_inherited": amplitude["amplitude_scalar_normalization_summary"]["scalar_normalization_witness_exported"],
        "damping_linear_to_nonlinear_gap_still_open": not damping["separation_summary"]["full_bridge_theorem_exported"],
        "component_gap_matrix_still_not_full_bridge": component_gap["completion_gap_summary"]["completion_map_partial_not_full"],
        "role_transfer_allowed": False,
        "strict_dynamic_source_exported": False,
        "raw_identity_claimed": False,
        "full_bridge_theorem_exported": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_FINITE_DIAGONAL_COMPLETION_MAP_CERTIFICATE__NO_SOURCE_OR_ROLE_TRANSFER",
        "status": "unique-finite-diagonal-completion-map-exported-on-z12__apd-factorization-inherited__no-source-theorem",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "diagonal completion map",
                "unique diagonal operator",
                "strict_over_legacy_quotient",
                "K_strict/K_legacy",
                "finite completion map",
                "pointwise quotient",
            ],
            "finding": "Existing necessity/sign reports contain pointwise quotients and sign separation, but no dedicated report packaged the unique finite diagonal operator Q:K_legacy_vector -> K_strict_vector with rank/determinant/no-scalar-only checks before this file.",
        },
        "finite_operator_definition": {
            "domain": [row["d"] for row in point_rows],
            "legacy_vector_name": "K_L(d)=K_legacy_ont(d)",
            "strict_vector_name": "K_S(d)=K_strict_gate(d)",
            "operator": "Q=diag(q_d), q_d=K_S(d)/K_L(d)",
            "factorization": "q_d=A(d)*P(d)*D(d)",
            "diagonal_matrix": diagonal_matrix_rows,
            "support_matrix_gf2": support_matrix_gf2,
        },
        "operator_rows": operator_rows,
        "finite_diagonal_completion_summary": summary,
        "cross_checks": {
            "guardrail_restores_legacy_as_intermediate": guardrail["legacy_kernel_intermediate_bridge_summary"]["legacy_kernel_restored_as_intermediate"],
            "operator_exists_and_reconstructs": summary["unique_diagonal_completion_map_exists"] and summary["max_abs_reconstruction_residual"] <= TOL,
            "operator_full_rank": summary["diagonal_operator_full_rank_on_finite_domain"],
            "apd_factorization_matches_necessity": summary["apd_factorization_inherited_exact"] and necessity["necessity_summary"]["exact_subsets_without_extra_scalar"] == ["alpha_normalization+phase_frequency_transport+damping_compression"],
            "scalar_only_rejected": summary["scalar_only_completion_fails"],
            "sign_bits_match_phase_certificate": summary["diagonal_sign_bits_match_positive_factor_phase_bits"],
            "no_source_role_or_bridge_claim": not summary["role_transfer_allowed"] and not summary["strict_dynamic_source_exported"] and not summary["raw_identity_claimed"] and not summary["full_bridge_theorem_exported"],
        },
        "proof_certificate": {
            "grep_step": "rg was used to distinguish this finite diagonal completion-map packaging from the existing necessity, sign-separation, amplitude, and damping reports.",
            "existence_step": "Because every audited legacy value K_L(d) is nonzero, q_d=K_S(d)/K_L(d) is defined for every d=0..11 and Q=diag(q_d) satisfies Q*K_L=K_S.",
            "uniqueness_step": "For a diagonal operator, each row equation q_d*K_L(d)=K_S(d) has the unique solution q_d=K_S(d)/K_L(d) because K_L(d) is nonzero.",
            "rank_step": "Every q_d is nonzero, so the diagonal support matrix has GF(2) rank 12 and the finite diagonal operator is full-rank on the audited coordinate domain.",
            "factorization_step": "The necessity report gives q_d-A(d)P(d)D(d) residuals with max absolute value below tolerance, so the diagonal bridge witness inherits the exact A/P/D factorization.",
            "non_scalar_step": "Using q_0 as a single scalar multiplier leaves a nonzero residual, so the finite bridge is not a scalar-only normalization map.",
            "theoretical_limit": "This is an exact finite-domain diagonal completion witness, not a strict dynamical derivation of A/P/D, not beta_tors->chi_11, not legacy physical-role transfer, and not ToE closure.",
        },
        "hard_limits": [
            "No raw identity K_legacy_ont == K_strict_gate is claimed.",
            "No strict dynamical source for A(d), P(d), D(d), omega/phi, beta/eta, or the transport cocycle is exported.",
            "No beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-to-strict finite diagonal completion map certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "On the audited finite `Z12` domain, the nonzero legacy vector admits a unique",
        "diagonal operator `Q=diag(K_strict(d)/K_legacy(d))` mapping it to the strict",
        "vector.  This is an exact finite bridge witness, not a source or role-transfer theorem.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["finite_diagonal_completion_summary"].items():
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
    print(json.dumps(payload["finite_diagonal_completion_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
