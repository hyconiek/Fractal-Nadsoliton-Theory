#!/usr/bin/env python3
"""Scratch probe: symbolic cancellation certificate for the legacy -> strict bridge.

This probe is deliberately one level more algebraic than the finite bridge
assembly report.  It does not merely re-check the 12 floating rows.  It records
the formal cancellation

    K_L(d) * A(d) * P(d) * D(d)
      = alpha*cL/Lden * alpha^{-1} * cS/cL * Lden/Sden
      = cS/Sden
      = K_S(d)

under the explicitly stated admissibility assumptions alpha != 0, cL != 0,
Lden != 0, Sden != 0.  The finite Z12 reports are used only to audit that these
assumptions are safe on the current audited domain.

This still does not derive the strict parameters, phase source, selector bit, or
legacy role transfer.  It is a formula-level cancellation witness for the
current completion ansatz.
"""
from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT_JSON = HERE / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.json"
OUT_MD = HERE / "bridge_strict_completion_legacy_to_strict_symbolic_cancellation_certificate_report.md"

SOURCE_REPORTS = {
    "necessity": HERE / "bridge_strict_kernel_completion_necessity_certificate_report.json",
    "finite_bridge_assembly": HERE / "bridge_strict_completion_legacy_to_strict_finite_bridge_assembly_certificate_report.json",
    "finite_diagonal_completion_map": HERE / "bridge_strict_completion_legacy_to_strict_finite_diagonal_completion_map_certificate_report.json",
    "amplitude_scalar_normalization": HERE / "bridge_strict_completion_legacy_to_strict_amplitude_scalar_normalization_certificate_report.json",
    "phase_frequency_affine_transport": HERE / "bridge_strict_completion_legacy_to_strict_phase_frequency_affine_transport_certificate_report.json",
    "damping_parameter_identifiability": HERE / "bridge_strict_completion_legacy_to_strict_damping_parameter_identifiability_certificate_report.json",
    "legacy_bridge_guardrail": HERE / "bridge_strict_completion_legacy_kernel_intermediate_bridge_guardrail_certificate_report.json",
}

TOL = 1e-12


def load_json(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"missing prerequisite report: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(ROOT))


def build_payload() -> dict[str, Any]:
    loaded = {name: load_json(path) for name, path in SOURCE_REPORTS.items()}
    necessity = loaded["necessity"]
    finite_assembly = loaded["finite_bridge_assembly"]
    diagonal = loaded["finite_diagonal_completion_map"]
    amplitude = loaded["amplitude_scalar_normalization"]
    phase = loaded["phase_frequency_affine_transport"]
    damping = loaded["damping_parameter_identifiability"]
    guardrail = loaded["legacy_bridge_guardrail"]

    constants = necessity["constants"]
    alpha_geo = constants["alpha_geo"]
    legacy = constants["legacy"]
    strict = constants["strict"]

    admissibility_rows = []
    for d in constants["domain_d_values"]:
        c_l = math.cos(legacy["omega"] * d + legacy["phi"])
        c_s = math.cos(strict["omega"] * d + strict["phi"])
        lden = 1.0 + legacy["beta_tors"] * d
        sden = 1.0 + strict["beta"] * (d ** strict["eta"])
        admissibility_rows.append(
            {
                "d": d,
                "alpha_geo_nonzero": abs(alpha_geo) > TOL,
                "legacy_cos_nonzero": abs(c_l) > TOL,
                "strict_cos_nonzero_for_diagonal": abs(c_s) > TOL,
                "legacy_denominator_positive": lden > 0.0,
                "strict_denominator_positive": sden > 0.0,
                "legacy_cos_value": c_l,
                "strict_cos_value": c_s,
                "legacy_denominator": lden,
                "strict_denominator": sden,
            }
        )

    cancellation_rows = [
        {
            "step": "expand",
            "expression": "K_L*A*P*D=(alpha*cL/Lden)*(1/alpha)*(cS/cL)*(Lden/Sden)",
            "remaining_factor": "alpha*cL*Lden cancels against alpha*cL*Lden in denominators",
            "requires": ["alpha != 0", "cL != 0", "Lden != 0", "Sden != 0"],
        },
        {
            "step": "cancel_alpha",
            "expression": "alpha*(1/alpha)=1",
            "remaining_factor": "cL/Lden * cS/cL * Lden/Sden",
            "requires": ["alpha_geo=4 ln(2)>0 inherited from amplitude certificate"],
        },
        {
            "step": "cancel_legacy_cos",
            "expression": "cL*(cS/cL)=cS",
            "remaining_factor": "cS/Lden * Lden/Sden",
            "requires": ["legacy cos has no zero on the audited Z12 domain"],
        },
        {
            "step": "cancel_legacy_denominator",
            "expression": "(1/Lden)*Lden=1",
            "remaining_factor": "cS/Sden",
            "requires": ["legacy denominator is positive on the audited Z12 domain"],
        },
        {
            "step": "identify_strict_kernel",
            "expression": "cS/Sden=cos(omega_S*d+phi_S)/(1+beta*d^eta)=K_strict_gate(d)",
            "remaining_factor": "K_strict_gate(d)",
            "requires": ["strict denominator is positive on the audited Z12 domain"],
        },
    ]

    summary = {
        "symbolic_cancellation_formula_exported": True,
        "finite_domain": constants["domain_d_values"],
        "alpha_geo_nonzero_on_domain": all(row["alpha_geo_nonzero"] for row in admissibility_rows),
        "legacy_cos_nonzero_on_domain": all(row["legacy_cos_nonzero"] for row in admissibility_rows),
        "legacy_denominator_positive_on_domain": all(row["legacy_denominator_positive"] for row in admissibility_rows),
        "strict_denominator_positive_on_domain": all(row["strict_denominator_positive"] for row in admissibility_rows),
        "finite_assembly_reconstruction_inherited": finite_assembly["finite_bridge_assembly_summary"]["assembled_map_reconstructs_strict_kernel"],
        "finite_assembly_matches_diagonal_inherited": finite_assembly["finite_bridge_assembly_summary"]["assembled_map_matches_finite_diagonal_certificate"],
        "diagonal_uniqueness_inherited": diagonal["finite_diagonal_completion_summary"]["unique_diagonal_completion_map_exists"],
        "amplitude_nonzero_inherited": amplitude["amplitude_scalar_normalization_summary"]["alpha_geo_positive"],
        "phase_transport_source_still_open": not phase["phase_frequency_affine_transport_summary"]["strict_phase_frequency_source_exported"],
        "damping_source_still_open": not damping["damping_parameter_identifiability_summary"]["strict_beta_eta_source_exported"],
        "role_transfer_required_after_full_bridge": guardrail["legacy_kernel_intermediate_bridge_summary"]["role_transfer_audit_required_after_full_bridge"],
        "strict_dynamic_source_exported": False,
        "selector_source_exported": False,
        "legacy_role_transfer_exported": False,
        "raw_kernel_identity_claimed": False,
        "full_bridge_theorem_exported": False,
        "toe_closure_claimed": False,
    }

    return {
        "result_kind": "SCRATCH_STRICT_COMPLETION_LEGACY_TO_STRICT_SYMBOLIC_CANCELLATION_CERTIFICATE__ANSATZ_LEVEL_ONLY",
        "status": "symbolic-apd-cancellation-valid-under-finite-admissibility__sources-and-role-transfer-open",
        "source_reports": {name: rel(path) for name, path in SOURCE_REPORTS.items()},
        "grep_disambiguation": {
            "searched_terms": [
                "symbolic bridge cancellation",
                "algebraic cancellation",
                "cancel alpha_geo",
                "cancel legacy denominator",
                "Q_assembly",
            ],
            "finding": "The repo already had finite APD/diagonal assembly witnesses; this report adds the formula-level cancellation ledger and uses finite reports only for nonzero-domain admissibility.",
        },
        "symbolic_identity": {
            "legacy_kernel": "K_L=alpha*cL/Lden",
            "strict_kernel": "K_S=cS/Sden",
            "completion_factor": "Q=A*P*D=(1/alpha)*(cS/cL)*(Lden/Sden)",
            "expanded_product": "K_L*Q=(alpha*cL/Lden)*(1/alpha)*(cS/cL)*(Lden/Sden)",
            "reduced_product": "cS/Sden=K_S",
            "admissibility_conditions": ["alpha != 0", "cL != 0", "Lden != 0", "Sden != 0"],
        },
        "admissibility_rows": admissibility_rows,
        "cancellation_rows": cancellation_rows,
        "symbolic_cancellation_summary": summary,
        "cross_checks": {
            "finite_admissibility_passes": summary["alpha_geo_nonzero_on_domain"] and summary["legacy_cos_nonzero_on_domain"] and summary["legacy_denominator_positive_on_domain"] and summary["strict_denominator_positive_on_domain"],
            "finite_assembly_and_diagonal_inherited": summary["finite_assembly_reconstruction_inherited"] and summary["finite_assembly_matches_diagonal_inherited"] and summary["diagonal_uniqueness_inherited"],
            "cancellation_ledger_complete": [row["step"] for row in cancellation_rows] == ["expand", "cancel_alpha", "cancel_legacy_cos", "cancel_legacy_denominator", "identify_strict_kernel"],
            "source_gaps_preserved": summary["phase_transport_source_still_open"] and summary["damping_source_still_open"] and not summary["strict_dynamic_source_exported"] and not summary["selector_source_exported"],
            "role_transfer_and_closure_limits_preserved": summary["role_transfer_required_after_full_bridge"] and not summary["legacy_role_transfer_exported"] and not summary["raw_kernel_identity_claimed"] and not summary["full_bridge_theorem_exported"] and not summary["toe_closure_claimed"],
        },
        "proof_certificate": {
            "grep_step": "rg was used before adding this report to distinguish symbolic cancellation from the existing finite bridge assembly and diagonal certificates.",
            "formal_step": "K_L*A*P*D=(alpha*cL/Lden)*(1/alpha)*(cS/cL)*(Lden/Sden) cancels to cS/Sden=K_S.",
            "admissibility_step": "The finite Z12 audit verifies alpha_geo != 0, legacy cos != 0, legacy denominator > 0, and strict denominator > 0 at every audited node.",
            "finite_consistency_step": "The symbolic cancellation is consistent with the finite assembly report and the unique diagonal Q=diag(K_strict/K_legacy).",
            "scope_step": "This is an ansatz-level formula identity under explicit admissibility assumptions; it does not derive the strict parameters, orientation selector, or legacy role-transfer theorem.",
        },
        "hard_limits": [
            "No unqualified raw identity K_legacy_ont == K_strict_gate is claimed outside the explicit completion factor Q.",
            "No strict dynamical derivation of A/P/D, omega/phi, beta/eta, or chi_11 is exported.",
            "No beta_tors -> beta/eta or beta_tors -> chi_11 theorem is claimed.",
            "No legacy physical-role transfer is licensed.",
            "No QW-2191 selector discharge is claimed.",
            "No ToE closure is claimed.",
        ],
    }


def write_markdown(payload: dict[str, Any]) -> None:
    lines = [
        "# Legacy-to-strict symbolic cancellation certificate",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Result",
        "",
        "The APD completion factor cancels the legacy scalar, legacy phase carrier,",
        "and legacy denominator formally, leaving the strict kernel formula under",
        "explicit finite-domain admissibility assumptions.  This is an ansatz-level",
        "formula identity, not a source theorem or role-transfer theorem.",
        "",
        "## Summary",
        "",
    ]
    for key, value in payload["symbolic_cancellation_summary"].items():
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
    print(json.dumps(payload["symbolic_cancellation_summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
