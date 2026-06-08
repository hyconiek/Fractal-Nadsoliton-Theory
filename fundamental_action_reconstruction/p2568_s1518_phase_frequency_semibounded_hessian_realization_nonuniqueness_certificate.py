#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from typing import Any

import numpy as np

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import DOC_FILES, REPO, ROOT, load_json, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

GEN = ROOT / "generated"
OUT = GEN / "p2568_s1518_phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate.json"
MD = GEN / "p2568_s1518_phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate.md"

SOURCE_FILES = {
    "P2567_MINIMAL_STATIONARY_HESSIAN_SADDLE": GEN / "p2567_s1517_phase_frequency_minimal_stationary_hessian_saddle_audit.json",
    "P2566_STATIONARITY_WEIGHT_CONE": GEN / "p2566_s1516_phase_frequency_selector_stationarity_weight_cone_audit.json",
    "P2561_POST_DAMPING_RESIDUAL_BRIDGE": GEN / "p2561_s1511_strict_completion_post_damping_residual_bridge_two_key_certificate.json",
}

HESSIAN_TARGETS = {
    "negative_identity_local_max_candidate": [[-1.0, 0.0], [0.0, -1.0]],
    "positive_identity_local_min_candidate": [[1.0, 0.0], [0.0, 1.0]],
    "anisotropic_negative_definite_candidate": [[-1.0, 0.0], [0.0, -2.0]],
}
NEGATIVE_EXPORT_FLAGS = [
    "semibounded_hessian_weight_source_exported",
    "strict_phase_frequency_source_exported",
    "strict_dynamical_source_for_A_P_D_exported",
    "strict_damping_beta_eta_source_exported",
    "bridge_theorem_exported",
    "legacy_to_strict_completion_bridge_exported",
    "role_transfer_theorem_exported",
    "role_bearing_ltotal_exported",
    "qw2191_discharged_by_this_certificate",
    "toe_closure_claimed",
]


def sha256_json(obj: Any) -> str:
    return hashlib.sha256(json.dumps(obj, sort_keys=True, ensure_ascii=False, separators=(",", ":")).encode()).hexdigest()


def rg_count(pattern: str) -> dict[str, Any]:
    proc = subprocess.run([
        "rg", "-n", pattern, "DIAGRAMS_KERNEL_TRANSFORMATION.md", "fundamental_action_reconstruction", "material_dowodowy",
        "-g", "*.py", "-g", "*.md", "-g", "*.tex", "-g", "!fundamental_action_reconstruction/generated/**",
    ], cwd=REPO, check=False, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    lines = sorted(line for line in proc.stdout.splitlines() if line)
    return {"count": len(lines), "samples": lines[:80]}


def rg_audit() -> dict[str, Any]:
    patterns = {
        "new_packet": "P2568|S1518|semibounded Hessian realization|Hessian realization nonuniqueness|phase frequency semibounded Hessian|phase/frequency semibounded Hessian",
        "intended_research_nonduplication": "higher-support|higher support|semibounded|semi-bounded|variational law|Hessian realization|positive.*Hessian.*phase|Hessian.*positive.*phase",
        "immediate_precursors": "P2567|S1517|minimal stationary Hessian|P2566|S1516|stationarity weight cone|P2561|S1511|strict_phase_frequency_source",
        "guardrails": "QW-2191|selector closure|source theorem|role-transfer theorem|role-bearing L_total|ToE closure",
    }
    return {"tool": "rg", "patterns": {name: rg_count(pattern) for name, pattern in patterns.items()}}


def theorem(payload: dict[str, Any], container: str) -> dict[str, Any]:
    return payload.get(container, {}).get("theorem_export", {})


def constraint_matrix(rows: list[dict[str, Any]]) -> np.ndarray:
    matrix_rows = []
    for row_name in ["gradient_phi", "gradient_omega", "hessian_omega_omega", "hessian_omega_phi", "hessian_phi_phi"]:
        values = []
        for row in rows:
            d = int(row["d"])
            sign = int(row["sign"])
            cos_theta = math.cos(float(row["theta"]))
            if row_name == "gradient_phi":
                values.append(float(row["gradient_phi_coefficient"]))
            elif row_name == "gradient_omega":
                values.append(float(row["gradient_omega_coefficient"]))
            else:
                coeff = -sign * cos_theta
                if row_name == "hessian_omega_omega":
                    values.append(coeff * d * d)
                elif row_name == "hessian_omega_phi":
                    values.append(coeff * d)
                elif row_name == "hessian_phi_phi":
                    values.append(coeff)
        matrix_rows.append(values)
    return np.array(matrix_rows, dtype=float)


def solve_target(matrix: np.ndarray, target_hessian: list[list[float]]) -> dict[str, Any]:
    target = np.array([0.0, 0.0, target_hessian[0][0], target_hessian[0][1], target_hessian[1][1]], dtype=float)
    weights, *_ = np.linalg.lstsq(matrix, target, rcond=None)
    residual = matrix @ weights - target
    hessian = [[float(target_hessian[0][0]), float(target_hessian[0][1])], [float(target_hessian[0][1]), float(target_hessian[1][1])]]
    eigenvalues = np.linalg.eigvalsh(np.array(hessian, dtype=float))
    if np.all(eigenvalues < 0.0):
        definiteness = "negative_definite_local_max_candidate"
    elif np.all(eigenvalues > 0.0):
        definiteness = "positive_definite_local_min_candidate"
    else:
        definiteness = "not_definite"
    return {
        "target_hessian": hessian,
        "target_eigenvalues": [float(value) for value in eigenvalues],
        "target_definiteness": definiteness,
        "solution_weight_count": int(weights.size),
        "solution_weights": [float(value) for value in weights],
        "has_positive_weights": bool(np.any(weights > 1e-12)),
        "has_negative_weights": bool(np.any(weights < -1e-12)),
        "weight_l2_norm": float(np.linalg.norm(weights)),
        "max_abs_residual": float(np.max(np.abs(residual))),
        "residual_vector": [float(value) for value in residual],
    }


def realization_audit(rows: list[dict[str, Any]]) -> dict[str, Any]:
    matrix = constraint_matrix(rows)
    rank = int(np.linalg.matrix_rank(matrix, tol=1e-10))
    target_rows = {name: solve_target(matrix, target) for name, target in HESSIAN_TARGETS.items()}
    return {
        "backend": "numpy.linalg.lstsq",
        "numpy_version": np.__version__,
        "constraint_matrix_shape": [int(value) for value in matrix.shape],
        "constraint_matrix_rank": rank,
        "solution_affine_nullity_for_fixed_hessian": int(matrix.shape[1] - rank),
        "target_rows": target_rows,
        "all_targets_realized_to_tolerance": all(row["max_abs_residual"] < 1e-10 for row in target_rows.values()),
        "both_local_max_and_local_min_hessians_realizable": (
            target_rows["negative_identity_local_max_candidate"]["max_abs_residual"] < 1e-10
            and target_rows["positive_identity_local_min_candidate"]["max_abs_residual"] < 1e-10
        ),
        "all_realizations_use_signed_weights": all(row["has_positive_weights"] and row["has_negative_weights"] for row in target_rows.values()),
    }


def build_payload() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    grep = rg_audit()
    p2567_payload = load_json(SOURCE_FILES["P2567_MINIMAL_STATIONARY_HESSIAN_SADDLE"])
    p2566_payload = load_json(SOURCE_FILES["P2566_STATIONARITY_WEIGHT_CONE"])
    p2561_payload = load_json(SOURCE_FILES["P2561_POST_DAMPING_RESIDUAL_BRIDGE"])
    p2567 = theorem(p2567_payload, "phase_frequency_minimal_stationary_hessian_saddle_audit")
    p2566 = theorem(p2566_payload, "phase_frequency_selector_stationarity_weight_cone_audit")
    p2561 = theorem(p2561_payload, "strict_completion_post_damping_residual_bridge_two_key_certificate")
    audit = realization_audit(p2566["stationarity_rows"])
    theorem_export: dict[str, Any] = {
        "theorem_name": "P2568_T1_phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate",
        "audited_chain": ["P2567/S1517", "P2566/S1516", "P2561/S1511"],
        "frontier_atom_under_attack": "strict_phase_frequency_source",
        "p2567_minimal_saddle_obstruction_inherited": p2567.get("minimal_signed_stationary_witnesses_fail_second_order_selector_test") is True,
        "p2566_stationarity_weight_cone_obstruction_inherited": p2566.get("nonnegative_natural_weight_stationarity_rejected") is True,
        "p2561_phase_frequency_residual_atom_inherited": "strict_phase_frequency_source" in p2561.get("residual_atoms_after_hypothetical_damping_source", []),
        "semibounded_hessian_realization_audit": audit,
        "semibounded_hessian_can_be_realized_but_not_sourced": audit["both_local_max_and_local_min_hessians_realizable"] and audit["all_realizations_use_signed_weights"],
        "hessian_sign_choice_is_extra_source_obligation": True,
        "recommended_next_honest_step": (
            "Do not treat higher-support signed semibounded Hessian realization as phase/frequency closure. P2568 shows local max and local min Hessians can both be manufactured with signed weights and a 7-dimensional affine freedom. The next honest step is to derive a positivity/measure/weight law from strict nadsoliton dynamics, or prove no such law can coexist with the strict phase/frequency target and QW-2191 guardrail."
        ),
        **{flag: False for flag in NEGATIVE_EXPORT_FLAGS},
    }
    gatekeepers = {
        "rg_non_duplication_audit_recorded": grep["tool"] == "rg",
        "p2567_inherited": theorem_export["p2567_minimal_saddle_obstruction_inherited"],
        "p2566_inherited": theorem_export["p2566_stationarity_weight_cone_obstruction_inherited"],
        "p2561_phase_frequency_atom_inherited": theorem_export["p2561_phase_frequency_residual_atom_inherited"],
        "constraint_rank_five": audit["constraint_matrix_rank"] == 5,
        "affine_nullity_seven": audit["solution_affine_nullity_for_fixed_hessian"] == 7,
        "all_targets_realized": audit["all_targets_realized_to_tolerance"],
        "local_max_and_min_both_realizable": audit["both_local_max_and_local_min_hessians_realizable"],
        "all_realizations_signed": audit["all_realizations_use_signed_weights"],
        "no_phase_source_exported": theorem_export["strict_phase_frequency_source_exported"] is False,
        "no_bridge_exported": theorem_export["bridge_theorem_exported"] is False,
        "no_role_transfer_exported": theorem_export["role_transfer_theorem_exported"] is False,
        "no_qw2191_discharge": theorem_export["qw2191_discharged_by_this_certificate"] is False,
        "no_toe_closure_claimed": theorem_export["toe_closure_claimed"] is False,
    }
    return {
        "packet_id": "P2568",
        "stage_id": "S1518",
        "status": "P2568_PHASE_FREQUENCY_SEMIBOUNDED_HESSIAN_REALIZATION_NONUNIQUENESS_NO_SOURCE_EXPORT_NO_BRIDGE_NO_ROLE_TRANSFER_NO_QW2191_NO_TOE",
        "source_files": {name: rel(path) for name, path in SOURCE_FILES.items()},
        "rg_non_duplication_audit": grep,
        "phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate": {
            "theorem_export": theorem_export,
            "source_fingerprints_sha256": {
                "P2567_MINIMAL_STATIONARY_HESSIAN_SADDLE": sha256_json(p2567_payload),
                "P2566_STATIONARITY_WEIGHT_CONE": sha256_json(p2566_payload),
                "P2561_POST_DAMPING_RESIDUAL_BRIDGE": sha256_json(p2561_payload),
            },
            "theorem_fingerprint_sha256": sha256_json(theorem_export),
        },
        "gatekeeper_checks": gatekeepers,
    }


def write_markdown(payload: dict[str, Any]) -> None:
    t = payload["phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate"]["theorem_export"]
    audit = t["semibounded_hessian_realization_audit"]
    lines = [
        "# P2568/S1518 phase/frequency semibounded Hessian realization nonuniqueness certificate", "",
        f"Status: `{payload['status']}`", "", "## Result", "",
        f"- Frontier atom under attack: `{t['frontier_atom_under_attack']}`.",
        f"- P2567 minimal-saddle obstruction inherited: `{t['p2567_minimal_saddle_obstruction_inherited']}`.",
        f"- Constraint matrix rank/nullity: `{audit['constraint_matrix_rank']}/{audit['solution_affine_nullity_for_fixed_hessian']}`.",
        f"- All target Hessians realized to tolerance: `{audit['all_targets_realized_to_tolerance']}`.",
        f"- Local max and local min Hessians both realizable: `{audit['both_local_max_and_local_min_hessians_realizable']}`.",
        f"- All realizations use signed weights: `{audit['all_realizations_use_signed_weights']}`.",
        f"- Hessian sign choice is extra source obligation: `{t['hessian_sign_choice_is_extra_source_obligation']}`.", "",
        "## Interpretation", "",
        "P2568 shows that moving beyond minimal three-node supports does not close the phase/frequency source.  With signed higher-support weights one can realize both negative-definite and positive-definite Hessians at the strict tuple, with seven residual affine degrees of freedom.  Semiboundedness can therefore be manufactured unless a strict weight/measure law is derived.", "",
        "## Recommended next honest step", "", t["recommended_next_honest_step"], "", "## Negative controls", "",
        "No semibounded-Hessian weight source, strict phase/frequency source, bridge theorem, role-transfer theorem, QW-2191 discharge, role-bearing `L_total`, or ToE closure is exported.", "", "## Fingerprint", "",
        f"`{payload['phase_frequency_semibounded_hessian_realization_nonuniqueness_certificate']['theorem_fingerprint_sha256']}`",
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def append_doc_sections() -> None:
    eq_section = """
`P2568/S1518` audits higher-support signed Hessian realization after the P2567 minimal-saddle obstruction.  The five linear constraints consisting of two stationarity equations plus a target Hessian have rank `5` on `12` weights, leaving a `7`-dimensional affine freedom.  The audit explicitly realizes negative-identity, positive-identity, and anisotropic negative-definite Hessians at the strict tuple, all using signed weights.  Therefore semibounded Hessian realization is not itself a phase/frequency source; the sign/measure/weight law remains the missing strict obligation.
""".strip()
    lag_section = """
`P2568/S1518` blocks promotion of higher-support signed semibounded Hessian fits into role-bearing phase/frequency `L_total` terms.  Both local max and local min Hessians can be manufactured post hoc with signed weights, so a strict source must provide the positivity/measure law, bridge role discipline, and QW-2191 handling.
""".strip()
    append_once(DOC_FILES["equation_sheet"], "## P2568/S1518 phase/frequency semibounded Hessian realization guard", "## P2568/S1518 phase/frequency semibounded Hessian realization guard\n\n" + eq_section)
    append_once(DOC_FILES["lagrangian_eom_draft"], "## P2568/S1518 phase/frequency semibounded Hessian realization Ltotal guard", "## P2568/S1518 phase/frequency semibounded Hessian realization Ltotal guard\n\n" + lag_section)


def main() -> None:
    payload = build_payload()
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
    write_markdown(payload)
    append_doc_sections()
    print(json.dumps(payload, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
