#!/usr/bin/env python3
from __future__ import annotations

import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
import scipy.linalg as la
import sympy as sp

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_ALPHA = GEN / "alpha_geo_strict_derived_v1.json"
IN_P1535 = GEN / "p1535_s485_strict_selector_tie_break_candidate_summary.json"
IN_P753 = GEN / "p753_current_strict_t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_audit_probe_summary.json"
IN_P2282 = GEN / "p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json"
IN_P2297 = GEN / "p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json"
IN_P2298 = GEN / "p2298_s1248_strict_non_gb_provider_source_nonadmissibility_theorem_probe.json"
OUT = GEN / "p2299_s1249_strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe.json"
MD = GEN / "p2299_s1249_strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe.md"


def load(path: Path) -> dict[str, Any]:
    if not path.exists():
        return {"_missing": str(path), "status": "OPEN_MISSING_ARTIFACT"}
    return json.loads(path.read_text(encoding="utf-8"))


def sha256_file(path: Path) -> str:
    if not path.exists():
        return ""
    return hashlib.sha256(path.read_bytes()).hexdigest()


def sha256_json(payload: Any) -> str:
    blob = json.dumps(payload, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(blob.encode("utf-8")).hexdigest()


def boolish(value: Any) -> bool:
    return bool(value) is True


def main() -> None:
    GEN.mkdir(exist_ok=True)
    alpha_packet = load(IN_ALPHA)
    p1535 = load(IN_P1535)
    p753 = load(IN_P753)
    p2282 = load(IN_P2282)
    p2297 = load(IN_P2297)
    p2298 = load(IN_P2298)

    alpha_value_expr = sp.sympify("4*log(2)")
    alpha_value_float = float(sp.N(alpha_value_expr, 30))
    shannon_weights_exact = [sp.exp(-alpha_value_expr * i).simplify() for i in range(8)]
    shannon_weights_numeric = np.array([float(sp.N(w, 30)) for w in shannon_weights_exact], dtype=float)
    shannon_weight_ratios = [float(shannon_weights_numeric[i + 1] / shannon_weights_numeric[i]) for i in range(7)]
    monotone_strict = bool(np.all(np.diff(shannon_weights_numeric) < 0.0))
    normalized = shannon_weights_numeric / float(np.sum(shannon_weights_numeric))
    entropy = float(-np.sum(normalized * np.log(normalized)))
    anisotropy_norm = float(la.norm(normalized - np.ones_like(normalized) / len(normalized), ord=2))

    p2297_probe = p2297.get("strict_non_gb_spatial_eom_provider_matrix_obstruction_probe", {})
    p2297_results = p2297_probe.get("provider_matrix_results", {})
    strict_core = p2297_results.get("strict_core_minimal_provider", {})
    augmented = p2297_results.get("p1990_augmented_provider_non_strict", {})
    formal_full = p2297_results.get("formal_full_residual_basis_provider", {})
    strict_core_report = strict_core.get("matrix_report", {})
    augmented_report = augmented.get("matrix_report", {})
    formal_full_report = formal_full.get("matrix_report", {})

    p2298_probe = p2298.get("strict_non_gb_provider_source_nonadmissibility_theorem_probe", {})
    p2298_nonadmissibility = p2298_probe.get("nonadmissibility_result", {})
    p2282_probe = p2282.get("strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe", {})
    gap_rows = p2282_probe.get("gap_rows", [])
    gap_status = {row.get("id", "UNKNOWN"): row.get("status", "UNKNOWN") for row in gap_rows}

    alternatives = p1535.get("alternatives", [])
    branch_ids = [item.get("branch_id") for item in alternatives]
    branch_b_noncyclic_available = "branch_B_noncyclic" in branch_ids
    selected_branch = (p1535.get("tie_break_result", {}) or {}).get("selected_branch_after_tie_break")
    qw2191_closed_by_tie_break = boolish(p1535.get("qw2191_closed"))
    shannon_lane_candidate_only = boolish(p753.get("current_strict_source_shannon_minimal_designated_pair12_entry_lane_remains_candidate_only"))
    provider_shift_required = boolish(p753.get("next_honest_move_requires_new_entry_object_or_provider_shift"))

    # Exact screening lemma: multiplying an inconsistent provider matrix A by any
    # nonzero scalar/source value alpha does not change rank(A) nor rank([A|b]).
    # Adding a selector/tie-break label without a new spatial-EOM operator column
    # also cannot change the P2297 provider matrix.  Therefore the Shannon value
    # is a real strict-side source for candidate construction, but not by itself a
    # full-family provider source.
    a_rank = int(strict_core_report.get("rank_A", 0))
    aug_rank = int(strict_core_report.get("rank_augmented", 0))
    alpha_nonzero = sp.simplify(alpha_value_expr) != 0
    scaled_rank_invariance_certificate = {
        "lemma": "rank(alpha*A)=rank(A) and rank([alpha*A|b])=rank([A|b]) for alpha != 0 after reparametrizing unknowns by alpha*c",
        "alpha_expr": str(alpha_value_expr),
        "alpha_nonzero": bool(alpha_nonzero),
        "source_rank_A": a_rank,
        "source_rank_augmented": aug_rank,
        "scaled_rank_A": a_rank,
        "scaled_rank_augmented": aug_rank,
        "inconsistency_preserved": bool(a_rank < aug_rank),
    }

    strict_shannon_attempt = {
        "attempt_id": "STRICT_SHANNON_NONCYCLIC_PROVIDER_SOURCE_ATTEMPT_V1",
        "uses_alpha_geo_strict_derived_v1": alpha_packet.get("status") == "actual_exported_strict_derived_source_upgrade_value",
        "alpha_geo_value": str(alpha_value_expr),
        "shannon_weights_exact_first_8": [str(w) for w in shannon_weights_exact],
        "shannon_weights_numeric_first_8": shannon_weights_numeric.tolist(),
        "successive_weight_ratios": shannon_weight_ratios,
        "weights_monotone_strict": monotone_strict,
        "normalized_entropy": entropy,
        "anisotropy_l2_from_uniform": anisotropy_norm,
        "branch_B_noncyclic_available_as_candidate": branch_b_noncyclic_available,
        "selected_tie_break_branch": selected_branch,
        "qw2191_closed_by_tie_break": qw2191_closed_by_tie_break,
        "shannon_lane_candidate_only": shannon_lane_candidate_only,
        "provider_shift_required_by_P753": provider_shift_required,
        "does_not_copy_residual_basis": True,
        "new_spatial_eom_operator_columns_exported": False,
        "scaled_rank_invariance_certificate": scaled_rank_invariance_certificate,
        "full_residual_family_cancellation_status": "FAILS_CURRENT_STRICT_PROVIDER_TEST",
        "failure_reason": (
            "4*ln(2) is accepted here as a strict-side Shannon source, but it supplies a value/orientation weight, "
            "not a new ADM/Bianchi-I spatial-EOM provider operator.  Without new operator columns, the P2297 strict-core "
            "matrix inconsistency is unchanged; with an explicit selector premise the route becomes non-strict."
        ),
    }

    non_strict_selector_branch = {
        "branch_id": "NON_STRICT_SHANNON_SELECTOR_PREMISE_BRANCH_V1",
        "status": "FORMALIZED_AS_NON_STRICT_BRANCH_ONLY",
        "premise": "Use alpha_geo_strict_derived_v1 = 4 ln 2 plus an explicit selector/tie-break premise to orient the provider lane.",
        "why_non_strict": [
            "P1535 tie-break candidate leaves qw2191_closed=false.",
            "P753 says the strict-source Shannon entry lane remains candidate-only and requires a provider shift.",
            "P2297/P2298 show the only full-family algebraic cancellation is the residual-copying formal basis, not a legal strict provider.",
            "An explicit selector/tie-break premise is exactly the extra symmetry-breaking datum forbidden from silent strict-core promotion by QW-2191 guardrails.",
        ],
        "may_be_used_for_control_experiments": True,
        "may_close_strict_G1_G2_G3": False,
        "may_discharge_QW2191": False,
        "may_claim_ToE_closure": False,
    }

    task3_impact = {
        "source_packet": "P2282",
        "gap_status_before_p2299": gap_status,
        "gap_status_after_p2299": gap_status,
        "closure_score_after_p2299": p2282_probe.get("closure_score"),
        "reason": "P2299 does not produce a legal strict full-family provider, so it must not update G1/G2/G3 closure rows.",
    }

    gatekeeper_checks = {
        "alpha_geo_strict_source_loaded": strict_shannon_attempt["uses_alpha_geo_strict_derived_v1"],
        "alpha_geo_is_four_ln2_not_legacy_import": alpha_packet.get("value") == "4 ln(2)",
        "shannon_weights_nonuniform": anisotropy_norm > 0.0 and monotone_strict,
        "branch_B_noncyclic_candidate_seen": branch_b_noncyclic_available,
        "candidate_does_not_copy_residual_basis": strict_shannon_attempt["does_not_copy_residual_basis"],
        "no_new_spatial_eom_operator_columns_exported": not strict_shannon_attempt["new_spatial_eom_operator_columns_exported"],
        "scaled_strict_core_inconsistency_preserved": scaled_rank_invariance_certificate["inconsistency_preserved"],
        "p2298_no_current_legal_strict_provider_confirmed": p2298_nonadmissibility.get("legal_strict_capable_sources") == [],
        "strict_shannon_provider_attempt_fails_full_residual": strict_shannon_attempt["full_residual_family_cancellation_status"] == "FAILS_CURRENT_STRICT_PROVIDER_TEST",
        "non_strict_selector_branch_formalized": non_strict_selector_branch["status"] == "FORMALIZED_AS_NON_STRICT_BRANCH_ONLY",
        "task3_g1_g2_g3_not_closed_by_p2299": all(status == "OPEN" for status in gap_status.values()),
        "qw2191_not_discharged": not qw2191_closed_by_tie_break,
        "no_selector_closure_claimed": True,
        "no_legacy_kernel_role_transfer": True,
        "no_toe_closure_claimed": True,
    }

    theorem_export = {
        "statement_id": "P2299_STRICT_SHANNON_PROVIDER_ATTEMPT_FAILURE_AND_NON_STRICT_BRANCH_THEOREM",
        "formal_statement": (
            "The strict-derived Shannon value alpha_geo_strict_derived_v1 = 4 ln 2 is admissible as a strict-side "
            "candidate source, but the current artifacts export no new ADM/Bianchi-I spatial-EOM provider columns from it. "
            "Thus the Shannon-scaled strict-core provider is rank-equivalent to the P2297 strict-core provider and remains "
            "inconsistent for the full non-GB residual family.  Any use of the Shannon/tie-break datum as an actual selector "
            "premise is therefore a separate non-strict selector-premise branch and cannot close strict G1/G2/G3 or QW-2191."
        ),
        "proof_bits": {
            "alpha_strict_source_exists": gatekeeper_checks["alpha_geo_strict_source_loaded"],
            "alpha_rescaling_preserves_strict_core_inconsistency": gatekeeper_checks["scaled_strict_core_inconsistency_preserved"],
            "no_new_eom_provider_columns": gatekeeper_checks["no_new_spatial_eom_operator_columns_exported"],
            "p2298_no_legal_provider": gatekeeper_checks["p2298_no_current_legal_strict_provider_confirmed"],
            "task3_gaps_remain_open": gatekeeper_checks["task3_g1_g2_g3_not_closed_by_p2299"],
        },
        "not_claimed": [
            "strict provider closure",
            "strict G1/G2/G3 closure",
            "QW-2191 discharge",
            "selector closure",
            "legacy-kernel role transfer",
            "ToE closure",
        ],
    }
    theorem_fingerprint = sha256_json(theorem_export)

    payload = {
        "schema_version": "p2299_s1249_v1",
        "packet_id": "P2299",
        "stage_id": "S1249",
        "produced_by": Path(__file__).name,
        "timestamp_utc": "2026-05-28T00:00:00+00:00",
        "status": "OPEN_OBSTRUCTION_WITH_STRICT_SHANNON_ATTEMPT_AND_NON_STRICT_BRANCH_TRACE",
        "result_kind": "STRICT_SHANNON_PROVIDER_SOURCE_ATTEMPT_FAILS_AND_NON_STRICT_SELECTOR_BRANCH_FORMALIZED",
        "strict_shannon_provider_source_attempt_and_non_strict_selector_branch_probe": {
            "probe_id": "P2299_S1249_STRICT_SHANNON_PROVIDER_SOURCE_ATTEMPT_AND_NON_STRICT_SELECTOR_BRANCH",
            "source_packets": {
                "alpha_geo_strict_derived_v1": "generated/alpha_geo_strict_derived_v1.json",
                "p1535": "generated/p1535_s485_strict_selector_tie_break_candidate_summary.json",
                "p753": "generated/p753_current_strict_t207_strict_source_shannon_minimal_designated_pair12_entry_lane_split_exhaustion_boundary_audit_probe_summary.json",
                "p2282": "generated/p2282_s1232_strict_task3_global_bianchi_i_g1_g2_g3_closure_matrix_probe.json",
                "p2297": "generated/p2297_s1247_strict_non_gb_spatial_eom_provider_matrix_obstruction_probe.json",
                "p2298": "generated/p2298_s1248_strict_non_gb_provider_source_nonadmissibility_theorem_probe.json",
            },
            "source_hashes": {
                "alpha_geo_strict_derived_v1_sha256": sha256_file(IN_ALPHA),
                "p1535_sha256": sha256_file(IN_P1535),
                "p753_sha256": sha256_file(IN_P753),
                "p2282_sha256": sha256_file(IN_P2282),
                "p2297_sha256": sha256_file(IN_P2297),
                "p2298_sha256": sha256_file(IN_P2298),
            },
            "strict_shannon_provider_attempt": strict_shannon_attempt,
            "non_strict_selector_premise_branch": non_strict_selector_branch,
            "task3_g1_g2_g3_impact": task3_impact,
            "upstream_provider_reports_used": {
                "p2297_strict_core_matrix_report": strict_core_report,
                "p2297_augmented_matrix_report": augmented_report,
                "p2297_formal_full_basis_matrix_report": formal_full_report,
            },
            "theorem_export": theorem_export,
            "theorem_fingerprint_sha256": theorem_fingerprint,
        },
        "recommended_next_honest_step": {
            "id": "P2300_candidate",
            "goal": "If strict closure is still desired, derive an actual ADM/Bianchi-I spatial-EOM provider operator from the Shannon/nad12-sigma source and re-run the provider matrix; otherwise keep the selector-premise branch explicitly non-strict and do not update G1/G2/G3.",
        },
        "gatekeeper_checks": gatekeeper_checks,
        "global_status": "OPEN_OBSTRUCTION_WITH_TRACE",
    }

    OUT.write_text(json.dumps(payload, indent=2, sort_keys=False) + "\n", encoding="utf-8")
    md = f"""# P2299/S1249 — strict Shannon provider-source attempt and non-strict selector branch\n\n- Status: `{payload['status']}`\n- Strict Shannon source: `alpha_geo_strict_derived_v1 = {alpha_packet.get('value')}`; numeric value `{alpha_value_float:.15f}`.\n- Candidate weights are non-uniform: `{gatekeeper_checks['shannon_weights_nonuniform']}`; first ratio `{shannon_weight_ratios[0]:.15f}`.\n- Strict provider attempt verdict: `{strict_shannon_attempt['full_residual_family_cancellation_status']}`.\n- Reason: Shannon supplies a strict-side value/orientation, but no new ADM/Bianchi-I spatial-EOM provider columns are exported.\n- Non-strict branch: `{non_strict_selector_branch['status']}`.\n- Task-3 G1/G2/G3 after P2299: `{task3_impact['gap_status_after_p2299']}`.\n- Theorem fingerprint: `{theorem_fingerprint}`\n\n## Guardrail statement\nP2299 accepts `4 ln 2` only as the strict-derived Shannon source already exported by `alpha_geo_strict_derived_v1`; it does not import the retired legacy kernel role, does not discharge QW-2191, and does not close G1/G2/G3.\n\n## Next honest step\n{payload['recommended_next_honest_step']['goal']}\n"""
    MD.write_text(md, encoding="utf-8")


if __name__ == "__main__":
    main()
