#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"
IN_P983 = GENERATED / "p983_current_nad12_sigma_pair_realization_side_provider_support_witness_actual_nonexport_probe_summary.json"
IN_P990 = GENERATED / "p990_current_nad12_sigma_pair_side_open_support_refinement_target_attempt_candidate_factor_integration_target_probe_summary.json"
IN_P998 = GENERATED / "p998_current_nad12_sigma_pair_side_candidate_factor_coherence_witness_target_probe_summary.json"
IN_P1008 = GENERATED / "p1008_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_actual_attempt_probe_summary.json"
IN_P1009 = GENERATED / "p1009_current_nad12_sigma_pair_side_sharper_same_lane_witness_attempt_verdict_or_sharper_same_lane_witness_refinement_nonexport_probe_summary.json"
IN_P1010 = GENERATED / "p1010_current_nad12_sigma_pair_side_sharper_same_lane_witness_target_probe_summary.json"
IN_AX9 = ROOT / "AX9_INFORMATIONAL_NADSOLITON_PRIMACY_AXIOM_PACKET.md"

OUT = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe.json"
OUT_SUMMARY = GENERATED / "p1011_current_strict_qw2191_lambda_branch_info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admission_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def contains_all(text: str, needles: list[str]) -> bool:
    lower = text.lower()
    return all(needle.lower() in lower for needle in needles)


def write_artifact(path: Path, payload: dict[str, Any]) -> None:
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_F949, IN_P983, IN_P990, IN_P998, IN_P1008, IN_P1009, IN_P1010, IN_AX9]
    missing = [rel(path) for path in prereq if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1011",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_artifact(OUT, artifact)
        write_artifact(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f949 = load_json(IN_F949)
    p983 = load_json(IN_P983)
    p990 = load_json(IN_P990)
    p998 = load_json(IN_P998)
    p1008 = load_json(IN_P1008)
    p1009 = load_json(IN_P1009)
    p1010 = load_json(IN_P1010)
    ax9_text = IN_AX9.read_text(encoding="utf-8")

    candidate_id = (
        "lambda_nad12_sigma_residual_shannon_information_primary_"
        "sparse_competitive_predictive_coding_like_selector_provider_"
        "shift_candidate_reference_lane_v1"
    )

    checks = [
        {
            "id": "f949_exports_noncyclic_pivot_to_lambda_branch",
            "pass": (
                f949.get("status")
                == "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
                and f949.get("same_lane_exhaustion_boundary_reached") is True
                and f949.get("preferred_noncyclic_pivot_family")
                == "Xi_nad12_sigma_residual_shannon_noncyclic_provider_split_target_v1"
                and f949.get("preferred_first_pivot_branch")
                == "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1"
            ),
            "details": "F949 already exports the honest noncyclic pivot into the Xi -> Lambda family.",
        },
        {
            "id": "ax9_keeps_information_primary_ontology_and_no_lower_information_layer",
            "pass": contains_all(
                ax9_text,
                [
                    "primordial information",
                    "no separate informational layer underneath it",
                ],
            ),
            "details": "AX9 keeps the ontology explicitly information-primary and forbids a lower informational substrate.",
        },
        {
            "id": "p983_keeps_lambda_provider_support_witness_future_only",
            "pass": (
                p983.get("status")
                == "PASS_CURRENT_NAD12_SIGMA_PAIR_REALIZATION_SIDE_PROVIDER_SUPPORT_WITNESS_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
                and p983.get("witness_name")
                == "Lambda_nad12_sigma_residual_shannon_pair_realization_side_provider_support_witness_v1"
                and p983.get("current_repo_has_exported_actual_realization_of_witness") is False
                and p983.get("witness_still_remains_future_only_not_actual_export") is True
            ),
            "details": "P983 keeps actual pair-realization-side provider support explicitly below export grade.",
        },
        {
            "id": "p990_exports_candidate_factor_integration_target",
            "pass": (
                p990.get("status")
                == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_OPEN_SUPPORT_REFINEMENT_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_INTEGRATION_REFINEMENT_TARGET_EXPORTED"
                and p990.get("t276_target_exported_on_current_repo_state") is True
            ),
            "details": "P990 exports one exact candidate-factor integration target on the active lane.",
        },
        {
            "id": "p998_exports_candidate_factor_coherence_witness_refinement_target",
            "pass": (
                p998.get("status")
                == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_CANDIDATE_FACTOR_COHERENCE_TARGET_ATTEMPT_EXACT_CANDIDATE_FACTOR_COHERENCE_WITNESS_REFINEMENT_TARGET_EXPORTED"
                and p998.get("t280_target_exported_on_current_repo_state") is True
            ),
            "details": "P998 exports one exact candidate-factor coherence-witness refinement target on the same lane.",
        },
        {
            "id": "p1008_places_the_lane_into_one_active_actual_attempt",
            "pass": (
                p1008.get("status")
                == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXPORTED"
                and p1008.get("t285_attempt_exported_on_current_repo_state") is True
                and p1008.get("t285_attempt_keeps_open_bundle_below_actual_support") is True
            ),
            "details": "P1008 places the lane into one active actual attempt without promoting success.",
        },
        {
            "id": "p1009_keeps_verdict_and_deeper_export_blocked",
            "pass": (
                p1009.get("status")
                == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_VERDICT_OR_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_NONEXPORT_AUDITED"
                and p1009.get("current_repo_has_lawful_verdict_for_exact_t285_attempt") is False
                and p1009.get("current_repo_has_exact_sharper_same_lane_witness_refinement_below_t285_attempt")
                is False
                and p1009.get("next_honest_move_is_freeze_exact_sharper_same_lane_witness_refinement_target_below_same_attempt")
                is True
            ),
            "details": "P1009 keeps lawful verdict and deeper same-lane export explicitly blocked.",
        },
        {
            "id": "p1010_exports_one_sharper_same_lane_target_below_the_same_attempt",
            "pass": (
                p1010.get("status")
                == "PASS_CURRENT_NAD12_SIGMA_PAIR_SIDE_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_ACTUAL_REALIZATION_ATTEMPT_EXACT_SHARPER_SAME_LANE_WITNESS_REFINEMENT_TARGET_EXPORTED"
                and p1010.get("t286_target_exported_on_current_repo_state") is True
                and p1010.get("t286_target_keeps_open_bundle_below_actual_support") is True
            ),
            "details": "P1010 exports the next sharper same-lane target below the same attempt.",
        },
    ]

    failed_checks = [item["id"] for item in checks if not item["pass"]]
    all_pass = all(item["pass"] for item in checks)

    theorem_result = {
        "info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted": all_pass,
        "candidate_reference_lane_grade": "reference_context_candidate_only" if all_pass else None,
        "candidate_reference_lane_id": candidate_id if all_pass else None,
        "information_primary_ontology_constraint_active": True if all_pass else None,
        "strict_selector_interface_exported": False if all_pass else None,
        "strict_selector_source_exported": False if all_pass else None,
        "provider_class_shift_realized": False if all_pass else None,
        "energy_based_identity_claim_admitted": False if all_pass else None,
        "neural_network_identity_claim_admitted": False if all_pass else None,
        "predictive_coding_like_reading_admitted_only_as_reference_context": all_pass,
        "next_honest_move_requires_exact_selector_interface_audit_for_candidate_lane": all_pass,
        "no_false_pass": True,
    }

    status = (
        "P1011_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_LANE_ADMITTED_SELECTOR_INTERFACE_BLOCKED"
        if all_pass
        else "P1011_REQUIRES_REVIEW"
    )

    artifact = {
        "stage": "P1011",
        "status": status,
        "as_of": AS_OF,
        "inputs": {
            "f949_noncyclic_pivot_packet_summary": rel(IN_F949),
            "ax9_information_primary_ontology_packet": rel(IN_AX9),
            "p983_lambda_provider_support_witness_nonexport_summary": rel(IN_P983),
            "p990_candidate_factor_integration_target_summary": rel(IN_P990),
            "p998_candidate_factor_coherence_witness_target_summary": rel(IN_P998),
            "p1008_lambda_actual_attempt_summary": rel(IN_P1008),
            "p1009_lambda_attempt_verdict_boundary_summary": rel(IN_P1009),
            "p1010_lambda_sharper_same_lane_target_summary": rel(IN_P1010),
        },
        "theorem_result": theorem_result,
        "candidate_reference_lane": {
            "candidate_id": theorem_result["candidate_reference_lane_id"],
            "pivot_family": f949.get("preferred_noncyclic_pivot_family"),
            "pivot_branch": f949.get("preferred_first_pivot_branch"),
            "candidate_grade": theorem_result["candidate_reference_lane_grade"],
            "active_attempt_stage": "t285_attempt_exported",
            "current_sharper_target_stage": "t286_target_exported",
            "reading_contract": (
                "information_primary_sparse_competitive_predictive_coding_like_"
                "selector_provider_shift_candidate_reference_only"
            ),
        },
        "checks": checks,
        "failed_checks": failed_checks,
        "current_honest_reading": [
            "The active Lambda noncyclic branch already has a strong informational/coherence scaffold, but still no strict selector interface.",
            "Because AX9 keeps information primary, the honest interpretation constraint is informational first, not energy-based first.",
            "The current branch may therefore be admitted only as one predictive-coding-like selector-provider shift candidate reference lane.",
        ],
        "recommended_next_packet": {
            "id": "F950_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_PROVIDER_SHIFT_CANDIDATE_REFERENCE_PACKET",
            "goal": "Export the active Lambda branch only as an info-primary SCPC-like selector-provider shift candidate reference lane.",
            "export_object_id": candidate_id,
        },
        "recommended_followup_probe": {
            "id": "P1012_CURRENT_STRICT_QW2191_LAMBDA_BRANCH_INFO_PRIMARY_SCPC_LIKE_SELECTOR_INTERFACE_NONEXPORT_AUDIT_PROBE",
            "goal": "Audit whether any exact strict selector interface already exists for the candidate lane before freezing one new interface target.",
        },
        "no_false_pass": True,
    }

    summary = {
        "stage": "P1011",
        "status": status,
        "as_of": AS_OF,
        "info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted": theorem_result[
            "info_primary_scpc_like_selector_provider_shift_candidate_reference_lane_admitted"
        ],
        "candidate_reference_lane_grade": theorem_result["candidate_reference_lane_grade"],
        "strict_selector_interface_exported": theorem_result["strict_selector_interface_exported"],
        "provider_class_shift_realized": theorem_result["provider_class_shift_realized"],
        "recommended_next_packet": artifact["recommended_next_packet"]["id"],
        "recommended_followup_probe": artifact["recommended_followup_probe"]["id"],
        "no_false_pass": True,
    }

    write_artifact(OUT, artifact)
    write_artifact(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
