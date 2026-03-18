#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"
TOL = 1.0e-15

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P729 = GENERATED / "p729_current_strict_t183_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_P730 = GENERATED / "p730_current_strict_t184_direction_free_shannon_residual_datum_pair12_orbit_direction_selection_bridge_nonexport_audit_probe_summary.json"
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_F647 = GENERATED / "f647_strict_witness_provider_export_packet_for_seed_v1_realization_attempt.json"
IN_P684 = GENERATED / "p684_current_strict_t173_w_break_rooted_directed_state_lift_audit_probe.json"

OUT_JSON = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def branch_score(weight_by_x: list[float], branch: dict[str, Any], n: int) -> float:
    total = 0.0
    for key, value in branch.items():
        idx = int(key) % n
        total += weight_by_x[idx] * float(value)
    return total


def sign_label(x: float, tol: float) -> str:
    if x > tol:
        return "positive"
    if x < -tol:
        return "negative"
    return "zero"


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_P730, IN_F301, IN_F647, IN_P684]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P731",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p729 = load_json(IN_P729)
    p730 = load_json(IN_P730)
    f301 = load_json(IN_F301)
    f647 = load_json(IN_F647)
    p684 = load_json(IN_P684)

    payload = (f647.get("constructed_source_object") or {}).get("exported_payload") or {}
    w_break_raw = payload.get("w_break_by_x")
    pair1_branch = ((f301.get("carrier_vectors") or {}).get("pair1")) or {}
    pair2_branch = ((f301.get("carrier_vectors") or {}).get("pair2")) or {}

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        passed = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": passed,
                "meaning": meaning,
            }
        )
        if not passed:
            blocking.append(check_id)

    valid_w_break = (
        isinstance(w_break_raw, list)
        and len(w_break_raw) == 12
        and all(isinstance(v, (int, float)) for v in w_break_raw)
    )
    if not valid_w_break:
        artifact = {
            "stage": "P731",
            "status": "INVALID_W_BREAK_PAYLOAD",
            "as_of": AS_OF,
            "expected": "F647 constructed_source_object.exported_payload.w_break_by_x is a length-12 numeric list",
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    w_break = [float(v) for v in w_break_raw]

    pair1_score = branch_score(w_break, pair1_branch, n=12)
    pair2_score = branch_score(w_break, pair2_branch, n=12)
    score_gap_abs = abs(pair1_score - pair2_score)
    score_sum_abs = abs(pair1_score + pair2_score)
    pair1_sign = sign_label(pair1_score, TOL)
    pair2_sign = sign_label(pair2_score, TOL)
    opposite_nonzero = (
        abs(pair1_score) > TOL and abs(pair2_score) > TOL and pair1_score * pair2_score < 0.0
    )

    add_check(
        "p729_pair12_orbit_direction_split_already_localized",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving pair1/pair2 ambiguity as opposite orbit-direction branches on the same residual-datum carrier.",
    )
    add_check(
        "p730_direction_free_shannon_lane_still_does_not_select_branch",
        bool(p730.get("current_direction_free_shannon_lane_selects_pair12_orbit_direction_branch")),
        False,
        "P730 already proves the current direction-free Shannon lane does not pick one of the surviving pair1/pair2 branches.",
    )
    add_check(
        "f647_w_break_payload_exported_and_fixed",
        {
            "constructed_source_object_exported": bool(f647.get("constructed_source_object_exported")),
            "admissible_S_sel_int": bool(f647.get("admissible_S_sel_int")),
            "kernel_split_safe": bool(((f647.get("strict_core_export_properties") or {}).get("kernel_split_safe"))),
        },
        {
            "constructed_source_object_exported": True,
            "admissible_S_sel_int": False,
            "kernel_split_safe": True,
        },
        "The current inversion-sensitive witness payload w_break is already exported on a kernel-split-safe strict-core witness lane, but remains below admissible S_sel_int.",
    )
    add_check(
        "p684_w_break_already_used_as_rooted_sign_source_only",
        {
            "root_sign_source": ((p684.get("root") or {}).get("pair") == "pair1"),
            "counts_as_strict_physical_orientation_datum": bool(((p684.get("sign_lift") or {}).get("counts_as_strict_physical_orientation_datum"))),
        },
        {
            "root_sign_source": True,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "P684 already uses w_break as a rooted sign source on the atlas lane, but explicitly not as a strict physical orientation datum.",
    )
    add_check(
        "w_break_witness_payload_pair12_branch_scores_nonzero_and_opposite",
        opposite_nonzero,
        True,
        "On the surviving F301 residual-datum carrier, the exported inversion-sensitive witness payload w_break gives opposite nonzero scores to the two remaining pair1/pair2 branches.",
    )
    add_check(
        "w_break_pair12_branch_scores_are_antisymmetric_under_inversion",
        score_sum_abs <= TOL,
        True,
        "Those two witness scores are antisymmetric under the current pair1/pair2 inversion pairing: score(pair1)+score(pair2)=0 within tolerance.",
    )
    add_check(
        "t185_w_break_witness_payload_pair12_orbit_direction_promotion_bridge_exported",
        False,
        False,
        "The repo still does not export the typed promotion bridge turning this witness-side branch separation into a selected strict source-side branch on full C_v1.",
    )

    status = (
        "PASS_W_BREAK_WITNESS_PAYLOAD_PAIR12_ORBIT_DIRECTION_PROMOTION_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P731_REQUIRES_REVIEW_CHANGED_W_BREAK_WITNESS_PAYLOAD_PAIR12_STATE"
    )

    artifact = {
        "stage": "P731",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_boundary_only",
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "P730": str(IN_P730.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
            "F647_artifact": str(IN_F647.relative_to(REPO)),
            "P684_artifact": str(IN_P684.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t185_target_name": "WBreakWitnessPayloadResidualDatumPair12OrbitDirectionPromotionBridge_global_C_v1_strict_v1",
        "t185_target_exported_on_current_repo_state": False,
        "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": opposite_nonzero,
        "pair12_w_break_branch_scores": {
            "pair1_delta_k_branch_score": pair1_score,
            "pair2_delta_minus_k_branch_score": pair2_score,
            "pair1_score_sign": pair1_sign,
            "pair2_score_sign": pair2_sign,
            "score_gap_abs": score_gap_abs,
            "score_sum_abs": score_sum_abs,
        },
        "audit_conclusion": {
            "current_repo_already_exports_inversion_sensitive_witness_side_pair12_branch_separation": opposite_nonzero,
            "current_repo_already_exports_typed_pair12_branch_promotion_bridge": False,
            "next_honest_move": (
                "attempt_a_typed_promotion_bridge_from_the_current_w_break_witness_payload_branch_separation_to_a_selected_source_side_pair12_branch_without_promoting_F647_to_admissible_S_sel_int_by_fiat"
            ),
        },
        "hard_limits": [
            "No T185 discharge claim.",
            "No claim that F647 is admissible S_sel_int.",
            "No identification of tau_src_candidate_v1 with the current constructed selector carrier by fiat.",
            "No promotion of the current witness-side score sign into a strict physical orientation datum.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P731",
        "status": status,
        "as_of": AS_OF,
        "t185_target_name": "WBreakWitnessPayloadResidualDatumPair12OrbitDirectionPromotionBridge_global_C_v1_strict_v1",
        "t185_target_exported_on_current_repo_state": False,
        "current_w_break_witness_payload_separates_pair12_orbit_direction_branches": opposite_nonzero,
        "pair1_w_break_branch_score_sign": pair1_sign,
        "pair2_w_break_branch_score_sign": pair2_sign,
        "w_break_pair12_branch_scores_are_antisymmetric": score_sum_abs <= TOL,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
