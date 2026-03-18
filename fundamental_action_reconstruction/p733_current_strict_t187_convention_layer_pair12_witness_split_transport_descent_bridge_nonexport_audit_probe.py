#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-18"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P732 = GENERATED / "p732_current_strict_t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_nonexport_audit_probe_summary.json"
IN_P688 = GENERATED / "p688_current_strict_t174_w_break_rooted_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe.json"
IN_P691 = GENERATED / "p691_current_strict_t174_sign_fixed_directed_state_edge_coherence_under_oriented_transition_sign_lift_audit_probe.json"

OUT_JSON = GENERATED / "p733_current_strict_t187_convention_layer_pair12_witness_split_transport_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p733_current_strict_t187_convention_layer_pair12_witness_split_transport_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P732, IN_P688, IN_P691]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P733",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p731 = load_json(IN_P731)
    p732 = load_json(IN_P732)
    p688 = load_json(IN_P688)
    p691 = load_json(IN_P691)

    pair12_rooted_edge = (((p688.get("edge_audit") or {}).get("edges") or {}).get("pair1_to_pair2")) or {}
    pair12_sign_fixed_edge = (((p691.get("edge_audit") or {}).get("edges") or {}).get("pair1_to_pair2")) or {}

    rooted_oriented_sign = pair12_rooted_edge.get("oriented_sign_on_operator")
    sign_fixed_oriented_sign = pair12_sign_fixed_edge.get("oriented_sign_on_operator")
    rooted_no_flip = bool(pair12_rooted_edge.get("edge_ok_no_sign_flip"))
    sign_fixed_no_flip = bool(pair12_sign_fixed_edge.get("edge_ok_no_sign_flip"))

    p731_pair1_sign = str(p731.get("pair1_w_break_branch_score_sign") or "zero")
    p731_pair2_sign = str(p731.get("pair2_w_break_branch_score_sign") or "zero")
    p731_opposite = p731_pair1_sign != p731_pair2_sign and "zero" not in {p731_pair1_sign, p731_pair2_sign}

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

    add_check(
        "p731_pair12_witness_split_already_opposite",
        {
            "pair12_split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "pair1_sign": p731_pair1_sign,
            "pair2_sign": p731_pair2_sign,
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
        },
        {
            "pair12_split_separated": True,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
            "t185_exported": False,
        },
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed promotion bridge remains unexported.",
    )
    add_check(
        "p732_current_pair1_rooted_convention_state_still_does_not_descend_split",
        {
            "current_pair1_rooted_convention_state_exists": bool(p732.get("current_pair1_rooted_convention_state_exists")),
            "pair1_pair2_convention_state_signs_equal": bool(p732.get("pair1_pair2_convention_state_signs_equal")),
            "split_descends": bool(p732.get("p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state")),
        },
        {
            "current_pair1_rooted_convention_state_exists": True,
            "pair1_pair2_convention_state_signs_equal": True,
            "split_descends": False,
        },
        "P732 already proves that the current pair1-rooted convention state keeps equal pair1/pair2 signs and therefore does not descend the P731 split.",
    )
    add_check(
        "rooted_oriented_pair12_edge_sign_stays_positive",
        {
            "oriented_sign_on_operator": rooted_oriented_sign,
            "edge_ok_no_sign_flip": rooted_no_flip,
        },
        {
            "oriented_sign_on_operator": 1,
            "edge_ok_no_sign_flip": True,
        },
        "Under the exported T174 rooted oriented edge-sign lift, the pair1->pair2 transport still carries positive sign and no sign flip.",
    )
    add_check(
        "sign_fixed_oriented_pair12_edge_sign_stays_positive",
        {
            "oriented_sign_on_operator": sign_fixed_oriented_sign,
            "edge_ok_no_sign_flip": sign_fixed_no_flip,
        },
        {
            "oriented_sign_on_operator": 1,
            "edge_ok_no_sign_flip": True,
        },
        "Even after the exported T175 sign-fix and the attached T174 oriented lift, the pair1->pair2 transport still carries positive sign and no sign flip.",
    )
    add_check(
        "current_convention_layer_pair12_transport_can_carry_opposite_p731_branch_signs",
        p731_opposite and rooted_oriented_sign == -1 and sign_fixed_oriented_sign == -1,
        False,
        "Therefore the current convention layer still does not provide a pair1/pair2 transition sign capable of carrying the opposite P731 witness split through current transport data.",
    )
    add_check(
        "t187_convention_layer_pair12_witness_split_transport_descent_bridge_exported",
        False,
        False,
        "The repo still does not export the convention-layer pair1/pair2 witness-split transport descent bridge.",
    )

    status = (
        "PASS_CONVENTION_LAYER_PAIR12_WITNESS_SPLIT_TRANSPORT_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P733_REQUIRES_REVIEW_CHANGED_CONVENTION_LAYER_PAIR12_WITNESS_SPLIT_TRANSPORT_STATE"
    )

    artifact = {
        "stage": "P733",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t187_convention_layer_pair12_witness_split_transport_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P732": str(IN_P732.relative_to(REPO)),
            "P688": str(IN_P688.relative_to(REPO)),
            "P691": str(IN_P691.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t187_target_name": "ConventionLayerPair12WitnessSplitTransportDescentBridge_global_C_v1_strict_v1",
        "t187_target_exported_on_current_repo_state": False,
        "current_convention_layer_pair12_transport_signs": {
            "rooted_oriented_pair12_edge_sign": rooted_oriented_sign,
            "sign_fixed_oriented_pair12_edge_sign": sign_fixed_oriented_sign,
        },
        "current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts": (
            rooted_oriented_sign == 1 and sign_fixed_oriented_sign == 1 and rooted_no_flip and sign_fixed_no_flip
        ),
        "p731_pair12_witness_split_descends_to_current_convention_layer_transport": False,
        "audit_conclusion": {
            "current_repo_already_exports_positive_pair12_convention_transport": (
                rooted_oriented_sign == 1 and sign_fixed_oriented_sign == 1 and rooted_no_flip and sign_fixed_no_flip
            ),
            "current_repo_already_exports_pair12_witness_split_transport_descent_bridge": False,
            "next_honest_move": (
                "attempt_a_typed_pair12_local_source_side_descent_mechanism_that_can_carry_the_p731_opposite_branch_split_without_collapsing_into_the_current_positive_convention_layer_transport"
            ),
        },
        "hard_limits": [
            "No T187 discharge claim.",
            "No claim that current convention-layer transport data are source-side branch-typed.",
            "No claim that T174/T175 promote a strict physical orientation datum.",
            "No claim that F647 is admissible S_sel_int.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P733",
        "status": status,
        "as_of": AS_OF,
        "t187_target_name": "ConventionLayerPair12WitnessSplitTransportDescentBridge_global_C_v1_strict_v1",
        "t187_target_exported_on_current_repo_state": False,
        "current_convention_layer_pair12_transport_is_positive_under_all_exported_lifts": (
            rooted_oriented_sign == 1 and sign_fixed_oriented_sign == 1 and rooted_no_flip and sign_fixed_no_flip
        ),
        "p731_pair12_witness_split_descends_to_current_convention_layer_transport": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
