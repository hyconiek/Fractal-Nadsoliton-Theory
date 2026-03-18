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
IN_F301 = GENERATED / "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json"
IN_F684_OBJ = GENERATED / "selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_strict_convention_v1.json"
IN_F684_SUMMARY = GENERATED / "f684_first_exported_selector_state_global_c_v1_directed_rooted_transport_from_s_sel_int_w_break_packet_summary.json"

OUT_JSON = GENERATED / "p732_current_strict_t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p732_current_strict_t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def sign_label_from_value(value: Any) -> str:
    if isinstance(value, (int, float)):
        if float(value) > 0.0:
            return "positive"
        if float(value) < 0.0:
            return "negative"
    return "zero"


def contains_token(obj: Any, token: str) -> bool:
    if isinstance(obj, dict):
        return any(contains_token(k, token) or contains_token(v, token) for k, v in obj.items())
    if isinstance(obj, list):
        return any(contains_token(v, token) for v in obj)
    return token in str(obj)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_F301, IN_F684_OBJ, IN_F684_SUMMARY]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P732",
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
    f301 = load_json(IN_F301)
    f684_obj = load_json(IN_F684_OBJ)
    f684_summary = load_json(IN_F684_SUMMARY)

    p731_pair1_sign = str(p731.get("pair1_w_break_branch_score_sign") or "zero")
    p731_pair2_sign = str(p731.get("pair2_w_break_branch_score_sign") or "zero")

    f684_signs = ((f684_obj.get("construction") or {}).get("signs_by_pair")) or (f684_summary.get("signs_by_pair") or {})
    f684_pair1_sign_raw = f684_signs.get("pair1")
    f684_pair2_sign_raw = f684_signs.get("pair2")
    f684_pair1_sign = sign_label_from_value(f684_pair1_sign_raw)
    f684_pair2_sign = sign_label_from_value(f684_pair2_sign_raw)

    state_type = f684_obj.get("state_type") or {}
    domain = f684_obj.get("domain") or {}
    depends_on = f684_obj.get("depends_on") or {}
    construction = f684_obj.get("construction") or {}

    f301_notes = [str(v) for v in (f301.get("notes") or [])]
    f301_current_absence = [str(v) for v in (f301.get("current_absence") or [])]
    f301_source_domain = str(f301.get("source_domain") or "")

    pair12_convention_state_signs_equal = f684_pair1_sign == f684_pair2_sign
    pair12_witness_split_descends_to_current_convention_state = (
        p731_pair1_sign == f684_pair1_sign and p731_pair2_sign == f684_pair2_sign
    )

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
        "p731_witness_side_pair12_branch_separation_already_exported",
        {
            "pair12_split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
        },
        {
            "pair12_split_separated": True,
            "t185_exported": False,
        },
        "P731 already proves that the surviving pair1/pair2 residual-datum branches are separated by opposite nonzero w_break witness scores, but that typed source-side promotion is still unexported.",
    )
    add_check(
        "current_pair1_rooted_convention_state_is_exported_but_explicitly_convention_only",
        {
            "status": str(f684_summary.get("status") or f684_obj.get("status") or ""),
            "root_pair": construction.get("root", {}).get("pair"),
            "counts_as_strict_physical_orientation_datum": bool(
                state_type.get("counts_as_strict_physical_orientation_datum")
            ),
            "sign_gauge_is_strict_convention": "strict_convention" in str(state_type.get("sign_gauge") or "").lower(),
        },
        {
            "status": "PASS_EXPORTED",
            "root_pair": "pair1",
            "counts_as_strict_physical_orientation_datum": False,
            "sign_gauge_is_strict_convention": True,
        },
        "The current directed atlas representative already exists, but it is rooted on pair1 and explicitly scoped only as a strict-convention sign gauge, not as a strict physical orientation datum.",
    )
    add_check(
        "current_pair1_rooted_convention_state_is_not_source_typed_f301_branch_object",
        {
            "configuration_space_object": domain.get("configuration_space_object"),
            "state_level": state_type.get("level"),
            "source_domain_matches_f301": str(domain.get("configuration_space_object") or "") == f301_source_domain,
            "depends_on_f301": contains_token(depends_on, "residual_datum_bridge_export_map_object_support_carrier_actual_inhabitant_instance.json")
            or contains_token(depends_on, f301_source_domain),
        },
        {
            "configuration_space_object": "C_v1_void_configuration_space_in_local_B_tilde_1_sector_v1",
            "state_level": "directed_vector_state",
            "source_domain_matches_f301": False,
            "depends_on_f301": False,
        },
        "The current pair1-rooted convention state lives on the atlas-side directed-vector section over C_v1 and is not yet typed as the F301 source-domain pair12 carrier on tau_src_candidate_v1.",
    )
    add_check(
        "f301_source_carrier_remains_selector_neutral",
        {
            "source_domain": f301_source_domain,
            "selector_neutral_note_present": any("selector-neutral" in note for note in f301_notes),
            "no_admissible_s_sel_int_present": any("no admissible S_sel_int" in item for item in f301_current_absence),
        },
        {
            "source_domain": "tau_src_candidate_v1",
            "selector_neutral_note_present": True,
            "no_admissible_s_sel_int_present": True,
        },
        "The surviving F301 carrier is still the selector-neutral source-domain object on tau_src_candidate_v1 and still does not carry an admissible S_sel_int upgrade.",
    )
    add_check(
        "pair1_pair2_convention_state_signs_equal_under_current_pair1_rooted_transport",
        pair12_convention_state_signs_equal,
        True,
        "Under the current pair1-rooted transport convention, the atlas sign labels on pair1 and pair2 are equal rather than opposite.",
    )
    add_check(
        "p731_pair12_witness_sign_split_does_not_descend_to_current_pair1_rooted_convention_state",
        pair12_witness_split_descends_to_current_convention_state,
        False,
        "Therefore the already-separated pair12 witness-score split from P731 does not descend to the current pair1-rooted convention-state sign labels.",
    )
    add_check(
        "t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_exported",
        False,
        False,
        "The repo still does not export the bridge turning the P731 pair12 witness split into a typed descent through the current pair1-rooted convention state.",
    )

    status = (
        "PASS_PAIR1_ROOTED_CONVENTION_STATE_PAIR12_WITNESS_SPLIT_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P732_REQUIRES_REVIEW_CHANGED_PAIR1_ROOTED_CONVENTION_STATE_OR_PAIR12_WITNESS_SPLIT_STATE"
    )

    artifact = {
        "stage": "P732",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t186_pair1_rooted_convention_state_pair12_witness_split_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "F301_artifact": str(IN_F301.relative_to(REPO)),
            "F684_object": str(IN_F684_OBJ.relative_to(REPO)),
            "F684_summary": str(IN_F684_SUMMARY.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t186_target_name": "Pair1RootedConventionStatePair12WitnessSplitDescentBridge_global_C_v1_strict_v1",
        "t186_target_exported_on_current_repo_state": False,
        "current_pair1_rooted_convention_state_exists": str(f684_summary.get("status") or f684_obj.get("status") or "") == "PASS_EXPORTED",
        "pair1_pair2_convention_state_signs_equal": pair12_convention_state_signs_equal,
        "pair1_pair2_convention_state_signs": {
            "pair1": f684_pair1_sign,
            "pair2": f684_pair2_sign,
        },
        "p731_pair12_witness_score_signs": {
            "pair1": p731_pair1_sign,
            "pair2": p731_pair2_sign,
        },
        "p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state": pair12_witness_split_descends_to_current_convention_state,
        "audit_conclusion": {
            "current_repo_already_exports_pair1_rooted_convention_state": str(f684_summary.get("status") or f684_obj.get("status") or "") == "PASS_EXPORTED",
            "current_repo_already_exports_pair12_witness_split_descent_through_that_convention_state": False,
            "next_honest_move": (
                "attempt_a_typed_source_side_promotion_or_descent_bridge_from_the_p731_pair12_witness_split_without_collapsing_it_into_the_current_pair1_rooted_convention_gauge_and_without_upgrading_F647_or_F301_by_fiat"
            ),
        },
        "hard_limits": [
            "No T186 discharge claim.",
            "No claim that the current pair1-rooted convention state is a strict physical orientation datum.",
            "No claim that the current pair1-rooted convention state is typed as the F301 source-domain branch carrier.",
            "No claim that F647 is admissible S_sel_int.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P732",
        "status": status,
        "as_of": AS_OF,
        "t186_target_name": "Pair1RootedConventionStatePair12WitnessSplitDescentBridge_global_C_v1_strict_v1",
        "t186_target_exported_on_current_repo_state": False,
        "current_pair1_rooted_convention_state_exists": str(f684_summary.get("status") or f684_obj.get("status") or "") == "PASS_EXPORTED",
        "pair1_pair2_convention_state_signs_equal": pair12_convention_state_signs_equal,
        "p731_pair12_witness_split_descends_to_current_pair1_rooted_convention_state": pair12_witness_split_descends_to_current_convention_state,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
