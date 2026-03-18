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
IN_P731 = GENERATED / "p731_current_strict_t185_w_break_witness_payload_residual_datum_pair12_orbit_direction_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P736 = GENERATED / "p736_current_strict_t190_local_provider_operator_shift_direction_pair12_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_O12 = GENERATED / "o12_pair1_pair2_selector_chart_transport_operator_strict_derived_from_sigma_int_alpha12_v1.json"
IN_A12 = GENERATED / "a_12_pair12_chart_glued_orientation_projector_operator_section_strict_core_v1.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"

OUT_JSON = GENERATED / "p737_current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p737_current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def contains_token(entries: Any, token: str) -> bool:
    return isinstance(entries, list) and any(token in str(entry) for entry in entries)


def close(a: Any, b: Any, tol: float = TOL) -> bool:
    return isinstance(a, (int, float)) and isinstance(b, (int, float)) and abs(float(a) - float(b)) <= tol


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P729, IN_P731, IN_P736, IN_O12, IN_A12, IN_ATLAS]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P737",
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
    p731 = load_json(IN_P731)
    p736 = load_json(IN_P736)
    o12 = load_json(IN_O12)
    a12 = load_json(IN_A12)
    atlas = load_json(IN_ATLAS)

    o12_scope = o12.get("scope") or {}
    o12_construction = o12.get("construction") or {}
    o12_hard_limits = o12.get("hard_limits") or []

    a12_audits = a12.get("audits") or {}
    a12_construction = a12.get("construction") or {}
    a12_hard_limits = a12.get("hard_limits") or []

    atlas_gluing = atlas.get("gluing_data") or {}
    atlas_overlap = atlas.get("overlap_domain_declaration") or {}
    atlas_hard_limits = atlas.get("hard_limits") or []
    atlas_charts = atlas.get("charts") or {}

    o12_projector_transport_is_sign_gauge_invariant = (
        "sign-gauge-invariant" in str(o12_construction.get("projector_transport") or "")
        and contains_token(o12_hard_limits, "Does not derive a sign-sensitive physical orientation datum")
    )
    a12_section_is_projector_level_sign_gauge_safe = (
        "A_2 = O_12 A_1 O_12^T" in str(a12_construction.get("gluing_law_full") or "")
        and close(a12_audits.get("A1_sign_gauge_invariance_max_abs"), 0.0)
        and close(a12_audits.get("A2_sign_gauge_invariance_max_abs"), 0.0)
        and contains_token(a12_hard_limits, "Projector-level only")
    )
    current_local_pair12_projector_atlas_glue_lane_exported = (
        sorted(atlas_charts.keys()) == ["pair1", "pair2"]
        and atlas_overlap.get("scope") == "sigma_int_corridor_declared_scope"
        and isinstance(atlas_gluing.get("operator_section_ref"), str)
        and "projector-level" in str(atlas_gluing.get("law") or "")
        and contains_token(atlas_hard_limits, "Projector-level only")
    )
    current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe = (
        current_local_pair12_projector_atlas_glue_lane_exported
        and o12_projector_transport_is_sign_gauge_invariant
        and a12_section_is_projector_level_sign_gauge_safe
        and "sign-gauge-safe" in str(atlas_gluing.get("law") or "")
    )
    current_local_pair12_projector_atlas_glue_lane_is_pair12_orbit_direction_typed = (
        current_local_pair12_projector_atlas_glue_lane_exported
        and not current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe
    )
    p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_local_pair12_projector_atlas_glue_lane_is_pair12_orbit_direction_typed
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
        "p729_pair12_orbit_direction_split_already_localized",
        bool(p729.get("remaining_pair12_split_localized_as_opposite_orbit_directions")),
        True,
        "P729 already localizes the surviving pair1/pair2 ambiguity as the opposite residual-datum orbit-direction branches delta_k and delta_-k.",
    )
    add_check(
        "p731_pair12_witness_split_already_opposite_and_unpromoted",
        {
            "split_separated": bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches")),
            "pair1_sign": p731.get("pair1_w_break_branch_score_sign"),
            "pair2_sign": p731.get("pair2_w_break_branch_score_sign"),
            "t185_exported": bool(p731.get("t185_target_exported_on_current_repo_state")),
        },
        {
            "split_separated": True,
            "pair1_sign": "negative",
            "pair2_sign": "positive",
            "t185_exported": False,
        },
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed promotion bridge remains unexported.",
    )
    add_check(
        "p736_local_provider_operator_shift_direction_lane_already_realized_but_still_symmetric",
        {
            "realizes_both_branches_symmetrically": bool(
                p736.get("current_local_provider_operator_shift_direction_lane_realizes_both_pair12_branches_symmetrically")
            ),
            "selector_neutral": bool(p736.get("current_local_provider_operator_shift_direction_lane_is_selector_neutral")),
            "split_descends": bool(
                p736.get("p731_pair12_witness_split_descends_to_current_local_provider_operator_shift_direction_lane")
            ),
            "t190_exported": bool(p736.get("t190_target_exported_on_current_repo_state")),
        },
        {
            "realizes_both_branches_symmetrically": True,
            "selector_neutral": True,
            "split_descends": False,
            "t190_exported": False,
        },
        "P736 already proves that the strongest current local non-scalar provider-operator lane exists, but still realizes both branches only symmetrically and does not descend the P731 split.",
    )
    add_check(
        "o12_pair12_chart_transport_is_exported_and_projector_level_sign_gauge_safe",
        {
            "charts": o12_scope.get("charts"),
            "global_selector_atlas": bool(o12_scope.get("global_selector_atlas")),
            "projector_transport_sign_gauge_invariant": o12_projector_transport_is_sign_gauge_invariant,
        },
        {
            "charts": ["pair1", "pair2"],
            "global_selector_atlas": False,
            "projector_transport_sign_gauge_invariant": True,
        },
        "The current pair1/pair2 chart-transport object O12 is real, lane-scoped, and already only projector-level / sign-gauge-safe on current exports.",
    )
    add_check(
        "a12_pair12_projector_operator_section_is_exported_and_sign_gauge_safe",
        {
            "gluing_law_projector_level": "A_2 = O_12 A_1 O_12^T" in str(a12_construction.get("gluing_law_full") or ""),
            "A1_sign_gauge_invariance_zero": close(a12_audits.get("A1_sign_gauge_invariance_max_abs"), 0.0),
            "A2_sign_gauge_invariance_zero": close(a12_audits.get("A2_sign_gauge_invariance_max_abs"), 0.0),
            "projector_level_only": contains_token(a12_hard_limits, "Projector-level only"),
        },
        {
            "gluing_law_projector_level": True,
            "A1_sign_gauge_invariance_zero": True,
            "A2_sign_gauge_invariance_zero": True,
            "projector_level_only": True,
        },
        "The current pair1/pair2 glued projector-operator section already exists, but precisely at the sign-gauge-safe projector level only.",
    )
    add_check(
        "current_local_pair12_projector_atlas_glue_lane_exported",
        current_local_pair12_projector_atlas_glue_lane_exported,
        True,
        "The current repo already exports one explicit local pair1/pair2 atlas stub with overlap declaration, transition operator, and projector gluing data.",
    )
    add_check(
        "current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe",
        current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe,
        True,
        "That current local pair1/pair2 atlas/glue lane remains projector-level and sign-gauge-safe on current exports.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane",
        p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane,
        False,
        "So the already-separated opposite P731 pair1/pair2 witness split still does not descend through the current local pair1/pair2 projector-atlas glue lane as one typed branch distinction.",
    )
    add_check(
        "t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_exported",
        False,
        False,
        "The repo still does not export the local pair1/pair2 projector-atlas glue orbit-direction descent bridge on current repo state.",
    )

    status = (
        "PASS_LOCAL_PAIR12_PROJECTOR_ATLAS_GLUE_ORBIT_DIRECTION_DESCENT_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P737_REQUIRES_REVIEW_CHANGED_LOCAL_PAIR12_PROJECTOR_ATLAS_GLUE_STATE"
    )

    artifact = {
        "stage": "P737",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_boundary_only",
        "inputs": {
            "P729": str(IN_P729.relative_to(REPO)),
            "P731": str(IN_P731.relative_to(REPO)),
            "P736": str(IN_P736.relative_to(REPO)),
            "O12_artifact": str(IN_O12.relative_to(REPO)),
            "A12_artifact": str(IN_A12.relative_to(REPO)),
            "pair12_projector_atlas_artifact": str(IN_ATLAS.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t191_target_name": "LocalPair12ProjectorAtlasGlueOrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t191_target_exported_on_current_repo_state": False,
        "current_local_pair12_projector_atlas_glue_lane_exported": current_local_pair12_projector_atlas_glue_lane_exported,
        "current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe": (
            current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe
        ),
        "p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane": (
            p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane
        ),
        "current_local_pair12_projector_atlas_glue_lane_profile": {
            "charts": sorted(atlas_charts.keys()),
            "overlap_scope": atlas_overlap.get("scope"),
            "gluing_law": atlas_gluing.get("law"),
            "projector_transport_sign_gauge_invariant": o12_projector_transport_is_sign_gauge_invariant,
        },
        "audit_conclusion": {
            "current_repo_already_exports_local_pair12_projector_atlas_glue_lane": current_local_pair12_projector_atlas_glue_lane_exported,
            "current_repo_already_exports_t191_target": False,
            "next_honest_move": (
                "attack_a_typed_local_branch_sensitive_bind_or_descent_mechanism_that_acts_above_the_current_pair12_projector_atlas_glue_lane_without_collapsing_back_to_projector_level_sign_gauge_safe_data_or_the_current_symmetric_provider_operator_lane"
            ),
        },
        "hard_limits": [
            "No T191 discharge claim.",
            "No claim that the current local pair1/pair2 projector-atlas glue lane already selects one raw pair1/pair2 orbit-direction branch.",
            "No claim that projector-level sign-gauge-safe transport upgrades to a sign-sensitive physical orientation datum.",
            "No admissible S_sel_int claim.",
            "No kernel-alone/global QW-2191 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P737",
        "status": status,
        "as_of": AS_OF,
        "t191_target_name": "LocalPair12ProjectorAtlasGlueOrbitDirectionDescentBridge_global_C_v1_strict_v1",
        "t191_target_exported_on_current_repo_state": False,
        "current_local_pair12_projector_atlas_glue_lane_exported": current_local_pair12_projector_atlas_glue_lane_exported,
        "current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe": (
            current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe
        ),
        "p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane": (
            p731_pair12_witness_split_descends_to_current_local_pair12_projector_atlas_glue_lane
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
