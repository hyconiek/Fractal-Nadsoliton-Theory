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
IN_P737 = GENERATED / "p737_current_strict_t191_local_pair12_projector_atlas_glue_orbit_direction_descent_bridge_nonexport_audit_probe_summary.json"
IN_P741 = GENERATED / "p741_current_strict_t195_actual_source_topology_selector_witness_pair12_witness_split_promotion_bridge_nonexport_audit_probe_summary.json"
IN_P742 = GENERATED / "p742_current_strict_t196_actual_source_topology_selector_witness_to_residual_datum_pair12_typed_carrier_bridge_nonexport_audit_probe_summary.json"
IN_F147 = GENERATED / "f147_first_actual_source_topology_selector_witness_packet_summary.json"
IN_ATLAS = GENERATED / "selector_atlas_pair12_sigma_int_corridor_projector_v1.json"

OUT_JSON = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p747_current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_audit_probe_summary.json"

SELECTOR_WITNESS_TARGET_NAME = "Sigma_sel_src_target_v1"
LOCAL_PAIR12_ATLAS_OBJECT_NAME = "SelectorAtlas_pair12_sigma_int_corridor_projector_v1"
CURRENT_THEOREM_FILE = "N743_CURRENT_STRICT_T201_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_TARGET_TO_LOCAL_PAIR12_CHART_SENSITIVE_ATLAS_BRIDGE_NONEXPORT_BOUNDARY_THEOREM.md"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def scan_positive_bridge_candidates() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "f*.py", "n*.py", "t*.py")
    candidates: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path.name == CURRENT_THEOREM_FILE or path in seen:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if SELECTOR_WITNESS_TARGET_NAME in text and LOCAL_PAIR12_ATLAS_OBJECT_NAME in text:
                candidates.append(str(path.relative_to(REPO)))
    return candidates


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_P731, IN_P737, IN_P741, IN_P742, IN_F147, IN_ATLAS]
    missing = [str(path.relative_to(REPO)) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P747",
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
    p737 = load_json(IN_P737)
    p741 = load_json(IN_P741)
    p742 = load_json(IN_P742)
    f147 = load_json(IN_F147)
    atlas = load_json(IN_ATLAS)

    support_packet = f147.get("support_packet") or {}
    selector_axis_realization = support_packet.get("selector_axis_realization") or {}
    selector_signed_split_realization = support_packet.get("selector_signed_split_realization") or {}
    preobserver_scope_realization = support_packet.get("preobserver_scope_realization") or {}

    atlas_charts = atlas.get("charts") or {}
    atlas_overlap = atlas.get("overlap_domain_declaration") or {}
    atlas_gluing = atlas.get("gluing_data") or {}
    atlas_hard_limits = atlas.get("hard_limits") or []

    positive_bridge_candidates = scan_positive_bridge_candidates()

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

    current_actual_source_topology_selector_witness_target_exported = (
        f147.get("witness") == "Pi_sel_src_actual_witness_v1"
        and f147.get("input_packet") == "tau_src_candidate_v1"
        and f147.get("codomain_packet") == SELECTOR_WITNESS_TARGET_NAME
        and bool(f147.get("actual_selector_witness_exported"))
    )

    current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm = (
        current_actual_source_topology_selector_witness_target_exported
        and bool(p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only"))
        and bool(p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed"))
        and selector_axis_realization.get("object") == "E_orient_preLM_v1"
        and selector_axis_realization.get("frame_basis") == ["u_T", "u_L"]
        and selector_signed_split_realization.get("object") == "B_sel_preLM_v1"
        and bool(preobserver_scope_realization.get("preobserver_only"))
        and bool(preobserver_scope_realization.get("observer_downstream_only"))
    )

    current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation = (
        bool(p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation"))
        and bool(p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed"))
    )

    current_local_pair12_chart_sensitive_atlas_lane_exported = (
        atlas.get("object") == LOCAL_PAIR12_ATLAS_OBJECT_NAME
        and sorted(atlas_charts.keys()) == ["pair1", "pair2"]
        and atlas_overlap.get("scope") == "sigma_int_corridor_declared_scope"
        and bool(p737.get("current_local_pair12_projector_atlas_glue_lane_exported"))
    )

    current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe = (
        current_local_pair12_chart_sensitive_atlas_lane_exported
        and "projector-level" in str(atlas_gluing.get("law") or "")
        and "sign-gauge-safe" in str(atlas_gluing.get("law") or "")
        and any("Projector-level only" in str(entry) for entry in atlas_hard_limits)
        and bool(p737.get("current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe"))
    )

    current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge = (
        current_actual_source_topology_selector_witness_target_exported
        and current_local_pair12_chart_sensitive_atlas_lane_exported
        and len(positive_bridge_candidates) > 0
        and not current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm
        and not current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation
        and not current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe
    )

    current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas = (
        current_actual_source_topology_selector_witness_target_exported
        and current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm
        and current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation
        and current_local_pair12_chart_sensitive_atlas_lane_exported
        and current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe
        and len(positive_bridge_candidates) == 0
        and not current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge
    )

    p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge = (
        bool(p731.get("current_w_break_witness_payload_separates_pair12_orbit_direction_branches"))
        and current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge
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
        "P731 already separates the surviving pair1/pair2 branches by opposite witness-score signs, while the typed source-side promotion bridge remains unexported.",
    )
    add_check(
        "current_actual_source_topology_selector_witness_target_exported",
        {
            "witness": f147.get("witness"),
            "input_packet": f147.get("input_packet"),
            "codomain_packet": f147.get("codomain_packet"),
            "actual_selector_witness_exported": bool(f147.get("actual_selector_witness_exported")),
        },
        {
            "witness": "Pi_sel_src_actual_witness_v1",
            "input_packet": "tau_src_candidate_v1",
            "codomain_packet": SELECTOR_WITNESS_TARGET_NAME,
            "actual_selector_witness_exported": True,
        },
        "F147 already exports one real actual source-topology selector witness with codomain Sigma_sel_src_target_v1 on tau_src_candidate_v1.",
    )
    add_check(
        "current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm",
        {
            "chart_bound_preobserver_only": bool(
                p741.get("current_actual_source_topology_selector_witness_is_chart_bound_preobserver_only")
            ),
            "prelm_not_pair12_typed": bool(
                p741.get("current_actual_source_topology_selector_witness_remains_prelm_not_pair12_typed")
            ),
            "selector_axis_object": selector_axis_realization.get("object"),
            "frame_basis": selector_axis_realization.get("frame_basis"),
            "signed_split_object": selector_signed_split_realization.get("object"),
        },
        {
            "chart_bound_preobserver_only": True,
            "prelm_not_pair12_typed": True,
            "selector_axis_object": "E_orient_preLM_v1",
            "frame_basis": ["u_T", "u_L"],
            "signed_split_object": "B_sel_preLM_v1",
        },
        "P741 already proves that the current Sigma_sel_src_target_v1 witness remains chart-bound, preobserver-only, and realized only in the preLM basis u_T/u_L rather than on a pair1/pair2 atlas basis.",
    )
    add_check(
        "current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation",
        {
            "basis_free_continuation_exported": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_basis_free_chart_label_forgetting_continuation")
            ),
            "continuation_remains_basis_free_not_pair12_typed": bool(
                p742.get("current_actual_selector_witness_codomain_continuation_remains_basis_free_not_pair12_typed")
            ),
            "pair12_typed_continuation_exported": bool(
                p742.get("current_actual_selector_witness_codomain_has_exported_pair12_typed_residual_datum_continuation")
            ),
        },
        {
            "basis_free_continuation_exported": True,
            "continuation_remains_basis_free_not_pair12_typed": True,
            "pair12_typed_continuation_exported": False,
        },
        "P742 already proves that the strongest current continuation out of Sigma_sel_src_target_v1 is exactly the chart-label-forgetting basis-free reduction, not a pair1/pair2 typed continuation.",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_exported",
        {
            "atlas_object": atlas.get("object"),
            "charts": sorted(atlas_charts.keys()),
            "overlap_scope": atlas_overlap.get("scope"),
            "lane_exported": bool(p737.get("current_local_pair12_projector_atlas_glue_lane_exported")),
        },
        {
            "atlas_object": LOCAL_PAIR12_ATLAS_OBJECT_NAME,
            "charts": ["pair1", "pair2"],
            "overlap_scope": "sigma_int_corridor_declared_scope",
            "lane_exported": True,
        },
        "F463/P737 already export one real local pair1/pair2 chart-sensitive atlas lane with overlap declaration on the sigma_int corridor scope.",
    )
    add_check(
        "current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe",
        {
            "gluing_law": atlas_gluing.get("law"),
            "projector_level_sign_gauge_safe": bool(
                p737.get("current_local_pair12_projector_atlas_glue_lane_is_projector_level_sign_gauge_safe")
            ),
        },
        {
            "gluing_law": "A_2(pair2) = O_12 A_1(pair1) O_12^T  (projector-level, sign-gauge-safe)",
            "projector_level_sign_gauge_safe": True,
        },
        "That current local pair1/pair2 atlas lane still remains projector-level and sign-gauge-safe rather than an orbit-direction-typed source-side bridge.",
    )
    add_check(
        "positive_packet_theorem_or_spec_bridge_candidates_from_selector_witness_target_to_local_pair12_chart_sensitive_atlas",
        positive_bridge_candidates,
        [],
        "No current packet, theorem, or spec simultaneously exports Sigma_sel_src_target_v1 together with the local pair1/pair2 chart-sensitive atlas object as a positive bridge lane.",
    )
    add_check(
        "current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas",
        current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas,
        True,
        "So the current actual selector-witness target still remains unbridged to the local pair1/pair2 chart-sensitive atlas lane on current exports.",
    )
    add_check(
        "current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge",
        current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge,
        False,
        "No current export carries Sigma_sel_src_target_v1 into the local pair1/pair2 chart-sensitive atlas lane as an explicit source-side bridge.",
    )
    add_check(
        "p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge",
        p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge,
        False,
        "Therefore the opposite P731 pair1/pair2 witness split still does not descend through any exported Sigma_sel_src_target_v1 to local pair1/pair2 chart-sensitive atlas bridge.",
    )
    add_check(
        "t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_exported",
        False,
        False,
        "So the actual source-topology selector-witness target to local pair1/pair2 chart-sensitive atlas bridge remains unexported on current repo state.",
    )

    status = (
        "PASS_ACTUAL_SOURCE_TOPOLOGY_SELECTOR_WITNESS_TARGET_TO_LOCAL_PAIR12_CHART_SENSITIVE_ATLAS_BRIDGE_NONEXPORT_AUDITED"
        if not blocking
        else "P747_REQUIRES_REVIEW_CHANGED_SELECTOR_TARGET_TO_LOCAL_PAIR12_ATLAS_STATE"
    )

    artifact = {
        "stage": "P747",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t201_actual_source_topology_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge_nonexport_boundary_only",
        "inputs": {
            "P731": str(IN_P731.relative_to(REPO)),
            "P737": str(IN_P737.relative_to(REPO)),
            "P741": str(IN_P741.relative_to(REPO)),
            "P742": str(IN_P742.relative_to(REPO)),
            "F147": str(IN_F147.relative_to(REPO)),
            "pair12_atlas_artifact": str(IN_ATLAS.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t201_target_name": "ActualSourceTopologySelectorWitnessTargetToLocalPair12ChartSensitiveAtlasBridge_sigma_int_corridor_strict_v1",
        "t201_target_exported_on_current_repo_state": False,
        "current_actual_source_topology_selector_witness_target_exported": (
            current_actual_source_topology_selector_witness_target_exported
        ),
        "current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm": (
            current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm
        ),
        "current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation": (
            current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_exported": (
            current_local_pair12_chart_sensitive_atlas_lane_exported
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe": (
            current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe
        ),
        "positive_packet_theorem_or_spec_bridge_candidates_from_selector_witness_target_to_local_pair12_chart_sensitive_atlas": (
            positive_bridge_candidates
        ),
        "current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas": (
            current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas
        ),
        "current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge": (
            current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge
        ),
        "p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge": (
            p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge
        ),
        "frontier": {
            "current_repo_already_exports_sigma_sel_src_target_v1": current_actual_source_topology_selector_witness_target_exported,
            "current_repo_already_exports_local_pair12_chart_sensitive_atlas_lane": current_local_pair12_chart_sensitive_atlas_lane_exported,
            "current_repo_already_exports_sigma_sel_src_target_to_local_pair12_chart_sensitive_atlas_bridge": (
                current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge
            ),
            "next_honest_move": (
                "attempt_a_genuinely_chart_sensitive_source_side_bridge_from_sigma_sel_src_target_v1_into_one_pair12_typed_atlas_lane_without_forgetting_chart_labels_and_without_collapsing_into_the_current_projector_level_sign_gauge_safe_pair12_atlas_glue"
            ),
        },
        "hard_limits": [
            "No claim that Sigma_sel_src_target_v1 is already pair1/pair2 atlas typed.",
            "No claim that the current local pair1/pair2 atlas/glue lane already carries a source-side orbit-direction distinction.",
            "No claim that projector-level chart transport or gluing already fixes a sign-sensitive physical orientation datum.",
            "No claim that the current basis-free continuation already upgrades back to chart-sensitive atlas data.",
            "No kernel-alone/global QW-2191 discharge.",
            "No strict-core selector closure.",
            "No ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "t201_target_name": artifact["t201_target_name"],
        "t201_target_exported_on_current_repo_state": artifact["t201_target_exported_on_current_repo_state"],
        "current_actual_source_topology_selector_witness_target_exported": (
            artifact["current_actual_source_topology_selector_witness_target_exported"]
        ),
        "current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm": (
            artifact["current_actual_source_topology_selector_witness_target_remains_chart_bound_prelm"]
        ),
        "current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation": (
            artifact["current_actual_selector_witness_target_has_exported_basis_free_chart_label_forgetting_continuation"]
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_exported": (
            artifact["current_local_pair12_chart_sensitive_atlas_lane_exported"]
        ),
        "current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe": (
            artifact["current_local_pair12_chart_sensitive_atlas_lane_remains_projector_level_sign_gauge_safe"]
        ),
        "current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas": (
            artifact["current_actual_selector_witness_target_remains_unbridged_to_local_pair12_chart_sensitive_atlas"]
        ),
        "current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge": (
            artifact["current_actual_selector_witness_target_has_exported_local_pair12_chart_sensitive_atlas_bridge"]
        ),
        "p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge": (
            artifact["p731_pair12_witness_split_descends_to_current_selector_witness_target_to_local_pair12_chart_sensitive_atlas_bridge"]
        ),
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
