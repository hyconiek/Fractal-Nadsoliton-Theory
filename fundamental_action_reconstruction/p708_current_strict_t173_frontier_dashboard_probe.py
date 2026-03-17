#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

# Core T173 frontier summaries (strict).
IN_N679 = GENERATED / "n679_current_strict_t172_strict_core_selector_closure_frontier_boundary_theorem_summary.json"
IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N681 = GENERATED / "n681_current_strict_t173_directed_output_sign_lift_obstruction_boundary_theorem_summary.json"
IN_N686 = GENERATED / "n686_current_strict_t173_global_axis_only_transition_edge_sign_flip_boundary_theorem_summary.json"
IN_N687 = GENERATED / "n687_current_strict_t173_global_edge_sign_coherence_obstruction_boundary_theorem_summary.json"
IN_P711 = GENERATED / "p711_current_strict_t173_previous_methodology_survival_and_global_gap_audit_probe_summary.json"
IN_P712 = GENERATED / "p712_current_strict_t176_existing_global_directed_sign_coherence_provider_nonexport_audit_probe_summary.json"
IN_P713 = GENERATED / "p713_current_strict_t176_multiroot_rooted_sign_lift_root_independence_audit_probe_summary.json"

# Convention-layer continuations (still below physical sign datum).
IN_N688 = GENERATED / "n688_current_strict_t174_global_oriented_transition_edge_sign_lift_discharge_theorem_summary.json"
IN_N690 = GENERATED / "n690_current_strict_t175_global_chart_sign_fixing_discharge_theorem_summary.json"
IN_N691 = GENERATED / "n691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_discharge_theorem_summary.json"

# Optional global dashboards for convenience.
IN_P441 = GENERATED / "p441_current_strict_global_closure_next_move_dashboard_probe_summary.json"
IN_P706 = GENERATED / "p706_current_release_7_strict_projective_operational_toe_os_closure_dashboard_probe_summary.json"

OUT_JSON = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe.json"
OUT_SUMMARY = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def tr(obj: dict[str, Any]) -> dict[str, Any]:
    val = obj.get("theorem_result")
    return val if isinstance(val, dict) else {}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    core = {
        "N679": IN_N679,
        "N680": IN_N680,
        "N681": IN_N681,
        "N686": IN_N686,
        "N687": IN_N687,
        "P711": IN_P711,
        "P712": IN_P712,
    }
    missing_core = [str(p.relative_to(REPO)) for p in core.values() if not p.exists()]
    if missing_core:
        artifact: dict[str, Any] = {
            "stage": "P708",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing_core,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    n679 = load_json(IN_N679)
    n680 = load_json(IN_N680)
    n681 = load_json(IN_N681)
    n686 = load_json(IN_N686)
    n687 = load_json(IN_N687)
    p711 = load_json(IN_P711)
    p712 = load_json(IN_P712)
    p713 = load_json(IN_P713) if IN_P713.exists() else None

    n688 = load_json(IN_N688) if IN_N688.exists() else None
    n690 = load_json(IN_N690) if IN_N690.exists() else None
    n691 = load_json(IN_N691) if IN_N691.exists() else None

    p441 = load_json(IN_P441) if IN_P441.exists() else None
    p706 = load_json(IN_P706) if IN_P706.exists() else None

    n679_tr = tr(n679)
    n680_tr = tr(n680)
    n681_tr = tr(n681)
    n686_tr = tr(n686)
    n687_tr = tr(n687)
    n688_tr = tr(n688) if isinstance(n688, dict) else {}
    n690_tr = tr(n690) if isinstance(n690, dict) else {}
    n691_tr = tr(n691) if isinstance(n691, dict) else {}

    recommended_next = "T173"
    if isinstance(p441, dict):
        recommended_next = str(p441.get("recommended_next_strict_target") or recommended_next)

    checks: list[dict[str, Any]] = []
    blocking: list[str] = []

    def add_check(check_id: str, actual: Any, expected: Any, meaning: str) -> None:
        ok = actual == expected
        checks.append(
            {
                "id": check_id,
                "actual": actual,
                "expected": expected,
                "pass": ok,
                "meaning": meaning,
            }
        )
        if not ok:
            blocking.append(check_id)

    # Post-T172 frontier facts.
    add_check(
        "projective_strict_core_selector_closure_discharged",
        bool(n680_tr.get("strict_core_selector_closure")),
        True,
        "Projective strict-core selector closure is discharged (N680).",
    )
    add_check(
        "kernel_alone_global_qw2191_discharge_false",
        bool(n679_tr.get("QW2191_kernel_alone_discharge") or n680_tr.get("QW2191_kernel_alone_discharge")),
        False,
        "Kernel-alone/global QW-2191 discharge remains false/unclaimed (N679/N680 hard limits).",
    )
    add_check(
        "directed_output_sign_lift_determined_in_strict_core_false",
        bool(n681_tr.get("directed_output_sign_lift_determined_in_strict_core")),
        False,
        "Directed output sign-lift is not determined in strict core (N681).",
    )
    add_check(
        "axis_only_transition_edge_sign_flips_present_true",
        bool(n686_tr.get("global_axis_only_transition_edge_sign_flips_present")),
        True,
        "Under fixed axis-only (α mod π) transitions, some overlap edges force a sign flip (N686).",
    )
    add_check(
        "global_edge_sign_coherence_solvable_by_chart_sign_relift_false",
        bool(n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")),
        False,
        "No per-chart Z2 relift solves full-edge sign coherence under fixed axis-only transitions (N687).",
    )
    add_check(
        "previous_methodology_local_survival_audited",
        bool(p711.get("previous_methodology_contains_reusable_strict_ingredients")),
        True,
        "The previous sigma_int/topological methodology still contributes reusable strict ingredients (P711).",
    )
    add_check(
        "previous_methodology_not_global_t173_discharge",
        bool(p711.get("previous_methodology_suffices_for_global_t173_discharge")),
        False,
        "That surviving previous methodology is still insufficient for a full global T173 discharge (P711).",
    )
    add_check(
        "t176_global_provider_not_yet_exported",
        bool(p712.get("t176_target_exported_on_current_repo_state")),
        False,
        "The exact next global provider target (T176-class) is not yet exported on current repo state (P712).",
    )
    # Convention-layer continuations (optional but expected on Release 7 state).
    if IN_N688.exists():
        add_check(
            "t174_oriented_edge_sign_lift_exported",
            bool(n688_tr.get("oriented_edge_sign_lift_exported")),
            True,
            "A convention-layer oriented edge sign-lift is exported (T174 via N688).",
        )
    if IN_N690.exists():
        add_check(
            "t175_chart_sign_fix_exported",
            bool(n690_tr.get("sign_fixed_directed_representative_exported")),
            True,
            "A convention-layer chart sign-fix (0-cochain) and sign-fixed directed representative are exported (T175 via N690).",
        )
    if IN_N691.exists():
        add_check(
            "t174_oriented_edge_sign_lift_from_sign_fixed_state_exported",
            bool(n691_tr.get("oriented_edge_sign_lift_exported")),
            True,
            "A convention-layer oriented edge sign-lift anchored to the sign-fixed state is exported (T174 via N691).",
        )

    # Hard-limit consistency: nothing here upgrades to ToE closure.
    add_check(
        "no_ToE_closure_claim",
        bool(
            n679_tr.get("ToE_closure")
            or n680_tr.get("ToE_closure")
            or n681_tr.get("ToE_closure")
            or n686_tr.get("ToE_closure")
            or n687_tr.get("ToE_closure")
            or n688_tr.get("ToE_closure")
            or n690_tr.get("ToE_closure")
            or n691_tr.get("ToE_closure")
        ),
        False,
        "All referenced theorem summaries keep ToE closure false (hard limits).",
    )

    ok = len(blocking) == 0
    status = "PASS_T173_FRONTIER_DASHBOARD_READY" if ok else "P708_REQUIRES_REVIEW_CHANGED_OR_INCOMPLETE_T173_FRONTIER_STATE"

    artifact = {
        "stage": "P708",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t173_frontier_dashboard_only",
        "inputs": {
            "N679": str(IN_N679.relative_to(REPO)),
            "N680": str(IN_N680.relative_to(REPO)),
            "N681": str(IN_N681.relative_to(REPO)),
            "N686": str(IN_N686.relative_to(REPO)),
            "N687": str(IN_N687.relative_to(REPO)),
            "P711": str(IN_P711.relative_to(REPO)),
            "P712": str(IN_P712.relative_to(REPO)),
            "P713": str(IN_P713.relative_to(REPO)) if IN_P713.exists() else None,
            "N688": str(IN_N688.relative_to(REPO)) if IN_N688.exists() else None,
            "N690": str(IN_N690.relative_to(REPO)) if IN_N690.exists() else None,
            "N691": str(IN_N691.relative_to(REPO)) if IN_N691.exists() else None,
            "P441": str(IN_P441.relative_to(REPO)) if IN_P441.exists() else None,
            "P706": str(IN_P706.relative_to(REPO)) if IN_P706.exists() else None,
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "recommended_next_strict_target": recommended_next,
        "t173_frontier_state": {
            "strict_core_selector_closure_projective": bool(n680_tr.get("strict_core_selector_closure")),
            "QW2191_kernel_alone_discharge": bool(
                n679_tr.get("QW2191_kernel_alone_discharge") or n680_tr.get("QW2191_kernel_alone_discharge")
            ),
            "directed_output_sign_lift_determined_in_strict_core": bool(n681_tr.get("directed_output_sign_lift_determined_in_strict_core")),
            "global_edge_sign_coherence_solvable_by_chart_sign_relift_under_axis_only_transitions": bool(
                n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")
            ),
            "previous_methodology_contains_reusable_strict_ingredients": bool(
                p711.get("previous_methodology_contains_reusable_strict_ingredients")
            ),
            "previous_methodology_suffices_for_global_t173_discharge": bool(
                p711.get("previous_methodology_suffices_for_global_t173_discharge")
            ),
            "t176_global_provider_exported": bool(p712.get("t176_target_exported_on_current_repo_state")),
            "all_root_independent_convention_section_exists": (
                bool(p713.get("root_independent_sign_vector")) and bool(p713.get("root_independent_output_vectors"))
                if isinstance(p713, dict)
                else None
            ),
            "supported_root_corridor_with_matching_convention_section": (
                bool(p713.get("supported_roots_sign_vector_agree")) and bool(p713.get("supported_roots_output_vectors_agree"))
                if isinstance(p713, dict)
                else None
            ),
            "supported_roots_for_current_w_break_candidate": p713.get("supported_roots") if isinstance(p713, dict) else None,
            "convention_layer_oriented_edge_sign_lift_exported": bool(n688_tr.get("oriented_edge_sign_lift_exported") or n691_tr.get("oriented_edge_sign_lift_exported")),
            "convention_layer_sign_fixed_directed_representative_exported": bool(n690_tr.get("sign_fixed_directed_representative_exported")),
            "operational_release_7_projective_os_closure_dashboard_status": (p706 or {}).get("status") if isinstance(p706, dict) else None,
        },
        "hard_limits": [
            "No kernel-alone/global QW-2191 discharge.",
            "No directed/sign-sensitive physical orientation datum promotion into strict core.",
            "No Standard Model host-matching claim in strict scope.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P708",
        "status": status,
        "as_of": AS_OF,
        "recommended_next_strict_target": recommended_next,
        "strict_core_selector_closure_projective": bool(n680_tr.get("strict_core_selector_closure")),
        "QW2191_kernel_alone_discharge": False,
        "directed_sign_sensitive_physical_orientation_in_strict_core": False,
        "previous_methodology_contains_reusable_strict_ingredients": bool(
            p711.get("previous_methodology_contains_reusable_strict_ingredients")
        ),
        "previous_methodology_suffices_for_global_t173_discharge": bool(
            p711.get("previous_methodology_suffices_for_global_t173_discharge")
        ),
        "t176_global_provider_exported": bool(p712.get("t176_target_exported_on_current_repo_state")),
        "all_root_independent_convention_section_exists": (
            bool(p713.get("root_independent_sign_vector")) and bool(p713.get("root_independent_output_vectors"))
            if isinstance(p713, dict)
            else None
        ),
        "supported_root_corridor_with_matching_convention_section": (
            bool(p713.get("supported_roots_sign_vector_agree")) and bool(p713.get("supported_roots_output_vectors_agree"))
            if isinstance(p713, dict)
            else None
        ),
        "supported_roots_for_current_w_break_candidate": p713.get("supported_roots") if isinstance(p713, dict) else None,
        "convention_layer_sign_tools_exported": {
            "T174_oriented_edge_sign_lift": bool(n688_tr.get("oriented_edge_sign_lift_exported") or n691_tr.get("oriented_edge_sign_lift_exported")),
            "T175_chart_sign_fix": bool(n690_tr.get("sign_fixed_directed_representative_exported")),
        },
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
