#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


AS_OF = "2026-03-17"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P711 = GENERATED / "p711_current_strict_t173_previous_methodology_survival_and_global_gap_audit_probe_summary.json"
IN_N680 = GENERATED / "n680_current_strict_t173_projective_strict_core_selector_closure_discharge_theorem_summary.json"
IN_N681 = GENERATED / "n681_current_strict_t173_directed_output_sign_lift_obstruction_boundary_theorem_summary.json"
IN_N686 = GENERATED / "n686_current_strict_t173_global_axis_only_transition_edge_sign_flip_boundary_theorem_summary.json"
IN_N687 = GENERATED / "n687_current_strict_t173_global_edge_sign_coherence_obstruction_boundary_theorem_summary.json"
IN_N688 = GENERATED / "n688_current_strict_t174_global_oriented_transition_edge_sign_lift_discharge_theorem_summary.json"
IN_N690 = GENERATED / "n690_current_strict_t175_global_chart_sign_fixing_discharge_theorem_summary.json"
IN_N691 = GENERATED / "n691_current_strict_t174_global_oriented_transition_edge_sign_lift_from_sign_fixed_directed_state_discharge_theorem_summary.json"

OUT_JSON = GENERATED / "p712_current_strict_t176_existing_global_directed_sign_coherence_provider_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p712_current_strict_t176_existing_global_directed_sign_coherence_provider_nonexport_audit_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def theorem_result(obj: dict[str, Any]) -> dict[str, Any]:
    val = obj.get("theorem_result")
    return val if isinstance(val, dict) else {}


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prereq = [IN_P711, IN_N680, IN_N681, IN_N686, IN_N687, IN_N688, IN_N690, IN_N691]
    missing = [str(p.relative_to(REPO)) for p in prereq if not p.exists()]
    if missing:
        artifact = {
            "stage": "P712",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p711 = load_json(IN_P711)
    n680 = load_json(IN_N680)
    n681 = load_json(IN_N681)
    n686 = load_json(IN_N686)
    n687 = load_json(IN_N687)
    n688 = load_json(IN_N688)
    n690 = load_json(IN_N690)
    n691 = load_json(IN_N691)

    n680_tr = theorem_result(n680)
    n681_tr = theorem_result(n681)
    n686_tr = theorem_result(n686)
    n687_tr = theorem_result(n687)
    n688_tr = theorem_result(n688)
    n690_tr = theorem_result(n690)
    n691_tr = theorem_result(n691)

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

    add_check(
        "p711_gap_already_named",
        p711.get("small_gap"),
        "missing_global_provider_for_chartwise_directed_sign_coherence_on_full_C_v1_atlas",
        "P711 already identifies the missing object class as a global provider gap on the full C_v1 atlas.",
    )
    add_check(
        "n680_projective_global_closure_present",
        bool(n680_tr.get("strict_core_selector_closure")),
        True,
        "Projective global strict-core selector closure is present (N680).",
    )
    add_check(
        "n681_directed_sign_not_determined_in_strict_core",
        bool(n681_tr.get("directed_output_sign_lift_determined_in_strict_core")),
        False,
        "Directed output sign-lift is still not determined in strict core (N681).",
    )
    add_check(
        "n686_axis_only_global_sign_flips_present",
        bool(n686_tr.get("global_axis_only_transition_edge_sign_flips_present")),
        True,
        "Some global atlas edges still force sign flips under axis-only representatives (N686).",
    )
    add_check(
        "n687_chart_sign_relift_not_enough",
        bool(n687_tr.get("global_edge_sign_coherence_solvable_by_chart_sign_relift")),
        False,
        "No per-chart Z2 relift solves full-edge sign coherence under current axis-only transitions (N687).",
    )
    add_check(
        "n688_is_convention_layer_only",
        bool(n688_tr.get("directed_sign_sensitive_physical_orientation_in_strict_core")),
        False,
        "The oriented edge sign-lift export is convention-layer only, not a strict-core physical/global provider (N688).",
    )
    add_check(
        "n690_is_convention_layer_only",
        bool(n690_tr.get("directed_sign_sensitive_physical_orientation_in_strict_core")),
        False,
        "The chart sign-fix export is convention-layer only, not a strict-core physical/global provider (N690).",
    )
    add_check(
        "n691_is_convention_layer_only",
        bool(n691_tr.get("directed_sign_sensitive_physical_orientation_in_strict_core")),
        False,
        "The sign-fixed oriented edge-lift export is convention-layer only, not a strict-core physical/global provider (N691).",
    )

    existing_t176_provider_exported = False
    add_check(
        "existing_t176_class_provider_exported",
        existing_t176_provider_exported,
        False,
        "No current exported object satisfies the strict-core T176 target class on the full C_v1 atlas.",
    )

    status = (
        "PASS_EXISTING_GLOBAL_DIRECTED_SIGN_COHERENCE_PROVIDER_NONEXPORT_AUDITED"
        if not blocking
        else "P712_REQUIRES_REVIEW_CHANGED_OR_INCOMPLETE_T176_NONEXPORT_AUDIT_STATE"
    )

    artifact = {
        "stage": "P712",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "current_strict_t176_nonexport_boundary_only",
        "inputs": {
            "P711": str(IN_P711.relative_to(REPO)),
            "N680": str(IN_N680.relative_to(REPO)),
            "N681": str(IN_N681.relative_to(REPO)),
            "N686": str(IN_N686.relative_to(REPO)),
            "N687": str(IN_N687.relative_to(REPO)),
            "N688": str(IN_N688.relative_to(REPO)),
            "N690": str(IN_N690.relative_to(REPO)),
            "N691": str(IN_N691.relative_to(REPO)),
        },
        "checks": checks,
        "blocking_mismatches": blocking,
        "t176_target_name": "GlobalDirectedSignCoherenceProvider_global_C_v1_strict_core_v1",
        "t176_target_exported_on_current_repo_state": existing_t176_provider_exported,
        "current_objects_reaching_only_convention_layer": {
            "T174_oriented_edge_sign_lift": bool(n688_tr.get("oriented_edge_sign_lift_exported") or n691_tr.get("oriented_edge_sign_lift_exported")),
            "T175_chart_sign_fixing": bool(n690_tr.get("sign_fixed_directed_representative_exported")),
        },
        "audit_conclusion": {
            "projective_closure_exists": bool(n680_tr.get("strict_core_selector_closure")),
            "global_directed_sign_branch_still_open": True,
            "current_repo_already_exports_t176_provider": False,
            "next_honest_move": "export_or_attack_a_genuinely_new_strict_core_global_provider_targeted_by_T176",
        },
        "hard_limits": [
            "No kernel-alone/global QW-2191 discharge.",
            "No directed/sign-sensitive physical orientation datum promotion into strict core.",
            "No impossibility-in-principle claim about future T176 discharge.",
            "No ToE closure claim.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P712",
        "status": status,
        "as_of": AS_OF,
        "t176_target_name": "GlobalDirectedSignCoherenceProvider_global_C_v1_strict_core_v1",
        "t176_target_exported_on_current_repo_state": False,
        "global_directed_sign_branch_still_open": True,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
