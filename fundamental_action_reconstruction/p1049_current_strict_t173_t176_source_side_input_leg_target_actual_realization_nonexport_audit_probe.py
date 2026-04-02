#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any


AS_OF = "2026-03-23"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_F947 = GENERATED / "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet.json"
IN_P947 = GENERATED / "p947_current_strict_t173_t176_t220_t222_chart_label_retaining_pair12_typed_seed_subinterface_source_side_input_leg_sufficiency_or_nonexport_audit_probe.json"
IN_P708 = GENERATED / "p708_current_strict_t173_frontier_dashboard_probe_summary.json"
IN_F949 = GENERATED / "f949_first_current_strict_qw2191_pair12_entry_point_same_lane_exhaustion_and_noncyclic_pivot_packet_summary.json"

OUT_JSON = GENERATED / "p1049_current_strict_t173_t176_source_side_input_leg_target_actual_realization_nonexport_audit_probe.json"
OUT_SUMMARY = GENERATED / "p1049_current_strict_t173_t176_source_side_input_leg_target_actual_realization_nonexport_audit_probe_summary.json"

F947_STATUS = (
    "F947_EXECUTED_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_SOURCE_SIDE_INPUT_LEG_TARGET_PACKET_NO_FALSE_PASS"
)
P947_STATUS = (
    "PASS_T173_T176_T220_T222_CHART_LABEL_RETAINING_PAIR12_TYPED_SEED_SUBINTERFACE_SOURCE_SIDE_INPUT_LEG_SUFFICIENCY_OR_NONEXPORT_AUDITED"
)
P708_STATUS = "PASS_T173_FRONTIER_DASHBOARD_READY"
F949_STATUS = (
    "PASS_CURRENT_STRICT_QW2191_PAIR12_ENTRY_POINT_SAME_LANE_EXHAUSTION_AND_NONCYCLIC_PIVOT_PACKET_EXPORTED"
)

BRIDGE_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_global_C_v1_strict_v1"
)
BRIDGE_OUTPUT_SCHEMA_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_output_schema_target_v1"
)
SOURCE_SIDE_INPUT_LEG_TARGET_ID = (
    "inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_v1"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def scan_actual_source_side_input_leg_suppliers() -> list[str]:
    patterns = ("F*.md", "N*.md", "T*.md", "P*.md", "f*.py", "n*.py", "t*.py", "p*.py")
    excluded = {
        Path(__file__).name,
        "P1049_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT_PROBE.md",
        "T302_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_ATTEMPT_SPEC.md",
        "p1050_current_strict_t173_t176_source_side_input_leg_target_actual_realization_attempt_probe.py",
        "F947_CURRENT_STRICT_T173_T176_INVERSION_SENSITIVE_PAIR12_BRANCH_SEPARATION_TO_CHART_SENSITIVE_TRANSPORTED_FLUX_CURRENT_LIKE_SECTION_BRIDGE_SOURCE_SIDE_INPUT_LEG_TARGET_PACKET.md",
        "P947_CURRENT_STRICT_T173_T176_T220_T222_CHART_LABEL_RETAINING_PAIR12_TYPED_SEED_SUBINTERFACE_SOURCE_SIDE_INPUT_LEG_SUFFICIENCY_OR_NONEXPORT_AUDIT_PROBE.md",
        "f947_current_strict_t173_t176_inversion_sensitive_pair12_branch_separation_to_chart_sensitive_transported_flux_current_like_section_bridge_source_side_input_leg_target_packet.py",
        "p947_current_strict_t173_t176_t220_t222_chart_label_retaining_pair12_typed_seed_subinterface_source_side_input_leg_sufficiency_or_nonexport_audit_probe.py",
    }
    hits: list[str] = []
    seen: set[Path] = set()
    for pattern in patterns:
        for path in sorted(ROOT.glob(pattern)):
            if path in seen or path.name in excluded:
                continue
            seen.add(path)
            text = path.read_text(encoding="utf-8")
            if (
                "source_side_input_leg" in text
                and BRIDGE_TARGET_ID in text
                and BRIDGE_OUTPUT_SCHEMA_TARGET_ID in text
                and SOURCE_SIDE_INPUT_LEG_TARGET_ID not in text
            ):
                hits.append(rel(path))
    return hits


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    prerequisites = [IN_F947, IN_P947, IN_P708, IN_F949]
    missing = [rel(path) for path in prerequisites if not path.exists()]
    if missing:
        artifact = {
            "stage": "P1049",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": missing,
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    f947 = load_json(IN_F947)
    p947 = load_json(IN_P947)
    p708 = load_json(IN_P708)
    f949 = load_json(IN_F949)
    named_actual_supplier_hits = scan_actual_source_side_input_leg_suppliers()

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

    f947_exact_target_is_frozen = (
        f947.get("status") == F947_STATUS
        and ((f947.get("target_object") or {}).get("object_id") == SOURCE_SIDE_INPUT_LEG_TARGET_ID)
        and ((f947.get("target_refs") or {}).get("bridge_global_target_ref") == BRIDGE_TARGET_ID)
        and ((f947.get("target_refs") or {}).get("bridge_output_schema_target_ref") == BRIDGE_OUTPUT_SCHEMA_TARGET_ID)
    )

    p947_no_actual_source_side_input_leg_supplier_is_exported = (
        p947.get("status") == P947_STATUS
        and ((p947.get("audit_conclusion") or {}).get("current_repo_already_exports_actual_source_side_input_leg_supplier") is False)
        and ((p947.get("audit_conclusion") or {}).get("route_local_seed_subinterface_even_if_actual_would_still_require_additional_full_c_v1_transported_section_lift") is True)
    )

    p708_t176_global_provider_still_unexported = (
        p708.get("status") == P708_STATUS
        and p708.get("t176_global_provider_exported") is False
    )

    f949_pair12_same_lane_is_exhausted_as_primary_route = (
        f949.get("status") == F949_STATUS
        and f949.get("same_lane_exhaustion_boundary_reached") is True
    )

    no_current_actual_source_side_input_leg_supplier_export_found = len(named_actual_supplier_hits) == 0

    target_actual_realization_still_nonexport = (
        f947_exact_target_is_frozen
        and p947_no_actual_source_side_input_leg_supplier_is_exported
        and p708_t176_global_provider_still_unexported
        and f949_pair12_same_lane_is_exhausted_as_primary_route
        and no_current_actual_source_side_input_leg_supplier_export_found
    )

    add_check(
        "f947_exact_target_is_frozen",
        f947_exact_target_is_frozen,
        True,
        "F947 already freezes the exact source-side input-leg target for the active T173/T176 bridge family.",
    )
    add_check(
        "p947_no_actual_source_side_input_leg_supplier_is_exported",
        p947_no_actual_source_side_input_leg_supplier_is_exported,
        True,
        "P947 already freezes that no current actual supplier of that leg is exported.",
    )
    add_check(
        "p708_t176_global_provider_still_unexported",
        p708_t176_global_provider_still_unexported,
        True,
        "The global T176 provider is still unexported on the current T173 frontier.",
    )
    add_check(
        "f949_pair12_same_lane_is_exhausted_as_primary_route",
        f949_pair12_same_lane_is_exhausted_as_primary_route,
        True,
        "The old pair12 same-lane descent is already exhausted as a primary route.",
    )
    add_check(
        "no_current_actual_source_side_input_leg_supplier_export_found",
        no_current_actual_source_side_input_leg_supplier_export_found,
        True,
        "No additional current export lawfully upgrades the frozen source-side input-leg target into an actual supplier.",
    )
    add_check(
        "target_actual_realization_still_nonexport",
        target_actual_realization_still_nonexport,
        True,
        "Therefore the exact source-side input-leg target still remains future-only on the current repo state.",
    )

    status = (
        "PASS_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDITED"
        if not blocking and target_actual_realization_still_nonexport
        else "FAIL_CURRENT_STRICT_T173_T176_SOURCE_SIDE_INPUT_LEG_TARGET_ACTUAL_REALIZATION_NONEXPORT_AUDIT"
    )

    artifact = {
        "stage": "P1049",
        "status": status,
        "as_of": AS_OF,
        "target_object_id": SOURCE_SIDE_INPUT_LEG_TARGET_ID,
        "target_actual_realization_exported_on_current_repo_state": not target_actual_realization_still_nonexport,
        "current_repo_already_exports_actual_source_side_input_leg_supplier": not no_current_actual_source_side_input_leg_supplier_export_found,
        "pair12_same_lane_primary_route_exhausted": f949_pair12_same_lane_is_exhausted_as_primary_route,
        "named_actual_source_side_input_leg_supplier_hits": named_actual_supplier_hits,
        "next_honest_move_is_freeze_one_exact_noncyclic_actual_realization_attempt": target_actual_realization_still_nonexport,
        "checks": checks,
        "blocking_checks": blocking,
        "no_false_pass": True,
    }
    write_json(OUT_JSON, artifact)

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_object_id": artifact["target_object_id"],
        "target_actual_realization_exported_on_current_repo_state": artifact[
            "target_actual_realization_exported_on_current_repo_state"
        ],
        "current_repo_already_exports_actual_source_side_input_leg_supplier": artifact[
            "current_repo_already_exports_actual_source_side_input_leg_supplier"
        ],
        "pair12_same_lane_primary_route_exhausted": artifact["pair12_same_lane_primary_route_exhausted"],
        "next_honest_move_is_freeze_one_exact_noncyclic_actual_realization_attempt": artifact[
            "next_honest_move_is_freeze_one_exact_noncyclic_actual_realization_attempt"
        ],
        "no_false_pass": True,
    }
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
