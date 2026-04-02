#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

AS_OF = "2026-03-29"
EXPECTED_TARGET = "MinimalONRDBoundaryToActiveBridgeExactReductionTarget_v1"
EXPECTED_BOUNDARY = "MinimalOrientedNonreciprocalDephasingNewImportBoundary_v1"
EXPECTED_BRIDGE = "ResidualDatumPair12OrbitDirectionSelectionBridge_global_C_v1_strict_v1"

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P1095 = GENERATED / "p1095_current_strict_t173_t176_post_f969_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_admission_probe_summary.json"
OUT_JSON = GENERATED / "p1096_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_probe.json"
OUT_SUMMARY = GENERATED / "p1096_current_strict_t173_t176_minimal_onrd_boundary_to_active_bridge_exact_reduction_target_probe_summary.json"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def write_json(path: Path, obj: dict[str, Any]) -> None:
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=True) + "\n", encoding="ascii")


def rel(path: Path) -> str:
    return str(path.relative_to(REPO))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)
    if not IN_P1095.exists():
        artifact = {
            "stage": "P1096",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [rel(IN_P1095)],
            "no_false_pass": True,
        }
        write_json(OUT_JSON, artifact)
        write_json(OUT_SUMMARY, artifact)
        print(OUT_SUMMARY)
        return

    p1095 = load_json(IN_P1095)
    admitted = p1095.get("status") == "PASS_CURRENT_STRICT_T173_T176_POST_F969_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_ADMITTED"
    status = (
        "PASS_CURRENT_STRICT_T173_T176_MINIMAL_ONRD_BOUNDARY_TO_ACTIVE_BRIDGE_EXACT_REDUCTION_TARGET_EXPORTED"
        if admitted
        else "P1096_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_TARGET_STATE"
    )

    artifact = {
        "stage": "P1096",
        "status": status,
        "as_of": AS_OF,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "depends_on": {
            "p1095_target_admission_summary": rel(IN_P1095),
        },
        "exported_target": {
            "object_id": EXPECTED_TARGET,
            "source_boundary_object_id": EXPECTED_BOUNDARY,
            "target_bridge_object_id": EXPECTED_BRIDGE,
            "search_mode_selected": "genuinely_new_inversion_sensitive_source_side_provider_class_route",
            "within_exported_noncyclic_provider_split_family": False,
            "typed_meaning": "one exact reduction target deriving the active bridge from the minimal ONRD boundary without claiming that the derivation already exists",
        },
        "admission_properties": {
            "admissible_as_exact_reduction_target": True,
            "counts_as_lawful_supplier": False,
            "counts_as_solution": False,
            "counts_as_strict_physical_orientation_datum": False,
        },
        "hard_limits": [
            "Does not export the reduction.",
            "Does not export a lawful supplier.",
            "Does not export a solution.",
            "Does not export a strict physical orientation datum.",
            "Does not discharge T183.",
            "Does not discharge T176.",
            "Does not discharge QW-2191.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": artifact["stage"],
        "status": artifact["status"],
        "as_of": artifact["as_of"],
        "target_object_id": artifact["exported_target"]["object_id"],
        "source_boundary_object_id": artifact["exported_target"]["source_boundary_object_id"],
        "target_bridge_object_id": artifact["exported_target"]["target_bridge_object_id"],
        "search_mode_selected": artifact["exported_target"]["search_mode_selected"],
        "within_exported_noncyclic_provider_split_family": artifact["exported_target"]["within_exported_noncyclic_provider_split_family"],
        "admissible_as_exact_reduction_target": True,
        "counts_as_lawful_supplier": False,
        "counts_as_solution": False,
        "counts_as_strict_physical_orientation_datum": False,
        "no_false_pass": True,
    }

    write_json(OUT_JSON, artifact)
    write_json(OUT_SUMMARY, summary)
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
