#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
GENERATED = ROOT / "generated"
OUT_JSON = (
    GENERATED
    / "p646_current_strict_witness_provider_scan_probe_for_seed_v1_realization_attempt.json"
)
OUT_SUMMARY = (
    GENERATED
    / "p646_current_strict_witness_provider_scan_probe_for_seed_v1_realization_attempt_summary.json"
)


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def matches_contract_signature(obj: dict[str, Any]) -> bool:
    if not isinstance(obj, dict):
        return False
    if obj.get("base_realization_attempt") != "S_sel_int_new_object_constructed_realization_attempt_v1":
        return False
    if obj.get("constructed_source_object_exported") is not True:
        return False
    if obj.get("no_false_pass") is not True:
        return False
    if not isinstance(obj.get("exported_constructed_source_object"), str):
        return False
    return True


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    f646 = load_json(
        GENERATED
        / "f646_strict_witness_provider_contract_packet_for_seed_v1_realization_attempt_summary.json"
    )

    candidates: list[dict[str, Any]] = []
    parse_errors: list[str] = []
    for path in sorted(GENERATED.glob("*.json")):
        try:
            obj = load_json(path)
        except Exception:
            parse_errors.append(path.name)
            continue
        if matches_contract_signature(obj):
            candidates.append(
                {
                    "path": str(path),
                    "exported_constructed_source_object": obj.get("exported_constructed_source_object"),
                }
            )

    checks_spec = [
        {
            "id": "f646_contract_present",
            "actual": f646["contract"]["id"],
            "expected": "strict_witness_provider_contract_for_seed_v1_realization_attempt_v1",
            "meaning": "F646 freezes the minimal witness-provider contract and scan signature",
        },
        {
            "id": "candidate_count",
            "actual": len(candidates),
            "expected": 0,
            "meaning": "conservative: no strict witness provider is exported yet on current repo state",
        },
    ]

    checks: list[dict[str, Any]] = []
    for item in checks_spec:
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": item["actual"] == item["expected"],
                "meaning": item["meaning"],
            }
        )

    artifact = {
        "stage": "P646",
        "lane": "current_strict_witness_provider_scan_only",
        "goal": "scan_generated_exports_for_a_strict_constructed_source_object_export_witness_provider_matching_F646_signature",
        "status": "CURRENT_REPO_EXPORTS_NO_STRICT_WITNESS_PROVIDER_MATCHING_F646_SIGNATURE_FOR_SEED_V1_AFTER_P646",
        "scan_result": {
            "candidates_count": len(candidates),
            "candidates": candidates,
            "parse_errors_count": len(parse_errors),
        },
        "checks": checks,
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    summary = {
        "stage": "P646",
        "status": artifact["status"],
        "lane": artifact["lane"],
        "scan_result": artifact["scan_result"],
        "strict_core_promotion": False,
        "full_closure_pass": False,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_JSON)


if __name__ == "__main__":
    main()

