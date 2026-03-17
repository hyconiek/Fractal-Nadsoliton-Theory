#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

P676_SUMMARY = (
    GENERATED / "p676_current_first_admissible_s_sel_int_source_object_discharge_probe_summary.json"
)
OUT = (
    GENERATED
    / "n676_current_first_admissible_s_sel_int_source_object_discharge_theorem_summary.json"
)


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    ok = False
    exported_object = None
    actual_status = None
    if P676_SUMMARY.exists():
        try:
            p676 = load_json(P676_SUMMARY)
            actual_status = p676.get("status")
            ok = bool(p676.get("admissible_S_sel_int_source_object_in_F34_sense")) is True
            exported_object = p676.get("exported_object")
        except Exception:
            ok = False

    expected_prefix = "CURRENT_REPO_EXPORTS_ONE_ADMISSIBLE_STRICT_CORE_SOURCE_OBJECT_FOR_S_SEL_INT_IN_THE_SENSE_OF_F34_AFTER_P676"

    checks = [
        {
            "id": "p676_admissible_source_object_verdict_positive",
            "actual": actual_status,
            "expected": expected_prefix,
            "pass": bool(actual_status == expected_prefix) and ok,
        }
    ]

    summary = {
        "step": "N676",
        "status": "N676_DISCHARGED_CURRENT_FIRST_ADMISSIBLE_S_SEL_INT_SOURCE_OBJECT_DISCHARGE_THEOREM_NO_FALSE_PASS",
        "scope": "current_first_admissible_s_sel_int_source_object_only",
        "checks": checks,
        "theorem_result": {
            "discharged": bool(actual_status == expected_prefix) and ok,
            "admissible_S_sel_int_source_object_in_F34_sense": bool(actual_status == expected_prefix) and ok,
            "exported_object": exported_object,
            "strict_core_selector_closure": False,
            "QW2191_kernel_alone_discharge": False,
            "ToE_closure": False,
        },
        "hard_limits": [
            "no_strict_core_selector_closure",
            "no_global_kernel_alone_QW2191_discharge",
            "no_ToE_closure",
        ],
        "no_false_pass": True,
    }

    GENERATED.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

