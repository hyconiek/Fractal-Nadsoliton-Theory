#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_F695 = GENERATED / "f695_first_actual_strict_projective_operational_toe_preclosure_support_packet_summary.json"
OUT = GENERATED / "p695_current_actual_strict_projective_operational_toe_preclosure_support_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F695.exists():
        summary = {
            "stage": "P695",
            "lane": "strict_projective_operational_toe_preclosure_only",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_F695.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f695 = load_json(IN_F695)

    checks_spec = [
        {
            "id": "support_packet_name_matches",
            "actual": f695.get("support_packet_name"),
            "expected": "Lambda_strict_projective_operational_toe_preclosure_support_v1",
        },
        {
            "id": "projective_strict_core_selector_closure_discharged",
            "actual": f695.get("projective_strict_core_selector_closure_discharged"),
            "expected": True,
        },
        {
            "id": "projective_output_ray_invariant_across_directed_lifts",
            "actual": f695.get("projective_output_ray_invariant_across_directed_lifts"),
            "expected": True,
        },
        {
            "id": "output_sign_lift_gauge_covariant",
            "actual": f695.get("output_sign_lift_gauge_covariant"),
            "expected": True,
        },
        {
            "id": "physical_mass_spectrum_proxy_computable_from_projective_closure",
            "actual": f695.get("physical_mass_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "strict_projective_operational_toe_preclosure_support_packet_exported",
            "actual": f695.get("strict_projective_operational_toe_preclosure_support_packet_exported"),
            "expected": True,
        },
        {
            "id": "kernel_alone_qw2191_discharge_false",
            "actual": f695.get("kernel_alone_qw2191_discharge"),
            "expected": False,
        },
        {
            "id": "actual_strict_core_toe_closure_not_discharged",
            "actual": f695.get("actual_strict_core_toe_closure_discharged"),
            "expected": False,
        },
        {
            "id": "actual_global_toe_closure_not_discharged",
            "actual": f695.get("actual_global_toe_closure_discharged"),
            "expected": False,
        },
    ]

    checks = []
    mismatches = []
    for item in checks_spec:
        ok = item["actual"] == item["expected"]
        checks.append(
            {
                "id": item["id"],
                "actual": item["actual"],
                "expected": item["expected"],
                "pass": ok,
            }
        )
        if not ok:
            mismatches.append(item["id"])

    if mismatches:
        status = "P695_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_PROJECTIVE_OPERATIONAL_TOE_PRECLOSURE_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_PRECLOSURE_SUPPORT_PACKET_AFTER_P695"

    summary = {
        "stage": "P695",
        "lane": "strict_projective_operational_toe_preclosure_only",
        "status": status,
        "as_of": AS_OF,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "strict_projective_operational_toe_preclosure_support_packet_exported": f695.get(
            "strict_projective_operational_toe_preclosure_support_packet_exported"
        ),
        "actual_strict_core_toe_closure_discharged": f695.get("actual_strict_core_toe_closure_discharged"),
        "actual_global_toe_closure_discharged": f695.get("actual_global_toe_closure_discharged"),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

