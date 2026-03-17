#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_F698 = GENERATED / "f698_first_actual_strict_projective_operational_toe_os_support_packet_summary.json"
OUT = GENERATED / "p698_current_actual_strict_projective_operational_toe_os_support_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F698.exists():
        summary = {
            "stage": "P698",
            "lane": "strict_projective_operational_toe_os_support_only",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_F698.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f698 = load_json(IN_F698)

    checks_spec = [
        {
            "id": "support_packet_name_matches",
            "actual": f698.get("support_packet_name"),
            "expected": "Lambda_strict_projective_operational_toe_os_support_v1",
        },
        {
            "id": "projective_strict_core_selector_closure_discharged",
            "actual": f698.get("projective_strict_core_selector_closure_discharged"),
            "expected": True,
        },
        {
            "id": "projective_output_ray_invariant_across_directed_lifts",
            "actual": f698.get("projective_output_ray_invariant_across_directed_lifts"),
            "expected": True,
        },
        {
            "id": "output_sign_lift_gauge_covariant",
            "actual": f698.get("output_sign_lift_gauge_covariant"),
            "expected": True,
        },
        {
            "id": "physical_mass_spectrum_proxy_computable_from_projective_closure",
            "actual": f698.get("physical_mass_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "selector_aligned_channel_spectrum_proxy_computable_from_projective_closure",
            "actual": f698.get("selector_aligned_channel_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "observer_limit_readout_computable_from_projective_closure_output_projector",
            "actual": f698.get("observer_limit_readout_computable_from_projective_closure_output_projector"),
            "expected": True,
        },
        {
            "id": "strict_projective_operational_toe_os_support_packet_exported",
            "actual": f698.get("strict_projective_operational_toe_os_support_packet_exported"),
            "expected": True,
        },
        {
            "id": "kernel_alone_qw2191_discharge_false",
            "actual": f698.get("kernel_alone_qw2191_discharge"),
            "expected": False,
        },
        {
            "id": "actual_strict_core_toe_closure_not_discharged",
            "actual": f698.get("actual_strict_core_toe_closure_discharged"),
            "expected": False,
        },
        {
            "id": "actual_global_toe_closure_not_discharged",
            "actual": f698.get("actual_global_toe_closure_discharged"),
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
        status = "P698_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_PACKET_AFTER_P698"

    summary = {
        "stage": "P698",
        "lane": "strict_projective_operational_toe_os_support_only",
        "status": status,
        "as_of": AS_OF,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "strict_projective_operational_toe_os_support_packet_exported": f698.get(
            "strict_projective_operational_toe_os_support_packet_exported"
        ),
        "actual_strict_core_toe_closure_discharged": f698.get("actual_strict_core_toe_closure_discharged"),
        "actual_global_toe_closure_discharged": f698.get("actual_global_toe_closure_discharged"),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

