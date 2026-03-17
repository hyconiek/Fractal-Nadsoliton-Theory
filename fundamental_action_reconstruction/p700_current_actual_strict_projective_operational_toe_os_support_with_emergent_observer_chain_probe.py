#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_F700 = GENERATED / "f700_first_actual_strict_projective_operational_toe_os_support_with_emergent_observer_chain_packet_summary.json"
OUT = GENERATED / "p700_current_actual_strict_projective_operational_toe_os_support_with_emergent_observer_chain_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F700.exists():
        summary = {
            "stage": "P700",
            "lane": "strict_projective_operational_toe_os_support_v2_only",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_F700.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f700 = load_json(IN_F700)

    checks_spec = [
        {
            "id": "support_packet_name_matches",
            "actual": f700.get("support_packet_name"),
            "expected": "Lambda_strict_projective_operational_toe_os_support_v2",
        },
        {
            "id": "projective_strict_core_selector_closure_discharged",
            "actual": f700.get("projective_strict_core_selector_closure_discharged"),
            "expected": True,
        },
        {
            "id": "projective_output_ray_invariant_across_directed_lifts",
            "actual": f700.get("projective_output_ray_invariant_across_directed_lifts"),
            "expected": True,
        },
        {
            "id": "output_sign_lift_gauge_covariant",
            "actual": f700.get("output_sign_lift_gauge_covariant"),
            "expected": True,
        },
        {
            "id": "physical_mass_spectrum_proxy_computable_from_projective_closure",
            "actual": f700.get("physical_mass_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "selector_aligned_channel_spectrum_proxy_computable_from_projective_closure",
            "actual": f700.get("selector_aligned_channel_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "observer_limit_readout_computable_from_projective_closure_output_projector",
            "actual": f700.get("observer_limit_readout_computable_from_projective_closure_output_projector"),
            "expected": True,
        },
        {
            "id": "projective_emergent_observer_chain_computable_from_projective_closure_output_projector",
            "actual": f700.get("projective_emergent_observer_chain_computable_from_projective_closure_output_projector"),
            "expected": True,
        },
        {
            "id": "strict_projective_operational_toe_os_support_packet_v2_exported",
            "actual": f700.get("strict_projective_operational_toe_os_support_packet_v2_exported"),
            "expected": True,
        },
        {
            "id": "kernel_alone_qw2191_discharge_false",
            "actual": f700.get("kernel_alone_qw2191_discharge"),
            "expected": False,
        },
        {
            "id": "actual_strict_core_toe_closure_not_discharged",
            "actual": f700.get("actual_strict_core_toe_closure_discharged"),
            "expected": False,
        },
        {
            "id": "actual_global_toe_closure_not_discharged",
            "actual": f700.get("actual_global_toe_closure_discharged"),
            "expected": False,
        },
        {
            "id": "actual_emergent_observer_closure_not_discharged",
            "actual": f700.get("actual_emergent_observer_closure"),
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
        status = "P700_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V2_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V2_PACKET_AFTER_P700"

    summary = {
        "stage": "P700",
        "lane": "strict_projective_operational_toe_os_support_v2_only",
        "status": status,
        "as_of": AS_OF,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "strict_projective_operational_toe_os_support_packet_v2_exported": f700.get(
            "strict_projective_operational_toe_os_support_packet_v2_exported"
        ),
        "actual_strict_core_toe_closure_discharged": f700.get("actual_strict_core_toe_closure_discharged"),
        "actual_global_toe_closure_discharged": f700.get("actual_global_toe_closure_discharged"),
        "actual_emergent_observer_closure": f700.get("actual_emergent_observer_closure"),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

