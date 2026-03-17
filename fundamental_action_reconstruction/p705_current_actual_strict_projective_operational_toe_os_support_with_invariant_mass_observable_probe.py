#!/usr/bin/env python3

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

AS_OF = "2026-03-17"

IN_F705 = GENERATED / "f705_first_actual_strict_projective_operational_toe_os_support_with_invariant_mass_observable_packet_summary.json"
OUT = GENERATED / "p705_current_actual_strict_projective_operational_toe_os_support_with_invariant_mass_observable_probe_summary.json"


def load_json(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_F705.exists():
        summary = {
            "stage": "P705",
            "lane": "strict_projective_operational_toe_os_support_v3_only",
            "status": "NOT_COMPUTABLE_MISSING_PREREQUISITES",
            "as_of": AS_OF,
            "missing": [str(IN_F705.relative_to(REPO))],
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    f705 = load_json(IN_F705)

    checks_spec = [
        {
            "id": "support_packet_name_matches",
            "actual": f705.get("support_packet_name"),
            "expected": "Lambda_strict_projective_operational_toe_os_support_v3",
        },
        {
            "id": "projective_strict_core_selector_closure_discharged",
            "actual": f705.get("projective_strict_core_selector_closure_discharged"),
            "expected": True,
        },
        {
            "id": "projective_output_ray_invariant_across_directed_lifts",
            "actual": f705.get("projective_output_ray_invariant_across_directed_lifts"),
            "expected": True,
        },
        {
            "id": "output_sign_lift_gauge_covariant",
            "actual": f705.get("output_sign_lift_gauge_covariant"),
            "expected": True,
        },
        {
            "id": "physical_mass_spectrum_proxy_computable_from_projective_closure",
            "actual": f705.get("physical_mass_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "selector_aligned_channel_spectrum_proxy_computable_from_projective_closure",
            "actual": f705.get("selector_aligned_channel_spectrum_proxy_computable_from_projective_closure"),
            "expected": True,
        },
        {
            "id": "basis_invariant_mass_observable_exported",
            "actual": f705.get("basis_invariant_mass_observable_exported"),
            "expected": True,
        },
        {
            "id": "observer_limit_readout_computable_from_projective_closure_output_projector",
            "actual": f705.get("observer_limit_readout_computable_from_projective_closure_output_projector"),
            "expected": True,
        },
        {
            "id": "projective_emergent_observer_chain_computable_from_projective_closure_output_projector",
            "actual": f705.get("projective_emergent_observer_chain_computable_from_projective_closure_output_projector"),
            "expected": True,
        },
        {
            "id": "strict_projective_operational_toe_os_support_packet_v3_exported",
            "actual": f705.get("strict_projective_operational_toe_os_support_packet_v3_exported"),
            "expected": True,
        },
        {
            "id": "kernel_alone_qw2191_discharge_false",
            "actual": f705.get("kernel_alone_qw2191_discharge"),
            "expected": False,
        },
        {
            "id": "actual_strict_core_toe_closure_not_discharged",
            "actual": f705.get("actual_strict_core_toe_closure_discharged"),
            "expected": False,
        },
        {
            "id": "actual_global_toe_closure_not_discharged",
            "actual": f705.get("actual_global_toe_closure_discharged"),
            "expected": False,
        },
        {
            "id": "actual_emergent_observer_closure_not_discharged",
            "actual": f705.get("actual_emergent_observer_closure"),
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
        status = "P705_REQUIRES_REVIEW_CHANGED_OR_INSUFFICIENT_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V3_STATE"
    else:
        status = "CURRENT_REPO_EXPORTS_ONE_ACTUAL_STRICT_PROJECTIVE_OPERATIONAL_TOE_OS_SUPPORT_V3_PACKET_AFTER_P705"

    summary = {
        "stage": "P705",
        "lane": "strict_projective_operational_toe_os_support_v3_only",
        "status": status,
        "as_of": AS_OF,
        "checks": checks,
        "blocking_mismatches": mismatches,
        "strict_projective_operational_toe_os_support_packet_v3_exported": f705.get(
            "strict_projective_operational_toe_os_support_packet_v3_exported"
        ),
        "actual_strict_core_toe_closure_discharged": f705.get("actual_strict_core_toe_closure_discharged"),
        "actual_global_toe_closure_discharged": f705.get("actual_global_toe_closure_discharged"),
        "actual_emergent_observer_closure": f705.get("actual_emergent_observer_closure"),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

