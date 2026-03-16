#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"
OUT = (
    GENERATED
    / "n520_current_first_strict_r18_declared_pair1_residual_zero_equations_value_instance_obstruction_theorem_summary.json"
)

EXPECTED_R18_STATUS = "PASS_PARTIAL_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET_READY"
EXPECTED_P459_STATUS = "PASS_COMPUTED_FROM_STRICT_DERIVED_INPUTS"
EXPECTED_P477_STATUSES = {
    "PASS_EVALUATION_ZERO_EQUATIONS_SATISFIED_UNDER_CURRENT_VALUE_INSTANCE",
    "PASS_EVALUATION_ZERO_EQUATIONS_VIOLATED_UNDER_CURRENT_VALUE_INSTANCE",
}


def _p(path: Path) -> str:
    try:
        return str(path.relative_to(REPO))
    except ValueError:
        return str(path)


def load_json(repo_relative_path: str) -> dict[str, Any]:
    return json.loads((REPO / repo_relative_path).read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    required = {
        "R18_summary": GENERATED
        / "r18_pair1_residual_declared_pullback_coefficient_class_reduction_packet_for_host_matching_route_summary.json",
        "P459_summary": GENERATED / "p459_current_strict_r16_declared_control_pullback_value_instantiation_probe_summary.json",
        "P477_summary": GENERATED
        / "p477_current_strict_r18_pair1_residual_zero_equations_value_instantiation_probe_summary.json",
    }
    missing = [k for k, p in required.items() if not p.is_file()]
    if missing:
        summary = {
            "stage": "N520",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "goal": "Package a value-instance-only obstruction: the current P459 strict-derived declared residual pullback violates the R18 declared pair1 residual zero equations (as evaluated by P477).",
            "missing_required_inputs": missing,
            "required_paths": {k: _p(v) for k, v in required.items()},
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "theorem_result": {"discharged": False, "reason": "missing required inputs"},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    sources = {
        "R18": json.loads(required["R18_summary"].read_text(encoding="utf-8")),
        "P459": json.loads(required["P459_summary"].read_text(encoding="utf-8")),
        "P477": json.loads(required["P477_summary"].read_text(encoding="utf-8")),
    }

    checks_spec = [
        {
            "id": "r18_pair1_reduction_packet_present",
            "actual": sources["R18"].get("status"),
            "expected": EXPECTED_R18_STATUS,
            "meaning": "R18 summary indicates the declared pair1 coefficient-class reduction packet is ready",
        },
        {
            "id": "p459_value_instance_present",
            "actual": sources["P459"].get("status"),
            "expected": EXPECTED_P459_STATUS,
            "meaning": "P459 exports a strict-derived value-instantiated declared control pullback object",
        },
        {
            "id": "p477_value_instance_evaluation_present",
            "actual": sources["P477"].get("status"),
            "expected": "one_of:" + ",".join(sorted(EXPECTED_P477_STATUSES)),
            "meaning": "P477 evaluates the R18 declared pair1 residual zero equations on the P459 value instance",
        },
        {
            "id": "p477_no_false_pass_flag_present",
            "actual": bool(sources["P477"].get("no_false_pass")),
            "expected": True,
            "meaning": "P477 summary explicitly asserts no_false_pass",
        },
    ]

    checks: list[dict[str, Any]] = []
    blocking_mismatches: list[str] = []
    for item in checks_spec:
        if item["id"] == "p477_value_instance_evaluation_present":
            ok = item["actual"] in EXPECTED_P477_STATUSES
        else:
            ok = item["actual"] == item["expected"]
        checks.append({**item, "pass": ok})
        if not ok:
            blocking_mismatches.append(item["id"])

    p477_status = str(sources["P477"].get("status", ""))
    violated = sources["P477"].get("violated_equations")
    if not isinstance(violated, list) or not all(isinstance(x, str) for x in violated):
        violated = []

    obstruction_present = p477_status == "PASS_EVALUATION_ZERO_EQUATIONS_VIOLATED_UNDER_CURRENT_VALUE_INSTANCE"
    obstruction_claim_ok = obstruction_present and bool(violated)

    if not obstruction_claim_ok:
        status = "N520_REQUIRES_REVIEW_OR_OBSTRUCTION_REMOVED"
        theorem_result = {
            "discharged": False,
            "reason": "P477 does not report a violation under the current value instance (or violated_equations is empty); N520 obstruction claim may no longer apply.",
        }
        required_next_step = "REVIEW_P477_AND_UPDATE_OR_RETIRE_N520_IF_THE_VALUE_INSTANCE_NO_LONGER_VIOLATES_R18_DECLARED_PAIR1_ZERO_EQUATIONS"
        no_false_pass = False if blocking_mismatches else True
    else:
        status = (
            "N520_DISCHARGED_CURRENT_FIRST_STRICT_R18_DECLARED_PAIR1_RESIDUAL_ZERO_EQUATIONS_VALUE_INSTANCE_OBSTRUCTION_THEOREM_NO_FALSE_PASS"
        )
        theorem_result = {
            "discharged": True,
            "obstruction_scope": "value_instance_only",
            "p459_value_instance_object": sources["P459"].get("constructed_object"),
            "r18_declared_pair1_zero_equations_satisfied": False,
            "violated_equations": violated,
        }
        required_next_step = "EXPORT_A_STRICT_DERIVED_PROVIDER_CLASS_OR_STATIONARITY_WITNESS_THAT_MAKES_THE_R18_DECLARED_PAIR1_ZERO_EQUATIONS_HOLD_OR_KEEP_THE_HOST_MATCHING_ROUTE_NEGATIVE"
        no_false_pass = True if not blocking_mismatches else False

    summary = {
        "stage": "N520",
        "status": status,
        "goal": "Package a value-instance-only obstruction: the current P459 strict-derived declared residual pullback violates the R18 declared pair1 residual zero equations (as evaluated by P477).",
        "scope": "current_strict_derived_value_instance_only",
        "checks": checks,
        "blocking_mismatches": blocking_mismatches,
        "theorem_result": theorem_result,
        "depends_on": [
            "R18_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET",
            "P459_CURRENT_STRICT_R16_DECLARED_CONTROL_PULLBACK_VALUE_INSTANTIATION_PROBE",
            "P477_CURRENT_STRICT_R18_PAIR1_RESIDUAL_ZERO_EQUATIONS_VALUE_INSTANTIATION_PROBE",
        ],
        "required_next_step": required_next_step,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "no_false_pass": bool(no_false_pass),
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

