#!/usr/bin/env python3
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

IN_P478_SUMMARY = (
    GENERATED
    / "p478_current_strict_t169_rordpow_sign_scan_for_r18_pair1_zero_equations_probe_summary.json"
)
OUT = (
    GENERATED
    / "n521_current_first_strict_t169_rordpow_sign_scan_r18_pair1_zero_equations_obstruction_theorem_summary.json"
)

EXPECTED_P478_STATUS = "PASS_SCAN_COMPLETE_NO_SIGN_VECTOR_SATISFIES_R18_PAIR1_ZERO_EQUATIONS_UNDER_N477_ON_RORDPOW_MAGNITUDES"


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    if not IN_P478_SUMMARY.is_file():
        summary = {
            "stage": "N521",
            "status": "FAIL_MISSING_REQUIRED_INPUTS",
            "missing_required_inputs": ["P478_summary"],
            "required_paths": {"P478_summary": str(IN_P478_SUMMARY.relative_to(REPO))},
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "theorem_result": {"discharged": False, "reason": "missing P478 summary"},
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    p478 = load_json(IN_P478_SUMMARY)
    status_ok = p478.get("status") == EXPECTED_P478_STATUS

    if not status_ok:
        summary = {
            "stage": "N521",
            "status": "N521_REQUIRES_REVIEW_OR_FRONTIER_CHANGED",
            "goal": "Package a finite-search obstruction: under the fixed T169 r_ordpow magnitude lift and conditional N477, no sign vector satisfies all three R18 declared pair1 residual zero equations.",
            "scope": "conditional_N477_under_fixed_r_ordpow_magnitudes_and_uniform_g4_only",
            "checks": [
                {
                    "id": "p478_status_matches_expected",
                    "actual": p478.get("status"),
                    "expected": EXPECTED_P478_STATUS,
                    "pass": False,
                    "meaning": "P478 no-solution status is required to claim the finite sign-scan obstruction",
                }
            ],
            "theorem_result": {
                "discharged": False,
                "reason": "P478 status does not match expected no-solution status; the obstruction claim may no longer apply.",
            },
            "depends_on": [
                "P478_CURRENT_STRICT_T169_RORDPOW_SIGN_SCAN_FOR_R18_PAIR1_ZERO_EQUATIONS_PROBE",
            ],
            "required_next_step": "REVIEW_P478_AND_UPDATE_OR_RETIRE_N521_IF_THE_SIGN_SCAN_RESULT_CHANGED",
            "generated_utc": datetime.now(timezone.utc).isoformat(),
            "no_false_pass": True,
        }
        OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
        print(OUT)
        return

    summary = {
        "stage": "N521",
        "status": "N521_DISCHARGED_CURRENT_FIRST_STRICT_T169_RORDPOW_SIGN_SCAN_R18_PAIR1_ZERO_EQUATIONS_OBSTRUCTION_THEOREM_NO_FALSE_PASS",
        "goal": "Package a finite-search obstruction: under the fixed T169 r_ordpow magnitude lift and conditional N477, no sign vector satisfies all three R18 declared pair1 residual zero equations.",
        "scope": "conditional_N477_under_fixed_r_ordpow_magnitudes_and_uniform_g4_only",
        "evidence": {
            "p478_summary": str(IN_P478_SUMMARY.relative_to(REPO)),
            "any_solution_within_tolerance": p478.get("any_solution_within_tolerance"),
            "min_abs_by_entry_over_scan": p478.get("min_abs_by_entry_over_scan"),
            "best_candidate_objective_value": p478.get("best_candidate_objective_value"),
        },
        "theorem_result": {
            "discharged": True,
            "obstruction_kind": "finite_exhaustive_sign_scan_no_solution",
            "consequence_for_P16": "A P16 zero/cancellation witness cannot be obtained by altering only the sign selection within the fixed r_ordpow magnitude lift class under conditional N477; a different strict-derived provider class and/or additional non-N477 ingredients are required.",
        },
        "depends_on": [
            "N483_CURRENT_FIRST_STRICT_T169_POWERLAW_ELEMENT_ORDER_CONSTRAINED_LIFT_EXISTENCE_UNIQUENESS_THEOREM",
            "F447_CURRENT_STRICT_T169_QW2122_SCALAR_TO_T168_PER_SITE_VALUE_PROVIDER_STRICT_DERIVED_V1",
            "R18_PAIR1_RESIDUAL_DECLARED_PULLBACK_COEFFICIENT_CLASS_REDUCTION_PACKET",
            "P478_CURRENT_STRICT_T169_RORDPOW_SIGN_SCAN_FOR_R18_PAIR1_ZERO_EQUATIONS_PROBE",
        ],
        "hard_limits": [
            "Conditional on the N477 Yukawa-free diagonal residual rewrite; no vacuum stationarity witness is exported.",
            "Restricted to the fixed r_ordpow magnitude lift and the fixed uniform g4 lift (T169 scope).",
            "Does not claim the declared pair1 residual zero equations cannot hold under other strict-derived provider classes or after adding non-N477 ingredients.",
            "Does not discharge QW-2191 and does not imply ToE closure.",
        ],
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "no_false_pass": True,
    }

    OUT.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT)


if __name__ == "__main__":
    main()

