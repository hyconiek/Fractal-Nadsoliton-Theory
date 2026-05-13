#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_P1376_SUMMARY = GEN / "p1376_full_ci_export_for_ci_coefficient_class_summary.json"
OUT = GEN / "p1377_ci_class_population_and_first_transport_drift_run_summary.json"


def rel_drift(a: float, b: float) -> float:
    return abs(a - b) / max(1.0, abs(a))


def main() -> None:
    p1376 = json.loads(IN_P1376_SUMMARY.read_text(encoding="utf-8"))

    coeffs = {
        "c_g3": {"A": 1.00, "B": 1.03},
        "c_g2": {"A": 1.00, "B": 0.98},
        "c_g1": {"A": 1.00, "B": 1.01},
    }

    drifts = {k: rel_drift(v["A"], v["B"]) for k, v in coeffs.items()}
    transport_drift = max(drifts.values())

    out = {
        "artifact": "p1377_ci_class_population_and_first_transport_drift_run_summary",
        "as_of": "2026-05-12",
        "input_dependency": IN_P1376_SUMMARY.name,
        "input_l_b1_02_readiness": p1376.get("l_b1_02_readiness"),
        "trial_coefficients": coeffs,
        "missing_coefficients": ["c_mix"],
        "drift_per_coefficient": drifts,
        "transport_drift": transport_drift,
        "epsilon_transport": "UNSET",
        "l_b1_02_status": "INCOMPLETE_CLASS_TRIAL_ONLY",
        "forbidden_claim": "do_not_issue_formal_pass_fail_before_full_class_and_epsilon_freeze",
        "next_packet": "P1378_CMIX_POPULATION_AND_EPSILON_TRANSPORT_FREEZE"
    }

    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[p1377] wrote {OUT}")


if __name__ == "__main__":
    main()
