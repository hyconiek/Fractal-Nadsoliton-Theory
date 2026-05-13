#!/usr/bin/env python3
"""P1414 SST-1C transport robustness checkpoint (strict-only, no legacy bridge)."""

from __future__ import annotations

import json
from pathlib import Path


def main() -> None:
    base = Path(__file__).resolve().parent
    gen = base / "generated"
    gen.mkdir(parents=True, exist_ok=True)

    # Strict preregistered transport cases around the converged selector candidate.
    transport_cases = [
        {"case": "T0_baseline", "selector_score": 0.4095},
        {"case": "T1_scale_shift", "selector_score": 0.4048},
        {"case": "T2_scheme_shift", "selector_score": 0.4011},
        {"case": "T3_boundary_shift", "selector_score": 0.3979},
    ]

    baseline = transport_cases[0]["selector_score"]
    eps_tol = 0.0100

    drifts = []
    robust_pass = True
    for c in transport_cases[1:]:
        drift = abs(c["selector_score"] - baseline)
        drifts.append({"case": c["case"], "abs_drift": round(drift, 6), "eps_tol": eps_tol, "pass": drift <= eps_tol})
        robust_pass = robust_pass and drift <= eps_tol

    verdict = "PASS_STRICT_TRANSPORT" if robust_pass else "FAIL_STRICT_TRANSPORT_DRIFT"

    summary = {
        "checkpoint_id": "P1414-SST1C",
        "as_of": "2026-05-13",
        "target": "F_nadsoliton => L_SM + L_GR",
        "mode": "strict_only_no_legacy_bridge",
        "source_dependency": "P1413 converged selector candidate",
        "transport_cases": transport_cases,
        "transport_drift_checks": drifts,
        "transport_robustness_pass": robust_pass,
        "verdict": verdict,
        "status": "NO_FALSE_PASS",
        "next_action": (
            "export selector_obstruction_v1 and pivot noncyclic provider class"
            if not robust_pass
            else "advance to selector uniqueness discharge attempt"
        ),
    }

    out = gen / "p1414_sst1c_transport_robustness_summary.json"
    out.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")

    if not robust_pass:
        obstruction = {
            "obstruction_id": "OBSTR-SST1C-TRANSPORT-v1",
            "reason": "transport drift beyond epsilon tolerance",
            "failed_cases": [d for d in drifts if not d["pass"]],
            "recommended_pivot": "new_noncyclic_provider_class",
        }
        (gen / "p1414_selector_obstruction_v1.json").write_text(json.dumps(obstruction, indent=2) + "\n", encoding="utf-8")

    print(json.dumps({"written": str(out), "verdict": verdict}, indent=2))


if __name__ == "__main__":
    main()
