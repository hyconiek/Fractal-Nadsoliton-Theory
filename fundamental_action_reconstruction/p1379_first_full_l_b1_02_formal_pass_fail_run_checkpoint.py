#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_P1377 = GEN / "p1377_ci_class_population_and_first_transport_drift_run_summary.json"
IN_P1378 = GEN / "p1378_cmix_population_and_epsilon_transport_freeze_summary.json"
OUT = GEN / "p1379_first_full_l_b1_02_formal_pass_fail_run_summary.json"


def main() -> None:
    p1377 = json.loads(IN_P1377.read_text(encoding="utf-8"))
    p1378 = json.loads(IN_P1378.read_text(encoding="utf-8"))

    drift_per = dict(p1377["drift_per_coefficient"])
    drift_per["c_mix"] = p1378["drift_mix"]

    epsilon = float(p1378["epsilon_transport_v1"])
    max_drift = max(drift_per.values())
    verdict = "PASS_V1" if max_drift <= epsilon else "FAIL_V1"

    out = {
        "artifact": "p1379_first_full_l_b1_02_formal_pass_fail_run_summary",
        "as_of": "2026-05-12",
        "input_dependencies": [IN_P1377.name, IN_P1378.name],
        "drift_per_coefficient": drift_per,
        "max_drift": max_drift,
        "epsilon_transport_v1": epsilon,
        "l_b1_02_formal_verdict": verdict,
        "scientific_scope_limit": "v1_transport_invariance_only_not_global_sm_gr_closure",
        "next_packet": "P1380_L_B1_01_SU3_SU2_U1_IMAGE_CLOSURE_THEOREM_ATTEMPT"
    }

    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[p1379] wrote {OUT}")


if __name__ == "__main__":
    main()
