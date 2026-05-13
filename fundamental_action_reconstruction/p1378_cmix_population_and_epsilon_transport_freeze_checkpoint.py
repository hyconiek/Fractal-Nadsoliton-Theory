#!/usr/bin/env python3
from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
IN_P1377 = GEN / "p1377_ci_class_population_and_first_transport_drift_run_summary.json"
OUT = GEN / "p1378_cmix_population_and_epsilon_transport_freeze_summary.json"


def rel_drift(a: float, b: float) -> float:
    return abs(a - b) / max(1.0, abs(a))


def main() -> None:
    p1377 = json.loads(IN_P1377.read_text(encoding="utf-8"))

    cmix = {"A": 1.00, "B": 0.97}
    drift_mix = rel_drift(cmix["A"], cmix["B"])
    epsilon_transport_v1 = 0.035

    out = {
        "artifact": "p1378_cmix_population_and_epsilon_transport_freeze_summary",
        "as_of": "2026-05-12",
        "input_dependency": IN_P1377.name,
        "input_l_b1_02_status": p1377.get("l_b1_02_status"),
        "cmix_population": cmix,
        "drift_mix": drift_mix,
        "epsilon_transport_v1": epsilon_transport_v1,
        "readiness": "FULL_CLASS_AND_EPSILON_FROZEN",
        "forbidden_claim": "do_not_claim_invariance_before_full_formal_run",
        "next_packet": "P1379_FIRST_FULL_L_B1_02_FORMAL_PASS_FAIL_RUN"
    }

    OUT.write_text(json.dumps(out, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(f"[p1378] wrote {OUT}")


if __name__ == "__main__":
    main()
