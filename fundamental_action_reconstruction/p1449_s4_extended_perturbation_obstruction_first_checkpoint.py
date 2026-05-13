#!/usr/bin/env python3
"""P1449 checkpoint: extended perturbation stress test with obstruction-first policy."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1449_s4_stress_scenarios_input.json"
SUM_PATH = ROOT / "generated" / "p1449_s4_stress_summary.json"
OBS_PATH = ROOT / "generated" / "p1449_s4_first_obstruction.json"


def main() -> None:
    payload = json.loads(IN_PATH.read_text(encoding="utf-8"))
    scenarios = payload["scenarios"]
    replay_tol = float(payload["replay_tol"])
    min_gain = float(payload["min_gain"])

    scope_ok = payload.get("scope_of_pass") == "LOCAL_ONLY_NON_GLOBAL_CLAIM"
    qw_ok = payload.get("strict_core_qw2191_closed") is False

    status = "PASS_LOCAL_STRESS_ONLY"
    first_obstruction = None

    if not (scope_ok and qw_ok):
        status = "FAIL_STRESS_SCOPE"
        first_obstruction = {"reason": "scope", "scope_ok": scope_ok, "qw_ok": qw_ok}
    else:
        for s in scenarios:
            gain = s["margin_after"] - s["margin_before"]
            if gain < min_gain:
                status = "FAIL_STRESS_MARGIN"
                first_obstruction = {"reason": "margin", "scenario": s, "gain": gain, "min_gain": min_gain}
                break
            if s["replay_gap"] > replay_tol:
                status = "FAIL_STRESS_REPLAY"
                first_obstruction = {
                    "reason": "replay",
                    "scenario": s,
                    "replay_gap": s["replay_gap"],
                    "replay_tol": replay_tol,
                }
                break

    summary = {
        "packet": "P1449",
        "status": status,
        "scope_of_pass": payload.get("scope_of_pass"),
        "strict_core_qw2191_closed": payload.get("strict_core_qw2191_closed"),
        "legacy_bridge_used": False,
        "scenario_count": len(scenarios),
        "checks": {
            "scope_ok": scope_ok,
            "qw2191_not_closed": qw_ok,
        },
    }

    SUM_PATH.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if first_obstruction is not None:
        OBS_PATH.write_text(json.dumps(first_obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBS_PATH.exists():
        OBS_PATH.unlink()

    print(f"[P1449] status={status} -> {SUM_PATH}")
    if first_obstruction is not None:
        print(f"[P1449] obstruction exported -> {OBS_PATH}")


if __name__ == "__main__":
    main()
