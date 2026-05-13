#!/usr/bin/env python3
"""P1447 checkpoint: strict extended-local transport robustness.

Consumes preregistered scenarios and exports local-only verdict.
"""

from __future__ import annotations

import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
IN_PATH = ROOT / "generated" / "p1447_s2_transport_scenarios_input.json"
OUT_PATH = ROOT / "generated" / "p1447_s2_transport_robustness_summary.json"


def load_input() -> dict:
    return json.loads(IN_PATH.read_text(encoding="utf-8"))


def evaluate(payload: dict) -> dict:
    replay_tol = float(payload["replay_tol"])
    scenarios = payload["scenarios"]

    margin_ok = all(s["margin_after"] >= s["margin_before"] for s in scenarios)
    replay_ok = all(s["replay_gap"] <= replay_tol for s in scenarios)

    scope = payload.get("scope_of_pass", "")
    qw_closed = bool(payload.get("strict_core_qw2191_closed", False))
    scope_ok = (scope == "LOCAL_ONLY_NON_GLOBAL_CLAIM") and (not qw_closed)

    if not scope_ok:
        status = "FAIL_SCOPE"
    elif not margin_ok:
        status = "FAIL_TRANSPORT_MARGIN"
    elif not replay_ok:
        status = "FAIL_TRANSPORT_REPLAY"
    else:
        status = "PASS_LOCAL_TRANSPORT_ONLY"

    return {
        "packet": "P1447",
        "status": status,
        "scope_of_pass": scope,
        "strict_core_qw2191_closed": qw_closed,
        "legacy_bridge_used": False,
        "replay_tol": replay_tol,
        "checks": {
            "margin_non_decreasing": margin_ok,
            "replay_within_tolerance": replay_ok,
            "scope_semantics_enforced": scope_ok,
        },
        "scenario_count": len(scenarios),
        "scenarios": scenarios,
    }


def main() -> None:
    payload = load_input()
    result = evaluate(payload)
    OUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUT_PATH.write_text(json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1447] status={result['status']} -> {OUT_PATH}")


if __name__ == "__main__":
    main()
