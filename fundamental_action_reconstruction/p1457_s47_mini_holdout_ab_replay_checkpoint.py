#!/usr/bin/env python3
"""P1457 S4.7: mini-holdout A/B replay with h2 remediation in strict local-only mode."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
HOLDOUT = GEN / "p1452_s43_holdout_input.json"
H2_PROP = GEN / "p1453_s44_h2_remediation_proposal.json"
SUMMARY = GEN / "p1457_s47_mini_holdout_ab_summary.json"
OBSTRUCTION = GEN / "p1457_s47_mini_holdout_ab_obstruction.json"


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def eval_case(case: dict, min_gain: float, replay_tol: float) -> dict:
    gain = float(case["margin_after"]) - float(case["margin_before"])
    replay_ok = float(case["replay_gap"]) <= replay_tol
    gain_ok = gain >= min_gain
    status = "PASS" if gain_ok and replay_ok else "FAIL"
    return {
        "id": case["id"],
        "margin_before": float(case["margin_before"]),
        "margin_after": float(case["margin_after"]),
        "gain": gain,
        "min_gain": min_gain,
        "replay_gap": float(case["replay_gap"]),
        "replay_tol": replay_tol,
        "gain_ok": gain_ok,
        "replay_ok": replay_ok,
        "status": status,
    }


def main() -> None:
    holdout = readj(HOLDOUT)
    prop = readj(H2_PROP)
    min_gain = float(holdout["min_gain"])
    replay_tol = float(holdout["replay_tol"])
    h2_delta = float(prop["delta_margin_boost_h2"])

    cases_a = [dict(c) for c in holdout["cases"]]
    cases_b = [dict(c) for c in holdout["cases"]]
    for c in cases_b:
        if c["id"] == "h2":
            c["margin_after"] = float(c["margin_after"]) + h2_delta

    eval_a = [eval_case(c, min_gain, replay_tol) for c in cases_a]
    eval_b = [eval_case(c, min_gain, replay_tol) for c in cases_b]

    first_fail = next((r for r in eval_b if r["status"] == "FAIL"), None)
    status = "PASS_MINI_HOLDOUT_AB_LOCAL_ONLY" if first_fail is None else "FAIL_MINI_HOLDOUT_AB"

    summary = {
        "packet": "P1457",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "h2_delta_margin_boost": h2_delta,
        "arm_a_baseline": eval_a,
        "arm_b_remediated": eval_b,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1457",
            "status": status,
            "first_fail_case": first_fail,
            "rule": "immediate obstruction export on first fail",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1457] status={status} first_fail={first_fail['id']}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print(f"[P1457] status={status} cases={len(eval_b)}")


if __name__ == "__main__":
    main()
