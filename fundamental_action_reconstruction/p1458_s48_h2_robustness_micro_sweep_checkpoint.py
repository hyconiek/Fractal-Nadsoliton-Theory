#!/usr/bin/env python3
"""P1458 S4.8: local-only robustness micro-sweep around h2 remediation delta."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
HOLDOUT = GEN / "p1452_s43_holdout_input.json"
PROP = GEN / "p1453_s44_h2_remediation_proposal.json"
SUMMARY = GEN / "p1458_s48_h2_robustness_sweep_summary.json"
OBSTRUCTION = GEN / "p1458_s48_h2_robustness_sweep_obstruction.json"


EPS = 0.0005


def readj(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def evaluate(cases: list[dict], min_gain: float, replay_tol: float) -> list[dict]:
    out = []
    for case in cases:
        gain = float(case["margin_after"]) - float(case["margin_before"])
        gain_ok = gain >= min_gain
        replay_ok = float(case["replay_gap"]) <= replay_tol
        out.append(
            {
                "id": case["id"],
                "gain": gain,
                "gain_ok": gain_ok,
                "replay_gap": float(case["replay_gap"]),
                "replay_ok": replay_ok,
                "status": "PASS" if gain_ok and replay_ok else "FAIL",
            }
        )
    return out


def main() -> None:
    holdout = readj(HOLDOUT)
    prop = readj(PROP)

    min_gain = float(holdout["min_gain"])
    replay_tol = float(holdout["replay_tol"])
    base_delta = float(prop["delta_margin_boost_h2"])
    deltas = [round(base_delta - EPS, 6), round(base_delta, 6), round(base_delta + EPS, 6)]

    variants = []
    first_fail = None
    for delta in deltas:
        cases = [dict(c) for c in holdout["cases"]]
        for case in cases:
            if case["id"] == "h2":
                case["margin_after"] = float(case["margin_after"]) + delta
        eval_rows = evaluate(cases, min_gain, replay_tol)
        all_pass = all(r["status"] == "PASS" for r in eval_rows)
        variant = {"delta_margin_boost_h2": delta, "all_pass": all_pass, "cases": eval_rows}
        variants.append(variant)
        if (not all_pass) and first_fail is None:
            failed_case = next(r for r in eval_rows if r["status"] == "FAIL")
            first_fail = {"delta_margin_boost_h2": delta, "case": failed_case}

    status = "PASS_H2_MICRO_SWEEP_LOCAL_ONLY" if first_fail is None else "FAIL_H2_MICRO_SWEEP"

    summary = {
        "packet": "P1458",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "strict_core_qw2191_closed": False,
        "legacy_bridge_used": False,
        "base_delta_margin_boost_h2": base_delta,
        "epsilon": EPS,
        "variants": variants,
    }
    SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    if first_fail is not None:
        obstruction = {
            "packet": "P1458",
            "status": status,
            "first_fail": first_fail,
            "rule": "immediate obstruction export on first fail",
        }
        OBSTRUCTION.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(f"[P1458] status={status} first_fail={first_fail['case']['id']} delta={first_fail['delta_margin_boost_h2']}")
    else:
        if OBSTRUCTION.exists():
            OBSTRUCTION.unlink()
        print(f"[P1458] status={status} variants={len(variants)}")


if __name__ == "__main__":
    main()
