#!/usr/bin/env python3
"""P1452 checkpoint: local holdout replay after targeted rerun."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
INP = ROOT / "generated" / "p1452_s43_holdout_input.json"
SUM = ROOT / "generated" / "p1452_s43_holdout_summary.json"
OBS = ROOT / "generated" / "p1452_s43_holdout_obstruction.json"


def main() -> None:
    payload = json.loads(INP.read_text(encoding="utf-8"))
    scope_ok = payload.get("scope_of_pass") == "LOCAL_ONLY_NON_GLOBAL_CLAIM"
    qw_ok = payload.get("strict_core_qw2191_closed") is False
    min_gain = float(payload["min_gain"])
    replay_tol = float(payload["replay_tol"])

    status = "PASS_HOLDOUT_LOCAL_ONLY"
    obstruction = None

    if not (scope_ok and qw_ok):
        status = "FAIL_HOLDOUT_SCOPE"
        obstruction = {"reason": "scope", "scope_ok": scope_ok, "qw_ok": qw_ok}
    else:
        for s in payload["cases"]:
            gain = s["margin_after"] - s["margin_before"]
            if gain < min_gain:
                status = "FAIL_HOLDOUT_MARGIN"
                obstruction = {"reason": "margin", "case": s, "gain": gain, "min_gain": min_gain}
                break
            if s["replay_gap"] > replay_tol:
                status = "FAIL_HOLDOUT_REPLAY"
                obstruction = {
                    "reason": "replay",
                    "case": s,
                    "replay_gap": s["replay_gap"],
                    "replay_tol": replay_tol,
                }
                break

    summary = {
        "packet": "P1452",
        "status": status,
        "scope_of_pass": payload.get("scope_of_pass"),
        "strict_core_qw2191_closed": payload.get("strict_core_qw2191_closed"),
        "legacy_bridge_used": False,
        "case_count": len(payload["cases"]),
    }

    SUM.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    if obstruction is not None:
        OBS.write_text(json.dumps(obstruction, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    elif OBS.exists():
        OBS.unlink()

    print(f"[P1452] status={status} -> {SUM}")
    if obstruction is not None:
        print(f"[P1452] obstruction exported -> {OBS}")


if __name__ == "__main__":
    main()
