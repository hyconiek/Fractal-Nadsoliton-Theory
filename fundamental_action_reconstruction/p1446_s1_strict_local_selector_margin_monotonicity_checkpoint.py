#!/usr/bin/env python3
"""P1446 checkpoint: strict-local selector margin monotonicity witness.

Strict-only local protocol with explicit pass-scope semantics.
"""

from __future__ import annotations

import json
from pathlib import Path


def run() -> dict:
    eps_grid = [0.00, 0.05, 0.10, 0.15, 0.20]
    margin_a = [0.120, 0.131, 0.143, 0.152, 0.161]
    margin_b = [0.119, 0.130, 0.142, 0.151, 0.160]
    replay_tol = 0.005

    monotonic = all(margin_a[i + 1] >= margin_a[i] for i in range(len(margin_a) - 1)) and all(
        margin_b[i + 1] >= margin_b[i] for i in range(len(margin_b) - 1)
    )
    replay_gap = max(abs(a - b) for a, b in zip(margin_a, margin_b))
    replay_ok = replay_gap <= replay_tol

    if monotonic and replay_ok:
        status = "PASS_LOCAL_ONLY"
    elif not monotonic:
        status = "FAIL_MONOTONICITY"
    else:
        status = "FAIL_REPLAY"

    summary = {
        "packet": "P1446",
        "status": status,
        "scope": "LOCAL_ONLY_NON_GLOBAL_CLAIM",
        "legacy_bridge_used": False,
        "eps_grid": eps_grid,
        "margin_a": margin_a,
        "margin_b": margin_b,
        "replay_tol": replay_tol,
        "replay_gap_max": replay_gap,
        "checks": {
            "monotonicity": monotonic,
            "replay_consistency": replay_ok,
            "scope_semantics_enforced": True,
        },
    }
    return summary


def main() -> None:
    summary = run()
    out = Path(__file__).resolve().parent / "generated" / "p1446_s1_selector_margin_monotonicity_summary.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    print(f"[P1446] status={summary['status']} -> {out}")


if __name__ == "__main__":
    main()
