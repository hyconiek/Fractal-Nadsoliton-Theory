#!/usr/bin/env python3
"""P1146: strict Shannon selector candidate grid nonclosure probe.

Constructs the strict-only Shannon-weighted candidate channel
S_cand(d) = exp(-alpha*d) * cos(omega*d+phi)/(1+beta*d**eta)
on a finite nonnegative grid and audits whether it yields a unique-selector-like
single-sign profile. Multiple sign sectors are treated as a nonclosure witness
proxy for QW-2191 compatibility discipline.
"""

from __future__ import annotations

import json
import math
from dataclasses import asdict, dataclass
from pathlib import Path


@dataclass(frozen=True)
class Params:
    omega: float = 0.18575
    phi: float = 0.16250
    beta: float = 1.0
    eta: float = 1.8
    alpha_geo_strict_derived_v1: float = 4.0 * math.log(2.0)


def k_strict_gate(d: float, p: Params) -> float:
    return math.cos(p.omega * d + p.phi) / (1.0 + p.beta * (d ** p.eta))


def s_cand(d: float, p: Params) -> float:
    w_sh = math.exp(-p.alpha_geo_strict_derived_v1 * d)
    return w_sh * k_strict_gate(d, p)


def main() -> None:
    p = Params()

    d_min = 0.0
    d_max = 24.0
    step = 0.1
    n = int(round((d_max - d_min) / step)) + 1

    grid = [d_min + i * step for i in range(n)]
    vals = [s_cand(d, p) for d in grid]

    pos = sum(1 for v in vals if v > 0.0)
    neg = sum(1 for v in vals if v < 0.0)
    zero = n - pos - neg

    sign_changes = 0
    prev = 0
    for v in vals:
        s = 1 if v > 0 else (-1 if v < 0 else 0)
        if s != 0 and prev != 0 and s != prev:
            sign_changes += 1
        if s != 0:
            prev = s

    uniqueness_proxy_pass = neg == 0 and sign_changes == 0
    verdict = "BLOCKED" if not uniqueness_proxy_pass else "CANDIDATE_PASS_ONLY"

    out = {
        "packet": "P1146",
        "status": "executed",
        "as_of": "2026-05-10",
        "strict_side_params": asdict(p),
        "domain": {"d_min": d_min, "d_max": d_max, "step": step, "points": n},
        "observables": {
            "min_s_cand": min(vals),
            "max_s_cand": max(vals),
            "positive_count": pos,
            "negative_count": neg,
            "zero_count": zero,
            "sign_change_count": sign_changes,
        },
        "audit": {
            "criterion": "single-sign/non-oscillatory proxy required for selector-like uniqueness hint",
            "uniqueness_proxy_pass": uniqueness_proxy_pass,
            "qw_2191_nonclosure_verdict": verdict,
            "interpretation": (
                "Grid-level proxy only. Not a theorem."
                "If oscillatory/sign-changing, candidate alone does not justify strict-core selector closure."
            ),
        },
    }

    out_path = Path(__file__).resolve().parent / "generated" / (
        "p1146_strict_shannon_selector_candidate_grid_nonclosure_probe_summary.json"
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(out, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(f"[P1146] wrote {out_path}")


if __name__ == "__main__":
    main()
