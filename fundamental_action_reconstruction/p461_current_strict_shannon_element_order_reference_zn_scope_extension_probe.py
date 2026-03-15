#!/usr/bin/env python3
from __future__ import annotations

import json
import math
from datetime import datetime, timezone
from pathlib import Path

ROOT = Path(__file__).resolve().parent
REPO = ROOT.parent
GENERATED = ROOT / "generated"

OUT_JSON = GENERATED / "p461_current_strict_shannon_element_order_reference_zn_scope_extension_probe.json"
OUT_SUMMARY = GENERATED / "p461_current_strict_shannon_element_order_reference_zn_scope_extension_probe_summary.json"


def zn_element_order(n: int, k: int) -> int:
    kk = int(k) % int(n)
    if kk == 0:
        return 1
    return int(n) // math.gcd(kk, int(n))


def f_2m_defect(ord_by_x: list[int], m: int) -> complex:
    n = len(ord_by_x)
    acc = 0.0 + 0.0j
    for x, v in enumerate(ord_by_x):
        ang = 4.0 * math.pi * float(m) * float(x) / float(n)
        acc += float(v) * complex(math.cos(ang), math.sin(ang))
    return acc


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n_list = [6, 8, 10, 12, 14, 16, 18, 20, 24]
    tol = 1e-12

    results = []
    n_with_all_pairs_cut: list[int] = []
    n_with_any_zero_defect: list[int] = []

    for n in n_list:
        ord_by_x = [zn_element_order(n, x) for x in range(n)]
        pair_m_max = (n // 2 - 1) if (n % 2 == 0) else ((n - 1) // 2)
        pairs = []
        all_cut = True
        any_zero = False

        for m in range(1, pair_m_max + 1):
            f = f_2m_defect(ord_by_x, m=m)
            re = float(f.real)
            im = float(f.imag)
            absf = float(abs(f))
            cut = bool(absf > tol)
            if not cut:
                all_cut = False
                any_zero = True
            theta_star = 0.5 * math.atan2(im, re) if cut else None
            pairs.append(
                {
                    "m": m,
                    "F_2m": {"re": re, "im": im, "abs": absf},
                    "cut_O2_to_Z2_axis_only": cut,
                    "theta_star_half_arg_if_cut": theta_star,
                }
            )

        if all_cut:
            n_with_all_pairs_cut.append(n)
        if any_zero:
            n_with_any_zero_defect.append(n)

        results.append(
            {
                "n": n,
                "pair_m_range": [1, pair_m_max],
                "ord_by_x": ord_by_x,
                "pairs": pairs,
                "all_pairs_cut": all_cut,
                "any_zero_defect": any_zero,
            }
        )

    recommendation = {
        "recommended_next_strict_target": None,
        "reason": None,
    }
    if any(n != 12 for n in n_with_all_pairs_cut):
        recommendation = {
            "recommended_next_strict_target": "SCOPE_EXTENSION_EXPORT_CANDIDATE",
            "reason": (
                "At least one n≠12 in the scanned list has all pair_m defects nonzero for ord_Zn, so a cautious next step "
                "is to export a probe-level mode-index assignment candidate for that n, keeping lane-scoped and no-false-pass."
            ),
        }
    else:
        recommendation = {
            "recommended_next_strict_target": "SCOPE_EXTENSION_BLOCKED_BY_ZERO_DEFECTS_ON_SCANNED_N",
            "reason": (
                "In the scanned n-list, n=12 is the only carrier where all pair_m defects are nonzero. "
                "Scope extension of the Shannon ord-reference cut beyond n=12 therefore appears blocked on the simple "
                "ord_Zn defect criterion, unless a new strict internal ingredient is added."
            ),
        }

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P461",
        "goal": "probe_scope_extension_of_shannon_element_order_reference_o2_cut_by_scanning_F_2m(ord_Zn)_nonzero_conditions",
        "inputs": {"n_list": n_list, "tolerance": tol},
        "results": results,
        "aggregate": {
            "n_with_all_pairs_cut": n_with_all_pairs_cut,
            "n_with_any_zero_defect": n_with_any_zero_defect,
        },
        "recommendation": recommendation,
        "hard_limits": [
            "Probe-only: does not promote any n≠12 results to theorem level.",
            "Does not export a strict mode-index assignment object for n≠12.",
            "Does not claim any QW-2191 discharge beyond declared n=12 lanes.",
            "Does not claim ToE closure.",
        ],
        "no_false_pass": True,
    }

    summary = {
        "stage": "P461",
        "status": "PASS_PROBE_READY",
        "n_scanned": n_list,
        "n_with_all_pairs_cut": n_with_all_pairs_cut,
        "recommendation": recommendation,
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()

