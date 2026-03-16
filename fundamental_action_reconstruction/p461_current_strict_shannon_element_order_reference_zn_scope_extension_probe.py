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

N514_THEOREM = (
    ROOT / "N514_CURRENT_FIRST_STRICT_ZN_ELEMENT_ORDER_REFERENCE_FOURIER_DEFECT_NONVANISHING_THEOREM.md"
)


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


def factorize(n: int) -> dict[int, int]:
    nn = int(n)
    if nn < 2:
        return {}
    out: dict[int, int] = {}
    d = 2
    while d * d <= nn:
        while nn % d == 0:
            out[d] = out.get(d, 0) + 1
            nn //= d
        d = 3 if d == 2 else d + 2
    if nn > 1:
        out[nn] = out.get(nn, 0) + 1
    return out


def v_p(k: int, p: int) -> int:
    kk = int(k)
    pp = int(p)
    v = 0
    while kk % pp == 0:
        kk //= pp
        v += 1
    return v


def prime_power_defect_F_k_ord(p: int, a: int, k: int) -> int:
    # N514 Lemma 3
    pp = int(p)
    aa = int(a)
    v = v_p(k, pp)
    if v < aa:
        num = pp ** (2 * v + 2) - 1
        den = pp + 1
        if num % den != 0:
            raise ValueError("non-integer prime-power defect formula (unexpected)")
        return -(num // den)
    # v >= a
    # sum_{t=0}^{a-1} p^{2t} = (p^{2a}-1)/(p^2-1)
    num = pp ** (2 * aa) - 1
    den = pp * pp - 1
    if den == 0 or num % den != 0:
        raise ValueError("non-integer geometric sum (unexpected)")
    s = num // den
    return 1 + (pp - 1) * pp * s


def defect_F_k_ord_exact(n: int, k: int) -> int:
    # N514 Theorem + Lemma 4: multiplicative in n for fixed k
    if n < 2:
        raise ValueError("n must be >=2")
    if k == 0:
        raise ValueError("k must be nonzero for the Fourier-degenerate pair defects")
    fac = factorize(n)
    out = 1
    for p, a in fac.items():
        out *= prime_power_defect_F_k_ord(p, a, k)
    return int(out)


def main() -> None:
    GENERATED.mkdir(exist_ok=True)

    n_list = [6, 8, 10, 12, 14, 16, 18, 20, 24]
    tol = 1e-12  # retained for float cross-check only

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
            k = 2 * int(m)
            re_int = defect_F_k_ord_exact(n, k=k)
            cut = bool(re_int != 0)
            if not cut:
                all_cut = False
                any_zero = True

            # float cross-check (no-false-pass hygiene): validate exact formula against direct sum
            f = f_2m_defect(ord_by_x, m=m)
            if abs(float(f.imag)) > tol:
                raise SystemExit(
                    json.dumps(
                        {
                            "stage": "P461",
                            "status": "FAIL_UNEXPECTED_NONREAL_DEFECT",
                            "n": n,
                            "m": m,
                            "F_2m_float": {"re": float(f.real), "im": float(f.imag), "abs": float(abs(f))},
                            "tolerance": tol,
                            "no_false_pass": True,
                        },
                        ensure_ascii=True,
                    )
                )
            if abs(float(f.real) - float(re_int)) > 1e-6:
                raise SystemExit(
                    json.dumps(
                        {
                            "stage": "P461",
                            "status": "FAIL_DEFECT_MISMATCH_EXACT_VS_FLOAT",
                            "n": n,
                            "m": m,
                            "k": k,
                            "F_2m_exact_re": int(re_int),
                            "F_2m_float_re": float(f.real),
                            "tolerance": 1e-6,
                            "no_false_pass": True,
                        },
                        ensure_ascii=True,
                    )
                )

            re = float(re_int)
            im = 0.0
            absf = float(abs(re_int))
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
        "recommended_next_strict_target": "SCOPE_EXTENSION_INFRA_ONLY",
        "reason": (
            "N514 now proves the nonvanishing defect condition F_k(ord_{Z_n})≠0 for all n≥2 and k∈{1,…,n-1}. "
            "Therefore this P461 scan is retained only as a computational regression/sanity harness. "
            "Any real scope-extension move must still separately export a typed carrier + an assignment object for the chosen n "
            "(e.g. Z_24 via F458/F468), without any physical promotion."
        ),
    }

    artifact = {
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "stage": "P461",
        "goal": "probe_scope_extension_of_shannon_element_order_reference_o2_cut_by_scanning_F_2m(ord_Zn)_nonzero_conditions",
        "inputs": {"n_list": n_list, "tolerance": tol},
        "method": {
            "defect_computation": "exact_integer_formula_from_N514 (prime-power formula + multiplicativity)",
            "float_cross_check": "direct_sum_float_cos_sin (must match exact within 1e-6; imag must be <=1e-12)",
            "theorem_ref": (str(N514_THEOREM.relative_to(REPO)) if N514_THEOREM.exists() else "N514"),
        },
        "results": results,
        "aggregate": {
            "n_with_all_pairs_cut": n_with_all_pairs_cut,
            "n_with_any_zero_defect": n_with_any_zero_defect,
        },
        "recommendation": recommendation,
        "hard_limits": [
            "Probe-only: regression/sanity scan; defect nonvanishing is theorem-level (N514).",
            "Does not export any typed mode-index assignment object for n≠12 (except where separately exported).",
            "Does not claim any physical promotion of n≠12 into QW-2190.",
            "Does not claim any global QW-2191 discharge beyond declared lanes.",
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
        "theorem_level_defect_nonvanishing_ref": (str(N514_THEOREM.relative_to(REPO)) if N514_THEOREM.exists() else "N514"),
        "no_false_pass": True,
    }

    OUT_JSON.write_text(json.dumps(artifact, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    OUT_SUMMARY.write_text(json.dumps(summary, indent=2, ensure_ascii=True) + "\n", encoding="ascii")
    print(OUT_SUMMARY)


if __name__ == "__main__":
    main()
