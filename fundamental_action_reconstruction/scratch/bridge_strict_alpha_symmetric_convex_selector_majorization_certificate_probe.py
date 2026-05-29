#!/usr/bin/env python3
"""Scratch probe: symmetric-convex / majorization selector certificate.

The previous Fourier minimum-ripple packet proved the integer convexity part for
one quadratic selector.  This packet takes the next honest step: isolate the
larger mathematical class of selectors for which the balanced alpha ledger

    (2, 2, 2, 1, 1)

is forced once the branch model has already supplied m=5 and n=8.

Model-internal certificate proved here:

1. Among positive integer ledgers with fixed branch count m and total n, the
   floor/ceil balanced ledger is the unique majorization-minimal canonical
   ledger.
2. Any symmetric strictly convex separable branch action, e.g.
   E_p(e)=sum_i e_i^p for p>1, is Schur-convex and therefore has its unique
   canonical minimum at the balanced ledger.
3. Pairwise smoothing gives a constructive certificate: if a>=b+2, replacing
   (a,b) by (a-1,b+1) strictly decreases each checked E_p.

This is a broader theorem-shaped reduction than the Fourier-specific Parseval
certificate, but it is still conditional.  It does not derive a strict
nadsoliton action whose branch-sector cost is symmetric strictly convex, nor
branch count, total exponent, denominator normalization, binary-rescale quotient,
or phase-labelled branches.  It therefore does not discharge QW-2191.
"""
from __future__ import annotations

import json
import math
from fractions import Fraction
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
MIN_RIPPLE = HERE / "bridge_strict_alpha_fourier_min_ripple_convexity_certificate_report.json"
OUT_JSON = HERE / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.json"
OUT_MD = HERE / "bridge_strict_alpha_symmetric_convex_selector_majorization_certificate_report.md"

TARGET_BRANCH_COUNT = 5
TARGET_BINARY_EXPONENT = 8
DENOMINATOR = 3
STRICT_TARGET_ETA = 9.0 / 5.0
NAD12_SUPPORT_SIZE = 12
ALPHA_SCALE = 4.0
POWERS_CHECKED = [2, 3, 4, 5, 6]
SCAN_MAX_BRANCH_COUNT = 9
SCAN_MAX_TOTAL_EXPONENT = 24


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def eta_from_product(product: Fraction, branch_count: int) -> float:
    correction = float(product) ** (1.0 / branch_count)
    return math.log(NAD12_SUPPORT_SIZE * correction) / math.log(ALPHA_SCALE)


def canonical_desc(values: tuple[int, ...]) -> tuple[int, ...]:
    return tuple(sorted(values, reverse=True))


def balanced_ledger(branch_count: int, total: int) -> tuple[int, ...]:
    base, remainder = divmod(total, branch_count)
    if base <= 0:
        raise ValueError("positive ledgers require total >= branch_count")
    return tuple([base + 1] * remainder + [base] * (branch_count - remainder))


def positive_partitions(total: int, parts: int, max_part: int | None = None) -> list[tuple[int, ...]]:
    if parts == 0:
        return [()] if total == 0 else []
    if max_part is None:
        max_part = total
    out: list[tuple[int, ...]] = []
    min_remaining = parts - 1
    for first in range(min(max_part, total - min_remaining), 0, -1):
        for rest in positive_partitions(total - first, parts - 1, min(first, total - first)):
            if len(rest) == parts - 1:
                out.append((first,) + rest)
    return out


def prefix_sums(ledger: tuple[int, ...]) -> list[int]:
    running = 0
    out = []
    for value in ledger:
        running += value
        out.append(running)
    return out


def majorizes(left: tuple[int, ...], right: tuple[int, ...]) -> bool:
    """Return True if canonical ``left`` majorizes canonical ``right``."""
    if sum(left) != sum(right) or len(left) != len(right):
        return False
    return all(a >= b for a, b in zip(prefix_sums(left)[:-1], prefix_sums(right)[:-1]))


def power_energy(ledger: tuple[int, ...], power: int) -> int:
    return sum(value**power for value in ledger)


def smoothing_step(ledger: tuple[int, ...]) -> dict[str, Any] | None:
    values = list(ledger)
    high_index = max(range(len(values)), key=lambda index: values[index])
    low_index = min(range(len(values)), key=lambda index: values[index])
    high = values[high_index]
    low = values[low_index]
    if high < low + 2:
        return None
    before = canonical_desc(tuple(values))
    values[high_index] -= 1
    values[low_index] += 1
    after = canonical_desc(tuple(values))
    return {
        "before": list(before),
        "after": list(after),
        "high": high,
        "low": low,
        "gap": high - low,
        "energy_drops": {
            f"p_{power}": power_energy(before, power) - power_energy(after, power)
            for power in POWERS_CHECKED
        },
    }


def smoothing_path(ledger: tuple[int, ...]) -> list[dict[str, Any]]:
    current = canonical_desc(ledger)
    path = []
    while True:
        step = smoothing_step(current)
        if step is None:
            break
        path.append(step)
        current = tuple(step["after"])
    return path


def target_majorization_certificate() -> dict[str, Any]:
    target = balanced_ledger(TARGET_BRANCH_COUNT, TARGET_BINARY_EXPONENT)
    ledgers = positive_partitions(TARGET_BINARY_EXPONENT, TARGET_BRANCH_COUNT)
    rows = []
    for ledger in ledgers:
        product = Fraction(2 ** sum(ledger), DENOMINATOR ** len(ledger))
        rows.append(
            {
                "ledger": list(ledger),
                "prefix_sums": prefix_sums(ledger),
                "majorizes_balanced": majorizes(ledger, target),
                "balanced_majorizes_this": majorizes(target, ledger),
                "power_energies": {f"p_{power}": power_energy(ledger, power) for power in POWERS_CHECKED},
                "smoothing_terminal": list(smoothing_path(ledger)[-1]["after"] if smoothing_path(ledger) else ledger),
                "product_fraction": f"{product.numerator}/{product.denominator}",
                "eta_residual_vs_9_5": eta_from_product(product, len(ledger)) - STRICT_TARGET_ETA,
            }
        )
    minimizers_by_power = {}
    for power in POWERS_CHECKED:
        key = f"p_{power}"
        min_value = min(row["power_energies"][key] for row in rows)
        minimizers_by_power[key] = {
            "min_value": min_value,
            "minimizers": [row["ledger"] for row in rows if row["power_energies"][key] == min_value],
        }
    return {
        "branch_count": TARGET_BRANCH_COUNT,
        "total_binary_exponent": TARGET_BINARY_EXPONENT,
        "balanced_ledger": list(target),
        "canonical_positive_ledgers": [list(ledger) for ledger in ledgers],
        "rows": rows,
        "all_ledgers_majorize_balanced": all(row["majorizes_balanced"] for row in rows),
        "balanced_majorizes_only_itself": [row["ledger"] for row in rows if row["balanced_majorizes_this"]] == [list(target)],
        "unique_minimizers_by_power": minimizers_by_power,
        "all_checked_powers_select_balanced": all(value["minimizers"] == [list(target)] for value in minimizers_by_power.values()),
        "sample_smoothing_paths": {
            "from_4_1_1_1_1": smoothing_path((4, 1, 1, 1, 1)),
            "from_3_2_1_1_1": smoothing_path((3, 2, 1, 1, 1)),
            "from_2_2_2_1_1": smoothing_path((2, 2, 2, 1, 1)),
        },
    }


def bounded_majorization_scan() -> dict[str, Any]:
    cases_checked = 0
    ledgers_checked = 0
    failures = []
    max_partition_count = 0
    max_partition_case: dict[str, int] | None = None
    sample_rows = []
    for branch_count in range(1, SCAN_MAX_BRANCH_COUNT + 1):
        for total in range(branch_count, SCAN_MAX_TOTAL_EXPONENT + 1):
            cases_checked += 1
            ledgers = positive_partitions(total, branch_count)
            ledgers_checked += len(ledgers)
            target = balanced_ledger(branch_count, total)
            if len(ledgers) > max_partition_count:
                max_partition_count = len(ledgers)
                max_partition_case = {"branch_count": branch_count, "total": total}
            all_majorize = all(majorizes(ledger, target) for ledger in ledgers)
            balanced_only_self = [ledger for ledger in ledgers if majorizes(target, ledger)] == [target]
            powers_ok = True
            for power in POWERS_CHECKED:
                scored = [(power_energy(ledger, power), ledger) for ledger in ledgers]
                min_score = min(score for score, _ledger in scored)
                minimizers = [ledger for score, ledger in scored if score == min_score]
                if minimizers != [target]:
                    powers_ok = False
                    failures.append({"branch_count": branch_count, "total": total, "power": power, "minimizers": [list(item) for item in minimizers], "expected": list(target)})
            if not all_majorize or not balanced_only_self:
                failures.append({"branch_count": branch_count, "total": total, "majorization_ok": all_majorize, "balanced_only_self": balanced_only_self, "expected": list(target)})
            if branch_count in {5, 7, 9} and total in {branch_count, branch_count + 3, 24}:
                sample_rows.append(
                    {
                        "branch_count": branch_count,
                        "total": total,
                        "partition_count": len(ledgers),
                        "balanced_ledger": list(target),
                        "all_majorize_balanced": all_majorize,
                        "checked_powers_select_balanced": powers_ok,
                    }
                )
    return {
        "scan_limits": {"max_branch_count": SCAN_MAX_BRANCH_COUNT, "max_total_exponent": SCAN_MAX_TOTAL_EXPONENT, "powers_checked": POWERS_CHECKED},
        "cases_checked": cases_checked,
        "canonical_ledgers_checked": ledgers_checked,
        "max_partition_count": max_partition_count,
        "max_partition_case": max_partition_case,
        "all_cases_passed": not failures,
        "failures": failures,
        "sample_rows": sample_rows,
    }


def main() -> None:
    min_ripple_report = load_json(MIN_RIPPLE)
    target = target_majorization_certificate()
    scan = bounded_majorization_scan()
    product = Fraction(2 ** TARGET_BINARY_EXPONENT, DENOMINATOR ** TARGET_BRANCH_COUNT)

    report = {
        "status": "OPEN_STRICT_ALPHA_SYMMETRIC_CONVEX_SELECTOR_MAJORIZATION_CERTIFICATE_NO_PHYSICAL_SELECTOR_THEOREM",
        "result_kind": "SCRATCH_STRICT_ALPHA_SYMMETRIC_CONVEX_SELECTOR_MAJORIZATION_CERTIFICATE_PROBE__NOT_A_THEOREM",
        "source_reports": {
            "fourier_min_ripple_convexity_certificate": str(MIN_RIPPLE.relative_to(ROOT)),
        },
        "previous_certificate_replay": {
            "result_kind": min_ripple_report["result_kind"],
            "balanced_ledger": min_ripple_report["target_eta_9_5_certificate"]["balanced_ledger"],
            "q_power_product": min_ripple_report["target_eta_9_5_certificate"]["q_power_product"],
        },
        "conditional_theorem_statement": {
            "premises": [
                "fixed positive integer branch count m",
                "fixed total binary exponent n",
                "canonical ledgers sorted in nonincreasing order",
                "symmetric strictly convex separable branch action, or any Schur-convex branch action",
            ],
            "majorization_claim": "Every positive integer ledger with fixed (m,n) majorizes the floor/ceil balanced ledger; the balanced ledger majorizes only itself canonically.",
            "smoothing_claim": "Repeated Robin-Hood/pairwise smoothing strictly decreases every checked convex power energy E_p=sum_i e_i^p for p>1 until the balanced ledger is reached.",
            "conclusion": "If a strict action exports such a symmetric strictly convex selector, the integer ledger selector is forced to the balanced floor/ceil ledger.",
            "physical_status": "majorization/convexity certificate only; strict nadsoliton action source is not derived here",
        },
        "target_eta_9_5_certificate": {
            "q_power_product": f"{product.numerator}/{product.denominator}",
            "eta_residual_vs_9_5": eta_from_product(product, TARGET_BRANCH_COUNT) - STRICT_TARGET_ETA,
            **target,
        },
        "bounded_majorization_scan": scan,
        "candidate_interpretation": {
            "supported_by_this_probe": bool(target["all_ledgers_majorize_balanced"] and target["all_checked_powers_select_balanced"] and scan["all_cases_passed"]),
            "content": "The selector obligation is sharpened: it is enough to derive any symmetric strictly convex / Schur-convex branch action; then (2,2,2,1,1) follows for m=5,n=8.",
            "why_this_is_more_proof_like": "The probe replaces a single quadratic Fourier functional by a majorization theorem class, pairwise smoothing certificates, and a bounded exhaustive scan across multiple convex powers.",
            "why_this_is_not_enough": "The probe does not derive the symmetric strictly convex branch action, branch count m=5, total exponent n=8, denominator 3, binary-rescale quotient, or phase-labelled branches from strict nadsoliton geometry.",
        },
        "hard_limits": [
            "No identity K_legacy_ont == K_strict_gate is exported.",
            "No legacy physical-role transfer is licensed.",
            "No theorem derives a symmetric strictly convex branch action from strict nadsoliton geometry.",
            "No theorem derives branch count m=5, total exponent n=8, denominator 3, binary-rescale quotient, or phase-labelled branches from strict core in this packet.",
            "This certificate is conditional on a Schur-convex selector/action premise; it is not a strict-core selector source.",
            "No Aut(Z_12)/N462/QW-2191 selector obstruction is discharged.",
            "No theorem derives eta=9/5 without the branch model, denominator/quotient convention, phase labelling, and convex selector/action premises.",
            "No QW-2191 selector discharge and no ToE closure are claimed.",
        ],
        "next_honest_step": "Search for or construct a strict nadsoliton action term whose branch-sector reduction is symmetric strictly convex / Schur-convex; otherwise record this as a clean conditional selector theorem schema, not a strict theorem.",
    }

    OUT_JSON.write_text(json.dumps(report, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    OUT_MD.write_text(
        "# Scratch strict-alpha symmetric-convex selector majorization certificate probe\n\n"
        "Status: majorization/convexity selector certificate for eta=9/5; no physical selector theorem.\n\n"
        "- Conditional theorem: for fixed positive integer `m,n`, every canonical ledger majorizes the floor/ceil balanced ledger.\n"
        "- Selector class: any symmetric strictly convex / Schur-convex branch action selects that balanced ledger.\n"
        f"- Target result: for `m={TARGET_BRANCH_COUNT}, n={TARGET_BINARY_EXPONENT}`, checked powers `{POWERS_CHECKED}` uniquely select `{target['balanced_ledger']}` and `q^5={product.numerator}/{product.denominator}`.\n"
        f"- Bounded scan: `{scan['cases_checked']}` `(m,n)` cases and `{scan['canonical_ledgers_checked']}` canonical ledgers checked; all majorization/convex-power checks passed: `{scan['all_cases_passed']}`.\n"
        "- Honest read: this reduces the selector source obligation to deriving a symmetric strictly convex branch action; it does not derive that action from strict geometry.\n"
        "- No false pass: no strict convex-action source theorem, no QW-2191 discharge, no ToE closure.\n",
        encoding="utf-8",
    )
    print(OUT_JSON)
    print(OUT_MD)


if __name__ == "__main__":
    main()
