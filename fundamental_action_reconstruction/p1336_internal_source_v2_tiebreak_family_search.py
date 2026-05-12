"""P1336: internal-source v2 search over constrained tie-break family."""

from __future__ import annotations

import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
GEN = ROOT / "generated"
OUT = GEN / "p1336_internal_source_v2_tiebreak_family_search_report_v1.json"

PHASE_GRID = [0.25 + 0.005 * i for i in range(61)]
AMP_GRID = [0.45 + 0.005 * j for j in range(51)]
DELTA = 0.03
ALPHAS = [i * 0.001 for i in range(0, 31)]  # [0.0, 0.03]


def score(phase: float, amp: float) -> float:
    return (phase - 0.40) + 0.9 * (amp - 0.58) - 0.6 * (phase - 0.40) ** 2


def dphase(phase: float) -> float:
    return 1.0 - 1.2 * (phase - 0.40)


def sign(x: float) -> int:
    return 1 if x >= 0 else -1


def selector_alpha(phase: float, amp: float, alpha: float) -> int:
    s = score(phase, amp)
    if abs(s) >= DELTA:
        return sign(s)
    return sign(s + alpha * dphase(phase))


def axiom_selector(phase: float, amp: float) -> int:
    return sign(score(phase, amp))


def main() -> None:
    results = []
    best = None

    for alpha in ALPHAS:
        mismatches = 0
        for p in PHASE_GRID:
            for a in AMP_GRID:
                if selector_alpha(p, a, alpha) != axiom_selector(p, a):
                    mismatches += 1
        results.append({"alpha": alpha, "mismatches": mismatches})
        if best is None or mismatches < best["mismatches"]:
            best = {"alpha": alpha, "mismatches": mismatches}

    exportable = best["mismatches"] == 0

    payload = {
        "packet_id": "P1336_INTERNAL_SOURCE_V2_TIEBREAK_FAMILY_SEARCH_REPORT_V1",
        "date_utc": "2026-05-12",
        "delta": DELTA,
        "alpha_grid": {"min": min(ALPHAS), "max": max(ALPHAS), "count": len(ALPHAS)},
        "best": best,
        "strict_internal_source_exportable": exportable,
        "status": "V2_EXPORTABLE" if exportable else "V2_NOT_EXPORTABLE",
        "notes": "Best alpha may collapse to boundary-neutral tie-break; requires theorem-level interpretation.",
        "qw2191_strict_status": "NOT_CLOSED",
    }

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    print(f"Wrote {OUT}")


if __name__ == "__main__":
    main()
