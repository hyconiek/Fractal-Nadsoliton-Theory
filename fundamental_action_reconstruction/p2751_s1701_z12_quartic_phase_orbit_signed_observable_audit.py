#!/usr/bin/env python3
"""P2751/S1701: Z12 quartic phase orbit signed-observable audit.

After P2750 found no existing artifact supplying the concrete odd source sign,
this pivots to one genuinely new finite typed observable: quartic Z12 phase
sums, quotiented by the same affine source reparametrisation used in P2745/P2747.
The goal is not to choose a polarity by hand, but to test whether a higher-degree
phase object exports an orbit-safe unpaired sign with a P2721-ready source law.
"""
from __future__ import annotations

import cmath
import hashlib
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

N = 12
UNITS = (1, 5, 7, 11)
TOL = 1e-9
GEN = ROOT / "generated"
P2750 = GEN / "p2750_s1700_concrete_odd_source_sign_value_inventory_no_go.json"
OUT = GEN / "p2751_s1701_z12_quartic_phase_orbit_signed_observable_audit.json"
MD = GEN / "p2751_s1701_z12_quartic_phase_orbit_signed_observable_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2750_pivot_obligation": r"P2750|genuinely new strict sign-source mechanism|different typed object|no-new-live-frontier",
    "phase_observable_boundary": r"quartic phase|cubic phase|Gauss-phase|coefficient orbit|orbit-safe signed",
    "p2721_boundary": r"P2721 coupling|P2721 polarity|lambda/P2721|coupling-polarity theorem",
    "closure_forbidden": r"QW-2191|selector closure|role transfer|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "quartic_phase_promoted_to_source_sign",
    "unique_quartic_orbit_polarity_exported",
    "p2721_coupling_theorem_exported",
    "lambda_p2721_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def run_rg(pattern: str) -> list[str]:
    cmd = ["rg", "-n", "--glob", "!generated/*.json", pattern, "AGENTS.md", "fundamental_action_reconstruction"]
    proc = subprocess.run(cmd, cwd=REPO, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, check=False)
    if proc.returncode not in (0, 1):
        raise RuntimeError(f"rg failed for {pattern!r}: {proc.stderr}")
    return [line for line in proc.stdout.splitlines() if line.strip()]


def evidence_scan() -> dict[str, Any]:
    rows = []
    for name, pattern in CONTENT_PATTERNS.items():
        hits = run_rg(pattern)
        rows.append({"content_lane": name, "pattern": pattern, "hit_count": len(hits), "sample_hits": hits[:8]})
    return {"content_pattern_count": len(rows), "rows": rows, "hit_counts": {r["content_lane"]: r["hit_count"] for r in rows}, "all_patterns_have_hits": all(r["hit_count"] > 0 for r in rows)}


def sign(x: float) -> int:
    return 1 if x > TOL else -1 if x < -TOL else 0


def quartic_phase_sum(a: int, b: int, c: int, d: int) -> complex:
    return sum(cmath.exp(2j * math.pi * ((a * n**4 + b * n**3 + c * n * n + d * n) % N) / N) for n in range(N))


def affine_image(a: int, b: int, c: int, d: int, u: int, t: int) -> tuple[int, int, int, int]:
    # Pull back q(n)=a*n^4+b*n^3+c*n^2+d*n by n -> u*n+t, dropping the constant term
    # exactly as the prior coefficient-orbit audits do.
    return (
        (a * u**4) % N,
        (u**3 * (b + 4 * a * t)) % N,
        (u**2 * (c + 3 * b * t + 6 * a * t * t)) % N,
        (u * (d + 2 * c * t + 3 * b * t * t + 4 * a * t**3)) % N,
    )


def quartic_phase_audit() -> dict[str, Any]:
    signs: dict[tuple[int, int, int, int], int] = {}
    sample_rows = []
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    z = quartic_phase_sum(a, b, c, d)
                    s = sign(z.imag)
                    key = (a, b, c, d)
                    signs[key] = s
                    if len(sample_rows) < 16 and s != 0:
                        sample_rows.append({"coefficients": list(key), "real": round(z.real, 12), "imag": round(z.imag, 12), "imag_sign": s})
    seen: set[tuple[int, int, int, int]] = set()
    orbits = []
    for quad in sorted(signs):
        if quad in seen:
            continue
        orbit = {affine_image(*quad, u, t) for u in UNITS for t in range(N)}
        seen |= orbit
        coeff = sum(signs[p] for p in orbit)
        row = {
            "orbit_id": len(orbits),
            "size": len(orbit),
            "signed_sum_coefficient": coeff,
            "positive_rows": sum(1 for p in orbit if signs[p] == 1),
            "negative_rows": sum(1 for p in orbit if signs[p] == -1),
            "zero_rows": sum(1 for p in orbit if signs[p] == 0),
            "representatives_sample": sorted([list(p) for p in orbit])[:10],
        }
        orbits.append(row)
    coeffs = sorted(o["signed_sum_coefficient"] for o in orbits if o["signed_sum_coefficient"] != 0)
    histogram: dict[str, int] = {}
    for coeff in coeffs:
        histogram[str(coeff)] = histogram.get(str(coeff), 0) + 1
    positive_hist = {str(k): v for k, v in sorted((abs(int(k)), v) for k, v in histogram.items() if int(k) > 0)}
    negative_hist = {str(abs(int(k))): v for k, v in sorted((int(k), v) for k, v in histogram.items() if int(k) < 0)}
    return {
        "typed_candidate": "imaginary-sign quartic phase sum Q(a,b,c,d)=sum_n exp(2*pi*i*(a*n^4+b*n^3+c*n^2+d*n)/12), quotiented by affine source reparametrisation",
        "modulus": N,
        "coefficient_quadruple_count": len(signs),
        "pointwise_positive_signs": sum(1 for s in signs.values() if s == 1),
        "pointwise_negative_signs": sum(1 for s in signs.values() if s == -1),
        "pointwise_zero_signs": sum(1 for s in signs.values() if s == 0),
        "global_signed_sum": sum(signs.values()),
        "sample_nonzero_pointwise_rows": sample_rows,
        "affine_orbit_count": len(orbits),
        "affine_orbit_sizes": sorted({o["size"] for o in orbits}),
        "nonzero_orbit_coefficient_count": len(coeffs),
        "nonzero_orbit_coefficients": coeffs,
        "nonzero_orbit_coefficient_histogram": histogram,
        "positive_nonzero_coefficients": sum(1 for c in coeffs if c > 0),
        "negative_nonzero_coefficients": sum(1 for c in coeffs if c < 0),
        "absolute_value_histogram_positive": positive_hist,
        "absolute_value_histogram_negative": negative_hist,
        "histogram_abs_paired": positive_hist == negative_hist,
        "nonzero_orbit_rows_sample": [o for o in orbits if o["signed_sum_coefficient"] != 0][:40],
        "finite_theorem": "Quartic Z12 phase sums form a genuine new typed observable beyond the P2750 inventory replay.  The full finite coefficient audit computes all 12^4 coefficient quadruples and affine coefficient orbits.  Nonzero orbit-safe signed coefficients exist, but the nonzero coefficient multiset remains paired by polarity, and no strict law selects one quartic orbit/polarity or couples it to P2721.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2750_boundary": scan["all_patterns_have_hits"],
        "new_quartic_phase_candidate_computed": audit["coefficient_quadruple_count"] == N**4,
        "has_nonzero_orbit_safe_signed_coefficients": audit["nonzero_orbit_coefficient_count"] > 0,
        "positive_negative_nonzero_counts_balanced": audit["positive_nonzero_coefficients"] == audit["negative_nonzero_coefficients"],
        "absolute_value_histogram_paired": audit["histogram_abs_paired"],
        "strict_quartic_orbit_polarity_source_exported": False,
        "p2721_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_p2749_p2750_completion": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The quartic phase observable supplies nonzero orbit-safe signed coefficients, but the coefficient family is polarity-paired and current artifacts export no strict law selecting one quartic orbit/polarity or an explicit P2721 coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["quartic_phase_audit"]
    lines = [
        "# P2751/S1701 Z12 quartic phase orbit signed observable audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite quartic audit",
        f"- coefficient_quadruple_count={a['coefficient_quadruple_count']}",
        f"- pointwise_signs=+{a['pointwise_positive_signs']} / -{a['pointwise_negative_signs']} / 0={a['pointwise_zero_signs']}",
        f"- global_signed_sum={a['global_signed_sum']}",
        f"- affine_orbit_count={a['affine_orbit_count']}",
        f"- affine_orbit_sizes={a['affine_orbit_sizes']}",
        f"- nonzero_orbit_coefficient_count={a['nonzero_orbit_coefficient_count']}",
        f"- positive_nonzero_coefficients={a['positive_nonzero_coefficients']}",
        f"- negative_nonzero_coefficients={a['negative_nonzero_coefficients']}",
        f"- histogram_abs_paired={a['histogram_abs_paired']}",
        "",
        "## Theorem statement",
        a["finite_theorem"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2750 = read_json(P2750)
    scan = evidence_scan()
    audit = quartic_phase_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2751_QUARTIC_PHASE_POLARITY_SOURCE_GAP",
        "input_hashes": {"P2750_CONCRETE_ODD_SOURCE_SIGN_VALUE_INVENTORY_NO_GO": sha(P2750)},
        "input_statuses": {"P2750_CONCRETE_ODD_SOURCE_SIGN_VALUE_INVENTORY_NO_GO": p2750.get("status")},
        "audited_candidate_class": "new quartic Z12 phase-sum signed observable after P2750 inventory no-go",
        "content_evidence_scan": scan,
        "quartic_phase_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the P2751 quartic phase observable to lambda/P2721 fixing.  It is a genuine new finite signed observable with nonzero affine-orbit coefficients, but its nonzero coefficient family remains polarity-paired and lacks a strict orbit/polarity source plus P2721 coupling theorem.  The next proof-grade move should either audit a new strict law selecting one quartic coefficient orbit/polarity with explicit P2721 coupling, or pivot to a different typed object; if neither is available, preserve the P2697-P2751 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2751/S1701 Z12 quartic phase orbit signed observable audit", "## P2751/S1701 Z12 quartic phase orbit signed observable audit\n\n`P2751/S1701` pivots beyond the P2750 current-artifact inventory to quartic `Z12` phase sums `Q(a,b,c,d)=sum_n exp(2*pi*i*(a*n^4+b*n^3+c*n^2+d*n)/12)` under affine source reparametrisation.  The finite audit computes all `12^4` coefficient quadruples and finds nonzero orbit-safe signed coefficients after affine quotienting, but the nonzero coefficient family remains polarity-paired.  No strict quartic coefficient-orbit source, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2751/S1701 quartic phase Ltotal guard", "## P2751/S1701 quartic phase Ltotal guard\n\n`P2751/S1701` adds no variational source term.  Although quartic `Z12` phase sums produce nonzero affine-orbit signed coefficients, the family remains polarity-paired and current artifacts do not export a strict coefficient-orbit/polarity selector law or `P2721` coupling theorem.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Z12 quartic phase orbit signed observable guardrail (P2751/S1701, 2026-06-15)", "## Current Z12 quartic phase orbit signed observable guardrail (P2751/S1701, 2026-06-15)\n\n- P2751 pivots beyond the P2750 current-artifact inventory to quartic `Z12` phase sums `Q(a,b,c,d)=sum_n exp(2*pi*i*(a*n^4+b*n^3+c*n^2+d*n)/12)` under affine source reparametrisation.\n- The finite audit computes all `12^4` coefficient quadruples and finds nonzero orbit-safe signed coefficients after affine quotienting, but the nonzero coefficient family remains polarity-paired and no strict quartic orbit/polarity source is exported.\n- Do not promote quartic phase coefficient choice to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a strict law selecting one nonzero quartic coefficient orbit/polarity and an explicit `P2721` coupling theorem.  Otherwise pivot to a different typed object or preserve the P2697-P2751 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
