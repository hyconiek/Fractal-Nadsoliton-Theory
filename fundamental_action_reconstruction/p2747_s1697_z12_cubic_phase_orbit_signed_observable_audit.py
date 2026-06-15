#!/usr/bin/env python3
"""P2747/S1697: Z12 cubic phase orbit signed observable audit.

P2746 closes the current Gauss coefficient-orbit selector attempt unless a new
strict sign law is supplied.  This audit pivots to a different typed object:
the imaginary sign of cubic Z12 phase sums, quotiented by affine source
reparametrisation.
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

GEN = ROOT / "generated"
P2746 = GEN / "p2746_s1696_gauss_phase_orbit_selector_source_law_no_go.json"
OUT = GEN / "p2747_s1697_z12_cubic_phase_orbit_signed_observable_audit.json"
MD = GEN / "p2747_s1697_z12_cubic_phase_orbit_signed_observable_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

N = 12
UNITS = (1, 5, 7, 11)
TOL = 1e-9
CONTENT_PATTERNS = {
    "post_p2746_pivot_obligation": r"P2746|pivot away from Gauss-phase selectors|genuinely new strict sign law|no-new-live-frontier",
    "cubic_phase_boundary": r"cubic phase|higher-phase|phase sum|signed observable|coefficient orbit",
    "affine_source_gauge_boundary": r"affine|source-gauge|orbit-safe signed|translation-gauge|polarity",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "cubic_phase_promoted_to_strict_source",
    "cubic_coefficient_orbit_selected_nonpremise",
    "cubic_polarity_selected_nonpremise",
    "p2721_polarity_coupling_exported",
    "lambda_fixed",
    "qw2191_discharged",
    "selector_closure_exported",
    "role_transfer_started",
    "ltotal_promoted",
    "toe_closure_exported",
]


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else {"missing": True, "path": rel(path)}


def sha(path: Path) -> str | None:
    return hashlib.sha256(path.read_bytes()).hexdigest() if path.exists() else None


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


def cubic_phase_sum(a: int, b: int, c: int) -> complex:
    return sum(cmath.exp(2j * math.pi * ((a * n**3 + b * n * n + c * n) % N) / N) for n in range(N))


def sign(x: float) -> int:
    return 1 if x > TOL else -1 if x < -TOL else 0


def affine_image(a: int, b: int, c: int, u: int, t: int) -> tuple[int, int, int]:
    # Pull back q(n)=a*n^3+b*n^2+c*n by n -> u*n+t, dropping the constant term.
    return (
        (a * u**3) % N,
        (u * u * (b + 3 * a * t)) % N,
        (u * (c + 2 * b * t + 3 * a * t * t)) % N,
    )


def cubic_phase_audit() -> dict[str, Any]:
    rows = []
    signs: dict[tuple[int, int, int], int] = {}
    for a in range(N):
        for b in range(N):
            for c in range(N):
                z = cubic_phase_sum(a, b, c)
                s = sign(z.imag)
                signs[(a, b, c)] = s
                rows.append({"a": a, "b": b, "c": c, "real": round(z.real, 12), "imag": round(z.imag, 12), "imag_sign": s})
    seen: set[tuple[int, int, int]] = set()
    orbits = []
    for triple in sorted(signs):
        if triple in seen:
            continue
        orbit = {affine_image(*triple, u, t) for u in UNITS for t in range(N)}
        seen |= orbit
        coeff = sum(signs[p] for p in orbit)
        orbits.append({
            "orbit_id": len(orbits),
            "size": len(orbit),
            "signed_sum_coefficient": coeff,
            "positive_rows": sum(1 for p in orbit if signs[p] == 1),
            "negative_rows": sum(1 for p in orbit if signs[p] == -1),
            "zero_rows": sum(1 for p in orbit if signs[p] == 0),
            "representatives_sample": sorted([list(p) for p in orbit])[:12],
        })
    coeffs = sorted(o["signed_sum_coefficient"] for o in orbits if o["signed_sum_coefficient"] != 0)
    histogram: dict[str, int] = {}
    for coeff in coeffs:
        histogram[str(coeff)] = histogram.get(str(coeff), 0) + 1
    return {
        "typed_candidate": "imaginary-sign cubic phase sum C(a,b,c)=sum_n exp(2*pi*i*(a*n^3+b*n^2+c*n)/12), quotiented by affine source reparametrisation",
        "modulus": N,
        "coefficient_triple_count": len(rows),
        "pointwise_positive_signs": sum(1 for r in rows if r["imag_sign"] == 1),
        "pointwise_negative_signs": sum(1 for r in rows if r["imag_sign"] == -1),
        "pointwise_zero_signs": sum(1 for r in rows if r["imag_sign"] == 0),
        "global_signed_sum": sum(signs.values()),
        "affine_orbit_count": len(orbits),
        "affine_orbit_sizes": sorted({o["size"] for o in orbits}),
        "nonzero_orbit_coefficient_count": len(coeffs),
        "nonzero_orbit_coefficients": coeffs,
        "nonzero_orbit_coefficient_histogram": histogram,
        "positive_nonzero_coefficients": sum(1 for c in coeffs if c > 0),
        "negative_nonzero_coefficients": sum(1 for c in coeffs if c < 0),
        "orbit_rows": orbits,
        "finite_theorem": "Cubic Z12 phase sums are a real pivot away from the P2746 Gauss-selector gap: 1728 coefficient triples produce 396 positive, 396 negative, and 936 zero imaginary signs.  Affine quotienting gives 180 coefficient orbits, with 44 nonzero signed-sum coefficients.  The nonzero coefficients are globally polarity-balanced by counts and values (-8/-4/-2/-1 paired with +1/+2/+4/+8), so the object supplies orbit-safe signed coefficients but still no strict source selecting one orbit/polarity or coupling it to P2721.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2746_pivot": scan["all_patterns_have_hits"],
        "new_cubic_phase_candidate_computed": audit["coefficient_triple_count"] == 1728,
        "has_nonzero_orbit_safe_signed_coefficients": audit["nonzero_orbit_coefficient_count"] > 0,
        "positive_negative_nonzero_counts_balanced": audit["positive_nonzero_coefficients"] == audit["negative_nonzero_coefficients"],
        "strict_cubic_orbit_source_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_lambda_p2721_source": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The cubic phase observable has many nonzero affine-orbit coefficients, but they remain polarity-balanced and current artifacts export no strict cubic coefficient-orbit source or P2721 coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["cubic_phase_audit"]
    lines = [
        "# P2747/S1697 Z12 cubic phase orbit signed observable audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite cubic-phase audit",
        f"- coefficient_triple_count={a['coefficient_triple_count']}",
        f"- pointwise signs: +{a['pointwise_positive_signs']} / -{a['pointwise_negative_signs']} / 0={a['pointwise_zero_signs']}",
        f"- global_signed_sum={a['global_signed_sum']}",
        f"- affine_orbit_count={a['affine_orbit_count']}",
        f"- affine_orbit_sizes={a['affine_orbit_sizes']}",
        f"- nonzero_orbit_coefficient_count={a['nonzero_orbit_coefficient_count']}",
        f"- positive_nonzero_coefficients={a['positive_nonzero_coefficients']}",
        f"- negative_nonzero_coefficients={a['negative_nonzero_coefficients']}",
        f"- nonzero_orbit_coefficient_histogram={a['nonzero_orbit_coefficient_histogram']}",
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
    p2746 = read_json(P2746)
    scan = evidence_scan()
    audit = cubic_phase_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2747_CUBIC_PHASE_POLARITY_SOURCE_GAP",
        "input_hashes": {"P2746_GAUSS_SELECTOR_NO_GO": sha(P2746)},
        "input_statuses": {"P2746_GAUSS_SELECTOR_NO_GO": p2746.get("status")},
        "audited_candidate_class": "Z12 cubic phase orbit signed observable after P2746 Gauss selector no-go",
        "content_evidence_scan": scan,
        "cubic_phase_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the P2747 cubic phase observable to lambda/P2721 fixing.  It gives many nonzero affine-orbit signed coefficients, but the coefficient family is still polarity-balanced and lacks a strict orbit/polarity source plus P2721 coupling theorem.  The next proof-grade move should audit exactly one missing premise: a strict cubic coefficient-orbit/polarity selector law with explicit P2721 coupling; if no such law is available, pivot to a different typed object or preserve the P2697-P2747 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2747/S1697 Z12 cubic phase orbit signed observable audit", "## P2747/S1697 Z12 cubic phase orbit signed observable audit\n\n`P2747/S1697` pivots away from the P2746 Gauss-selector gap to cubic `Z12` phase sums `C(a,b,c)=sum_n exp(2*pi*i*(a*n^3+b*n^2+c*n)/12)` under affine source reparametrisation.  The finite audit finds `1728` coefficient triples with `396` positive, `396` negative, and `936` zero imaginary signs; affine quotienting gives `180` coefficient orbits and `44` nonzero signed-sum coefficients.  The nonzero coefficient family remains polarity-balanced, and no strict cubic coefficient-orbit source, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2747/S1697 cubic phase Ltotal guard", "## P2747/S1697 cubic phase Ltotal guard\n\n`P2747/S1697` adds no variational source term.  Although cubic `Z12` phase sums produce nonzero affine-orbit signed coefficients, the family is polarity-balanced and current artifacts do not export a strict coefficient-orbit/polarity selector law or `P2721` coupling theorem.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Z12 cubic phase orbit signed observable guardrail (P2747/S1697, 2026-06-14)", "## Current Z12 cubic phase orbit signed observable guardrail (P2747/S1697, 2026-06-14)\n\n- P2747 pivots away from the P2746 Gauss-selector gap to cubic `Z12` phase sums `C(a,b,c)=sum_n exp(2*pi*i*(a*n^3+b*n^2+c*n)/12)` under affine source reparametrisation.\n- The finite computation finds `1728` coefficient triples with `396` positive, `396` negative, and `936` zero imaginary signs; affine quotienting gives `180` coefficient orbits and `44` nonzero signed-sum coefficients, but the nonzero family remains polarity-balanced.\n- Do not promote cubic phase coefficient choice to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a strict law selecting one nonzero cubic coefficient orbit/polarity and an explicit `P2721` coupling theorem.  Otherwise pivot to a different typed object or preserve the P2697-P2747 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
