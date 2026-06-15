#!/usr/bin/env python3
"""P2745/S1695: Z12 quadratic Gauss-phase signed observable audit.

P2744 leaves a pivot obligation outside finite Z12 sign-character/frame/cycle-spectrum
observables.  This audit tests a genuinely different finite object: the imaginary
phase sign of quadratic Gauss sums over Z12 characters, then quotients coefficient
pairs by affine source reparametrisation.
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
OUT = GEN / "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.json"
MD = GEN / "p2745_s1695_z12_quadratic_gauss_phase_signed_observable_audit.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"
INPUTS = {"P2744_SPECTRAL_ASYMMETRY_NO_GO": GEN / "p2744_s1694_z12_cycle_spectral_asymmetry_no_go.json"}

N = 12
UNITS = (1, 5, 7, 11)
TOL = 1e-9
CONTENT_PATTERNS = {
    "post_p2744_pivot_obligation": r"P2744|pivot outside finite `?Z12`? sign-character/frame/cycle-spectrum|nonzero spectral asymmetry|no-new-live-frontier",
    "gauss_phase_boundary": r"Gauss|quadratic phase|phase sign|character twist|signed observable",
    "affine_source_gauge_boundary": r"affine|source-gauge|orbit-safe signed value|translation-gauge|representative choice",
    "p2721_closure_boundary": r"P2721 polarity|lambda/P2721|QW-2191|selector closure|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "quadratic_gauss_phase_promoted_to_strict_source",
    "gauss_coefficient_orbit_selected_nonpremise",
    "gauss_polarity_selected_nonpremise",
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


def gauss_sum(a: int, b: int) -> complex:
    return sum(cmath.exp(2j * math.pi * ((a * n * n + b * n) % N) / N) for n in range(N))


def sign(x: float) -> int:
    return 1 if x > TOL else -1 if x < -TOL else 0


def affine_image(a: int, b: int, u: int, t: int) -> tuple[int, int]:
    # Pull back q(n)=a*n^2+b*n by n -> u*n+t, dropping the constant term.
    return ((a * u * u) % N, (u * (b + 2 * a * t)) % N)


def gauss_phase_audit() -> dict[str, Any]:
    rows = []
    signs: dict[tuple[int, int], int] = {}
    for a in range(N):
        for b in range(N):
            z = gauss_sum(a, b)
            s = sign(z.imag)
            signs[(a, b)] = s
            rows.append({"a": a, "b": b, "real": round(z.real, 12), "imag": round(z.imag, 12), "imag_sign": s})
    seen: set[tuple[int, int]] = set()
    orbits = []
    for pair in sorted(signs):
        if pair in seen:
            continue
        orbit = {affine_image(pair[0], pair[1], u, t) for u in UNITS for t in range(N)}
        seen |= orbit
        coeff = sum(signs[p] for p in orbit)
        orbits.append({"orbit_id": len(orbits), "size": len(orbit), "signed_sum_coefficient": coeff, "positive_rows": sum(1 for p in orbit if signs[p] == 1), "negative_rows": sum(1 for p in orbit if signs[p] == -1), "zero_rows": sum(1 for p in orbit if signs[p] == 0), "representatives": sorted([list(p) for p in orbit])[:12]})
    nonzero = [o for o in orbits if o["signed_sum_coefficient"] != 0]
    paired_coefficients = sorted(o["signed_sum_coefficient"] for o in nonzero)
    return {
        "typed_candidate": "imaginary-sign quadratic Gauss phase G(a,b)=sum_n exp(2*pi*i*(a*n^2+b*n)/12), quotiented by affine source reparametrisation",
        "modulus": N,
        "coefficient_pair_count": len(rows),
        "pointwise_positive_signs": sum(1 for r in rows if r["imag_sign"] == 1),
        "pointwise_negative_signs": sum(1 for r in rows if r["imag_sign"] == -1),
        "pointwise_zero_signs": sum(1 for r in rows if r["imag_sign"] == 0),
        "affine_orbit_count": len(orbits),
        "affine_orbit_sizes": sorted({o["size"] for o in orbits}),
        "nonzero_orbit_coefficient_count": len(nonzero),
        "nonzero_orbit_coefficients": paired_coefficients,
        "positive_nonzero_coefficients": sum(1 for c in paired_coefficients if c > 0),
        "negative_nonzero_coefficients": sum(1 for c in paired_coefficients if c < 0),
        "global_signed_sum": sum(signs.values()),
        "pointwise_rows_sample": rows[:16],
        "orbit_rows": orbits,
        "finite_theorem": "Quadratic Gauss phases over Z12 are a real pivot outside the P2744 cycle-spectrum test: 144 coefficient pairs produce 20 positive, 20 negative, and 104 zero imaginary signs.  Affine quotienting gives 40 coefficient orbits, with 8 nonzero signed-sum coefficients appearing in opposite polarities [-2,-2,-1,-1,1,1,2,2].  Thus the object supplies orbit-safe signed coefficients, but only as an unselected polarity family; current artifacts still export no strict law choosing a coefficient orbit/sign and no P2721 coupling theorem.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_pivot_boundaries": scan["all_patterns_have_hits"],
        "new_candidate_computed": audit["coefficient_pair_count"] == 144,
        "has_nonzero_orbit_safe_signed_coefficients": audit["nonzero_orbit_coefficient_count"] > 0,
        "polarity_family_unpaired_by_strict_source": False,
        "strict_coefficient_orbit_source_exported": False,
        "p2721_polarity_coupling_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_lambda_p2721_source": all(facts.values()),
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "P2745 finds nonzero affine-orbit Gauss-phase signed coefficients, but they occur in opposite polarities and no current strict artifact chooses a coefficient orbit/sign or couples it to P2721.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["gauss_phase_audit"]
    lines = [
        "# P2745/S1695 Z12 quadratic Gauss-phase signed observable audit",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite Gauss-phase audit",
        f"- coefficient_pair_count={a['coefficient_pair_count']}",
        f"- pointwise signs: +{a['pointwise_positive_signs']} / -{a['pointwise_negative_signs']} / 0={a['pointwise_zero_signs']}",
        f"- affine_orbit_count={a['affine_orbit_count']}",
        f"- affine_orbit_sizes={a['affine_orbit_sizes']}",
        f"- nonzero_orbit_coefficient_count={a['nonzero_orbit_coefficient_count']}",
        f"- nonzero_orbit_coefficients={a['nonzero_orbit_coefficients']}",
        f"- global_signed_sum={a['global_signed_sum']}",
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
    inputs = {key: read_json(path) for key, path in INPUTS.items()}
    scan = evidence_scan()
    audit = gauss_phase_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2745_GAUSS_PHASE_POLARITY_SOURCE_GAP" if not acceptance["accepted_as_lambda_p2721_source"] else "P2745_REQUIRES_MANUAL_REVIEW",
        "input_hashes": {key: sha(path) for key, path in INPUTS.items()},
        "input_statuses": {key: value.get("status") for key, value in inputs.items()},
        "audited_candidate_class": "Z12 quadratic Gauss-phase signed observable after P2744 spectral-asymmetry no-go",
        "content_evidence_scan": scan,
        "gauss_phase_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not promote the P2745 quadratic Gauss-phase observable to lambda/P2721 fixing yet.  It is stronger than P2744 in that affine quotienting leaves 8 nonzero signed orbit coefficients, but those coefficients are polarity-paired and current artifacts export no strict coefficient-orbit/sign source or P2721 coupling theorem.  The next proof-grade move should audit exactly one missing premise: a strict law selecting one nonzero Gauss coefficient orbit and polarity with an explicit P2721 coupling theorem; if no such law is available, preserve the P2697-P2745 no-new-live-frontier certificate or pivot to a different typed object.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2745/S1695 Z12 quadratic Gauss-phase signed observable audit", "## P2745/S1695 Z12 quadratic Gauss-phase signed observable audit\n\n`P2745/S1695` pivots outside finite `Z12` sign-character/frame/cycle-spectrum observables to the imaginary phase sign of quadratic Gauss sums `G(a,b)=sum_n exp(2*pi*i*(a*n^2+b*n)/12)`.  The finite audit finds `144` coefficient pairs with `20` positive, `20` negative, and `104` zero imaginary signs; affine quotienting gives `40` coefficient orbits and `8` nonzero signed-sum coefficients `[-2,-2,-1,-1,1,1,2,2]`.  This is a real orbit-safe signed observable family, but no strict coefficient-orbit/sign source, `P2721` polarity-coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2745/S1695 Gauss-phase Ltotal guard", "## P2745/S1695 Gauss-phase Ltotal guard\n\n`P2745/S1695` adds no variational source term.  Although the quadratic Gauss-phase audit finds nonzero affine-orbit signed coefficients, current artifacts do not export a strict law selecting one coefficient orbit/polarity or coupling it to `P2721`.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current Z12 quadratic Gauss-phase signed observable guardrail (P2745/S1695, 2026-06-14)", "## Current Z12 quadratic Gauss-phase signed observable guardrail (P2745/S1695, 2026-06-14)\n\n- P2745 pivots outside finite `Z12` sign-character/frame/cycle-spectrum observables to the imaginary phase sign of quadratic Gauss sums `G(a,b)=sum_n exp(2*pi*i*(a*n^2+b*n)/12)` under affine source reparametrisation.\n- The finite computation finds `144` coefficient pairs with `20` positive, `20` negative, and `104` zero imaginary signs; affine quotienting gives `40` coefficient orbits and `8` nonzero signed-sum coefficients `[-2,-2,-1,-1,1,1,2,2]`.\n- Do not promote this Gauss-phase observable family to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a strict law selecting one nonzero coefficient orbit/polarity and an explicit `P2721` coupling theorem.  A next admissible move is exactly that missing-premise audit, or preservation of the P2697-P2745 no-new-live-frontier certificate/pivot to a different typed object.\n")
    return payload


if __name__ == "__main__":
    main()
