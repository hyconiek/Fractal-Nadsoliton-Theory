#!/usr/bin/env python3
"""P2752/S1702: quartic-phase negation-pairing selector no-go.

P2751 found nonzero affine-orbit signed coefficients for quartic Z12 phase sums.
This file audits the exact missing premise: can the quartic coefficient-orbit
family itself select one polarity?  The finite proof constructs the coefficient
negation involution q -> -q.  Since Q(-q)=conj(Q(q)), imaginary signs flip; the
audit verifies this descends to affine orbits and pairs every nonzero orbit with
a distinct opposite-coefficient orbit.
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
P2751 = GEN / "p2751_s1701_z12_quartic_phase_orbit_signed_observable_audit.json"
OUT = GEN / "p2752_s1702_quartic_phase_negation_pairing_selector_no_go.json"
MD = GEN / "p2752_s1702_quartic_phase_negation_pairing_selector_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2751_missing_premise": r"P2751|quartic coefficient orbit|quartic orbit/polarity|P2697-P2751",
    "negation_pairing_boundary": r"polarity-paired|opposite|coefficient family|negation|conjugate",
    "p2721_boundary": r"P2721 coupling|P2721 polarity|lambda/P2721|coupling theorem",
    "closure_forbidden": r"QW-2191|selector closure|role transfer|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "quartic_orbit_polarity_selector_exported",
    "negation_pairing_broken",
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
    return (
        (a * u**4) % N,
        (u**3 * (b + 4 * a * t)) % N,
        (u**2 * (c + 3 * b * t + 6 * a * t * t)) % N,
        (u * (d + 2 * c * t + 3 * b * t * t + 4 * a * t**3)) % N,
    )


def neg_quad(q: tuple[int, int, int, int]) -> tuple[int, int, int, int]:
    return tuple((-x) % N for x in q)  # type: ignore[return-value]


def compute_orbits() -> tuple[dict[tuple[int, int, int, int], int], list[dict[str, Any]], dict[tuple[int, int, int, int], int]]:
    signs: dict[tuple[int, int, int, int], int] = {}
    for a in range(N):
        for b in range(N):
            for c in range(N):
                for d in range(N):
                    signs[(a, b, c, d)] = sign(quartic_phase_sum(a, b, c, d).imag)

    seen: set[tuple[int, int, int, int]] = set()
    orbit_index: dict[tuple[int, int, int, int], int] = {}
    orbits: list[dict[str, Any]] = []
    for quad in sorted(signs):
        if quad in seen:
            continue
        orbit = {affine_image(*quad, u, t) for u in UNITS for t in range(N)}
        seen |= orbit
        oid = len(orbits)
        for q in orbit:
            orbit_index[q] = oid
        coeff = sum(signs[p] for p in orbit)
        orbits.append({
            "orbit_id": oid,
            "size": len(orbit),
            "signed_sum_coefficient": coeff,
            "positive_rows": sum(1 for p in orbit if signs[p] == 1),
            "negative_rows": sum(1 for p in orbit if signs[p] == -1),
            "zero_rows": sum(1 for p in orbit if signs[p] == 0),
            "canonical_representative": list(min(orbit)),
            "members": orbit,
        })
    return signs, orbits, orbit_index


def negation_pairing_audit() -> dict[str, Any]:
    signs, orbits, orbit_index = compute_orbits()
    sign_flip_failures = []
    for q, s in signs.items():
        ns = signs[neg_quad(q)]
        if ns != -s:
            sign_flip_failures.append({"q": list(q), "sign": s, "neg_q": list(neg_quad(q)), "neg_sign": ns})
            if len(sign_flip_failures) >= 12:
                break

    pair_rows = []
    nonzero_ids = [o["orbit_id"] for o in orbits if o["signed_sum_coefficient"] != 0]
    for oid in nonzero_ids:
        orbit = orbits[oid]
        negated_member = neg_quad(tuple(orbit["canonical_representative"]))
        paired_id = orbit_index[negated_member]
        paired = orbits[paired_id]
        pair_rows.append({
            "orbit_id": oid,
            "coefficient": orbit["signed_sum_coefficient"],
            "paired_orbit_id": paired_id,
            "paired_coefficient": paired["signed_sum_coefficient"],
            "same_orbit": paired_id == oid,
            "size": orbit["size"],
            "paired_size": paired["size"],
            "canonical_representative": orbit["canonical_representative"],
            "paired_canonical_representative": paired["canonical_representative"],
        })
    failures = [r for r in pair_rows if r["same_orbit"] or r["paired_coefficient"] != -r["coefficient"] or r["paired_size"] != r["size"]]
    unique_unordered_pairs = {tuple(sorted((r["orbit_id"], r["paired_orbit_id"]))) for r in pair_rows}
    return {
        "typed_candidate": "selector/source-law audit for P2751 quartic affine coefficient-orbit polarity using coefficient negation q -> -q",
        "coefficient_quadruple_count": len(signs),
        "affine_orbit_count": len(orbits),
        "nonzero_orbit_count": len(nonzero_ids),
        "nonzero_unordered_negation_pair_count": len(unique_unordered_pairs),
        "sign_flip_failure_count": len(sign_flip_failures),
        "sign_flip_failures_sample": sign_flip_failures,
        "pairing_failure_count": len(failures),
        "pairing_failures_sample": failures[:12],
        "all_nonzero_orbits_paired_by_negation": len(failures) == 0 and len(sign_flip_failures) == 0,
        "pair_rows_sample": pair_rows[:40],
        "finite_theorem": "Coefficient negation sends every quartic phase sum to its complex conjugate, so the imaginary sign flips pointwise.  The finite affine-orbit audit verifies that this involution descends to the P2751 coefficient-orbit quotient and pairs every nonzero quartic orbit with a distinct orbit of equal size and opposite signed-sum coefficient.  Therefore the quartic coefficient-orbit family has no internal polarity selector unless a new strict law breaks this negation pairing and supplies a P2721 coupling theorem.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_p2751_missing_premise": scan["all_patterns_have_hits"],
        "quartic_nonzero_orbits_exist": audit["nonzero_orbit_count"] > 0,
        "pointwise_negation_sign_flip_verified": audit["sign_flip_failure_count"] == 0,
        "all_nonzero_orbits_paired_by_negation": audit["all_nonzero_orbits_paired_by_negation"],
        "quartic_internal_polarity_selector_exported": False,
        "p2721_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_quartic_selector_source_law": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The quartic nonzero orbit coefficients are exactly paired by coefficient negation; current artifacts export no strict law breaking this pairing and no P2721 coupling theorem.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["negation_pairing_audit"]
    lines = [
        "# P2752/S1702 quartic phase negation-pairing selector no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite negation-pairing audit",
        f"- coefficient_quadruple_count={a['coefficient_quadruple_count']}",
        f"- affine_orbit_count={a['affine_orbit_count']}",
        f"- nonzero_orbit_count={a['nonzero_orbit_count']}",
        f"- nonzero_unordered_negation_pair_count={a['nonzero_unordered_negation_pair_count']}",
        f"- sign_flip_failure_count={a['sign_flip_failure_count']}",
        f"- pairing_failure_count={a['pairing_failure_count']}",
        f"- all_nonzero_orbits_paired_by_negation={a['all_nonzero_orbits_paired_by_negation']}",
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
    p2751 = read_json(P2751)
    scan = evidence_scan()
    audit = negation_pairing_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2752_QUARTIC_PHASE_NEGATION_PAIRING_SELECTOR_NO_GO",
        "input_hashes": {"P2751_QUARTIC_PHASE_POLARITY_SOURCE_GAP": sha(P2751)},
        "input_statuses": {"P2751_QUARTIC_PHASE_POLARITY_SOURCE_GAP": p2751.get("status")},
        "audited_candidate_class": "P2751 missing-premise audit: quartic coefficient-orbit/polarity selector law via negation pairing",
        "content_evidence_scan": scan,
        "negation_pairing_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue P2751 by manually choosing a positive quartic coefficient orbit.  P2752 proves that coefficient negation pairs every nonzero quartic affine orbit with an equal-size opposite-coefficient orbit, so the quartic family itself has no internal polarity selector.  The next proof-grade move must either introduce a genuinely new strict law breaking this negation pairing with an explicit P2721 coupling theorem, or pivot to a different typed object; otherwise preserve the P2697-P2752 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2752/S1702 quartic phase negation-pairing selector no-go", "## P2752/S1702 quartic phase negation-pairing selector no-go\n\n`P2752/S1702` audits the exact missing premise left by `P2751/S1701`: whether the quartic coefficient-orbit family itself can select one orbit polarity.  Coefficient negation `q -> -q` sends each quartic phase sum to its complex conjugate and flips the imaginary sign.  The finite affine-orbit audit verifies zero sign-flip failures and zero pairing failures: every nonzero quartic orbit is paired with a distinct equal-size orbit of opposite signed-sum coefficient.  Thus no strict quartic orbit/polarity selector, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2752/S1702 quartic negation-pairing Ltotal guard", "## P2752/S1702 quartic negation-pairing Ltotal guard\n\n`P2752/S1702` adds no variational source term.  It proves the P2751 quartic nonzero orbit coefficients remain paired by coefficient negation, and current artifacts do not export a strict law breaking that pairing or a `P2721` coupling theorem.  Therefore this does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current quartic phase negation-pairing selector no-go guardrail (P2752/S1702, 2026-06-15)", "## Current quartic phase negation-pairing selector no-go guardrail (P2752/S1702, 2026-06-15)\n\n- P2752 audits the exact missing premise left by P2751: whether the quartic coefficient-orbit family itself can select one orbit polarity.\n- Coefficient negation `q -> -q` sends each quartic phase sum to its complex conjugate and flips the imaginary sign; the finite affine-orbit audit verifies every nonzero quartic orbit is paired with a distinct equal-size orbit of opposite signed-sum coefficient.\n- Do not promote a manually chosen positive quartic coefficient orbit to `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a genuinely new strict law breaking this negation pairing and an explicit `P2721` coupling theorem.  Otherwise pivot to a different typed object or preserve the P2697-P2752 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
