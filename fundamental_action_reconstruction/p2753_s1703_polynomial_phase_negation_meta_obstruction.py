#!/usr/bin/env python3
"""P2753/S1703: polynomial phase negation meta-obstruction.

P2752 closed the quartic coefficient-orbit selector attempt by an explicit
q -> -q pairing.  This audit prevents an unproductive loop over quintic/sextic
phase sums by proving and finitely checking the general obstruction for the
whole imaginary-sign polynomial phase family: coefficient negation conjugates
the phase sum, flips imaginary sign, and therefore pairs every nonzero sign
unless a new strict law breaks the negation involution and couples to P2721.
"""
from __future__ import annotations

import cmath
import hashlib
import itertools
import json
import math
import subprocess
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

N = 12
MAX_DEGREE = 5
TOL = 1e-9
GEN = ROOT / "generated"
P2752 = GEN / "p2752_s1702_quartic_phase_negation_pairing_selector_no_go.json"
OUT = GEN / "p2753_s1703_polynomial_phase_negation_meta_obstruction.json"
MD = GEN / "p2753_s1703_polynomial_phase_negation_meta_obstruction.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2752_loop_cut": r"P2752|negation pairing|different typed object|P2697-P2752",
    "polynomial_phase_lane": r"quartic phase|cubic phase|Gauss-phase|phase sum|coefficient orbit",
    "p2721_boundary": r"P2721 coupling|P2721 polarity|lambda/P2721|coupling theorem",
    "closure_forbidden": r"QW-2191|selector closure|role transfer|L_total|ToE closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "polynomial_phase_lane_promoted_to_selector",
    "negation_involution_broken",
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


def phase_sum(coeffs: tuple[int, ...]) -> complex:
    degree = len(coeffs)
    total = 0j
    for n in range(N):
        exponent = sum(coeffs[j - 1] * (n**j) for j in range(1, degree + 1)) % N
        total += cmath.exp(2j * math.pi * exponent / N)
    return total


def neg_coeffs(coeffs: tuple[int, ...]) -> tuple[int, ...]:
    return tuple((-c) % N for c in coeffs)


def finite_degree_audit() -> dict[str, Any]:
    degree_rows = []
    total_vectors = 0
    total_failures = 0
    for degree in range(1, MAX_DEGREE + 1):
        positive = negative = zero = failures = 0
        sample_failures = []
        sample_nonzero_pairs = []
        for coeffs in itertools.product(range(N), repeat=degree):
            z = phase_sum(coeffs)
            s = sign(z.imag)
            ns = sign(phase_sum(neg_coeffs(coeffs)).imag)
            if s == 1:
                positive += 1
            elif s == -1:
                negative += 1
            else:
                zero += 1
            if ns != -s:
                failures += 1
                if len(sample_failures) < 8:
                    sample_failures.append({"coefficients": list(coeffs), "sign": s, "neg_coefficients": list(neg_coeffs(coeffs)), "neg_sign": ns})
            if s != 0 and len(sample_nonzero_pairs) < 6:
                sample_nonzero_pairs.append({"coefficients": list(coeffs), "sign": s, "neg_coefficients": list(neg_coeffs(coeffs)), "neg_sign": ns})
        count = N**degree
        total_vectors += count
        total_failures += failures
        degree_rows.append({
            "degree": degree,
            "coefficient_vector_count": count,
            "positive_imag_signs": positive,
            "negative_imag_signs": negative,
            "zero_imag_signs": zero,
            "global_signed_sum": positive - negative,
            "sign_flip_failure_count": failures,
            "sample_failures": sample_failures,
            "sample_nonzero_pairs": sample_nonzero_pairs,
            "positive_negative_balanced": positive == negative,
        })
    return {
        "modulus": N,
        "max_degree_finitely_checked": MAX_DEGREE,
        "total_coefficient_vectors_checked": total_vectors,
        "total_sign_flip_failures": total_failures,
        "degree_rows": degree_rows,
        "all_degrees_positive_negative_balanced": all(r["positive_negative_balanced"] for r in degree_rows),
        "all_degrees_sign_flip_verified": total_failures == 0,
        "general_symbolic_obstruction": "For any coefficient vector q over Z12 and polynomial phase sum S(q)=sum_n exp(2*pi*i*q(n)/12), coefficient negation gives S(-q)=conj(S(q)).  Hence Im(S(-q))=-Im(S(q)).  Any phase-sum selector rule invariant under availability of q and -q cannot choose a nonzero polarity without an added strict law that breaks q <-> -q and an explicit P2721 coupling theorem.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2752_loop_cut": scan["all_patterns_have_hits"],
        "finite_degrees_1_to_5_checked": audit["max_degree_finitely_checked"] == 5,
        "all_checked_degrees_sign_flip_verified": audit["all_degrees_sign_flip_verified"],
        "all_checked_degrees_positive_negative_balanced": audit["all_degrees_positive_negative_balanced"],
        "new_strict_negation_breaking_law_exported": False,
        "p2721_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_polynomial_phase_selector_source": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The polynomial phase-sum imaginary-sign lane is closed under coefficient negation, which conjugates the sum and pairs every nonzero sign; no strict negation-breaking law or P2721 coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["polynomial_phase_negation_audit"]
    lines = [
        "# P2753/S1703 polynomial phase negation meta-obstruction",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite degree audit",
        f"- modulus={a['modulus']}",
        f"- max_degree_finitely_checked={a['max_degree_finitely_checked']}",
        f"- total_coefficient_vectors_checked={a['total_coefficient_vectors_checked']}",
        f"- total_sign_flip_failures={a['total_sign_flip_failures']}",
        f"- all_degrees_positive_negative_balanced={a['all_degrees_positive_negative_balanced']}",
        f"- all_degrees_sign_flip_verified={a['all_degrees_sign_flip_verified']}",
        "",
        "## Degree rows",
    ]
    for row in a["degree_rows"]:
        lines.append(
            f"- degree={row['degree']}: count={row['coefficient_vector_count']}, +={row['positive_imag_signs']}, -={row['negative_imag_signs']}, 0={row['zero_imag_signs']}, failures={row['sign_flip_failure_count']}"
        )
    lines.extend(["", "## Theorem statement", a["general_symbolic_obstruction"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2752 = read_json(P2752)
    scan = evidence_scan()
    audit = finite_degree_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2753_POLYNOMIAL_PHASE_NEGATION_META_OBSTRUCTION",
        "input_hashes": {"P2752_QUARTIC_PHASE_NEGATION_PAIRING_SELECTOR_NO_GO": sha(P2752)},
        "input_statuses": {"P2752_QUARTIC_PHASE_NEGATION_PAIRING_SELECTOR_NO_GO": p2752.get("status")},
        "audited_candidate_class": "meta-obstruction for imaginary-sign polynomial phase-sum selector attempts under coefficient negation",
        "content_evidence_scan": scan,
        "polynomial_phase_negation_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue the same polynomial phase-sum lane by merely increasing degree.  P2753 proves the coefficient-negation obstruction for the whole imaginary-sign polynomial phase family and finitely verifies degrees 1 through 5 with zero sign-flip failures.  The next proof-grade move must either introduce a genuinely new strict negation-breaking source law with explicit P2721 coupling, or pivot outside polynomial phase-sum imaginary-sign observables; otherwise preserve the P2697-P2753 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2753/S1703 polynomial phase negation meta-obstruction", "## P2753/S1703 polynomial phase negation meta-obstruction\n\n`P2753/S1703` cuts the post-P2752 loop of simply raising polynomial degree.  For any coefficient vector `q`, the imaginary-sign polynomial phase sum satisfies `S(-q)=conj(S(q))`, so coefficient negation flips the imaginary sign.  The finite audit checks all coefficient vectors for degrees `1` through `5` over `Z12` with zero sign-flip failures and balanced positive/negative counts in every degree.  Thus the polynomial phase-sum imaginary-sign lane exports no strict negation-breaking source law, no `P2721` coupling theorem, no `lambda/P2721` fixing, no `QW-2191` discharge, selector closure, role transfer, `L_total`, or ToE closure.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2753/S1703 polynomial phase negation Ltotal guard", "## P2753/S1703 polynomial phase negation Ltotal guard\n\n`P2753/S1703` adds no variational source term.  It proves the polynomial phase-sum imaginary-sign lane remains paired by coefficient negation unless a new strict negation-breaking law is exported with a `P2721` coupling theorem.  Therefore it does not promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure.\n")
    append_once(AGENTS, "Current polynomial phase negation meta-obstruction guardrail (P2753/S1703, 2026-06-15)", "## Current polynomial phase negation meta-obstruction guardrail (P2753/S1703, 2026-06-15)\n\n- P2753 cuts the post-P2752 loop of merely increasing polynomial degree in `Z12` phase-sum imaginary-sign observables.\n- For any coefficient vector `q`, coefficient negation gives `S(-q)=conj(S(q))` and flips the imaginary sign; the finite audit checks degrees `1` through `5` with zero sign-flip failures and balanced positive/negative counts in every degree.\n- Do not continue polynomial phase-sum degree escalation as `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure without a genuinely new strict negation-breaking source law and an explicit `P2721` coupling theorem.  Otherwise pivot outside polynomial phase-sum imaginary-sign observables or preserve the P2697-P2753 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
