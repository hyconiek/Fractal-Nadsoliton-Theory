#!/usr/bin/env python3
"""P2757/S1707: entropy pair-current polarity-selector invariant no-go.

P2756 proved Aut-cancellation for the whole local antisymmetric entropy
pair-current basis unless a strict law selects a directed step/polarity.  This
audits that exact missing premise inside the same data: can intrinsic,
polarity-safe signatures of the pair-current feature vector choose one of the
opposite steps?  No.  Opposite steps have feature vectors v and -v, so every
sign-blind/internal signature is identical on the pair.  Sign-sensitive choices
can choose, but only by importing the missing orientation/polarity premise.
"""
from __future__ import annotations

import hashlib
import json
import math
import subprocess
from functools import reduce
from pathlib import Path
from typing import Any

from p2468_s1418_strict_pointwise_interval_decimal_chunked_opposite_tail_replay_ledger_certificate import REPO, ROOT, rel
from p2506_s1456_strict_damping_rg_minimum_roughness_selector_candidate_certificate import append_once

N = 12
QUANTA = 4
UNITS = (1, 5, 7, 11)
OPPOSITE = {1: 11, 11: 1, 5: 7, 7: 5}
GEN = ROOT / "generated"
P2756 = GEN / "p2756_s1706_entropy_pair_current_basis_aut_cancellation_theorem.json"
OUT = GEN / "p2757_s1707_entropy_pair_current_polarity_selector_invariant_no_go.json"
MD = GEN / "p2757_s1707_entropy_pair_current_polarity_selector_invariant_no_go.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2756_step_polarity_obligation": r"P2756|strict step/polarity source|entropy pair-current|P2697-P2756",
    "selector_boundary": r"P2721|lambda/P2721|QW-2191|orientation torsor|selector closure",
    "polarity_boundary": r"polarity|opposite step|1/11|5/7|Aut\(Z12\)",
    "closure_forbidden": r"role transfer|L_total|ToE closure|bridge closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "entropy_pair_current_polarity_selector_exported",
    "sign_blind_signature_selects_opposite_step",
    "sign_sensitive_choice_accepted_without_premise",
    "p2721_pair_current_coupling_theorem_exported",
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


def compositions(total: int, parts: int):
    if parts == 1:
        yield (total,)
        return
    for first in range(total + 1):
        for rest in compositions(total - first, parts - 1):
            yield (first,) + rest


def pair_basis() -> list[tuple[int, int]]:
    return [(a, b) for a in range(QUANTA + 1) for b in range(a + 1, QUANTA + 1)]


def feature_vector(levels: tuple[int, ...], step: int, basis: list[tuple[int, int]]) -> tuple[int, ...]:
    index = {pair: k for k, pair in enumerate(basis)}
    vec = [0] * len(basis)
    for i, a in enumerate(levels):
        b = levels[(i + step) % N]
        if a == b:
            continue
        lo, hi = (a, b) if a < b else (b, a)
        vec[index[(lo, hi)]] += 1 if (a, b) == (lo, hi) else -1
    return tuple(vec)


def gcd_abs(values: tuple[int, ...]) -> int:
    vals = [abs(v) for v in values if v]
    return reduce(math.gcd, vals) if vals else 0


def sign_blind_signature(vec: tuple[int, ...]) -> tuple[Any, ...]:
    absvals = tuple(abs(v) for v in vec)
    return (
        sum(absvals),
        sum(v * v for v in vec),
        sum(1 for v in vec if v),
        max(absvals) if absvals else 0,
        gcd_abs(vec),
        tuple(sorted(absvals, reverse=True)),
    )


def selector_audit() -> dict[str, Any]:
    basis = pair_basis()
    pair_count = 0
    nonzero_pair_count = 0
    sign_blind_signature_failures = []
    negation_failures = []
    sign_sensitive_lex_choices = {u: 0 for u in UNITS}
    sample_sign_sensitive_choices = []
    for counts in compositions(QUANTA, N):
        levels = tuple(counts)
        vectors = {u: feature_vector(levels, u, basis) for u in UNITS}
        for u in (1, 5):
            pair_count += 1
            v = vectors[u]
            ov = vectors[OPPOSITE[u]]
            if any(v) or any(ov):
                nonzero_pair_count += 1
            if tuple(-x for x in v) != ov and len(negation_failures) < 8:
                negation_failures.append({"counts": list(counts), "step": u, "vector": list(v), "opposite_step": OPPOSITE[u], "opposite_vector": list(ov)})
            if sign_blind_signature(v) != sign_blind_signature(ov) and len(sign_blind_signature_failures) < 8:
                sign_blind_signature_failures.append({"counts": list(counts), "step": u, "signature": sign_blind_signature(v), "opposite_step": OPPOSITE[u], "opposite_signature": sign_blind_signature(ov)})
            if v != ov and (any(v) or any(ov)):
                chosen = u if v > ov else OPPOSITE[u]
                sign_sensitive_lex_choices[chosen] += 1
                if len(sample_sign_sensitive_choices) < 8:
                    sample_sign_sensitive_choices.append({"counts": list(counts), "chosen_by_lexicographic_signed_vector": chosen, "step": u, "vector": list(v), "opposite_step": OPPOSITE[u], "opposite_vector": list(ov)})
    return {
        "modulus": N,
        "quanta": QUANTA,
        "composition_count": math.comb(QUANTA + N - 1, N - 1),
        "opposite_step_pair_rows": pair_count,
        "nonzero_opposite_step_pair_rows": nonzero_pair_count,
        "basis_dimension": len(basis),
        "negation_failure_count": len(negation_failures),
        "sample_negation_failures": negation_failures,
        "sign_blind_signature_failure_count": len(sign_blind_signature_failures),
        "sample_sign_blind_signature_failures": sign_blind_signature_failures,
        "sign_sensitive_lexicographic_choice_counts": sign_sensitive_lex_choices,
        "sample_sign_sensitive_lexicographic_choices": sample_sign_sensitive_choices,
        "all_opposite_vectors_are_negatives": len(negation_failures) == 0,
        "all_sign_blind_signatures_match_opposite_steps": len(sign_blind_signature_failures) == 0,
        "sign_sensitive_rules_can_choose_only_by_using_vector_sign": any(sign_sensitive_lex_choices.values()),
        "theorem_statement": "Because opposite entropy pair-current feature vectors satisfy v(-u)=-v(u), every intrinsic sign-blind signature F(v)=F(-v) gives equal scores to opposite steps.  A sign-sensitive rule can choose one side only by using the very polarity that P2721 requires as a strict source, so it is a premise rather than a selector export.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2756_polarity_obligation": scan["all_patterns_have_hits"],
        "opposite_vectors_are_exact_negatives": audit["all_opposite_vectors_are_negatives"],
        "sign_blind_signatures_match_all_opposite_steps": audit["all_sign_blind_signatures_match_opposite_steps"],
        "sign_sensitive_rules_are_premise_dependent": audit["sign_sensitive_rules_can_choose_only_by_using_vector_sign"],
        "strict_step_or_polarity_selector_exported": False,
        "p2721_pair_current_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_entropy_pair_current_polarity_selector": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "All sign-blind/internal signatures give equal scores to opposite entropy-current steps; sign-sensitive rules can choose only by importing the missing polarity premise, and no P2721 coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["entropy_pair_current_polarity_selector_audit"]
    lines = [
        "# P2757/S1707 entropy pair-current polarity-selector invariant no-go",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite polarity-selector audit",
        f"- modulus={a['modulus']}",
        f"- quanta={a['quanta']}",
        f"- composition_count={a['composition_count']}",
        f"- opposite_step_pair_rows={a['opposite_step_pair_rows']}",
        f"- nonzero_opposite_step_pair_rows={a['nonzero_opposite_step_pair_rows']}",
        f"- basis_dimension={a['basis_dimension']}",
        f"- negation_failure_count={a['negation_failure_count']}",
        f"- sign_blind_signature_failure_count={a['sign_blind_signature_failure_count']}",
        f"- sign_sensitive_lexicographic_choice_counts={a['sign_sensitive_lexicographic_choice_counts']}",
        "",
        "## Theorem",
        a["theorem_statement"],
        "",
        "## Recommendation",
        payload["decision"]["next_honest_step"],
    ]
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2756 = read_json(P2756)
    scan = evidence_scan()
    audit = selector_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2757_ENTROPY_PAIR_CURRENT_POLARITY_SELECTOR_INVARIANT_NO_GO",
        "input_hashes": {"P2756_ENTROPY_PAIR_CURRENT_BASIS_AUT_CANCELLATION_THEOREM": sha(P2756)},
        "input_statuses": {"P2756_ENTROPY_PAIR_CURRENT_BASIS_AUT_CANCELLATION_THEOREM": p2756.get("status")},
        "audited_candidate_class": "intrinsic polarity selector signatures for the P2756 entropy pair-current basis",
        "content_evidence_scan": scan,
        "entropy_pair_current_polarity_selector_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not continue the local entropy pair-current lane by adding sign-blind scoring signatures or sign-sensitive tie-breaks.  P2757 proves that sign-blind/internal signatures cannot distinguish opposite steps, while sign-sensitive rules choose only by importing the missing polarity premise.  The next proof-grade move must pivot outside local entropy pair-current observables to a genuinely new strict typed object/source, or export an independent strict step/polarity source with explicit P2721 coupling; otherwise preserve the P2697-P2757 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2757/S1707 entropy pair-current polarity-selector invariant no-go", "## P2757/S1707 entropy pair-current polarity-selector invariant no-go\n\n`P2757/S1707` audits the exact missing premise after P2756: whether intrinsic signatures of the entropy pair-current feature vector can select a directed step/polarity.  For all `2730` opposite-step rows in the finite four-quanta `Z12` scan, opposite vectors are exact negatives and every sign-blind signature tested from absolute/support/norm data matches across the pair.  Sign-sensitive lexicographic rules can choose, but only by using the missing vector sign/polarity as a premise.  Thus no strict step/polarity source, `P2721` coupling theorem, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, bridge closure, or ToE closure is exported.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2757/S1707 entropy pair-current polarity-selector Ltotal guard", "## P2757/S1707 entropy pair-current polarity-selector Ltotal guard\n\n`P2757/S1707` adds no variational source term.  It shows that sign-blind entropy pair-current signatures cannot select a directed polarity, while sign-sensitive tie-breaks import the missing premise.  Therefore it cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure without an independent strict step/polarity source and explicit `P2721` coupling theorem.\n")
    append_once(AGENTS, "Current entropy pair-current polarity-selector invariant no-go guardrail (P2757/S1707, 2026-06-15)", "## Current entropy pair-current polarity-selector invariant no-go guardrail (P2757/S1707, 2026-06-15)\n\n- P2757 audits the exact missing premise after P2756: whether intrinsic signatures of the entropy pair-current feature vector can select a directed step/polarity.\n- Opposite step vectors are exact negatives on every finite row, every sign-blind/internal signature matches across opposite steps, and sign-sensitive rules choose only by importing the missing polarity premise.\n- Do not continue local entropy pair-current observables via sign-blind scoring or sign-sensitive tie-breaks as `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must pivot outside local entropy pair-current observables to a genuinely new strict typed object/source, export an independent strict step/polarity source with explicit `P2721` coupling, or preserve the P2697-P2757 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
