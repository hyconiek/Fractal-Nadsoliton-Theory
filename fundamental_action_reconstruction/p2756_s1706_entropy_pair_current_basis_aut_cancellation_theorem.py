#!/usr/bin/env python3
"""P2756/S1706: entropy pair-current basis Aut-cancellation theorem.

P2755 used one concrete directed entropy current.  This follow-up removes the
"toy-current" weakness by auditing the whole finite class of translation-summed
nearest-neighbour entropy pair currents on the four-quanta Z12 entropy alphabet.
Any antisymmetric local pair law A(h_i,h_j) is a linear combination of unordered
entropy-level pair basis functions.  The directed basis can be nonzero, but the
opposite directed step negates every basis coefficient, so Aut(Z12)-symmetric
handling over steps {1,5,7,11} cancels the entire vector space, not merely the
single P2755 formula.
"""
from __future__ import annotations

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
QUANTA = 4
UNITS = (1, 5, 7, 11)
OPPOSITE = {1: 11, 11: 1, 5: 7, 7: 5}
GEN = ROOT / "generated"
P2755 = GEN / "p2755_s1705_entropy_gradient_current_aut_cancellation_audit.json"
OUT = GEN / "p2756_s1706_entropy_pair_current_basis_aut_cancellation_theorem.json"
MD = GEN / "p2756_s1706_entropy_pair_current_basis_aut_cancellation_theorem.md"
STRICT_EQUATION_SHEET = ROOT / "STRICT_EQUATION_SHEET_KERNEL_TO_LTOTAL_CURRENT.md"
STRICT_LAGRANGIAN_DRAFT = ROOT / "STRICT_KERNEL_LAGRANGIAN_AND_EOM_DRAFT.md"
AGENTS = REPO / "AGENTS.md"

CONTENT_PATTERNS = {
    "post_p2755_toy_current_boundary": r"P2755|entropy-current|directed entropy current|Aut-averaged entropy-current|P2697-P2755",
    "selector_boundary": r"P2721|lambda/P2721|QW-2191|orientation torsor|selector closure",
    "aut_z12_boundary": r"Aut\(Z12\)|inversion|steps `1/11`|steps `1/11` and `5/7`|orientation-reversing",
    "closure_forbidden": r"role transfer|L_total|ToE closure|bridge closure",
}
NEGATIVE_EXPORT_FLAGS = [
    "entropy_pair_current_basis_promoted_to_selector_source",
    "aut_averaged_pair_current_basis_nonzero",
    "strict_pair_current_step_or_polarity_selector_exported",
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


def entropy_level_for_count(c: int) -> int:
    """Five-level entropy alphabet induced by counts 0..4.

    The actual values are e_c=-(c/4)log(c/4), with c=0 treated as 0.  For the
    finite pair-current basis only the level labels matter; c=1 and c=3 are
    distinct labels because the local probability values differ even if future
    degeneracies are inspected separately.
    """
    return c


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


def rank_int(rows: list[tuple[int, ...]]) -> int:
    """Small exact rational rank via fraction-free Gaussian elimination."""
    from fractions import Fraction

    mat = [[Fraction(x) for x in row] for row in rows if any(row)]
    if not mat:
        return 0
    m = len(mat)
    n = len(mat[0])
    rank = 0
    col = 0
    while rank < m and col < n:
        pivot = None
        for r in range(rank, m):
            if mat[r][col] != 0:
                pivot = r
                break
        if pivot is None:
            col += 1
            continue
        mat[rank], mat[pivot] = mat[pivot], mat[rank]
        pv = mat[rank][col]
        mat[rank] = [x / pv for x in mat[rank]]
        for r in range(m):
            if r != rank and mat[r][col] != 0:
                factor = mat[r][col]
                mat[r] = [x - factor * y for x, y in zip(mat[r], mat[rank])]
        rank += 1
        col += 1
    return rank


def basis_audit() -> dict[str, Any]:
    basis = pair_basis()
    all_directed_rows: list[tuple[int, ...]] = []
    all_aut_rows: list[tuple[int, ...]] = []
    nonzero_by_step = {u: 0 for u in UNITS}
    opposite_failures = []
    aut_failures = []
    sample_nonzero = []
    for counts in compositions(QUANTA, N):
        levels = tuple(entropy_level_for_count(c) for c in counts)
        vectors = {u: feature_vector(levels, u, basis) for u in UNITS}
        for u, vec in vectors.items():
            all_directed_rows.append(vec)
            if any(vec):
                nonzero_by_step[u] += 1
                if len(sample_nonzero) < 8:
                    sample_nonzero.append({"counts": list(counts), "step": u, "feature_vector": list(vec), "opposite_step": OPPOSITE[u], "opposite_feature_vector": list(vectors[OPPOSITE[u]])})
        for u in (1, 5):
            paired_sum = tuple(a + b for a, b in zip(vectors[u], vectors[OPPOSITE[u]]))
            if any(paired_sum) and len(opposite_failures) < 8:
                opposite_failures.append({"counts": list(counts), "step": u, "paired_sum": list(paired_sum)})
        aut_sum = tuple(sum(vectors[u][j] for u in UNITS) for j in range(len(basis)))
        all_aut_rows.append(aut_sum)
        if any(aut_sum) and len(aut_failures) < 8:
            aut_failures.append({"counts": list(counts), "aut_sum": list(aut_sum), "vectors": {str(k): list(v) for k, v in vectors.items()}})
    directed_rank = rank_int(all_directed_rows)
    aut_rank = rank_int(all_aut_rows)
    return {
        "modulus": N,
        "quanta": QUANTA,
        "composition_count": math.comb(QUANTA + N - 1, N - 1),
        "entropy_level_count": QUANTA + 1,
        "antisymmetric_pair_basis_dimension": len(basis),
        "basis_pairs": [list(pair) for pair in basis],
        "directed_feature_rank": directed_rank,
        "aut_averaged_feature_rank": aut_rank,
        "nonzero_directed_feature_counts_by_step": nonzero_by_step,
        "opposite_pair_failure_count": len(opposite_failures),
        "sample_opposite_failures": opposite_failures,
        "aut_sum_failure_count": len(aut_failures),
        "sample_aut_failures": aut_failures,
        "sample_nonzero_directed_features": sample_nonzero,
        "all_opposite_steps_cancel_as_vectors": len(opposite_failures) == 0,
        "aut_averaged_pair_current_basis_identically_zero": len(aut_failures) == 0 and aut_rank == 0,
        "theorem_statement": "For every translation-summed nearest-neighbour antisymmetric entropy pair law A(h_i,h_j), the finite basis vector for step -u is the negative of the basis vector for step u.  Since Aut(Z12) includes opposite unit pairs 1/11 and 5/7, the selector-free Aut-summed pair-current vector is zero for every coefficient choice in the full antisymmetric pair-current basis.",
    }


def acceptance_matrix(scan: dict[str, Any], audit: dict[str, Any]) -> dict[str, Any]:
    facts = {
        "content_grep_detects_post_p2755_boundary": scan["all_patterns_have_hits"],
        "full_antisymmetric_pair_basis_audited": audit["antisymmetric_pair_basis_dimension"] == math.comb(QUANTA + 1, 2),
        "directed_pair_current_basis_nontrivial": audit["directed_feature_rank"] > 0,
        "opposite_step_vector_pairing_verified": audit["all_opposite_steps_cancel_as_vectors"],
        "aut_averaged_pair_current_basis_identically_zero": audit["aut_averaged_pair_current_basis_identically_zero"],
        "strict_step_or_polarity_selector_exported": False,
        "p2721_pair_current_coupling_theorem_exported": False,
    }
    return {
        "facts": facts,
        "accepted_as_strict_entropy_pair_current_selector": False,
        "missing_criteria": [k for k, v in facts.items() if not v],
        "blocker": "The whole antisymmetric entropy pair-current basis is nontrivial only directionally.  Aut(Z12)-symmetric summing cancels every basis vector, and no strict step/polarity selector or P2721 coupling theorem is exported.",
    }


def write_markdown(payload: dict[str, Any]) -> None:
    a = payload["entropy_pair_current_basis_audit"]
    lines = [
        "# P2756/S1706 entropy pair-current basis Aut-cancellation theorem",
        "",
        f"Status: `{payload['status']}`",
        "",
        "## Finite basis audit",
        f"- modulus={a['modulus']}",
        f"- quanta={a['quanta']}",
        f"- composition_count={a['composition_count']}",
        f"- entropy_level_count={a['entropy_level_count']}",
        f"- antisymmetric_pair_basis_dimension={a['antisymmetric_pair_basis_dimension']}",
        f"- directed_feature_rank={a['directed_feature_rank']}",
        f"- aut_averaged_feature_rank={a['aut_averaged_feature_rank']}",
        f"- opposite_pair_failure_count={a['opposite_pair_failure_count']}",
        f"- aut_sum_failure_count={a['aut_sum_failure_count']}",
        "",
        "## Nonzero directed feature counts by step",
    ]
    for step, count in a["nonzero_directed_feature_counts_by_step"].items():
        lines.append(f"- step={step}: nonzero_feature_vectors={count}")
    lines.extend(["", "## Theorem", a["theorem_statement"], "", "## Recommendation", payload["decision"]["next_honest_step"]])
    MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main() -> dict[str, Any]:
    GEN.mkdir(exist_ok=True)
    p2755 = read_json(P2755)
    scan = evidence_scan()
    audit = basis_audit()
    acceptance = acceptance_matrix(scan, audit)
    payload = {
        "status": "P2756_ENTROPY_PAIR_CURRENT_BASIS_AUT_CANCELLATION_THEOREM",
        "input_hashes": {"P2755_ENTROPY_GRADIENT_CURRENT_AUT_CANCELLATION_AUDIT_NO_GO": sha(P2755)},
        "input_statuses": {"P2755_ENTROPY_GRADIENT_CURRENT_AUT_CANCELLATION_AUDIT_NO_GO": p2755.get("status")},
        "audited_candidate_class": "full translation-summed nearest-neighbour antisymmetric Shannon entropy pair-current basis on four-quanta Z12 entropy alphabet",
        "content_evidence_scan": scan,
        "entropy_pair_current_basis_audit": audit,
        "acceptance_matrix": acceptance,
        "decision": {
            "negative_export_flags": {flag: False for flag in NEGATIVE_EXPORT_FLAGS},
            "reason": acceptance["blocker"],
            "next_honest_step": "Do not search this same local entropy pair-current class by changing the antisymmetric function A: P2756 audits the full finite antisymmetric pair-current basis and proves Aut(Z12) cancellation for every coefficient choice.  The next proof-grade move must either export a strict law selecting a directed entropy-current step/polarity with explicit P2721 coupling, or pivot outside local entropy pair-current observables to a genuinely new typed object; otherwise preserve the P2697-P2756 no-new-live-frontier certificate.",
        },
    }
    OUT.write_text(json.dumps(payload, indent=2, sort_keys=True, ensure_ascii=False) + "\n", encoding="utf-8")
    write_markdown(payload)
    append_once(STRICT_EQUATION_SHEET, "P2756/S1706 entropy pair-current basis Aut-cancellation theorem", "## P2756/S1706 entropy pair-current basis Aut-cancellation theorem\n\n`P2756/S1706` strengthens P2755 from one toy entropy current to the full finite class of translation-summed nearest-neighbour antisymmetric entropy pair currents on the four-quanta `Z12` entropy alphabet.  The directed pair-current basis is nontrivial, but the finite audit verifies that the step `-u` vector is always the negative of the step `u` vector; because `Aut(Z12)` pairs `1/11` and `5/7`, the Aut-summed pair-current basis has rank `0` with zero opposite-pair and Aut-sum failures.  Thus changing the antisymmetric local pair law cannot export `P2721` coupling, `lambda/P2721` fixing, `QW-2191` discharge, selector closure, role transfer, `L_total`, bridge closure, or ToE closure without a new strict step/polarity source.\n")
    append_once(STRICT_LAGRANGIAN_DRAFT, "P2756/S1706 entropy pair-current basis Ltotal guard", "## P2756/S1706 entropy pair-current basis Ltotal guard\n\n`P2756/S1706` adds no variational source term.  It proves that the full local antisymmetric entropy pair-current basis cancels under selector-free `Aut(Z12)` handling, so varying the local entropy-pair function cannot promote role-bearing `L_total`, selector closure, bridge closure, role transfer, or ToE closure without a strict step/polarity source and explicit `P2721` coupling theorem.\n")
    append_once(AGENTS, "Current entropy pair-current basis Aut-cancellation guardrail (P2756/S1706, 2026-06-15)", "## Current entropy pair-current basis Aut-cancellation guardrail (P2756/S1706, 2026-06-15)\n\n- P2756 removes the P2755 toy-current limitation by auditing the full finite basis of translation-summed nearest-neighbour antisymmetric Shannon entropy pair currents on the four-quanta `Z12` entropy alphabet.\n- The directed basis is nontrivial, but every opposite step vector cancels (`1/11`, `5/7`), and the selector-free `Aut(Z12)`-summed pair-current basis has rank `0` with zero opposite-pair and Aut-sum failures.\n- Do not continue this lane by changing the antisymmetric local entropy-pair function as `lambda/P2721` fixing, `QW-2191` discharge, selector closure, bridge closure, role transfer, `L_total`, or ToE closure.  A next admissible move must export a strict step/polarity source with explicit `P2721` coupling, pivot outside local entropy pair-current observables, or preserve the P2697-P2756 no-new-live-frontier certificate.\n")
    return payload


if __name__ == "__main__":
    main()
